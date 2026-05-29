#!/usr/bin/env bash
# reconstruct_fusions.sh — reference-free reconstruction of fusion junctions from a
# GEAC fusion *evidence BAM* (records tagged FX:Z:GENEA::GENEB by
# `geac experimental fusions --reads-output`).
#
# For each fusion (distinct FX tag) it:
#   1. extracts that fusion's reads into a per-fusion BAM,
#   2. emits a FASTA (mates kept adjacent via samtools collate),
#   3. assembles them de novo with CAP3 — reference-free,
#   4. optionally aligns the contigs to a reference with minimap2, so the breakpoint
#      shows up as a split / supplementary alignment (the contig maps partly to one
#      partner gene and partly to the other).
#
# A good --reference (-r) is a small FASTA of just the two partner gene bodies; the
# full genome works too but the junction is easier to read off a 2-gene mini-ref.
#
# Usage:
#   scripts/reconstruct_fusions.sh -b evidence.bam [-o outdir] [-r ref.fa] \
#                                  [-f "GENEA::GENEB"] [-t threads]
#
#   -b  evidence BAM (required)
#   -o  output directory (default: fusion_reconstruction)
#   -r  reference FASTA to align contigs against (optional; needs minimap2)
#   -f  process only this one fusion label (default: every FX tag in the BAM)
#   -t  threads (default: 4)
#
# Requires: samtools (>=1.12, for `view -d TAG:VALUE`), cap3.
#           minimap2 is only needed when -r is given.
#
# Note (macOS/arm64): cap3 is x86-only; install via conda's osx-64 subdir:
#   CONDA_SUBDIR=osx-64 conda create -n cap3 -c bioconda -c conda-forge cap3
#   conda activate cap3   # then run this script

set -euo pipefail

usage() { sed -n '2,32p' "$0"; exit "${1:-1}"; }

BAM=""; OUTDIR="fusion_reconstruction"; REF=""; ONLY_FX=""; THREADS=4
while getopts ":b:o:r:f:t:h" opt; do
  case "$opt" in
    b) BAM=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    r) REF=$OPTARG ;;
    f) ONLY_FX=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    h) usage 0 ;;
    *) usage 1 ;;
  esac
done

[[ -n "$BAM" ]] || { echo "error: -b <evidence.bam> is required" >&2; usage 1; }
[[ -f "$BAM" ]] || { echo "error: BAM not found: $BAM" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || { echo "error: '$1' not found in PATH${2:-}" >&2; exit 1; }; }
need samtools ""
need cap3 " — install via conda (see the macOS/arm64 note in the header)"
[[ -n "$REF" ]] && need minimap2 " — required when -r is given"

mkdir -p "$OUTDIR"

# Sanitize an FX label (e.g. "EWSR1::FLI1") into a filesystem-safe token.
sanitize() { printf '%s' "$1" | tr -c 'A-Za-z0-9._-' '_'; }

# Discover the set of fusions to process (one per line).
fxlist=$(mktemp)
trap 'rm -f "$fxlist"' EXIT
if [[ -n "$ONLY_FX" ]]; then
  printf '%s\n' "$ONLY_FX" > "$fxlist"
else
  samtools view "$BAM" \
    | grep -oE 'FX:Z:[^[:space:]]+' \
    | sed 's/^FX:Z://' \
    | sort -u > "$fxlist"
fi

n_fusions=$(grep -c . "$fxlist" || true)
[[ "$n_fusions" -gt 0 ]] || {
  echo "no FX:Z: tags found in $BAM — was it produced by 'geac experimental fusions --reads-output'?" >&2
  exit 1
}
echo "fusions to process: $n_fusions"

summary="$OUTDIR/summary.tsv"
printf 'fusion\treads\tcontigs\tmax_contig_len\tbreakpoint\n' > "$summary"

while IFS= read -r fx; do
  [[ -n "$fx" ]] || continue
  tok=$(sanitize "$fx")
  d="$OUTDIR/$tok"
  mkdir -p "$d"
  echo ">> $fx"

  # 1. Exact FX-tag subset (samtools splits TAG:VALUE on the first ':', so the "::"
  #    inside the gene-pair label is preserved as part of the value).
  samtools view -b -d "FX:$fx" "$BAM" > "$d/reads.bam"
  nreads=$(samtools view -c "$d/reads.bam")

  # 2. FASTA for assembly (CAP3 takes FASTA; collate first so paired mates are
  #    adjacent in the file).
  samtools collate -@ "$THREADS" -Ou "$d/reads.bam" \
    | samtools fasta -@ "$THREADS" - > "$d/reads.fa" 2>/dev/null

  # 3. Reference-free de novo assembly. CAP3 writes <input>.cap.contigs next to
  #    the input FASTA; symlink it to a stable contigs.fa for downstream steps.
  ncontigs=0; maxlen=0
  if [[ "$nreads" -gt 0 ]]; then
    if cap3 "$d/reads.fa" > "$d/cap3.log" 2>&1; then
      ln -sf "reads.fa.cap.contigs" "$d/contigs.fa"
      ncontigs=$(grep -c '^>' "$d/contigs.fa" || true)
      maxlen=$(awk '/^>/{next}{print length($0)}' "$d/contigs.fa" | sort -rn | head -1)
      maxlen=${maxlen:-0}
    else
      echo "   cap3 failed (see $d/cap3.log)" >&2
    fi
  fi

  # 4. Optional: align contigs to reference; a contig spanning the junction aligns
  #    chimerically (primary + supplementary), and the split point is the breakpoint.
  #    asm5 = assembled contigs, high identity to the same-species reference.
  bp="n/a"
  if [[ -n "$REF" && "$ncontigs" -gt 0 ]]; then
    minimap2 -ax asm5 -t "$THREADS" "$REF" "$d/contigs.fa" 2>/dev/null \
      | samtools sort -@ "$THREADS" -o "$d/contigs.bam" -
    samtools index "$d/contigs.bam"

    # Read the breakpoint off each junction-spanning (chimeric) contig. We look
    # only at primary records (-F 2304 drops secondary+supplementary) that carry
    # an SA:Z: tag: the primary maps to one partner, the SA segment to the other,
    # and each partner's breakpoint is the reference coordinate at the clipped end
    # of its CIGAR (alignment start if the clip leads, alignment end if it trails).
    bpfile="$d/breakpoint.txt"
    printf 'contig\tpartner_A\tpartner_B\n' > "$bpfile"
    samtools view -F 2304 "$d/contigs.bam" \
      | awk -F'\t' 'BEGIN{OFS="\t"}
        # reference bases consumed by a CIGAR (M/D/N/=/X)
        function refspan(c,   i,ch,n,op,sp){ sp=0; n=0
          for(i=1;i<=length(c);i++){ ch=substr(c,i,1)
            if(ch>="0"&&ch<="9"){ n=n*10+ch }
            else{ op=ch; if(op=="M"||op=="D"||op=="N"||op=="="||op=="X") sp+=n; n=0 } }
          return sp }
        function leadclip(c){ return (c ~ /^[0-9]+[SH]/) ? 1 : 0 }
        {
          sa=""
          for(i=12;i<=NF;i++) if($i ~ /^SA:Z:/) sa=substr($i,6)
          if(sa=="") next
          pbrk = leadclip($6) ? $4 : $4 + refspan($6)            # primary partner
          split(sa,segs,";"); split(segs[1],f,",")               # first SA segment
          sbrk = leadclip(f[4]) ? f[2] : f[2] + refspan(f[4])    # other partner
          print $1, $3":"pbrk, f[1]":"sbrk
        }' >> "$bpfile"

    # Summarize: first chimeric contig as "A|B", or note none was found.
    if [[ $(grep -c . "$bpfile") -gt 1 ]]; then
      bp=$(awk -F'\t' 'NR==2{print $2"|"$3}' "$bpfile")
    else
      bp="none"
    fi
  fi

  printf '%s\t%s\t%s\t%s\t%s\n' "$fx" "$nreads" "$ncontigs" "$maxlen" "$bp" >> "$summary"
  echo "   reads=$nreads contigs=$ncontigs max_contig_len=$maxlen breakpoint=$bp -> $d"
done < "$fxlist"

echo "done. per-fusion outputs under $OUTDIR/, summary: $summary"
