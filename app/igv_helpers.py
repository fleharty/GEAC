"""Pure-Python / DuckDB helpers for IGV session generation.

Kept separate from geac_explorer.py so they can be unit-tested without
importing the Streamlit runtime or any visualisation libraries.
"""

from __future__ import annotations

import duckdb
import pandas as pd


def make_vcf(df: pd.DataFrame) -> str:
    """Build a multi-sample VCF string from a DataFrame of alt_bases rows.

    Coordinate convention:
        GEAC stores 0-based positions; VCF is 1-based, so POS = pos + 1.
        ref_allele is always the single anchor base.

    Allele encoding:
        SNV       alt_allele = "T"     → REF=ref_allele,           ALT=alt_allele
        deletion  alt_allele = "-ACG"  → REF=ref_allele+"ACG",     ALT=ref_allele
        insertion alt_allele = "+ACGT" → REF=ref_allele,           ALT=ref_allele+"ACGT"

    Output: one row per unique (CHROM, POS, REF, ALT), one FORMAT/sample column
    per distinct sample_id (sorted).  Samples absent at a locus are filled with
    the missing-data sentinel "./.:.:.:.".
    """
    _HEADER = "\n".join([
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt alleles">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fraction">',
    ])

    if df.empty:
        return _HEADER + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"

    samples = sorted(df["sample_id"].unique().tolist())

    def _vcf_alleles(ref: str, alt: str) -> tuple[str, str]:
        if alt.startswith("-"):
            return ref + alt[1:], ref
        if alt.startswith("+"):
            return ref, ref + alt[1:]
        return ref, alt

    # Build a format string for each row, then pivot to (variant_key × sample).
    tmp = df.copy()
    tmp["_fmt"] = (
        "./."
        + ":" + tmp["total_depth"].astype(int).astype(str)
        + ":" + tmp["ref_count"].astype(int).astype(str) + "," + tmp["alt_count"].astype(int).astype(str)
        + ":" + tmp["vaf"].astype(str)
    )

    # Deduplicate: if the same sample appears twice at the same variant key, keep first.
    tmp = tmp.drop_duplicates(subset=["chrom", "pos", "ref_allele", "alt_allele", "sample_id"])

    pivot = tmp.pivot(
        index=["chrom", "pos", "ref_allele", "alt_allele"],
        columns="sample_id",
        values="_fmt",
    ).reindex(columns=samples).fillna("./.:.:.:.").reset_index()

    alleles = pivot.apply(
        lambda r: _vcf_alleles(r["ref_allele"], r["alt_allele"]), axis=1
    )
    pivot["vcf_pos"] = pivot["pos"].astype(int) + 1
    pivot["vcf_ref"] = [a[0] for a in alleles]
    pivot["vcf_alt"] = [a[1] for a in alleles]
    pivot = pivot.sort_values(["chrom", "vcf_pos"])

    col_header = "\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
    )
    rows = []
    for _, r in pivot.iterrows():
        sample_cols = [str(r[s]) for s in samples]
        rows.append("\t".join(
            [r["chrom"], str(int(r["vcf_pos"])), ".", r["vcf_ref"], r["vcf_alt"], ".", ".", ".", "GT:DP:AD:VAF"]
            + sample_cols
        ))

    return _HEADER + "\n" + col_header + "\n" + "\n".join(rows) + "\n"


def insert_size_active_part(
    is_lo: int, is_hi: int, is_min: int, is_max: int
) -> str | None:
    """Return the insert-size banner fragment when the filter is non-default, else None.

    Kept separate from geac_explorer.py so it can be unit-tested without
    importing Streamlit.
    """
    if is_lo > is_min or is_hi < is_max:
        return f"insert size: {is_lo}–{is_hi}"
    return None


def per_read_warning_note(recompute_vaf: bool) -> str:
    """Return the mode-appropriate body text for the per-read filter warning banner.

    Kept separate from geac_explorer.py so it can be unit-tested without
    importing Streamlit.
    """
    if recompute_vaf:
        return (
            "alt_count is re-aggregated from reads passing the filter; "
            "original_vaf shows the unfiltered VAF for comparison. "
            "ref_count, total_depth, and strand/overlap columns still reflect "
            "unfiltered locus-level values."
        )
    else:
        return (
            "Loci with no alt reads passing these filters are hidden. "
            "alt_count and VAF reflect all reads."
        )


def query_distinct_samples(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where: str,
) -> list[str]:
    """Return sorted list of distinct sample_ids matching *where* in *table_expr*.

    Must query the database directly rather than inspecting a display-limited
    DataFrame: the IGV cap warning must reflect the full dataset regardless of
    how many rows are shown in the UI.
    """
    return (
        con.execute(
            f"SELECT DISTINCT sample_id FROM {table_expr} WHERE {where} ORDER BY sample_id"
        )
        .df()["sample_id"]
        .tolist()
    )


def resolve_index_uri(bam_uri: str, explicit_index: str | None) -> str | None:
    """Return the index URI for a BAM/CRAM.

    Uses *explicit_index* if non-empty; otherwise infers from *bam_uri* extension:
      *.bam  → *.bam.bai
      *.cram → *.cram.crai
    Returns None if the extension is unrecognised and no explicit index was given.
    """
    if explicit_index and explicit_index.strip():
        return explicit_index.strip()
    lower = bam_uri.lower()
    if lower.endswith(".bam"):
        return bam_uri + ".bai"
    if lower.endswith(".cram"):
        return bam_uri + ".crai"
    if lower.endswith(".vcf.gz") or lower.endswith(".vcf.bgz"):
        return bam_uri + ".tbi"
    if lower.endswith(".bcf"):
        return bam_uri + ".csi"
    return None


def gs_to_signed_url(gs_uri: str, expiry_minutes: int = 60) -> str:
    """Convert a gs:// URI to a signed HTTPS URL valid for *expiry_minutes*.

    Requires service-account credentials (ADC or GOOGLE_APPLICATION_CREDENTIALS).
    Falls back to a direct storage.googleapis.com URL for publicly readable objects.
    """
    from datetime import timedelta

    try:
        import google.auth
        from google.cloud import storage

        bucket_name, blob_path = gs_uri[5:].split("/", 1)
        credentials, _ = google.auth.default()
        client = storage.Client(credentials=credentials)
        blob = client.bucket(bucket_name).blob(blob_path)
        return blob.generate_signed_url(
            version="v4",
            expiration=timedelta(minutes=expiry_minutes),
            method="GET",
        )
    except Exception:
        # Fall back to direct HTTPS (works for public buckets)
        bucket_name, blob_path = gs_uri[5:].split("/", 1)
        return f"https://storage.googleapis.com/{bucket_name}/{blob_path}"


def load_manifest(path: str) -> dict:
    """Load a manifest TSV and return a dict keyed by sample_id.

    Expected columns: sample_id/bam_path or collaborator_sample_id/duplex_output_bam.
    Optional columns: bai_path/duplex_output_bam_index, variants_path/final_annotated_variants.

    The returned dict maps sample_id → {"bam": str, "bai": str|None, "variants_tsv": str|None}.
    """
    import pandas as pd

    def first_attr(row, names: list[str]) -> str | None:
        for name in names:
            if hasattr(row, name):
                val = getattr(row, name)
                if pd.notna(val) and str(val).strip():
                    return str(val)
        return None

    mdf = pd.read_csv(path.strip(), sep="\t")
    result = {}
    for row in mdf.itertuples(index=False):
        sid = first_attr(row, ["sample_id", "collaborator_sample_id"])
        bam = first_attr(row, ["bam_path", "duplex_output_bam"])
        if not sid or not bam:
            continue
        result[sid] = {
            "bam": bam,
            "bai": first_attr(row, ["bai_path", "duplex_output_bam_index"]),
            "variants_tsv": first_attr(row, ["variants_path", "final_annotated_variants"]),
        }
    return result
