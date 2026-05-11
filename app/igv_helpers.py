"""Pure-Python / DuckDB helpers for IGV session generation.

Kept separate from geac_explorer.py so they can be unit-tested without
importing the Streamlit runtime or any visualisation libraries.
"""

from __future__ import annotations

import duckdb
import pandas as pd


def make_bed(df: pd.DataFrame) -> str:
    """Build a BED string from a DataFrame of alt_bases rows.

    Coordinate convention (GEAC uses 0-based positions, BED is 0-based half-open):

    SNV / insertion:
        BED [pos, pos+1) — single-base interval at the variant position.

    Deletion (alt_allele starts with '-', e.g. "-ACG"):
        pos       = anchor position (the reference base before the deletion;
                    same VCF-convention anchor that GEAC stores for every indel).
        del_len   = number of deleted bases = len(alt_allele) - 1
        BED start = pos      (anchor — one base before the deletion gap)
        BED end   = pos + len(alt_allele)
                  = pos + 1 + del_len   (exclusive, covers anchor + all deleted bases)

    The '-' prefix in alt_allele has length 1, which happens to cancel the +1
    needed for the half-open BED end, so `pos + len(alt_allele)` is the correct
    exclusive end without any separate adjustment.

    Where multiple records share the same (chrom, pos), the largest end coordinate
    is used so that the longest deletion at a locus is fully covered.

    Target-boundary note: the Rust collector marks a deletion as on-target if any
    part of its deleted range overlaps a target interval (not just the anchor). So
    a deletion whose anchor is before the target but whose deleted bases span or
    overlap the target will be included when "On target" is selected. The BED
    interval reflects the full deletion extent — intentional, since truncating at
    the target boundary would misrepresent the variant in IGV.
    """
    def _end(row) -> int:
        alt = str(row["alt_allele"])
        if row["variant_type"] == "deletion" and alt.startswith("-"):
            # pos is the anchor; deleted bases span pos+1 .. pos+(len(alt)-1).
            # BED exclusive end to include the last deleted base:
            #   pos + (len(alt)-1) + 1 = pos + len(alt)
            return int(row["pos"]) + len(alt)
        return int(row["pos"]) + 1

    tmp = df.copy()
    tmp["bed_end"] = tmp.apply(_end, axis=1)
    positions = (
        tmp.groupby(["chrom", "pos"])["bed_end"]
        .max()
        .reset_index()
        .sort_values(["chrom", "pos"])
    )
    lines = [
        f"{row.chrom}\t{int(row.pos)}\t{int(row.bed_end)}"
        for row in positions.itertuples(index=False)
    ]
    return "\n".join(lines) + "\n"


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

    Expected columns: collaborator_sample_id, duplex_output_bam,
    duplex_output_bam_index (optional), final_annotated_variants (optional).

    The returned dict maps sample_id → {"bam": str, "bai": str|None, "variants_tsv": str|None}.
    """
    import pandas as pd

    mdf = pd.read_csv(path.strip(), sep="\t")
    result = {}
    for row in mdf.itertuples(index=False):
        bai = (
            str(row.duplex_output_bam_index)
            if hasattr(row, "duplex_output_bam_index") and pd.notna(row.duplex_output_bam_index)
            else None
        )
        variants = (
            str(row.final_annotated_variants)
            if hasattr(row, "final_annotated_variants") and pd.notna(row.final_annotated_variants)
            else None
        )
        result[str(row.collaborator_sample_id)] = {
            "bam": str(row.duplex_output_bam),
            "bai": bai,
            "variants_tsv": variants,
        }
    return result
