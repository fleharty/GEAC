use rust_htslib::bam::{self, Read, record::{Cigar, CigarString, Record}};
use std::path::{Path, PathBuf};

/// Write a minimal FASTA reference for chr1 (all A's, `len` bp) plus its .fai index.
pub fn write_reference(dir: &Path, len: usize) -> PathBuf {
    let fa = dir.join("ref.fa");
    let seq = "A".repeat(len);
    std::fs::write(&fa, format!(">chr1\n{}\n", seq)).unwrap();

    // FAI format: name, length, byte_offset_of_seq, bases_per_line, bytes_per_line
    // ">chr1\n" is 6 bytes, sequence is one line of `len` bases + newline.
    let fai = dir.join("ref.fa.fai");
    std::fs::write(&fai, format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1)).unwrap();

    fa
}

/// A single locus to write into a test BAM: (0-based position, alt base, n_alt_reads, n_ref_reads).
pub type Locus = (i64, u8, usize, usize);

/// Write a coordinate-sorted, indexed BAM with synthetic reads.
///
/// For each locus `(alt_pos, alt_base, n_alt, n_ref)`, writes:
///   - `n_alt` reads of length `read_len` starting 10 bp before `alt_pos`, with `alt_base`
///     at the position that maps to `alt_pos`.
///   - `n_ref` reads at the same start position, all matching the A-repeat reference.
///
/// Loci are written in ascending position order so the BAM is coordinate-sorted.
pub fn write_bam(
    dir: &Path,
    filename: &str,
    sample_id: &str,
    ref_len: usize,
    mut loci: Vec<Locus>,
    read_len: usize,
) -> PathBuf {
    loci.sort_by_key(|l| l.0);

    let bam_path = dir.join(filename);

    let mut header = bam::header::Header::new();

    let mut hd = bam::header::HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);

    let mut sq = bam::header::HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", ref_len as i32);
    header.push_record(&sq);

    let mut rg = bam::header::HeaderRecord::new(b"RG");
    rg.push_tag(b"ID", "rg1");
    rg.push_tag(b"SM", sample_id);
    header.push_record(&rg);

    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    let quals = vec![40u8; read_len];
    let cigar = CigarString(vec![Cigar::Match(read_len as u32)]);

    for (alt_pos, alt_base, n_alt, n_ref) in &loci {
        let read_start = (alt_pos - 10).max(0);
        let alt_offset = (alt_pos - read_start) as usize;

        let mut alt_seq = vec![b'A'; read_len];
        alt_seq[alt_offset] = *alt_base;

        for i in 0..*n_alt {
            let mut rec = Record::new();
            rec.set(format!("alt_{}_{i}", alt_pos).as_bytes(), Some(&cigar), &alt_seq, &quals);
            rec.set_flags(0); // explicitly mark as mapped, forward strand, primary
            rec.set_tid(0);
            rec.set_pos(read_start);
            rec.set_mapq(60);
            rec.push_aux(b"RG", bam::record::Aux::String("rg1")).unwrap();
            writer.write(&rec).unwrap();
        }

        let ref_seq = vec![b'A'; read_len];
        for i in 0..*n_ref {
            let mut rec = Record::new();
            rec.set(format!("ref_{}_{i}", alt_pos).as_bytes(), Some(&cigar), &ref_seq, &quals);
            rec.set_flags(0); // explicitly mark as mapped, forward strand, primary
            rec.set_tid(0);
            rec.set_pos(read_start);
            rec.set_mapq(60);
            rec.push_aux(b"RG", bam::record::Aux::String("rg1")).unwrap();
            writer.write(&rec).unwrap();
        }
    }

    drop(writer);

    // Use samtools to sort and index so the BAM and BAI are standard-compliant.
    let sorted = dir.join(format!("{}.sorted.bam", filename));
    let status = std::process::Command::new("samtools")
        .args(["sort", "-o", sorted.to_str().unwrap(), bam_path.to_str().unwrap()])
        .status()
        .expect("samtools sort failed");
    assert!(status.success(), "samtools sort failed");

    let status = std::process::Command::new("samtools")
        .args(["index", sorted.to_str().unwrap()])
        .status()
        .expect("samtools index failed");
    assert!(status.success(), "samtools index failed");

    sorted
}

/// Run `geac` with the given arguments and return the process output.
pub fn run_geac(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_geac"))
        .args(args)
        .output()
        .expect("failed to execute geac")
}

/// Run `geac` and panic with stdout/stderr if it exits non-zero.
pub fn assert_geac_success(args: &[&str]) {
    let out = run_geac(args);
    if !out.status.success() {
        panic!(
            "geac {:?} failed ({})\nstdout: {}\nstderr: {}",
            args,
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr),
        );
    }
}

/// Count rows in a Parquet file using an in-memory DuckDB connection.
pub fn parquet_count(path: &Path) -> i64 {
    let conn = duckdb::Connection::open_in_memory().unwrap();
    conn.query_row(
        &format!("SELECT COUNT(*) FROM read_parquet('{}')", path.display()),
        [],
        |row| row.get(0),
    )
    .unwrap()
}

/// Fetch a single integer column from a Parquet file with the given WHERE clause.
pub fn parquet_query_i32(path: &Path, column: &str, where_clause: &str) -> i32 {
    let conn = duckdb::Connection::open_in_memory().unwrap();
    conn.query_row(
        &format!(
            "SELECT {column} FROM read_parquet('{}') WHERE {where_clause} LIMIT 1",
            path.display()
        ),
        [],
        |row| row.get(0),
    )
    .unwrap()
}

/// Count records in a BAM file using a plain sequential reader.
pub fn count_bam_reads(bam_path: &Path) -> usize {
    let mut reader = bam::Reader::from_path(bam_path).unwrap();
    reader.records().count()
}

/// Fetch a single string column from a Parquet file with the given WHERE clause.
pub fn parquet_query_str(path: &Path, column: &str, where_clause: &str) -> String {
    let conn = duckdb::Connection::open_in_memory().unwrap();
    conn.query_row(
        &format!(
            "SELECT {column} FROM read_parquet('{}') WHERE {where_clause} LIMIT 1",
            path.display()
        ),
        [],
        |row| row.get(0),
    )
    .unwrap()
}
