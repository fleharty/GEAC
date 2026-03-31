use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString, Record},
    Read,
};
use std::path::{Path, PathBuf};

// ── Reference helpers ─────────────────────────────────────────────────────────

/// Write a FASTA reference for chr1 filled with the given ASCII base, plus its .fai.
/// The filename is derived from the base character (e.g. 'C' → "ref_C.fa").
pub fn write_reference_with_base(dir: &Path, len: usize, base: u8) -> PathBuf {
    let name = format!("ref_{}.fa", base as char);
    let fa = dir.join(&name);
    let seq = std::iter::repeat(base as char)
        .take(len)
        .collect::<String>();
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fai = dir.join(format!("{name}.fai"));
    std::fs::write(&fai, format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1)).unwrap();
    fa
}

/// Write a BED file on chr1 with the given 0-based half-open intervals.
pub fn write_bed(dir: &Path, filename: &str, intervals: &[(u32, u32)]) -> PathBuf {
    let path = dir.join(filename);
    let content: String = intervals
        .iter()
        .map(|(s, e)| format!("chr1\t{s}\t{e}\n"))
        .collect();
    std::fs::write(&path, content).unwrap();
    path
}

// ── Coverage BAM helper ───────────────────────────────────────────────────────

/// Per-read specification for `write_coverage_bam`.
pub struct CovRead {
    /// 0-based start position of the read on chr1
    pub pos: i64,
    /// SAM flags (e.g. 0 = normal, 0x400 = duplicate, 0x1|0x2|0x40 = proper-pair R1)
    pub flags: u16,
    /// Mapping quality
    pub mapq: u8,
    /// SAM TLEN (0 = not set; positive on R1, negative on R2)
    pub insert_size: i32,
    /// Mate position (only used when insert_size != 0)
    pub mate_pos: i64,
}

impl CovRead {
    pub fn regular(pos: i64) -> Self {
        Self {
            pos,
            flags: 0,
            mapq: 60,
            insert_size: 0,
            mate_pos: 0,
        }
    }
    pub fn duplicate(pos: i64) -> Self {
        Self {
            pos,
            flags: 0x400,
            mapq: 60,
            insert_size: 0,
            mate_pos: 0,
        }
    }
    pub fn mapq0(pos: i64) -> Self {
        Self {
            pos,
            flags: 0,
            mapq: 0,
            insert_size: 0,
            mate_pos: 0,
        }
    }
    /// Proper-pair R1 read with the given insert size (TLEN = +insert_size).
    pub fn r1_paired(pos: i64, insert_size: i32) -> Self {
        let mate_pos = pos + insert_size as i64 - 1; // approximate mate start
        Self {
            pos,
            flags: 0x1 | 0x2 | 0x40,
            mapq: 60,
            insert_size,
            mate_pos,
        }
    }
}

/// Write a coordinate-sorted, indexed BAM from explicit per-read specifications.
/// All reads have all-A sequence of `read_len` bases; CIGAR = `{read_len}M`.
pub fn write_coverage_bam(
    dir: &Path,
    filename: &str,
    sample_id: &str,
    ref_len: usize,
    reads: Vec<CovRead>,
    read_len: usize,
) -> PathBuf {
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

    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();

    let seq = vec![b'A'; read_len];
    let quals = vec![40u8; read_len];
    let cigar = CigarString(vec![Cigar::Match(read_len as u32)]);

    for (i, r) in reads.iter().enumerate() {
        let mut rec = Record::new();
        rec.set(format!("read_{i}").as_bytes(), Some(&cigar), &seq, &quals);
        rec.set_flags(r.flags);
        rec.set_tid(0);
        rec.set_pos(r.pos);
        rec.set_mapq(r.mapq);
        if r.insert_size != 0 {
            rec.set_mtid(0);
            rec.set_mpos(r.mate_pos);
            rec.set_insert_size(r.insert_size as i64);
        }
        rec.push_aux(b"RG", bam::record::Aux::String("rg1"))
            .unwrap();
        writer.write(&rec).unwrap();
    }

    drop(writer);

    let sorted = dir.join(format!("{filename}.sorted.bam"));
    let status = std::process::Command::new("samtools")
        .args([
            "sort",
            "-o",
            sorted.to_str().unwrap(),
            bam_path.to_str().unwrap(),
        ])
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

// ── Query helpers ─────────────────────────────────────────────────────────────

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

    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    let quals = vec![40u8; read_len];
    let cigar = CigarString(vec![Cigar::Match(read_len as u32)]);

    for (alt_pos, alt_base, n_alt, n_ref) in &loci {
        let read_start = (alt_pos - 10).max(0);
        let alt_offset = (alt_pos - read_start) as usize;

        let mut alt_seq = vec![b'A'; read_len];
        alt_seq[alt_offset] = *alt_base;

        for i in 0..*n_alt {
            let mut rec = Record::new();
            rec.set(
                format!("alt_{}_{i}", alt_pos).as_bytes(),
                Some(&cigar),
                &alt_seq,
                &quals,
            );
            rec.set_flags(0); // explicitly mark as mapped, forward strand, primary
            rec.set_tid(0);
            rec.set_pos(read_start);
            rec.set_mapq(60);
            rec.push_aux(b"RG", bam::record::Aux::String("rg1"))
                .unwrap();
            writer.write(&rec).unwrap();
        }

        let ref_seq = vec![b'A'; read_len];
        for i in 0..*n_ref {
            let mut rec = Record::new();
            rec.set(
                format!("ref_{}_{i}", alt_pos).as_bytes(),
                Some(&cigar),
                &ref_seq,
                &quals,
            );
            rec.set_flags(0); // explicitly mark as mapped, forward strand, primary
            rec.set_tid(0);
            rec.set_pos(read_start);
            rec.set_mapq(60);
            rec.push_aux(b"RG", bam::record::Aux::String("rg1"))
                .unwrap();
            writer.write(&rec).unwrap();
        }
    }

    drop(writer);

    // Use samtools to sort and index so the BAM and BAI are standard-compliant.
    let sorted = dir.join(format!("{}.sorted.bam", filename));
    let status = std::process::Command::new("samtools")
        .args([
            "sort",
            "-o",
            sorted.to_str().unwrap(),
            bam_path.to_str().unwrap(),
        ])
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

/// Per-locus specification for `write_paired_bam`.
///
/// Fields: `(alt_pos, alt_base, n_both_alt, n_both_ref, n_r2_only_alt)`
/// - `n_both_alt`: pairs where both R1 (forward) and R2 (reverse) carry the alt base
/// - `n_both_ref`: pairs where both reads carry the reference base
/// - `n_r2_only_alt`: pairs where R1=ref and R2=alt (R2-only artefact pattern)
pub type PairedLocus = (i64, u8, usize, usize, usize);

/// Write overlapping paired-end reads for testing fwd/rev strand count logic.
///
/// Each pair is fully overlapping: R1 and R2 both start at `alt_pos - 10` covering
/// the variant site. R1 is forward (flag 0x63), R2 is reverse (flag 0x93). Reads
/// in the same pair share a QNAME so `tally_pileup` treats them as overlapping.
pub fn write_paired_bam(
    dir: &Path,
    filename: &str,
    sample_id: &str,
    ref_len: usize,
    loci: Vec<PairedLocus>,
    read_len: usize,
) -> PathBuf {
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

    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    let quals = vec![40u8; read_len];
    let cigar = CigarString(vec![Cigar::Match(read_len as u32)]);

    // R1: paired|proper|mate_rev|first = 0x63
    // R2: paired|proper|rev|second    = 0x93
    let r1_flags: u16 = 0x1 | 0x2 | 0x20 | 0x40;
    let r2_flags: u16 = 0x1 | 0x2 | 0x10 | 0x80;

    for (alt_pos, alt_base, n_both_alt, n_both_ref, n_r2_only_alt) in &loci {
        let read_start = (alt_pos - 10).max(0);
        let alt_offset = (alt_pos - read_start) as usize;

        let mut alt_seq = vec![b'A'; read_len];
        alt_seq[alt_offset] = *alt_base;
        let ref_seq = vec![b'A'; read_len];

        macro_rules! write_pair {
            ($label:expr, $r1_seq:expr, $r2_seq:expr) => {{
                let mut r1 = Record::new();
                r1.set($label, Some(&cigar), $r1_seq, &quals);
                r1.set_flags(r1_flags);
                r1.set_tid(0);
                r1.set_pos(read_start);
                r1.set_mapq(60);
                r1.set_mtid(0);
                r1.set_mpos(read_start);
                r1.set_insert_size(read_len as i64);
                r1.push_aux(b"RG", bam::record::Aux::String("rg1")).unwrap();
                writer.write(&r1).unwrap();

                let mut r2 = Record::new();
                r2.set($label, Some(&cigar), $r2_seq, &quals);
                r2.set_flags(r2_flags);
                r2.set_tid(0);
                r2.set_pos(read_start);
                r2.set_mapq(60);
                r2.set_mtid(0);
                r2.set_mpos(read_start);
                r2.set_insert_size(-(read_len as i64));
                r2.push_aux(b"RG", bam::record::Aux::String("rg1")).unwrap();
                writer.write(&r2).unwrap();
            }};
        }

        for i in 0..*n_both_alt {
            let qname = format!("both_alt_{}_{i}", alt_pos);
            write_pair!(qname.as_bytes(), &alt_seq, &alt_seq);
        }
        for i in 0..*n_both_ref {
            let qname = format!("both_ref_{}_{i}", alt_pos);
            write_pair!(qname.as_bytes(), &ref_seq, &ref_seq);
        }
        for i in 0..*n_r2_only_alt {
            let qname = format!("r2_only_{}_{i}", alt_pos);
            write_pair!(qname.as_bytes(), &ref_seq, &alt_seq);
        }
    }

    drop(writer);

    let sorted = dir.join(format!("{}.sorted.bam", filename));
    let status = std::process::Command::new("samtools")
        .args([
            "sort",
            "-o",
            sorted.to_str().unwrap(),
            bam_path.to_str().unwrap(),
        ])
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

// ── Cycle-number test BAM ─────────────────────────────────────────────────────

/// Specification for a single read in `write_cycle_bam`.
pub struct CycleTestRead {
    /// 0-based leftmost reference position of the first *aligned* base.
    pub pos: i64,
    /// SAM flags (e.g. 0 = forward, 0x10 = reverse-strand).
    pub flags: u16,
    /// Hard clips at the *start* of the stored CIGAR (left side in reference coords).
    pub leading_hard_clips: u32,
    /// Hard clips at the *end* of the stored CIGAR (right side in reference coords).
    pub trailing_hard_clips: u32,
    /// Stored sequence (hard-clipped bases must be omitted; length = Match bases).
    pub seq: Vec<u8>,
    /// Base qualities, same length as `seq`.
    pub quals: Vec<u8>,
}

/// Write a coordinate-sorted, indexed BAM from explicit `CycleTestRead` specifications.
///
/// The CIGAR for each read is built as `[leading_hard_clips H,] {seq_len} M [, trailing_hard_clips H]`.
pub fn write_cycle_bam(
    dir: &Path,
    filename: &str,
    sample_id: &str,
    ref_len: usize,
    reads: Vec<CycleTestRead>,
) -> PathBuf {
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

    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();

    for (i, r) in reads.iter().enumerate() {
        let match_len = r.seq.len() as u32;
        let mut ops = Vec::new();
        if r.leading_hard_clips > 0 {
            ops.push(Cigar::HardClip(r.leading_hard_clips));
        }
        ops.push(Cigar::Match(match_len));
        if r.trailing_hard_clips > 0 {
            ops.push(Cigar::HardClip(r.trailing_hard_clips));
        }
        let cigar = CigarString(ops);

        let mut rec = Record::new();
        rec.set(
            format!("cycle_read_{i}").as_bytes(),
            Some(&cigar),
            &r.seq,
            &r.quals,
        );
        rec.set_flags(r.flags);
        rec.set_tid(0);
        rec.set_pos(r.pos);
        rec.set_mapq(60);
        rec.push_aux(b"RG", bam::record::Aux::String("rg1"))
            .unwrap();
        writer.write(&rec).unwrap();
    }

    drop(writer);

    let sorted = dir.join(format!("{filename}.sorted.bam"));
    let status = std::process::Command::new("samtools")
        .args([
            "sort",
            "-o",
            sorted.to_str().unwrap(),
            bam_path.to_str().unwrap(),
        ])
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

/// Check whether a table exists in a DuckDB file.
pub fn duckdb_table_exists(db: &Path, table: &str) -> bool {
    let conn = duckdb::Connection::open(db).unwrap();
    conn.execute_batch(&format!("SELECT 1 FROM {} LIMIT 0", table))
        .is_ok()
}

/// Count rows in a DuckDB table.
pub fn duckdb_count(db: &Path, table: &str) -> i64 {
    let conn = duckdb::Connection::open(db).unwrap();
    conn.query_row(&format!("SELECT COUNT(*) FROM {}", table), [], |r| r.get(0))
        .unwrap()
}

/// Fetch a single i64 column from a Parquet file with the given WHERE clause.
pub fn parquet_query_i64(path: &Path, column: &str, where_clause: &str) -> i64 {
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

/// Fetch a single f32 column from a Parquet file with the given WHERE clause.
pub fn parquet_query_f32(path: &Path, column: &str, where_clause: &str) -> f32 {
    let conn = duckdb::Connection::open_in_memory().unwrap();
    conn.query_row(
        &format!(
            "SELECT {column} FROM read_parquet('{}') WHERE {where_clause} LIMIT 1",
            path.display()
        ),
        [],
        |row| row.get::<_, f64>(0).map(|v| v as f32),
    )
    .unwrap()
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

/// Return the column names in a Parquet file.
pub fn parquet_columns(path: &Path) -> Vec<String> {
    let conn = duckdb::Connection::open_in_memory().unwrap();
    let mut stmt = conn
        .prepare(&format!(
            "DESCRIBE SELECT * FROM read_parquet('{}')",
            path.display()
        ))
        .unwrap();
    stmt.query_map([], |row| row.get(0))
        .unwrap()
        .collect::<Result<Vec<String>, _>>()
        .unwrap()
}

/// Return the column names in a DuckDB table.
pub fn duckdb_columns(db: &Path, table: &str) -> Vec<String> {
    let conn = duckdb::Connection::open(db).unwrap();
    let mut stmt = conn
        .prepare(&format!("DESCRIBE SELECT * FROM {}", table))
        .unwrap();
    stmt.query_map([], |row| row.get(0))
        .unwrap()
        .collect::<Result<Vec<String>, _>>()
        .unwrap()
}
