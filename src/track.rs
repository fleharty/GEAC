/// Pre-computed annotation tracks loaded from BEDGraph files.
///
/// Each track covers a set of genomic intervals with a numeric score (0.0–1.0 for
/// mappability tracks, arbitrary for others). The typical use-case is mappability
/// tracks such as ENCODE GEM, genmap, or Umap.
///
/// # Usage
///
/// ```text
/// geac coverage --track gem150:gem_150mer.bedgraph --track umap50:umap_k50.bedgraph
/// ```
///
/// This produces two extra nullable Float32 columns in the output Parquet named
/// `gem150` and `umap50`.
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};

// ── Single track ──────────────────────────────────────────────────────────────

/// One BEDGraph track for a single chromosome: sorted intervals with scores.
/// Invariant: intervals are sorted by `start` and non-overlapping (BEDGraph spec).
struct ChromTrack {
    /// Parallel arrays: start positions (0-based), end positions (exclusive), scores.
    starts: Vec<i64>,
    ends:   Vec<i64>,
    scores: Vec<f32>,
}

impl ChromTrack {
    fn new() -> Self {
        Self { starts: Vec::new(), ends: Vec::new(), scores: Vec::new() }
    }

    fn push(&mut self, start: i64, end: i64, score: f32) {
        self.starts.push(start);
        self.ends.push(end);
        self.scores.push(score);
    }

    /// Sort all intervals by start position.
    fn sort(&mut self) {
        if self.starts.is_empty() {
            return;
        }
        // Build an index-sorted permutation.
        let n = self.starts.len();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by_key(|&i| self.starts[i]);
        let starts  = order.iter().map(|&i| self.starts[i]).collect();
        let ends    = order.iter().map(|&i| self.ends[i]).collect();
        let scores  = order.iter().map(|&i| self.scores[i]).collect();
        self.starts = starts;
        self.ends   = ends;
        self.scores = scores;
    }

    /// Binary-search for the score at `pos` (0-based). Returns `None` if `pos` is
    /// not covered by any interval.
    fn get(&self, pos: i64) -> Option<f32> {
        // Find the last interval whose start <= pos.
        let idx = self.starts.partition_point(|&s| s <= pos);
        if idx == 0 {
            return None;
        }
        let i = idx - 1;
        if pos < self.ends[i] {
            Some(self.scores[i])
        } else {
            None
        }
    }
}

// ── Single named track (all chromosomes) ─────────────────────────────────────

/// A single named annotation track loaded from a BEDGraph file.
pub struct AnnotationTrack {
    chroms: HashMap<String, ChromTrack>,
}

impl AnnotationTrack {
    /// Load from a BEDGraph file (gzip-compressed files are NOT supported yet;
    /// use `bgzip -d` to decompress first).
    pub fn load(path: &Path) -> Result<Self> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open track file: {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut chroms: HashMap<String, ChromTrack> = HashMap::new();

        for (lineno, line) in reader.lines().enumerate() {
            let line = line.with_context(|| format!("read error at line {}", lineno + 1))?;
            let line = line.trim();

            // Skip blank lines and BEDGraph header/browser/track directives.
            if line.is_empty() || line.starts_with('#') || line.starts_with("browser")
                || line.starts_with("track")
            {
                continue;
            }

            let mut fields = line.splitn(4, '\t');
            let chrom = fields.next()
                .with_context(|| format!("missing chrom at line {}", lineno + 1))?;
            let start_str = fields.next()
                .with_context(|| format!("missing start at line {}", lineno + 1))?;
            let end_str = fields.next()
                .with_context(|| format!("missing end at line {}", lineno + 1))?;
            let score_str = fields.next()
                .with_context(|| format!("missing score at line {}", lineno + 1))?;

            let start: i64 = start_str.parse()
                .with_context(|| format!("invalid start '{}' at line {}", start_str, lineno + 1))?;
            let end: i64 = end_str.parse()
                .with_context(|| format!("invalid end '{}' at line {}", end_str, lineno + 1))?;
            let score: f32 = score_str.trim().parse()
                .with_context(|| format!("invalid score '{}' at line {}", score_str, lineno + 1))?;

            if start >= end {
                bail!("degenerate interval [{start}, {end}) at line {}", lineno + 1);
            }

            chroms.entry(chrom.to_string()).or_insert_with(ChromTrack::new).push(start, end, score);
        }

        // Sort each chromosome track so binary search is valid.
        for ct in chroms.values_mut() {
            ct.sort();
        }

        Ok(Self { chroms })
    }

    /// Look up the score at `pos` (0-based). Tries `chrom` as-is first, then with
    /// or without the `chr` prefix to bridge NCBI↔UCSC naming.
    pub fn get(&self, chrom: &str, pos: i64) -> Option<f32> {
        if let Some(ct) = self.chroms.get(chrom) {
            return ct.get(pos);
        }
        // Try bridging chr prefix.
        let alt = if chrom.starts_with("chr") {
            chrom[3..].to_string()
        } else {
            format!("chr{chrom}")
        };
        self.chroms.get(&alt)?.get(pos)
    }
}

// ── TrackSet ─────────────────────────────────────────────────────────────────

/// An ordered collection of named tracks.  The order determines the column order
/// in the output Parquet and must match `CoverageRecord::track_values`.
pub struct TrackSet {
    /// Track names in order.
    names: Vec<String>,
    tracks: Vec<AnnotationTrack>,
}

impl TrackSet {
    /// Parse `NAME:path` specifiers (as supplied via `--track`) and load each file.
    pub fn load(specs: &[String]) -> Result<Self> {
        let mut names = Vec::with_capacity(specs.len());
        let mut tracks = Vec::with_capacity(specs.len());

        for spec in specs {
            let (name, path_str) = spec.split_once(':')
                .with_context(|| format!(
                    "invalid --track spec '{spec}': expected NAME:FILE (e.g. gem150:gem.bedgraph)"
                ))?;
            if name.is_empty() {
                bail!("--track spec '{spec}': track name must not be empty");
            }
            let path = Path::new(path_str);
            let track = AnnotationTrack::load(path)
                .with_context(|| format!("loading track '{name}' from '{path_str}'"))?;
            names.push(name.to_string());
            tracks.push(track);
        }

        Ok(Self { names, tracks })
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Look up scores for all tracks at `(chrom, pos)`. Returns a `Vec` of length
    /// `self.len()` in the same order as `self.names()`.
    pub fn lookup(&self, chrom: &str, pos: i64) -> Vec<Option<f32>> {
        self.tracks.iter().map(|t| t.get(chrom, pos)).collect()
    }

    /// Return a `Vec<None>` of the correct length (for zero-depth positions).
    pub fn none_values(&self) -> Vec<Option<f32>> {
        vec![None; self.len()]
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_bedgraph(dir: &std::path::Path, content: &str) -> std::path::PathBuf {
        let path = dir.join("test.bedgraph");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn basic_lookup() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = write_bedgraph(dir.path(), "chr1\t0\t100\t0.9\nchr1\t100\t200\t0.5\n");
        let track = AnnotationTrack::load(&p).unwrap();
        assert_eq!(track.get("chr1", 0),   Some(0.9));
        assert_eq!(track.get("chr1", 99),  Some(0.9));
        assert_eq!(track.get("chr1", 100), Some(0.5));
        assert_eq!(track.get("chr1", 199), Some(0.5));
        assert_eq!(track.get("chr1", 200), None);
    }

    #[test]
    fn gap_between_intervals_returns_none() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = write_bedgraph(dir.path(), "chr1\t0\t50\t1.0\nchr1\t100\t150\t0.8\n");
        let track = AnnotationTrack::load(&p).unwrap();
        assert_eq!(track.get("chr1", 75), None);
    }

    #[test]
    fn chr_prefix_bridging() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = write_bedgraph(dir.path(), "1\t0\t100\t0.7\n");
        let track = AnnotationTrack::load(&p).unwrap();
        // Track has "1"; look up with "chr1" — should bridge.
        assert_eq!(track.get("chr1", 50), Some(0.7));
        // Reverse bridge: track has "1", look up with "1" directly.
        assert_eq!(track.get("1", 50), Some(0.7));
    }

    #[test]
    fn track_directives_and_comments_skipped() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = write_bedgraph(dir.path(),
            "track type=bedGraph name=test\nbrowser position chr1:1-100\n# comment\nchr1\t10\t20\t0.5\n"
        );
        let track = AnnotationTrack::load(&p).unwrap();
        assert_eq!(track.get("chr1", 15), Some(0.5));
        assert_eq!(track.get("chr1", 9),  None);
    }

    #[test]
    fn unsorted_input_is_sorted_on_load() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = write_bedgraph(dir.path(), "chr1\t100\t200\t0.5\nchr1\t0\t100\t0.9\n");
        let track = AnnotationTrack::load(&p).unwrap();
        assert_eq!(track.get("chr1", 50),  Some(0.9));
        assert_eq!(track.get("chr1", 150), Some(0.5));
    }

    #[test]
    fn trackset_lookup_ordered() {
        let dir = tempfile::TempDir::new().unwrap();
        let p1 = write_bedgraph(dir.path(), "chr1\t0\t100\t0.9\n");
        // Use a different filename for the second track.
        let p2 = dir.path().join("track2.bedgraph");
        std::fs::write(&p2, "chr1\t0\t100\t0.4\n").unwrap();

        let specs = vec![
            format!("gem150:{}", p1.display()),
            format!("umap50:{}", p2.display()),
        ];
        let ts = TrackSet::load(&specs).unwrap();
        assert_eq!(ts.names(), &["gem150", "umap50"]);
        let vals = ts.lookup("chr1", 50);
        assert_eq!(vals, vec![Some(0.9), Some(0.4)]);
    }

    #[test]
    fn trackset_none_values_correct_length() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = write_bedgraph(dir.path(), "chr1\t0\t10\t1.0\n");
        let specs = vec![
            format!("a:{}", p.display()),
            format!("b:{}", p.display()),
            format!("c:{}", p.display()),
        ];
        let ts = TrackSet::load(&specs).unwrap();
        assert_eq!(ts.none_values(), vec![None, None, None]);
    }
}
