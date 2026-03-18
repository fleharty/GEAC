use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// Shared state updated by the pileup loop and read by the reporter thread.
#[derive(Clone)]
pub struct ProgressState {
    pub positions_processed: Arc<AtomicU64>,
    pub reads_processed: Arc<AtomicU64>,
    pub alt_bases_found: Arc<AtomicU64>,
    pub current_locus: Arc<Mutex<String>>,
    stop: Arc<AtomicBool>,
}

impl ProgressState {
    pub fn new() -> Self {
        Self {
            positions_processed: Arc::new(AtomicU64::new(0)),
            reads_processed: Arc::new(AtomicU64::new(0)),
            alt_bases_found: Arc::new(AtomicU64::new(0)),
            current_locus: Arc::new(Mutex::new(String::from("—"))),
            stop: Arc::new(AtomicBool::new(false)),
        }
    }

    pub fn update_locus(&self, chrom: &str, pos: i64) {
        if let Ok(mut locus) = self.current_locus.lock() {
            locus.clear();
            locus.push_str(chrom);
            locus.push(':');
            // Display as 1-based for readability
            locus.push_str(&(pos + 1).to_string());
        }
    }
}

pub struct ProgressReporter {
    state: ProgressState,
    handle: Option<JoinHandle<()>>,
}

impl ProgressReporter {
    /// Start a background reporter thread that prints stats every `interval` seconds.
    /// Returns immediately if `interval` is 0.
    pub fn start(interval_secs: u64) -> (Self, ProgressState) {
        let state = ProgressState::new();
        let reporter_state = state.clone();

        let handle = if interval_secs == 0 {
            None
        } else {
            let s = reporter_state.clone();
            let interval = Duration::from_secs(interval_secs);
            Some(thread::spawn(move || {
                let start = Instant::now();
                let mut last_reads: u64 = 0;
                let mut last_tick = Instant::now();

                loop {
                    thread::sleep(Duration::from_millis(500));

                    if s.stop.load(Ordering::Relaxed) {
                        break;
                    }

                    if last_tick.elapsed() < interval {
                        continue;
                    }

                    let now = Instant::now();
                    let elapsed = start.elapsed().as_secs_f64();
                    let positions = s.positions_processed.load(Ordering::Relaxed);
                    let reads = s.reads_processed.load(Ordering::Relaxed);
                    let alts = s.alt_bases_found.load(Ordering::Relaxed);
                    let locus = s
                        .current_locus
                        .lock()
                        .map(|l| l.clone())
                        .unwrap_or_else(|_| "?".to_string());

                    let tick_secs = now.duration_since(last_tick).as_secs_f64();
                    let reads_per_sec = if tick_secs > 0.0 {
                        (reads - last_reads) as f64 / tick_secs
                    } else {
                        0.0
                    };

                    eprintln!(
                        "[{elapsed:>6.0}s] locus={locus:<25} positions={positions:<10} reads={reads:<12} reads/sec={reads_per_sec:<10.0} alt_bases={alts}"
                    );

                    last_reads = reads;
                    last_tick = now;
                }
            }))
        };

        (Self { state: state.clone(), handle }, state)
    }

    /// Print a final summary and stop the background thread.
    pub fn finish(mut self, start: Instant) {
        self.state.stop.store(true, Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }

        let elapsed = start.elapsed().as_secs_f64();
        let positions = self.state.positions_processed.load(Ordering::Relaxed);
        let reads = self.state.reads_processed.load(Ordering::Relaxed);
        let alts = self.state.alt_bases_found.load(Ordering::Relaxed);
        let reads_per_sec = if elapsed > 0.0 { reads as f64 / elapsed } else { 0.0 };

        eprintln!(
            "[done] elapsed={elapsed:.1}s positions={positions} reads={reads} reads/sec={reads_per_sec:.0} alt_bases={alts}"
        );
    }
}
