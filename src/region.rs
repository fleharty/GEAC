use std::path::Path;

use anyhow::Result;

use crate::targets::TargetIntervals;

pub enum RegionInput {
    Single(String),
    Intervals(TargetIntervals),
}

/// Parse a `--region` value as either a BED/interval-list file or a plain region string.
///
/// If the value is an existing file path it is loaded as a `TargetIntervals` (BED or Picard
/// interval list); otherwise it is treated as a region string (e.g. `"chr1:1000-2000"`).
pub fn parse_region_input(region: Option<&str>) -> Result<Option<RegionInput>> {
    match region {
        None => Ok(None),
        Some(s) => {
            let path = Path::new(s);
            if path.exists() {
                Ok(Some(RegionInput::Intervals(TargetIntervals::load(path)?)))
            } else {
                Ok(Some(RegionInput::Single(s.to_owned())))
            }
        }
    }
}
