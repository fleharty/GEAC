from . import ai_plot_builder
from . import cohort
from . import duplex_simplex
from . import error_spectrum
from . import fragmentomics
from . import mnv_candidates
from . import overlap_agreement
from . import panel_of_normals
from . import pipeline_comparison
from . import read_type_comparison
from . import reads
from . import strand_bias
from . import summary
from . import tumor_normal
from . import vaf_distribution

# Named aliases — use these instead of TAB_MODULES[N] so that reordering the
# tuple below does not silently break tab dispatch in geac_explorer.py.
TAB_SUMMARY = summary
TAB_VAF_DISTRIBUTION = vaf_distribution
TAB_ERROR_SPECTRUM = error_spectrum
TAB_STRAND_BIAS = strand_bias
TAB_COHORT = cohort
TAB_READS = reads
TAB_DUPLEX_SIMPLEX = duplex_simplex
TAB_TUMOR_NORMAL = tumor_normal
TAB_PANEL_OF_NORMALS = panel_of_normals
TAB_PIPELINE_COMPARISON = pipeline_comparison
TAB_READ_TYPE_COMPARISON = read_type_comparison
TAB_OVERLAP_AGREEMENT = overlap_agreement
TAB_MNV_CANDIDATES = mnv_candidates
TAB_AI_PLOT_BUILDER = ai_plot_builder
TAB_FRAGMENTOMICS = fragmentomics

TAB_MODULES = (
    TAB_SUMMARY,
    TAB_VAF_DISTRIBUTION,
    TAB_ERROR_SPECTRUM,
    TAB_STRAND_BIAS,
    TAB_COHORT,
    TAB_READS,
    TAB_DUPLEX_SIMPLEX,
    TAB_TUMOR_NORMAL,
    TAB_PANEL_OF_NORMALS,
    TAB_PIPELINE_COMPARISON,
    TAB_READ_TYPE_COMPARISON,
    TAB_OVERLAP_AGREEMENT,
    TAB_MNV_CANDIDATES,
    TAB_AI_PLOT_BUILDER,
    TAB_FRAGMENTOMICS,
)
