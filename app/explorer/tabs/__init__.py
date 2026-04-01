from . import cohort
from . import duplex_simplex
from . import error_spectrum
from . import overlap_agreement
from . import panel_of_normals
from . import pipeline_comparison
from . import read_type_comparison
from . import reads
from . import strand_bias
from . import tumor_normal
from . import vaf_distribution

TAB_MODULES = (
    vaf_distribution,
    error_spectrum,
    strand_bias,
    overlap_agreement,
    cohort,
    reads,
    duplex_simplex,
    tumor_normal,
    panel_of_normals,
    pipeline_comparison,
    read_type_comparison,
)
