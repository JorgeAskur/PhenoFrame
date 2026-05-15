"""PyMaize: procedural maize plant modeling and phenotyping trait computation."""

from ._version import __version__
from .wrapper import (
    Leaf,
    LeafDesc,
    LeafSplineTraits,
    Maize,
    Tiller,
    TillerDesc,
    build_example,
    from_xml_to_obj,
)
from .traits import compute_traits_from_descriptor, write_traits_xml

__all__ = [
    "__version__",
    "Leaf",
    "LeafDesc",
    "LeafSplineTraits",
    "Maize",
    "Tiller",
    "TillerDesc",
    "build_example",
    "compute_traits_from_descriptor",
    "from_xml_to_obj",
    "write_traits_xml",
]
