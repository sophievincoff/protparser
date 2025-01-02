# src/protparser/__init__.py

# Import key classes/functions from submodules
from .alphafold import AlphaFoldStructure
from .rcsb import RCSBStructure

# Optional: Define __all__ for wildcard imports
__all__ = [
    "AlphaFoldStructure",
    "RCSBStructure"
]
