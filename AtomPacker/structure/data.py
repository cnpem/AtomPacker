# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `get_lattice_constants`
function, and the `lattice_constants` dictionary.

The `lattice_constants` dictionary contains the lattice constants for
different crystal systems. The `get_lattice_constants` function is used to get
the lattice constants for a given atom type and lattice type.
"""

__all__ = ["get_lattice_constants", "lattice_constants"]

from typing import Any, Dict, Optional, Tuple
from warnings import warn


def _get_bcc(
    atom_type: str, a: Optional[float], b: Optional[float], c: Optional[float]
) -> float | None:
    """
    Get the lattice constant for BCC lattice.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.
    a : float, optional
        The lattice constant `a`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    b : float, optional
        The lattice constant `b`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    c : float, optional
        The lattice constant `c`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.

    Returns
    -------
    float | None
        The lattice constant `a` for BCC lattice. If None, use
        `ase.data.reference_states`.

    Raises
    ------
    ValueError
        If the lattice constant `b` and `c` are provided.

    Warns
    -----
    UserWarning
        If the lattice constant `a` is not found in `lattice_constants`.
    """
    # Check constant `b` and `c`
    if b is not None or c is not None:
        raise ValueError(
            "The lattice constant `b` and `c` are not required for BCC \
lattice."
        )

    # Check constant `a`
    if isinstance(a, type(None)):
        if isinstance(a, float):
            return a
    else:
        try:
            lattice_constants[atom_type]["bcc"]
        except KeyError:
            warn(
                f"Constant `a` for {atom_type} not found for `bcc` in \
lattice_constants. Using `ase.data.reference_states`."
            )
            return None


def _get_fcc(
    atom_type: str, a: Optional[float], b: Optional[float], c: Optional[float]
) -> float | None:
    """
    Get the lattice constant for FCC lattice.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.
    a : float, optional
        The lattice constant `a`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    b : float, optional
        The lattice constant `b`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    c : float, optional
        The lattice constant `c`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.

    Returns
    -------
    float | None
        The lattice constant `a` for FCC lattice. If None, use
        `ase.data.reference_states`.

    Raises
    ------
    ValueError
        If the lattice constant `b` and `c` are provided.

    Warns
    -----
    UserWarning
        If the lattice constant `a` is not found in `lattice_constants`.
    """
    # Check constant `b` and `c`
    if b is not None or c is not None:
        raise ValueError(
            "The lattice constant `b` and `c` are not required for FCC \
lattice."
        )

    # Check constant `a`
    if isinstance(a, type(None)):
        if isinstance(a, float):
            return a
    else:
        try:
            lattice_constants[atom_type]["fcc"]
        except KeyError:
            warn(
                f"Constant `a` for {atom_type} not found for `fcc` in \
lattice_constants. Using `ase.data.reference_states`."
            )
            return None


def _get_hcp(
    atom_type: str, a: Optional[float], b: Optional[float], c: Optional[float]
) -> Tuple[float, float] | None:
    """
    Get the lattice constant for HCP lattice.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.
    a : float, optional
        The lattice constant `a`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    b : float, optional
        The lattice constant `b`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    c : float, optional
        The lattice constant `c`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.

    Returns
    -------
    Tuple[float, float] | None
        The lattice constants `(a, c)` for HCP lattice. If None, use
        `ase.data.reference_states`.

    Raises
    ------
    ValueError
        If the lattice constant `b` is provided.
        If the lattice constant `a` is provided but not `c`.
        If the lattice constant `c` is provided but not `a`.

    Warns
    -----
    UserWarning
        If the lattice constants `a` and `c` are not found in
        `lattice_constants`.
    """
    # Check constant `a` and `c`
    if a is None and c is not None:
        raise ValueError(
            "As you set the lattice constant `c`, the lattice constant `a` is \
required."
        )
    if a is not None and c is None:
        raise ValueError(
            "As you set the lattice constant `a`, the lattice constant `c` is \
required."
        )

    # Check constant `b`
    if b is not None:
        raise ValueError(
            "The lattice constant `b` is not required for HCP \
lattice."
        )

    # Check constant `a` and `c`
    if a is not None and c is not None:
        if isinstance(a, float) and isinstance(c, float):
            return (a, c)
    else:
        try:
            lattice_constants[atom_type]["hcp"]
        except KeyError:
            warn(
                f"Constants `a` and `c` for {atom_type} not found for `hcp` \
in lattice_constants. Using `ase.data.reference_states`."
            )
            return None


def _get_sc(
    atom_type: str, a: Optional[float], b: Optional[float], c: Optional[float]
) -> float | None:
    """
     Get the lattice constant for SC lattice.

    Parameters
     ----------
     atom_type : str
        The atomic symbol of the atom.
     a : float, optional
        The lattice constant `a`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
     b : float, optional
        The lattice constant `b`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
     c : float, optional
        The lattice constant `c`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.

     Returns
     -------
     float | None
         The lattice constant `a` for SC lattice. If None, use
        `ase.data.reference_states`.
    """
    # Check constant `b` and `c`
    if b is not None or c is not None:
        raise ValueError(
            "The lattice constant `b` and `c` are not required for SC \
lattice."
        )

    # Check constant `a`
    if isinstance(a, type(None)):
        if isinstance(a, float):
            return a
    else:
        try:
            lattice_constants[atom_type]["sc"]
        except KeyError:
            warn(
                f"Constant `a` for {atom_type} not found for `sc` in \
lattice_constants. Using `ase.data.reference_states`."
            )
            return None


def get_lattice_constants(
    atom_type: str,
    lattice_type: str,
    a: Optional[float],
    b: Optional[float],
    c: Optional[float],
) -> Tuple[float, float, float] | Tuple[float, float] | float | None:
    """
    Get the lattice constants for a given atom type and lattice type.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.
    lattice_type : str
        The type of lattice in the cluster. The available lattice types are
        'bcc', 'fcc', 'hcp' and 'sc'.
    a : float, optional
        The lattice constant `a`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    b : float, optional
        The lattice constant `b`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.
    c : float, optional
        The lattice constant `c`. If provided, it will be used as is. If not,
        the value will be fetched from `AtomPacker.data.lattice_constants`.
        If it's not found there, the experimental value from `ase.data` will
        be used.

    Returns
    -------
    Tuple[float, float, float] | Tuple[float, float] | float | None
        The lattice constants for the given atom type and lattice type. If
        None, use `ase.data.reference_states`.

    Raises
    ------
    ValueError
        If the lattice type is not valid.
    """
    match lattice_type:
        case "bcc":
            return _get_bcc(atom_type, a, b, c)
        case "fcc":
            return _get_fcc(atom_type, a, b, c)
        case "hcp":
            return _get_hcp(atom_type, a, b, c)
        case "sc":
            return _get_sc(atom_type, a, b, c)
        case _:
            raise ValueError(f"Invalid lattice type: {lattice_type}.")


"""
The lattice constants from Martienssen & Warlimont [1]_.

References
----------
.. [1] Martienssen, W., & Warlimont, H. (Eds.). (2005). Springer Handbook of
    Condensed Matter and Materials Data. Springer Berlin Heidelberg.
    https://doi.org/10.1007/3-540-30437-1.
"""
lattice_constants: Dict[str, Dict[str, Any]] = {
    # H (1)
    "H": {
        "hexagonal": {"a": 3.771, "c": 6.156},  # crystallographic
        # "fcc": {"a": 5.334},  # allotropic and high-pressure modifications
        # "hcp": {"a": 3.771, "c": 6.156}, # allotropic and high-pressure modifications
    },
    # He (2)
    "He": {
        "hexagonal": {"a": 3.577, "c": 5.842},  # crystallographic
        # "fcc": {"a": 4.240},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.110},  # allotropic and high-pressure modifications
    },
    # Li (3)
    "Li": {
        "bcc": {"a": 3.5093},  # crystallographic
        # "bcc": {"a": 3.093},  # allotropic and high-pressure modifications
        # "hcp": {"a": 3.111, "c": 5.093},  # allotropic and high-pressure modifications
    },
    # Be (4)
    "Be": {
        "hexagonal": {"a": 2.2857, "c": 3.5839},  # crystallographic
        # "bcc": {"a": 2.5515},  # allotropic and high-pressure modifications
        # "hcp": {"a": 2.2857, "c": 3.5839},  # allotropic and high-pressure modifications
    },
    # B (5)
    # None
    # C (6)
    "C": {
        "diamond": {"a": 3.5671},  # crystallographic
        "graphite": {"a": 2.4612, "c": 6.7090},  # crystallographic
    },
    # N (7)
    "N": {
        "sc": {"a": 5.659},  # crystallographic
        # "cubic": {"a": 5.659},  # allotropic and high-pressure modifications
        # "hexagonal": {"a": 4.046, "c": 6.629},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 3.957, "c": 5.101},  # allotropic and high-pressure modifications
    },
    # O (8)
    "O": {
        "monoclinic": {
            "a": 5.403,
            "b": 3.429,
            "c": 5.086,
            "gamma": 132.53,
        },  # crystallographic
        # "monoclinic": {"a": 5.403, "b": 3.429, "c": 5.086, "gamma": 132.53},  # allotropic and high-pressure modifications
        # "rhombohedral": {"a": 4.210, "alpha": 46.27},  # allotropic and high-pressure modifications
        # "cubic": {"a": 6.830},  # allotropic and high-pressure modifications
    },
    # F (9)
    "F": {
        "monoclinic": {
            "a": 5.50,
            "b": 3.28,
            "c": 7.28,
            "beta": 102.17,
        },  # crystallographic
        # "monoclinic": {"a": 5.50, "b": 3.28, "c": 7.28, "beta": 102.17},  # allotropic and high-pressure modifications
        "cubic": {"a": 6.67},  # allotropic and high-pressure modifications
    },
    # Ne (10)
    "Ne": {
        "fcc": {"a": 4.4622},  # crystallographic
    },
    # Na (11)
    "Na": {
        "bcc": {"a": 4.2096},  # crystallographic
        # "bcc": {"a": 4.2096},  # allotropic and high-pressure modifications
        # "hcp": {"a": 3.767, "c": 6.154},  # allotropic and high-pressure modifications
    },
    # Mg (12)
    "Mg": {
        "hexagonal": {"a": 3.2093, "c": 5.2107},  # crystallographic
    },
    # Al (13)
    "Al": {
        "fcc": {"a": 3.6149},  # crystallographic
        # "fcc": {"a": 4.0496},  # allotropic and high-pressure modifications
        # "hcp": {"a": 2.693, "c": 4.398},  # allotropic and high-pressure modifications
    },
    # Si (14)
    "Si": {
        "diamond": {"a": 5.43102},  # crystallographic
        # "diamond": {"a": 5.4306},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 4.686, "c": 2.585},  # allotropic and high-pressure modifications
        # "sc": {"a": 6.3600},  # allotropic and high-pressure modifications
        # "hexagonal": {"a": 3.80, "c": 6.28},  # allotropic and high-pressure modifications
    },
    # P (15)
    "P": {
        "orthorhombic": {"a": 3.3136, "b": 10.478, "c": 4.3763},  # crystallographic
    },
    # S (16)
    "S": {
        "orthorhombic": {"a": 10.464, "b": 12.8660, "c": 24.4860},  # crystallographic
    },
    # Cl (17)
    "Cl": {
        "orthorhombic": {"a": 6.24, "b": 4.48, "c": 8.26},  # crystallographic
    },
    # Ar (18)
    "Ar": {
        "fcc": {"a": 5.312},  # crystallographic
        # "fcc": {"a": 5.312},  # allotropic and high-pressure modifications
        "hcp": {"a": 3.760, "c": 6.141},  # allotropic and high-pressure modifications
    },
    # K (19)
    "K": {
        "bcc": {"a": 5.321},  # crystallographic
    },
    # Ca (20)
    "Ca": {
        "fcc": {"a": 5.5884},  # crystallographic
        # "bcc": {"a": 4.480},  # allotropic and high-pressure modifications
        # "fcc": {"a": 5.5884},  # allotropic and high-pressure modifications
    },
    # Sc (21)
    "Sc": {
        "hexagonal": {"a": 3.3088, "c": 5.2680},  # crystallographic
        # "hcp": {"a": 3.3088, "c": 5.2680},  # allotropic and high-pressure modifications
    },
    # Ti (22)
    "Ti": {
        "hexagonal": {"a": 2.9503, "c": 4.6836},  # crystallographic
        # "hcp": {"a": 2.9503, "c": 4.6836},  # allotropic and high-pressure modifications
        # "bcc": {"a": 3.3065},  # allotropic and high-pressure modifications
    },
    # V (23)
    "V": {
        "bcc": {"a": 3.0238},  # crystallographic
    },
    # Cr (24)
    "Cr": {
        "bcc": {"a": 2.8847},  # crystallographic
        # "bcc": {"a": 2.8847},  # allotropic and high-pressure modifications
        # "bcc": {"a": 2.8820},  # allotropic and high-pressure modifications
    },
    # Mn (25)
    "Mn": {
        "bcc": {"a": 8.9219},  # crystallographic
        # "bcc": {"a": 8.9219},  # allotropic and high-pressure modifications
        # "cubic": {"a": 6.3152},  # allotropic and high-pressure modifications
        # "fcc": {"a": 3.8624},  # allotropic and high-pressure modifications
        # "bcc": {"a": 3.0806},  # allotropic and high-pressure modifications
    },
    # Fe (26)
    "Fe": {
        "bcc": {"a": 2.8665},  # crystallographic
        # "bcc": {"a": 2.8665},  # allotropic and high-pressure modifications
        # "fcc": {"a": 3.6467},  # allotropic and high-pressure modifications
        # "bcc": {"a": 2.9135},  # allotropic and high-pressure modifications
        # "hcp": {"a": 2.485, "c": 3.990},  # allotropic and high-pressure modifications
    },
    # Co (27)
    "Co": {
        "hcp": {"a": 2.5071, "c": 4.0694},  # crystallographic
        # "fcc": {"a": 3.5445}, # allotropic and high-pressure modifications
        # "hcp": {"a": 2.5071, "c": 4.0694},  # allotropic and high-pressure modifications
    },
    # Ni (28)
    "Ni": {
        "fcc": {"a": 3.5241},  # crystallographic
    },
    # Cu (29)
    "Cu": {
        "fcc": {"a": 3.6149},  # crystallographic
    },
    # Zn (30)
    "Zn": {
        "hexagonal": {"a": 2.6644, "c": 4.9494},  # crystallographic
    },
    # Ga (31)
    "Ga": {
        "orthorhombic": {"a": 4.5192, "b": 7.6586, "c": 4.5258},  # crystallographic
        # "orthorhombic": {"a": 4.5192, "b": 7.6586, "c": 4.5258},  # allotropic and high-pressure modifications
        # "orthorhombic": {"a": 10.593, "b": 13.523, "c": 5.203},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 2.808, "c": 4.458},  # allotropic and high-pressure modifications
    },
    # Ge (32)
    "Ge": {
        "diamond": {"a": 5.659},  # crystallographic
        # "fcc": {"a": 5.6574},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 4.884, "c": 2.692},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 5.93, "c": 6.98},  # allotropic and high-pressure modifications
        # "bcc": {"a": 6.92},  # allotropic and high-pressure modifications
    },
    # As (33)
    "As": {
        "rhombohedral": {"a": 4.1320, "alpha": 54.12},  # crystallographic
        # "rhombohedral": {"a": 4.1320, "alpha": 54.12},  # allotropic and high-pressure modifications
        "orthorhombic": {
            "a": 3.62,
            "b": 10.85,
            "c": 4.48,
        },  # allotropic and high-pressure modifications
    },
    # Se (34)
    "Se": {
        "hexagonal": {"a": 4.3655, "c": 4.9576},  # crystallographic
        # "hexagonal": {"a": 4.3655, "c": 4.9576},  # allotropic and high-pressure modifications
        # "monoclinic": {"a": 9.054, "b": 9.083, "c": 2.336, "gamma": 90.82},  # allotropic and high-pressure modifications
        # "monoclinic": {"a": 15.018 "b": 14.713, "c": 8.879, "gamma": 93.6},  # allotropic and high-pressure modifications
    },
    # Br (35)
    "Br": {
        "orthorhombic": {"a": 6.68, "b": 4.49, "c": 8.74},  # crystallographic
    },
    # Kr (36)
    "Kr": {
        "fcc": {"a": 5.6459},  # crystallographic
    },
    # Rb (37)
    "Rb": {
        "bcc": {"a": 5.703},  # crystallographic
    },
    # Sr (38)
    "Sr": {
        "fcc": {"a": 6.084},  # crystallographic
        # "fcc": {"a": 6.084},  # allotropic and high-pressure modifications
        # "hcp": {"a": 4.280, "c": 7.050},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.870},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.437},  # allotropic and high-pressure modifications
    },
    # Y (39)
    "Y": {
        "hexagonal": {"a": 3.6482, "c": 5.7318},  # crystallographic
        # "hcp": {"a": 3.6482, "c": 5.7318},  # allotropic and high-pressure modifications
    },
    # Zr (40)
    "Zr": {
        "hexagonal": {"a": 3.2317, "c": 5.1476},  # crystallographic
        # "hcp": {"a": 3.2317, "c": 5.1476},  # allotropic and high-pressure modifications
        # "bcc": {"a": 3.609},  # allotropic and high-pressure modifications
    },
    # Nb (41)
    "Nb": {
        "bcc": {"a": 3.3007},  # crystallographic
    },
    # Mo (42)
    "Mo": {
        "bcc": {"a": 3.1470},  # crystallographic
    },
    # Tc (43)
    "Tc": {
        "hexagonal": {"a": 2.738, "c": 4.394},  # crystallographic
    },
    # Ru (44)
    "Ru": {
        "hexagonal": {"a": 2.7053, "c": 4.2814},  # crystallographic
    },
    # Rh (45)
    "Rh": {
        "fcc": {"a": 2.8032},  # crystallographic
    },
    # Pd (46)
    "Pd": {
        "fcc": {"a": 3.8901},
    },
    # Ag (47)
    "Ag": {
        "fcc": {"a": 4.0861},  # crystallographic
    },
    # Cd (48)
    "Cd": {
        "hexagonal": {"a": 2.9788, "c": 5.6167},  # crystallographic
    },
    # In (49)
    "In": {
        "tetragonal": {"a": 4.5990, "c": 4.9470},  # crystallographic
    },
    # Sn (50)
    "Sn": {
        "tetragonal": {"a": 58.197},  # crystallographic
    },
    # Sb (51)
    "Sb": {
        "rhombohedral": {"a": 4.5065, "alpha": 57.11},  # crystallographic
        # "rhombohedral": {"a": 4.5065, "alpha": 57.11},  # allotropic and high-pressure modifications
        # "cubic": {"a": 2.992},  # allotropic and high-pressure modifications
        # "hcp": {"a": 3.376, "c": 5.341},  # allotropic and high-pressure modifications
        # "monoclinic": {"a": 5.56, "b": 4.04, "c": 4.22, "beta": 86.0},  # allotropic and high-pressure modifications
    },
    # Te (52)
    "Te": {
        "hexagonal": {"a": 4.4561, "c": 5.9271},  # crystallographic
        # "hexagonal": {"a": 4.4561, "c": 5.9271},  # allotropic and high-pressure modifications
        # "rhombohedral": {"a": 4.69, "alpha": 53.30},  # allotropic and high-pressure modifications
        # "rhombohedral": {"a": 3.002, "alpha": 103.3},  # allotropic and high-pressure modifications
    },
    # I (53)
    "I": {
        "orthorhombic": {"a": 7.268, "b": 4.797, "c": 9.797},  # crystallographic
    },
    # Xe (54)
    "Xe": {
        "fcc": {"a": 6.132},  # crystallographic
    },
    # Cs (55)
    "Cs": {
        "bcc": {"a": 6.141},  # crystallographic
        # "bcc": {"a": 6.141},  # allotropic and high-pressure modifications
        # "fcc": {"a": 5.984},  # allotropic and high-pressure modifications
        # "fcc": {"a": 5.800},  # allotropic and high-pressure modifications
    },
    # Ba (56)
    "Ba": {
        "bcc": {"a": 5.023},  # crystallographic
    },
    # La (57)
    "La": {
        "hexagonal": {"a": 3.7740, "c": 12.171},  # crystallographic
        # "hexagonal": {"a": 3.7740, "c": 12.171},  # allotropic and high-pressure modifications
        # "fcc": {"a": 5.3045},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.265},  # allotropic and high-pressure modifications
        # "fcc": {"a": 5.170},  # allotropic and high-pressure modifications
    },
    # Ce (58)
    "Ce": {
        "fcc": {"a": 5.1610},  # crystallographic
        # "fcc": {"a": 5.1610},  # allotropic and high-pressure modifications
        # "hexagonal": {"a": 3.673, "c": 11.802},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.82},  # allotropic and high-pressure modifications
    },
    # Pr (59)
    "Pr": {
        "hexagonal": {"a": 3.6721, "c": 11.8326},  # crystallographic
        # "hexagonal": {"a": 3.6721, "c": 11.8326},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.13},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.88},  # allotropic and high-pressure modifications
    },
    # Nd (60)
    "Nd": {
        "hexagonal": {"a": 3.6582, "c": 11.7966},  # crystallographic
        # "hexagonal": {"a": 3.6582, "c": 11.7966},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.13},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.80},  # allotropic and high-pressure modifications
    },
    # Pm (61)
    "Pm": {
        "hexagonal": {"a": 3.65, "c": 11.65},  # crystallographic
        # "hexagonal": {"a": 3.65, "c": 11.65},  # allotropic and high-pressure modifications
    },
    # Sm (62)
    "Sm": {
        "hexagonal": {"a": 3.6290, "c": 26.207},  # crystallographic
        # "trigonal": {"a": 3.629, "c": 26.207},  # crystallographic
        # "hexagonal": {"a": 3.6180, "c": 11.66},  # allotropic and high-pressure modifications
    },
    # Eu (63)
    "Eu": {
        "bcc": {"a": 4.5827},  # crystallographic
    },
    # Gd (64)
    "Gd": {
        "hexagonal": {"a": 3.6336, "c": 5.7810},  # crystallographic
        # "hcp": {"a": 3.6336, "c": 5.7810},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.06},  # allotropic and high-pressure modifications
        # "trigonal": {"a": 3.61, "c": 26.03},  # allotropic and high-pressure modifications
    },
    # Tb (65)
    "Tb": {
        "hexagonal": {"a": 3.6055, "c": 5.6966},  # crystallographic
        # "hcp": {"a": 3.6055, "c": 5.56966},  # allotropic and high-pressure modifications
        # "trigonal": {"a": 3.41, "c": 24.50},  # allotropic and high-pressure modifications
    },
    # Dy (66)
    "Dy": {
        "hexagonal": {"a": 3.5915, "c": 5.6501},  # crystallographic
        # "hcp": {"a": 3.5915, "c": 5.6501},  # allotropic and high-pressure modifications
        # "orthorhombic": {"a": 3.595, "b": 6.184, "c": 5.678},  # allotropic and high-pressure modifications
        # "trigonal": {"a": 3.436, "c": 24.83},  # allotropic and high-pressure modifications
    },
    # Ho (67)
    "Ho": {
        "hexagonal": {"a": 3.5778, "c": 5.6178},  # crystallographic
        # "hcp": {"a": 3.5778, "c": 5.6178},  # allotropic and high-pressure modifications
        # "trigonal": {"a": 3.34, "c": 24.50},  # allotropic and high-pressure modifications
    },
    # Er (68)
    "Er": {
        "hexagonal": {"a": 3.5592, "c": 5.5850},  # crystallographic
        # "hcp": {"a": 3.5592, "c": 5.5850},  # allotropic and high-pressure modifications
    },
    # Tm (69)
    "Tm": {
        "hexagonal": {"a": 3.5375, "c": 5.5540},  # crystallographic
        # "hcp": {"a": 3.5375, "c": 5.5540},  # allotropic and high-pressure modifications
    },
    # Yb (70)
    "Yb": {
        "fcc": {"a": 5.4848},  # crystallographic
        # "fcc": {"a": 5.4848},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.44},  # allotropic and high-pressure modifications
        # "hcp": {"a": 3.8799, "c": 6.3859},  # allotropic and high-pressure modifications
    },
    # Lu (71)
    "Lu": {
        "hexagonal": {"a": 3.5052, "c": 5.5494},  # crystallographic
        # "hcp": {"a": 3.5052, "c": 5.5494},  # allotropic and high-pressure modifications
    },
    # Hf (72)
    "Hf": {
        "hexagonal": {"a": 3.1946, "c": 5.0511},  # crystallographic
        # "hcp": {"a": 3.1946, "c": 5.0511},  # allotropic and high-pressure modifications
        # "bcc": {"a": 3.610},  # allotropic and high-pressure modifications
    },
    # Ta (73)
    "Ta": {
        "bcc": {"a": 3.3031},  # crystallographic
    },
    # W (74)
    "W": {
        "bcc": {"a": 3.1651},  # crystallographic
    },
    # Re (75)
    "Re": {
        "hexagonal": {"a": 2.7608, "c": 4.4580},  # crystallographic
    },
    # Os (76)
    "Os": {
        "hexagonal": {"a": 2.7348, "c": 4.3913},  # crystallographic
    },
    # Ir (77)
    "Ir": {
        "fcc": {"a": 3.8391},  # crystallographic
    },
    # Pt (78)
    "Pt": {
        "fcc": {"a": 3.9233},
    },
    # Au (79)
    "Au": {
        "fcc": {"a": 4.0784},  # crystallographic
    },
    # Hg (80)
    "Hg": {
        "rhombohedral": {
            "a": 3.005,
            "alpha": 70.53,
        },  # crystallographic
        # "tetragonal": {"a": 3.995, "c": 2.825},  # allotropic and high-pressure modifications
    },
    # Tl (81)
    "Tl": {
        "hexagonal": {"a": 3.4563, "c": 5.5263},  # crystallographic
        # "bcc": {"a": 3.4563},  # allotropic and high-pressure 879
        # "hcp": {"a": 3.4563, "c": 5.5263},  # allotropic and high-pressure modifications
    },
    # Pb (82)
    "Pb": {
        "fcc": {"a": 4.9502},  # crystallographic
    },
    # Bi (83)
    "Bi": {
        "rhombohedral": {"a": 4.7460, "alpha": 57.23},  # crystallographic
    },
    # Po (84)
    "Po": {
        "sc": {"a": 3.366},  # crystallographic
        # "cubic": {"a": 3.366},  # allotropic and high-pressure modifications
        # "rhombohedral": {"a": 3.373, "alpha": 98.08},  # allotropic and high-pressure modifications
    },
    # At (85)
    # None
    # Rn (86)
    # None
    # Fr (87)
    # None
    # Ra (88)
    "Ra": {
        "bcc": {"a": 5.148},  # crystallographic
    },
    # Ac (89)
    # None
    # Th (90)
    "Th": {
        "fcc": {"a": 5.0851},  # crystallographic
        # "fcc": {"a": 5.0851},  # allotropic and high-pressure modifications
        # "bcc": {"a": 4.11},  # allotropic and high-pressure modifications
    },
    # Pa (91)
    "Pa": {
        "tetragonal": {"a": 3.945, "c": 3.242},  # crystallographic
        # "tetragonal": {"a": 3.945, "c": 3.242},  # allotropic and high-pressure modifications
    },
    # U (92)
    "U": {
        "orthorhombic": {"a": 2.8538, "b": 5.8680, "c": 4.9557},  # crystallographic
        # "orthorhombic": {"a": 2.8538, "b": 5.8680, "c": 4.9557},  # allotropic and high-pressure modifications
        # "tetragnonal": {"a": 10.759, "c": 5.654},  # allotropic and high-pressure modifications
        # "bcc": {"a": 3.524},  # allotropic and high-pressure modifications
    },
    # Np (93)
    "Np": {
        "orthorhombic": {"a": 6.663, "b": 4.723, "c": 4.887},  # crystallographic
        # "orthorhombic": {"a": 6.663, "b": 4.723, "c": 4.887},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 4.896, "c": 3.387},  # allotropic and high-pressure modifications
        # "bcc": {"a": 3.52},  # allotropic and high-pressure modifications
    },
    # Pu (94)
    "Pu": {
        "monoclinic": {
            "a": 6.183,
            "b": 4.822,
            "c": 10.968,
            "gamma": 101.78,
        },  # crystallographic
        # "monoclinic": {"a": 6.183, "b": 4.822, "c": 10.968, "gamma": 101.78},  # allotropic and high-pressure modifications
        # "monoclinic": {"a": 9.284, "b": 10.463, "c": 7.859, "gamma": 92.13},  # allotropic and high-pressure modifications
        # "orthorhombic": {"a": 3.1587, "b": 5.7682, "c": 10.162},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.6371},  # allotropic and high-pressure modifications
        # "tetragonal": {"a": 3.3261, "c": 4.4630},  # allotropic and high-pressure modifications
        # "bcc": {"a": 5.703},  # allotropic and high-pressure modifications
    },
    # Am (95)
    "Am": {
        "hexagonal": {"a": 3.468, "c": 11.241},  # crystallographic
        # "hexagonal": {"a": 3.468, "c": 11.241},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.894},  # allotropic and high-pressure modifications
        # "orthorhombic": {"a": 3.063, "b": 5.968, "c": 5.169},  # allotropic and high-pressure modifications
    },
    # Cm (96)
    "Cm": {
        "hexagonal": {"a": 3.496, "c": 11.331},  # crystallographic
        # "hexagonal": {"a": 3.496, "c": 11.331},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.381},  # allotropic and high-pressure modifications
    },
    # Bk (97)
    "Bk": {
        "hexagonal": {"a": 3.416, "c": 11.069},  # crystallographic
        # "hexagonal": {"a": 3.416, "c": 11.069},  # allotropic and high-pressure modifications
        # "fcc": {"a": 4.997},  # allotropic and high-pressure modifications
    },
    # Cf (98)
    # None
    # Es (99)
    # None
    # Fm (100)
    # None
    # Md (101)
    # None
    # No (102)
    # None
    # Lr (103)
    # None
    # Rf (104)
    # None
    # Db (105)
    # None
    # Sg (106)
    # None
    # Bh (107)
    # None
    # Hs (108)
    # None
    # Mt (109)
    # None
    # Ds (110)
    # None
}
