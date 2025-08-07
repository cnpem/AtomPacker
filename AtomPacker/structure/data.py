# This source code is part of the AtomPacker package and is distributed
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

from warnings import warn


def _get_bcc(atom_type: str) -> float | None:
    """
    Get the lattice constant for BCC lattice.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.

    Returns
    -------
    float | None
        The lattice constant `a` for BCC lattice. If None, use
        `ase.data.reference_states`.

    Warns
    -----
    UserWarning
        If the lattice constant `a` is not found in `lattice_constants`.
    """
    try:
        return lattice_constants[atom_type]["bcc"]["a"]
    except KeyError:
        warn(
            f"Constant `a` for {atom_type} not found for `bcc` in \
lattice_constants. Using `ase.data.reference_states`."
        )
        return None


def _get_fcc(atom_type: str) -> float | None:
    """
    Get the lattice constant for FCC lattice.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.

    Returns
    -------
    float | None
        The lattice constant `a` for FCC lattice. If None, use
        `ase.data.reference_states`.

    Warns
    -----
    UserWarning
        If the lattice constant `a` is not found in `lattice_constants`.
    """
    try:
        return lattice_constants[atom_type]["fcc"]["a"]
    except KeyError:
        warn(
            f"Constant `a` for {atom_type} not found for `fcc` in \
lattice_constants. Using `ase.data.reference_states`."
        )
        return None


def _get_hcp(atom_type: str) -> tuple[float, float] | None:
    """
    Get the lattice constant for HCP lattice.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.

    Returns
    -------
    Tuple[float, float] | None
        The lattice constants `(a, c)` for HCP lattice. If None, use
        `ase.data.reference_states`.

    Warns
    -----
    UserWarning
        If the lattice constants `a` and `c` are not found in
        `lattice_constants`.
    """
    try:
        return (
            lattice_constants[atom_type]["hcp"]["a"],
            lattice_constants[atom_type]["hcp"]["c"],
        )
    except KeyError:
        warn(
            f"Constants `a` and `c` for {atom_type} not found for `hcp` \
in lattice_constants. Using `ase.data.reference_states`."
        )
        return None


def _get_sc(atom_type: str) -> float | None:
    """
     Get the lattice constant for SC lattice.

    Parameters
     ----------
     atom_type : str
        The atomic symbol of the atom.

     Returns
     -------
     float | None
         The lattice constant `a` for SC lattice. If None, use
        `ase.data.reference_states`.
    """
    try:
        return lattice_constants[atom_type]["sc"]["a"]
    except KeyError:
        warn(
            f"Constant `a` for {atom_type} not found for `sc` in \
lattice_constants. Using `ase.data.reference_states`."
        )
        return None


def get_lattice_constants(
    atom_type: str, lattice_type: str
) -> tuple[float, float, float] | tuple[float, float] | float | None:
    """
    Get the lattice constants for a given atom type and lattice type.

    Parameters
    ----------
    atom_type : str
        The atomic symbol of the atom.
    lattice_type : str
        The type of lattice in the cluster. The available lattice types are
        'bcc', 'fcc', 'hcp' and 'sc'.

    Returns
    -------
    tuple[float, float, float] | tuple[float, float] | float | None
        The lattice constants for the given atom type and lattice type. If
        None, use `ase.data.reference_states`.

    Note
    ----
        The lattice constants will be fetched from
        `AtomPacker.data.lattice_constants` if available. If not,
        the experimental values from `ase.data` will be used.

    Raises
    ------
    ValueError
        If the lattice type is not valid.
    """
    match lattice_type:
        case "bcc":
            return _get_bcc(atom_type)
        case "fcc":
            return _get_fcc(atom_type)
        case "hcp":
            return _get_hcp(atom_type)
        case "sc":
            return _get_sc(atom_type)
        case _:
            raise ValueError(f"Invalid lattice type: {lattice_type}.")


"""
The lattice constants (in Å) from Martienssen & Warlimont [1]_.

Note
----
The standard phase (X) is the one used in the `get_lattice_constants` function.
The different phases, polymorphs, or structural modifications (e.g., α-X,
β-X, γ-X, δ-X, ε-X, ζ-X, X-I, X-II, X-III, X-IV) are included in this
dictionary for reference purposes only.

References
----------
.. [1] Martienssen, W., & Warlimont, H. (Eds.). (2005). Springer Handbook of
    Condensed Matter and Materials Data. Springer Berlin Heidelberg.
    https://doi.org/10.1007/3-540-30437-1.
"""
lattice_constants: dict[str, dict[str, object]] = {
    # H (1)
    "H": {
        "hexagonal": {
            "a": 3.771,
            "c": 6.156,
        },  # 4.2 K
    },
    "α-H": {
        "fcc": {"a": 5.334},  # <1.25 K
    },
    "β-H": {
        "hcp": {
            "a": 3.771,
            "c": 6.156,
        },  # <13.81 K
    },
    # He (2)
    "He": {
        "hexagonal": {
            "a": 3.577,
            "c": 5.842,
        },  # 1.5 K
    },
    "α-He": {
        "hcp": {
            "a": 3.577,
            "c": 5.842,
        },  # <0.95 K
    },
    "β-He": {
        "fcc": {"a": 4.240},  # 1.6 K and 0.125 GPa
    },
    "γ-He": {
        "bcc": {"a": 4.110},  # 1.73 K and 0.030 GPa
    },
    # Li (3)
    "Li": {
        "bcc": {"a": 3.5093},  # crystallographic
    },
    "α-Li": {
        "hcp": {
            "a": 3.111,
            "c": 5.093,
        },  # <72 K
    },
    "β-Li": {
        "bcc": {"a": 3.093},  # room temperature
    },
    "γ-Li": {"fcc": {"a": None}},
    # Be (4)
    "Be": {
        "hexagonal": {
            "a": 2.2857,
            "c": 3.5839,
        },  # crystallographic
    },
    "α-Be": {
        "hcp": {
            "a": 2.2857,
            "c": 3.5839,
        },  # room temperature
    },
    "β-Be": {
        "bcc": {"a": 2.5515},
    },  # >1523 K
    # B (5)
    # None
    # C (6)
    "C": {
        "diamond": {"a": 3.5671},  # room temperature and >60 GPa
        "graphite": {
            "a": 2.4612,
            "c": 6.7090,
        },  # room temperature and standard pressure
    },
    # N (7)
    "N": {
        "sc": {"a": 5.659},  # 4.2 K
    },
    "α-N": {
        "sc": {"a": 5.659},
    },  # <20 K
    "β-N": {
        "hexagonal": {
            "a": 4.046,
            "c": 6.629,
        },  # > 35.6 K
    },
    "γ-N": {
        "tetragonal": {
            "a": 3.957,
            "c": 5.101,
        },  # 20 K and >3.3 GPa
    },
    # O (8)
    "O": {
        "monoclinic": {
            "a": 5.403,
            "b": 3.429,
            "c": 5.086,
            "gamma": 132.53,
        },  # <23 K
    },
    "α-O": {
        "monoclinic": {
            "a": 5.403,
            "b": 3.429,
            "c": 5.086,
            "gamma": 132.53,
        },  # <23 K
    },
    "β-O": {
        "rhombohedral": {
            "a": 4.210,
            "alpha": 46.27,
        },  # >23.9 K
    },
    "γ-O": {
        "sc": {"a": 6.830},  # >43.6 K
    },
    # F (9)
    "F": {
        "monoclinic": {
            "a": 5.50,
            "b": 3.28,
            "c": 7.28,
            "beta": 102.17,
        },  # <45.6 K
    },
    "α-F": {
        "monoclinic": {
            "a": 5.50,
            "b": 3.28,
            "c": 7.28,
            "beta": 102.17,
        },  # <45.6 K
    },
    "β-F": {
        "cubic": {"a": 6.67},  # >45.6 K
    },
    # Ne (10)
    "Ne": {
        "fcc": {"a": 4.4622},  # 4.2 K
    },
    # Na (11)
    "Na": {
        "bcc": {"a": 4.2096},  # crystallographic
    },
    "α-Na": {
        "hcp": {
            "a": 3.767,
            "c": 6.154,
        },  # <36 K
    },
    "β-Na": {
        "bcc": {"a": 4.2096},  # room temperature
    },
    # Mg (12)
    "Mg": {
        "hexagonal": {
            "a": 3.2093,
            "c": 5.2107,
        },  # crystallographic
    },
    # Al (13)
    "Al": {
        "fcc": {"a": 3.6149},  # crystallographic
    },
    "α-Al": {
        "fcc": {"a": 4.0496},  # room temperature and standard pressure
    },
    "β-Al": {
        "hcp": {
            "a": 2.693,
            "c": 4.398,
        },  # >20.5 GPa
    },
    # Si (14)
    "Si": {
        "diamond": {"a": 5.43102},  # 295.65 K
    },
    "α-Si": {
        "fcc": {"a": 5.4306},  # room temperature and standard pressure
    },
    "β-Si": {
        "tetragonal": {
            "a": 4.686,
            "c": 2.585,
        },  # >9.5 GPa
    },
    "γ-Si": {
        "sc": {"a": 6.3600},  # >16.0 GPa
    },
    "δ-Si": {
        "hexagonal": {
            "a": 3.80,
            "c": 6.28,
        },  # decompressed β-Si
    },
    # P (15)
    "P": {
        "orthorhombic": {
            "a": 3.3136,
            "b": 10.478,
            "c": 4.3763,
        },  # crystallographic
    },
    # S (16)
    "S": {
        "orthorhombic": {
            "a": 10.464,
            "b": 12.8660,
            "c": 24.4860,
        },  # crystallographic
    },
    # Cl (17)
    "Cl": {
        "orthorhombic": {
            "a": 6.24,
            "b": 4.48,
            "c": 8.26,
        },  # 113 K
    },
    # Ar (18)
    "Ar": {
        "fcc": {"a": 5.312},  # 4.2 K
    },
    "α-Ar": {
        "fcc": {"a": 5.312},  # <83.8 K
    },
    "β-Ar": {
        "hcp": {
            "a": 3.760,
            "c": 6.141,
        },  # >83.8 K
    },
    # K (19)
    "K": {
        "bcc": {"a": 5.321},  # crystallographic
    },
    # Ca (20)
    "Ca": {
        "fcc": {"a": 5.5884},  # crystallographic
    },
    "α-Ca": {
        "fcc": {"a": 5.5884},  # room temperature
    },
    "γ-Ca": {
        "bcc": {"a": 4.480},  # >1010 K
    },
    # Sc (21)
    "Sc": {
        "hexagonal": {
            "a": 3.3088,
            "c": 5.2680,
        },  # crystallographic
    },
    "α-Sc": {
        "hcp": {
            "a": 3.3088,
            "c": 5.2680,
        },  # room temperature
    },
    "β-Sc": {
        "bcc": {"a": None},  # >1607 K
    },
    # Ti (22)
    "Ti": {
        "hexagonal": {
            "a": 2.9503,
            "c": 4.6836,
        },  # crystallographic
    },
    "α-Ti": {
        "hcp": {
            "a": 2.9503,
            "c": 4.6836,
        },  # room temperature
    },
    "β-Ti": {
        "bcc": {"a": 3.3065},  # >1173 K
    },
    # V (23)
    "V": {
        "bcc": {"a": 3.0238},  # crystallographic
    },
    # Cr (24)
    "Cr": {
        "bcc": {"a": 2.8847},  # crystallographic
    },
    "α-Cr": {
        "bcc": {"a": 2.8847},  # room temperature and standard pressure
    },
    "β-Cr": {
        "bcc": {"a": 2.8820},  # high pressure
    },
    # Mn (25)
    "Mn": {
        "bcc": {"a": 8.9219},  # crystallographic
    },
    "α-Mn": {
        "bcc": {"a": 8.9219},  # room temperature and standard pressure
    },
    "β-Mn": {
        "sc": {"a": 6.3152},  # >1000 K
    },
    "γ-Mn": {
        "fcc": {"a": 3.8624},  # >1368 K
    },
    "δ-Mn": {
        "bcc": {"a": 3.0806},  # 1408 K
    },
    # Fe (26)
    "Fe": {
        "bcc": {"a": 2.8665},  # crystallographic
    },
    "α-Fe": {
        "bcc": {"a": 2.8665},  # room temperature and standard pressure
    },
    "β-Fe": {
        "fcc": {"a": 3.6467},  # >1183 K
    },
    "δ-Fe": {
        "bcc": {"a": 2.9135},  # 1663 K
    },
    "ε-Fe": {
        "hcp": {
            "a": 2.485,
            "c": 3.990,
        },  # 13.0 GPa
    },
    # Co (27)
    "Co": {
        "hcp": {
            "a": 2.5071,
            "c": 4.0694,
        },  # crystallographic
    },
    "α-Co": {
        "fcc": {"a": 3.5445},  # >661 K
    },
    "ε-Co": {
        "hcp": {
            "a": 2.5071,
            "c": 4.0694,
        },  # room temperature and standard pressure
    },
    # Ni (28)
    "Ni": {
        "fcc": {"a": 3.5241},  # crystallographic
    },
    # Cu (29)
    "Cu": {
        "fcc": {"a": 3.6149},  # 298 K
    },
    # Zn (30)
    "Zn": {
        "hexagonal": {
            "a": 2.6644,
            "c": 4.9494,
        },  # 298 K
    },
    # Ga (31)
    "Ga": {
        "orthorhombic": {
            "a": 4.5192,
            "b": 7.6586,
            "c": 4.5258,
        },  # crystallographic
    },
    "α-Ga": {
        "orthorhombic": {
            "a": 4.5192,
            "b": 7.6586,
            "c": 4.5258,
        },  # room temperature and standard pressure
    },
    "β-Ga": {
        "tetragonal": {
            "a": 2.808,
            "c": 4.458,
        },  # >1.2 GPa
    },
    "γ-Ga": {
        "orthorhombic": {
            "a": 10.593,
            "b": 13.523,
            "c": 5.203,
        },  # 220 K and >3.0 GPa
    },
    # Ge (32)
    "Ge": {
        "diamond": {"a": 5.659},  # room temperature
    },
    "α-Ge": {
        "fcc": {"a": 5.6574},  # room temperature and standard pressure
    },
    "β-Ge": {
        "tetragonal": {
            "a": 4.884,
            "c": 2.692,
        },  # >12.0 GPa
    },
    "γ-Ge": {
        "tetragonal": {
            "a": 5.93,
            "c": 6.98,
        },  # decompressed β-Ge
    },
    "δ-Ge": {
        "bcc": {"a": 6.92},  # >12.0 GPa
    },
    # As (33)
    "As": {
        "rhombohedral": {
            "a": 4.1320,
            "alpha": 54.12,
        },  # crystallographic
    },
    "α-As": {
        "rhombohedral": {
            "a": 4.1320,
            "alpha": 54.12,
        },
    },
    "ε-As": {
        "orthorhombic": {
            "a": 3.62,
            "b": 10.85,
            "c": 4.48,
        },  # >721 K
    },
    # Se (34)
    "Se": {
        "hexagonal": {
            "a": 4.3655,
            "c": 4.9576,
        },  # crystallographic
    },
    "α-Se": {
        "monoclinic": {
            "a": 9.054,
            "b": 9.083,
            "c": 2.336,
            "gamma": 90.82,
        },  # room temperature
    },
    "β-Se": {
        "monoclinic": {
            "a": 15.018,
            "b": 14.713,
            "c": 8.879,
            "gamma": 93.6,
        },  # room temperature
    },
    "γ-Se": {
        "hexagonal": {
            "a": 4.3655,
            "c": 4.9576,
        },  # room temperature
    },
    # Br (35)
    "Br": {
        "orthorhombic": {
            "a": 6.68,
            "b": 4.49,
            "c": 8.74,
        },  # 123 K
    },
    # Kr (36)
    "Kr": {
        "fcc": {"a": 5.6459},  # 4.2 K
    },
    # Rb (37)
    "Rb": {
        "bcc": {"a": 5.703},  # crystallographic
    },
    # Sr (38)
    "Sr": {
        "fcc": {"a": 6.084},  # crystallographic
    },
    "α-Sr": {
        "fcc": {"a": 6.084},  # room temperature
    },
    "β-Sr": {
        "hcp": {
            "a": 4.280,
            "c": 7.050,
        },  # >486 K
    },
    "γ-Sr": {
        "bcc": {"a": 4.870},  # >878 K
    },
    "Sr-II": {
        "bcc": {"a": 4.437},  # >3.5 GPa
    },
    # Y (39)
    "Y": {
        "hexagonal": {
            "a": 3.6482,
            "c": 5.7318,
        },  # crystallographic
    },
    "α-Y": {
        "hcp": {
            "a": 3.6482,
            "c": 5.7318,
        },  # room temperature
    },
    "β-Y": {
        "bcc": {"a": None},  # >1752 K
    },
    # Zr (40)
    "Zr": {
        "hexagonal": {
            "a": 3.2317,
            "c": 5.1476,
        },  # crystallographic
    },
    "α-Zr": {
        "hcp": {
            "a": 3.2317,
            "c": 5.1476,
        },  # room temperature
    },
    "β-Zr": {
        "bcc": {"a": 3.609},  # >1138 K
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
        "hexagonal": {
            "a": 2.738,
            "c": 4.394,
        },  # crystallographic
    },
    # Ru (44)
    "Ru": {
        "hexagonal": {
            "a": 2.7053,
            "c": 4.2814,
        },  # crystallographic
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
        "fcc": {"a": 4.0861},  # 298 K
    },
    # Cd (48)
    "Cd": {
        "hexagonal": {
            "a": 2.9788,
            "c": 5.6167,
        },  # 298 K
    },
    # In (49)
    "In": {
        "tetragonal": {
            "a": 4.5990,
            "c": 4.9470,
        },  # crystallographic
    },
    # Sn (50)
    "Sn": {
        "tetragonal": {"a": 58.197},  # 300 K
    },
    "α-Sn": {
        "fcc": {"a": 56.4892},  # <291 K
    },
    "β-Sn": {
        "tetragonal": {
            "a": 5.8316,
            "c": 3.1815,
        },  # room temperature
    },
    "γ-Sn": {
        "tetragonal": {
            "a": 3.70,
            "c": 3.37,
        },  # >9.0 GPa
    },
    # Sb (51)
    "Sb": {
        "rhombohedral": {
            "a": 4.5065,
            "alpha": 57.11,
        },  # crystallographic
    },
    "Sb-I": {
        "rhombohedral": {
            "a": 4.5065,
            "alpha": 57.11,
        },  # room temperature and standard pressure
    },
    "Sb-II": {
        "sc": {"a": 2.992},  # >5.0 GPa
    },
    "Sb-III": {
        "hcp": {
            "a": 3.376,
            "c": 5.341,
        },  # >7.5 GPa
    },
    "Sb-IV": {
        "monoclinic": {
            "a": 5.56,
            "b": 4.04,
            "c": 4.22,
            "beta": 86.0,
        },  # 14.0 GPa
    },
    # Te (52)
    "Te": {
        "hexagonal": {
            "a": 4.4561,
            "c": 5.9271,
        },  # crystallographic
    },
    "α-Te": {
        "hexagonal": {
            "a": 4.4561,
            "c": 5.9271,
        },  # room temperature and standard pressure
    },
    "β-Te": {
        "rhombohedral": {
            "a": 4.69,
            "alpha": 53.30,
        },  # >2.0 GPa
    },
    "γ-Te": {
        "rhombohedral": {
            "a": 3.002,
            "alpha": 103.3,
        },  # >7.0 GPa
    },
    # I (53)
    "I": {
        "orthorhombic": {
            "a": 7.268,
            "b": 4.797,
            "c": 9.797,
        },  # crystallographic
    },
    # Xe (54)
    "Xe": {
        "fcc": {"a": 6.132},  # 4.2 K
    },
    # Cs (55)
    "Cs": {
        "bcc": {"a": 6.141},  # crystallographic
    },
    "Cs-I": {
        "bcc": {"a": 6.141},  # room temperature and standard pressure
    },
    "Cs-II": {
        "fcc": {"a": 5.984},  # >2.37 GPa
    },
    "Cs-III": {
        "fcc": {"a": 5.800},  # >4.22 GPa
    },
    # Ba (56)
    "Ba": {
        "bcc": {"a": 5.023},  # crystallographic
    },
    # La (57)
    "La": {
        "hexagonal": {
            "a": 3.7740,
            "c": 12.171,
        },  # crystallographic
    },
    "α-La": {
        "hexagonal": {
            "a": 3.7740,
            "c": 12.171,
        },  # room temperature and standard pressure
    },
    "β-La": {
        "fcc": {"a": 5.3045},  # >613 K
    },
    "γ-La": {
        "bcc": {"a": 4.265},  # >1141 K
    },
    "β-2-La": {
        "fcc": {"a": 5.170},  # >2.0 GPa
    },
    # Ce (58)
    "Ce": {
        "fcc": {"a": 5.1610},  # crystallographic
    },
    "α-Ce": {
        "fcc": {"a": 5.1610},  # room temperature and standard pressure
    },
    "β-Ce": {
        "hexagonal": {
            "a": 3.673,
            "c": 11.802,
        },  # >263 K
    },
    "γ-Ce": {
        "fcc": {"a": None},  # <95 K
    },
    "α-2-Ce": {
        "fcc": {"a": 4.82},  # >1.5 GPa
    },
    "Ce-III": {
        "orthorhombic": {
            "a": None,
            "b": None,
            "c": None,
        },  # 5.1 GPa
    },
    # Pr (59)
    "Pr": {
        "hexagonal": {
            "a": 3.6721,
            "c": 11.8326,
        },  # crystallographic
    },
    "α-Pr": {
        "hexagonal": {
            "a": 3.6721,
            "c": 11.8326,
        },  # room temperature and standard pressure
    },
    "β-Pr": {
        "bcc": {"a": 4.13},  # >1094 K
    },
    "γ-Pr": {
        "fcc": {"a": 4.88},  # >4.0 GPa
    },
    # Nd (60)
    "Nd": {
        "hexagonal": {
            "a": 3.6582,
            "c": 11.7966,
        },  # crystallographic
    },
    "α-Nd": {
        "hexagonal": {
            "a": 3.6582,
            "c": 11.7966,
        },  # room temperature and standard pressure
    },
    "β-Nd": {
        "bcc": {"a": 4.13},  # >1135 K
    },
    "γ-Nd": {
        "fcc": {"a": 4.80},  # >5.0 GPa
    },
    # Pm (61)
    "Pm": {
        "hexagonal": {
            "a": 3.65,
            "c": 11.65,
        },  # crystallographic
    },
    "α-Pm": {
        "hexagonal": {
            "a": 3.65,
            "c": 11.65,
        },  # room temperature and standard pressure
    },
    "β-Pm": {
        "bcc": {"a": None},  # >1163 K
    },
    # Sm (62)
    "Sm": {
        "hexagonal": {
            "a": 3.6290,
            "c": 26.207,
        },  # crystallographic
    },
    "α-Sm": {
        "trigonal": {
            "a": 3.629,
            "c": 26.207,
        },  # room temperature and standard pressure
    },
    "β-Sm": {
        "bcc": {"a": None},  # >1190 K
    },
    "γ-Sm": {
        "hexagonal": {
            "a": 3.6180,
            "c": 11.66,
        },  # >4.0 GPa
    },
    # Eu (63)
    "Eu": {
        "bcc": {"a": 4.5827},  # crystallographic
    },
    # Gd (64)
    "Gd": {
        "hexagonal": {
            "a": 3.6336,
            "c": 5.7810,
        },  # crystallographic
    },
    "α-Gd": {
        "hcp": {
            "a": 3.6336,
            "c": 5.7810,
        },  # room temperature and standard pressure
    },
    "β-Gd": {
        "bcc": {"a": 4.06},  # >1535 K
    },
    "γ-Gd": {
        "trigonal": {
            "a": 3.61,
            "c": 26.03,
        },  # >3.0 GPa
    },
    # Tb (65)
    "Tb": {
        "hexagonal": {
            "a": 3.6055,
            "c": 5.6966,
        },  # crystallographic
    },
    "α-Tb": {
        "hcp": {
            "a": 3.6055,
            "c": 5.6966,
        },  # room temperature and standard pressure
    },
    "β-Tb": {
        "bcc": {"a": None},  # >1589 K
    },
    "γ-Tb": {
        "trigonal": {
            "a": 3.41,
            "c": 24.50,
        },  # >6.0 G
    },
    # Dy (66)
    "Dy": {
        "hexagonal": {
            "a": 3.5915,
            "c": 5.6501,
        },  # crystallographic
    },
    "α-Dy": {
        "hcp": {
            "a": 3.5915,
            "c": 5.6501,
        },  # room temperature and standard pressure
    },
    "β-Dy": {
        "bcc": {"a": None},  # >1243 K
    },
    "α-2-Dy": {
        "orthorhombic": {
            "a": 3.595,
            "b": 6.184,
            "c": 5.678,
        },  # <86 K
    },
    "γ-Dy": {
        "trigonal": {
            "a": 3.436,
            "c": 24.83,
        },  # >7.5 GPa
    },
    # Ho (67)
    "Ho": {
        "hexagonal": {
            "a": 3.5778,
            "c": 5.6178,
        },  # crystallographic
    },
    "α-Ho": {
        "hcp": {
            "a": 3.5778,
            "c": 5.6178,
        },  # room temperature and standard pressure
    },
    "β-Ho": {
        "bcc": {"a": None},  # high temperature
    },
    "γ-Ho": {
        "trigonal": {
            "a": 3.34,
            "c": 24.50,
        },  # >4.0 GPa
    },
    # Er (68)
    "Er": {
        "hexagonal": {
            "a": 3.5592,
            "c": 5.5850,
        },  # crystallographic
    },
    "α-Er": {
        "hcp": {
            "a": 3.5592,
            "c": 5.5850,
        },  # room temperature and standard pressure
    },
    "β-Er": {
        "bcc": {"a": None},  # high temperature
    },
    # Tm (69)
    "Tm": {
        "hexagonal": {
            "a": 3.5375,
            "c": 5.5540,
        },  # crystallographic
    },
    "α-Tm": {
        "hcp": {
            "a": 3.5375,
            "c": 5.5540,
        },  # room temperature and standard pressure
    },
    "β-Tm": {
        "bcc": {"a": None},  # high temperature
    },
    "Tm-II": {
        "trigonal": {
            "a": None,
            "c": None,
        },  # >6.0 GPa
    },
    # Yb (70)
    "Yb": {
        "fcc": {"a": 5.4848},  # crystallographic
    },
    "α-Yb": {
        "fcc": {"a": 5.4848},  # room temperature and standard pressure
    },
    "β-Yb": {
        "bcc": {"a": 4.44},  # >1005 K
    },
    "γ-Yb": {
        "hcp": {
            "a": 3.8799,
            "c": 6.3859,
        },  # <270 K
    },
    # Lu (71)
    "Lu": {
        "hexagonal": {
            "a": 3.5052,
            "c": 5.5494,
        },  # crystallographic
    },
    "α-Lu": {
        "hcp": {
            "a": 3.5052,
            "c": 5.5494,
        },  # room temperature and standard pressure
    },
    "β-Lu": {
        "bcc": {"a": None},  # >1005 K
    },
    "Lu-II": {
        "triangular": {
            "a": None,
            "c": None,
        },  # >23.0 GPa
    },
    # Hf (72)
    "Hf": {
        "hexagonal": {
            "a": 3.1946,
            "c": 5.0511,
        },  # crystallographic
    },
    "α-Hf": {
        "hcp": {
            "a": 3.1946,
            "c": 5.0511,
        },  # room temperature and standard pressure
    },
    "β-Hf": {
        "bcc": {"a": 3.610},  # >2268 K
    },
    "Hf-II": {
        "hexagonal": {
            "a": None,
            "c": None,
        },  # >38.8 GPa
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
        "hexagonal": {
            "a": 2.7608,
            "c": 4.4580,
        },  # crystallographic
    },
    # Os (76)
    "Os": {
        "hexagonal": {
            "a": 2.7348,
            "c": 4.3913,
        },  # crystallographic
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
        "fcc": {"a": 4.0784},  # 298 K
    },
    # Hg (80)
    "Hg": {
        "rhombohedral": {
            "a": 3.005,
            "alpha": 70.53,
        },  # 225 K
    },
    "α-Hg": {
        "rhombohedral": {
            "a": 3.005,
            "alpha": 70.53,
        },  # 225 K
    },
    "β-Hg": {
        "tetragonal": {
            "a": 3.995,
            "c": 2.825,
        },  # 77 K and high pressure
    },
    # Tl (81)
    "Tl": {
        "hexagonal": {
            "a": 3.4563,
            "c": 5.5263,
        },  # crystallographic
    },
    "α-Tl": {
        "hcp": {
            "a": 3.4563,
            "c": 5.5263,
        },  # room temperature
    },
    "β-Tl": {
        "bcc": {"a": 3.4563},  # >503 K
    },
    "γ-Tl": {
        "fcc": {"a": None},  # high pressure
    },
    # Pb (82)
    "Pb": {
        "fcc": {"a": 4.9502},  # crystallographic
    },
    "Pb-I": {
        "fcc": {"a": 4.9502},  # room temperature and standard pressure
    },
    "Pb-II": {
        "hcp": {
            "a": 3.265,
            "c": 5.387,
        },  # >10.3 GPa
    },
    # Bi (83)
    "Bi": {
        "rhombohedral": {
            "a": 4.7460,
            "alpha": 57.23,
        },  # crystallographic
    },
    "α-Bi": {
        "rhombohedral": {
            "a": 4.7460,
            "alpha": 57.23,
        },  # room temperature and standard pressure
    },
    "β-Bi": {
        "monoclinic": {
            "a": None,
            "b": None,
            "c": None,
        },  # >0.28 GPa
    },
    "γ-Bi": {
        "monoclinic": {
            "a": 6.05,
            "b": 4.20,
            "c": 4.65,
        },  # >0.28 GPa
    },
    "ζ-Bi": {
        "bcc": {"a": 38.00},  # >9.0 GPa
    },
    # Po (84)
    "Po": {
        "sc": {"a": 3.366},  # crystallographic
    },
    "α-Po": {
        "sc": {"a": 3.366},  # room temperature and standard pressure
    },
    "β-Po": {
        "rhombohedral": {
            "a": 3.373,
            "alpha": 98.08,
        },  # >327 K
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
    },
    "α-Th": {
        "fcc": {"a": 5.0851},  # room temperature and standard pressure
    },
    "β-Th": {
        "bcc": {"a": 4.11},  # >1673 K
    },
    # Pa (91)
    "Pa": {
        "tetragonal": {
            "a": 3.945,
            "c": 3.242,
        },  # crystallographic
    },
    "α-Pa": {
        "tetragonal": {
            "a": 3.945,
            "c": 3.242,
        },  # room temperature and standard pressure
    },
    "β-Pa": {
        "bcc": {"a": None},  # >1443 K
    },
    # U (92)
    "U": {
        "orthorhombic": {
            "a": 2.8538,
            "b": 5.8680,
            "c": 4.9557,
        },  # crystallographic
    },
    "α-U": {
        "orthorhombic": {
            "a": 2.8538,
            "b": 5.8680,
            "c": 4.9557,
        },  # room temperature and standard pressure
    },
    "β-U": {
        "tetragnonal": {
            "a": 10.759,
            "c": 5.654,
        },  # >935 K
    },
    "γ-U": {
        "bcc": {"a": 3.524},  # >1045 K
    },
    # Np (93)
    "Np": {
        "orthorhombic": {
            "a": 6.663,
            "b": 4.723,
            "c": 4.887,
        },  # crystallographic
    },
    "α-Np": {
        "orthorhombic": {
            "a": 6.663,
            "b": 4.723,
            "c": 4.887,
        },  # room temperature and standard pressure
    },
    "β-Np": {
        "tetragonal": {
            "a": 4.896,
            "c": 3.387,
        },  # >553 K
    },
    "γ-Np": {
        "bcc": {"a": 3.52},  # >850 K
    },
    # Pu (94)
    "Pu": {
        "monoclinic": {
            "a": 6.183,
            "b": 4.822,
            "c": 10.968,
            "gamma": 101.78,
        },  # crystallographic
    },
    "α-Pu": {
        "monoclinic": {
            "a": 6.183,
            "b": 4.822,
            "c": 10.968,
            "gamma": 101.78,
        },  # room temperature and standard pressure
    },
    "β-Pu": {
        "monoclinic": {
            "a": 9.284,
            "b": 10.463,
            "c": 7.859,
            "gamma": 92.13,
        },  # >395 K
    },
    "γ-Pu": {
        "orthorhombic": {
            "a": 3.1587,
            "b": 5.7682,
            "c": 10.162,
        },  # >508 K
    },
    "δ-Pu": {
        "fcc": {"a": 4.6371},  # >592 K
    },
    "δ-2-Pu": {
        "tetragonal": {
            "a": 3.3261,
            "c": 4.4630,
        },  # >726 K
    },
    "ε-Pu": {
        "bcc": {"a": 5.703},  # >744 K
    },
    # Am (95)
    "Am": {
        "hexagonal": {
            "a": 3.468,
            "c": 11.241,
        },  # crystallographic
    },
    "α-Am": {
        "hexagonal": {
            "a": 3.468,
            "c": 11.241,
        },  # room temperature and standard pressure
    },
    "β-Am": {
        "fcc": {"a": 4.894},  # >878 K
    },
    "γ-Am": {
        "orthorhombic": {
            "a": 3.063,
            "b": 5.968,
            "c": 5.169,
        },  # >15.0 GPa
    },
    # Cm (96)
    "Cm": {
        "hexagonal": {"a": 3.496, "c": 11.331},  # crystallographic
    },
    "α-Cm": {
        "hexagonal": {
            "a": 3.496,
            "c": 11.331,
        },  # room temperature and standard pressure
    },
    "β-Cm": {
        "fcc": {"a": 4.381},  # >1449 K
    },
    # Bk (97)
    "Bk": {
        "hexagonal": {"a": 3.416, "c": 11.069},  # crystallographic
    },
    "α-Bk": {
        "hexagonal": {
            "a": 3.416,
            "c": 11.069,
        },  # room temperature and standard pressure
    },
    "β-Bk": {
        "fcc": {"a": 4.997},  # >1183 K
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
