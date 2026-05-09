from pathlib import Path
from typing import Literal, Optional, Dict, Any, List, Tuple, Union
import os

from backend.mcp_server import mcp
from backend.mcp_server.util.comm import generate_work_path, run_command


@mcp.tool()
def generate_bulk_structure(
    element: str,
    crystal_structure: Literal[
        "sc", "fcc", "bcc", "hcp", "diamond", "zincblende", "rocksalt"
    ] = "fcc",
    a: float = None,
    c: float = None,
    cubic: bool = False,
    orthorhombic: bool = False,
    file_format: Literal["cif", "poscar"] = "cif",
) -> Dict[str, Any]:
    """
    Generate a bulk crystal structure using ASE's `bulk` function.

    Args:
        element (str): The chemical symbol of the element (e.g., 'Cu', 'Si', 'NaCl').
        crystal_structure (str): The type of crystal structure to generate. Options include:
            - 'sc' (simple cubic), a is needed
            - 'fcc' (face-centered cubic), a is needed
            - 'bcc' (body-centered cubic), a is needed
            - 'hcp' (hexagonal close-packed), a is needed, if c is None, c will be set to sqrt(8/3) * a.
            - 'diamond' (diamond cubic structure), a is needed
            - 'zincblende' (zinc blende structure), a is needed, two elements are needed, e.g., 'GaAs'
            - 'rocksalt' (rock salt structure), a is needed, two elements are needed, e.g., 'NaCl'
        a (float, optional): Lattice constant in Angstroms. Required for all structures.
        c (float, optional): Lattice constant for the c-axis in Angstroms. Required for 'hcp' structure.
        cubic (bool, optional): If constructing a cubic supercell for fcc, bcc, diamond, zincblende, or rocksalt structures.
        orthorhombic (bool, optional): If constructing orthorhombic cell for 'hcp' structure.
        file_format (str, optional): The format of the output file. Options are 'cif' or 'poscar'. Default is 'cif'.

    Notes: all crystal need the lattice constant a, which is the length of the unit cell (or conventional cell).

    Returns:
        structure_file: The path to generated structure file.
        cell: The cell parameters of the generated structure as a list of lists.
        coordinate: The atomic coordinates of the generated structure as a list of lists.

    Examples:
    >>> # FCC Cu
    >>> cu_fcc = generate_bulk_structure('Cu', 'fcc', a=3.6)
    >>>
    >>> # HCP Mg with custom c-axis
    >>> mg_hcp = generate_bulk_structure('Mg', 'hcp', a=3.2, c=5.2, orthorhombic=True)
    >>>
    >>> # Diamond Si
    >>> si_diamond = generate_bulk_structure('Si', 'diamond', a=5.43, cubic=True)
    >>> # Zincblende GaAs
    >>> gaas_zincblende = generate_bulk_structure('GaAs', 'zincblende', a=5.65, cubic=True)
    """
    if a is None:
        raise ValueError(
            "Lattice constant 'a' must be provided for all crystal structures."
        )

    from ase.build import bulk

    special_params = {}

    if crystal_structure == "hcp":
        if c is not None:
            special_params["c"] = c
        special_params["orthorhombic"] = orthorhombic

    if crystal_structure in ["fcc", "bcc", "diamond", "zincblende"]:
        special_params["cubic"] = cubic

    try:
        structure = bulk(
            name=element,
            crystalstructure=crystal_structure,
            a=a,
            **special_params,
        )
    except Exception as e:
        raise ValueError(f"Generate structure failed: {str(e)}") from e

    work_path = generate_work_path(create=True)

    if file_format == "cif":
        structure_file = f"{work_path}/{element}_{crystal_structure}.cif"
        structure.write(structure_file, format="cif")
    elif file_format == "poscar":
        structure_file = f"{work_path}/{element}_{crystal_structure}.vasp"
        structure.write(structure_file, format="vasp")
    else:
        raise ValueError("Unsupported file format. Use 'cif' or 'poscar'.")

    return {
        "structure_file": str(Path(structure_file).absolute()),
        "cell": structure.get_cell().tolist(),
        "coordinate": structure.get_positions().tolist(),
    }


@mcp.tool()
def generate_molecule_structure(
    molecule_name: Literal[
        "PH3", "P2", "CH3CHO", "H2COH", "CS", "OCHCHO", "C3H9C", "CH3COF",
        "CH3CH2OCH3", "HCOOH", "HCCl3", "HOCl", "H2", "SH2", "C2H2", "C4H4NH",
        "CH3SCH3", "SiH2_s3B1d", "CH3SH", "CH3CO", "CO", "ClF3", "SiH4",
        "C2H6CHOH", "CH2NHCH2", "isobutene", "HCO", "bicyclobutane", "LiF",
        "Si", "C2H6", "CN", "ClNO", "S", "SiF4", "H3CNH2", "methylenecyclopropane",
        "CH3CH2OH", "F", "NaCl", "CH3Cl", "CH3SiH3", "AlF3", "C2H3", "ClF", "PF3",
        "PH2", "CH3CN", "cyclobutene", "CH3ONO", "SiH3", "C3H6_D3h", "CO2", "NO",
        "trans-butane", "H2CCHCl", "LiH", "NH2", "CH", "CH2OCH2", "C6H6",
        "CH3CONH2", "cyclobutane", "H2CCHCN", "butadiene", "C", "H2CO", "CH3COOH",
        "HCF3", "CH3S", "CS2", "SiH2_s1A1d", "C4H4S", "N2H4", "OH", "CH3OCH3",
        "C5H5N", "H2O", "HCl", "CH2_s1A1d", "CH3CH2SH", "CH3NO2", "Cl", "Be", "BCl3",
        "C4H4O", "Al", "CH3O", "CH3OH", "C3H7Cl", "isobutane", "Na", "CCl4",
        "CH3CH2O", "H2CCHF", "C3H7", "CH3", "O3", "P", "C2H4", "NCCN", "S2", "AlCl3",
        "SiCl4", "SiO", "C3H4_D2d", "H", "COF2", "2-butyne", "C2H5", "BF3", "N2O",
        "F2O", "SO2", "H2CCl2", "CF3CN", "HCN", "C2H6NH", "OCS", "B", "ClO",
        "C3H8", "HF", "O2", "SO", "NH", "C2F4", "NF3", "CH2_s3B1d", "CH3CH2Cl",
        "CH3COCl", "NH3", "C3H9N", "CF4", "C3H6_Cs", "Si2H6", "HCOOCH3", "O", "CCH",
        "N", "Si2", "C2H6SO", "C5H8", "H2CF2", "Li2", "CH2SCH2", "C2Cl4", "C3H4_C3v",
        "CH3COCH3", "F2", "CH4", "SH", "H2CCO", "CH3CH2NH2", "Li", "N2", "Cl2", "H2O2",
        "Na2", "BeH", "C3H4_C2v", "NO2", "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
        "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
        "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
        "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
        "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl",
        "Mc", "Lv", "Ts", "Og",
    ] = "H2O",
    cell: Optional[List[List[float]]] = None,
    vacuum: Optional[float] = 5.0,
    output_file_format: Literal["cif", "poscar", "abacus"] = "abacus",
) -> Dict[str, Any]:
    """
    Generate molecule structure from ase's collection of molecules or single atoms.
    Args:
        molecule_name: The name of the molecule or atom to generate. It can be a chemical symbol (e.g., 'H', 'O', 'C') or
                       a molecule name in g2 collection contained in ASE's collections.
        cell: The cell parameters for the generated structure. Default is a 10x10x10 Angstrom cell. Units in angstrom.
        vacuum: The vacuum space to add around the molecule. Default is 5.0 Angstrom.
        output_file_format: The format of the output file. Default is 'abacus'. 'poscar' represents POSCAR format used by VASP.
    Returns:
        A dictionary containing:
        - structure_file: The absolute path to the generated structure file.
        - cell: The cell parameters of the generated structure as a list of lists.
        - coordinate: The atomic coordinates of the generated structure as a list of lists.
    """
    from ase import Atoms
    from ase.build import molecule
    from ase.data import chemical_symbols
    from ase.collections import g2

    if cell is None:
        cell = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]

    format_to_use = output_file_format
    if format_to_use == "poscar":
        format_to_use = "vasp"

    if molecule_name in g2.names:
        atoms = molecule(molecule_name)
        atoms.set_cell(cell)
        atoms.center(vacuum=vacuum)
    elif molecule_name in chemical_symbols and molecule_name != "X":
        atoms = Atoms(symbol=molecule_name, positions=[[0, 0, 0]], cell=cell)
    else:
        raise ValueError(f"Molecule '{molecule_name}' not found in ASE collections or chemical symbols.")

    if output_file_format == "abacus":
        stru_file_path = Path(f"{molecule_name}.stru").absolute()
    else:
        stru_file_path = Path(f"{molecule_name}.{format_to_use}").absolute()

    atoms.write(stru_file_path, format=format_to_use)

    return {
        "structure_file": str(stru_file_path),
        "cell": atoms.get_cell().tolist(),
        "coordinate": atoms.get_positions().tolist(),
    }


@mcp.tool()
def search_properties_from_database(
    gas: Literal[
        "c2h4", "c2h2", "c3h8", "ch4", "cl2", "co", "co2",
        "h2", "hcn", "n2", "nh3", "no2", "so2",
    ],
    cation: Literal[
        "allylmim", "allylmpy", "allylPy", "amim", "bbbu",
        "bmpy", "c1mim", "C1Py", "c2mim", "C2Py", "c3mim", "C3Py", "c4mim",
        "C4Py", "c6mim", "C6Py", "c8mim", "c12mim", "c18mim", "ch", "cynaolmim",
        "cynaolmpy", "cynaolPy", "dea", "dmea", "empy", "epp", "epr", "epy",
        "hbbu", "hdbu", "Hmim", "hmpy", "hobmim", "mbn", "mdbu", "mea", "meoim",
        "mmpy", "n444", "n1122", "N2222", "nbmim", "nh444", "obbu", "p1444", "p4444",
        "p44414", "p44416", "phenethylmim", "phenethylmpy", "phenethylPy", "py5", "py6",
        "pyr14", "pyr14OCH3", "s222", "sbmim", "tea", "tmg",
    ],
    anion: Literal[
        "Alanine", "alcl4", "arginine", "asparagine-", "asparticacid", "bcn4",
        "be", "bf4", "bh2cn2", "br", "BuNHC3SO3-", "C8H9SO4-", "cf3coo", "CF3SO3-",
        "CH2OHCH2COO-", "CH3CH2COO-", "CH3CHOHCOO-", "chc", "cl", "clo4", "cpc",
        "cysteine", "dca", "diethNC3SO3-", "dimeNC3SO3-", "fecl4", "glutamicacid",
        "glutamine-", "Glycine-", "HCOO-", "HeNHC3SO3-", "hso4", "im", "isobuNHC3SO3-",
        "Isoleucine", "isoproNHC3SO3-", "Lac-", "Leucine", "lysine", "meso4", "no3",
        "oac", "pf6", "Pheneylalanine", "Proline", "scn", "Serine", "threonine-",
        "tryptophan_1", "tyrosine", "Valine",
    ],
    properties: List[
        Literal[
            "solvation_free_energy",
            "free_volume",
            "total_free_volume",
            "free_volume_fraction",
            "self_diffusion_coefficient",
            "binding_energy",
        ]
    ] = None,
) -> Dict[str, Union[str, float]]:
    """
    Search properties of gas in ionic liquid systems from the database.

    Args:
        gas (str): The gas species to search for. Options include 'c2h4', 'c2h2', 'c3h8', 'ch4', 'cl2', 'co', 'co2',
                   'h2', 'hcn', 'n2', 'nh3', 'no2', 'so2'.
        cation (str): The cation species to search for. Options include various imidazolium, pyridinium, and other cations.
        anion (str): The anion species to search for. Options include various carboxylates, sulfonates, and other anions.
        properties (list): List of properties to retrieve. Below properties are available:
            - "solvation_free_energy": Solvation free energy of the gas in the ionic liquid system, the unit is kcal/mol.
            - "free_volume": Free volume of the system, the unit is Angstrom^3.
            - "total_free_volume": Total free volume of the system, the unit is Angstrom^3.
            - "free_volume_fraction": Free volume fraction of the system, the unit is percentage.
            - "self_diffusion_coefficient": Self-diffusion coefficient of the gas in the ionic liquid system, the unit is m^2/s.
            - "binding_energy": Binding energy beween gas and cation_anion in ionic liquid system, the unit is eV.

    Returns:
        dict: A dictionary containing the requested properties for the specified gas, cation, and anion.
    """
    import pickle

    if properties is None:
        properties = ["solvation_free_energy"]

    dataset_path = os.environ.get("ILS4GAS_TRAIN_EB_PATH", None)
    if dataset_path is None:
        raise ValueError(
            "Environment variable ILS4GAS_TRAIN_EB_PATH is not set. "
            "Please set it to the path of the dataset."
        )

    dataset = pickle.load(open(dataset_path, "rb"))

    key_name = f"{gas}_{cation}_{anion}"
    if key_name not in dataset:
        raise ValueError(f"Data for {key_name} not found in the dataset.")

    values = {}
    for prop in properties:
        if prop in dataset[key_name]:
            values[prop] = dataset[key_name][prop]
        else:
            values[prop] = "Property not available"
    return values


@mcp.tool()
def search_dft_feature(
    mol_name: Literal[
        "c2h2", "allylmim", "Alanine", "alcl4", "arginine", "asparagine-",
        "asparticacid", "bcn4", "be", "bf4", "bh2cn2", "BuNHC3SO3-", "C8H9SO4-",
        "cf3coo", "CF3SO3-", "CH2OHCH2COO-", "CH3BENNSO3-", "CH3CH2COO-",
        "CH3CHOHCOO-", "chc", "clo4", "cpc", "cysteine", "dca", "diethNC3SO3-",
        "dimeNC3SO3-", "glutamicacid", "glutamine-", "Glycine-", "HCOO-",
        "HeNHC3SO3-", "hso4", "im", "isobuNHC3SO3-", "Isoleucine", "isoproNHC3SO3-",
        "Lac-", "Leucine", "lysine", "meso4", "MeSO4-", "Methionine", "no3", "oac",
        "pf6", "Pheneylalanine", "phenylalanine", "Proline", "scn", "Serine",
        "threonine-", "tryptophan_1", "tyrosine", "Valine", "fecl4",
        "1_Alanine", "1_alcl4", "1_arginine", "1_asparagine-", "1_asparticacid",
        "1_bcn4", "1_be", "1_bf4", "1_bh2cn2", "1_BuNHC3SO3-", "1_C8H9SO4-",
        "1_cf3coo", "1_CF3SO3-", "1_CH2OHCH2COO-", "1_CH3BENNSO3-", "1_CH3CH2COO-",
        "1_CH3CHOHCOO-", "1_chc", "1_clo4", "1_cpc", "1_cysteine", "1_dca",
        "1_diethNC3SO3-", "1_dimeNC3SO3-", "1_glutamicacid", "1_glutamine-",
        "1_Glycine-", "1_HCOO-", "1_HeNHC3SO3-", "1_hso4", "1_im", "1_isobuNHC3SO3-",
        "1_Isoleucine", "1_isoproNHC3SO3-", "1_Lac-", "1_Leucine", "1_lysine",
        "1_meso4", "1_MeSO4-", "1_Methionine", "1_no3", "1_oac", "1_pf6",
        "1_Pheneylalanine", "1_phenylalanine", "1_Proline", "1_scn", "1_Serine",
        "1_threonine-", "1_tryptophan_1", "1_tyrosine", "1_Valine", "1_fecl4",
        "amim", "bbbu", "bmpy", "c12mim", "c18mim", "c1mim", "C1Py", "c2mim", "C2Py",
        "c3mim", "C3Py", "c4mim", "C4Py", "c6mim", "C6Py", "c8mim", "ch", "cynaolmim",
        "cynaolmpy", "cynaolPy", "dea", "dmea", "empy", "epp", "epr", "epy", "hbbu",
        "hdbu", "Hmim", "hmpy", "hobmim", "mbn", "mdbu", "mea", "meoim", "mmpy",
        "n1122", "N2222", "n444", "nbmim", "nh444", "obbu", "p1444", "p44414",
        "p44416", "p4444", "phenethylmim", "phenethylmpy", "phenethylPy", "py5", "py6",
        "pyr14", "pyr14OCH3", "s222", "sbmim", "tea", "tmg",
        "c2h4", "c3h8", "ch4", "cl2", "co", "co2", "hcn", "n2", "nh3", "so2",
    ],
) -> Dict[str, float]:
    """
    Search DFT features of a molecule from the database.

    Args:
        mol_name (str): The name of the molecule to search for. Options include various small molecules and atoms.

    Returns:
        dict: A dictionary containing the DFT features of the specified molecule.

    The properties contained in the DFT feature dataset include:
        'AtomNum': int, the number of atoms in the molecule.
        'Weight': float, the molecular weight in amu.
        'Mol_Radius': float, the radius of the molecule in Angstrom.
        'Mol_Size_Short': float, the short size of the molecule in Angstrom.
        'Mol_Size_2': float, the second size of the molecule in Angstrom.
        'Mol_Size_L': float, the long size of the molecule in Angstrom.
        'Length_Ratio': float, the ratio of the long size to the short size.
        'Len_Div_Diameter': float, the ratio of the long size to the diameter of the molecule.
        'MPP': float, the maximum polarizability of the molecule in Angstrom^3.
        'SDP': float, the static dipole polarizability of the molecule in Angstrom^3.
        'Dipole_Moment': float, the dipole moment of the molecule in Debye.
        'Volume': float, the volume of the molecule in Angstrom^3.
        'Density': float, the density of the molecule in g/cm^3.
        'ESPmin': float, the minimum electrostatic potential of the molecule in eV.
        'ESPmax': float, the maximum electrostatic potential of the molecule in eV.
        'Overall_Surface_Area': float, the overall surface area of the molecule in Angstrom^2.
        'Pos_Surface_Area': float, the positive surface area of the molecule in Angstrom^2.
        'Neg_Surface_Area': float, the negative surface area of the molecule in Angstrom^2.
        'Overall_Average': float, the overall average electrostatic potential of the molecule in eV.
        'Pos_Average': float, the positive average electrostatic potential of the molecule in eV.
        'Neg_Average': float, the negative average electrostatic potential of the molecule in eV.
        'Overall_Variance': float, the overall variance of the electrostatic potential of the molecule in eV^2.
        'Nu': float, the nucleophilicity of the molecule.
        'Pi': float, the electrophilicity of the molecule.
        'MPI': float, the maximum polarizability index of the molecule.
        'Nonpolar_Area': float, the nonpolar area of the molecule in Angstrom^2.
        'Polar_Area': float, the polar area of the molecule in Angstrom^2.
        'ALIEmin': float, the minimum absolute local ionization energy of the molecule in eV.
        'ALIEmax': float, the maximum absolute local ionization energy of the molecule in eV.
        'ALIE_Ave': float, the average absolute local ionization energy of the molecule in eV.
        'ALIE_Var': float, the variance of the absolute local ionization energy of the molecule in eV^2.
        'LEAmin': float, the minimum local electrophilicity of the molecule in eV.
        'LEAmax': float, the maximum local electrophilicity of the molecule in eV.
        'LEA_Ave': float, the average local electrophilicity of the molecule in eV.
    """
    import pickle

    dataset_path = os.environ.get("ILS4GAS_QM_FEATURE_PATH", None)
    if dataset_path is None:
        raise ValueError(
            "Environment variable ILS4GAS_QM_FEATURE_PATH is not set. "
            "Please set it to the path of the QM feature dataset."
        )

    dataset = pickle.load(open(dataset_path, "rb"))
    if mol_name not in dataset:
        raise ValueError(f"Data for {mol_name} not found in the QM feature dataset.")

    return dataset[mol_name]


@mcp.tool()
def generate_novel_ions(
    base_molecule_smile: str,
    molecule_type: Literal["cation", "anion"],
    num_molecules: int = 10,
    temperature: float = 0.7,
) -> str:
    """Generate novel ionic liquid molecules using a machine learning pipeline.

    Runs an external bash script that uses FlashFormer to generate new ion
    candidates, computes SA Score and similarity, and outputs ranked results.

    Args:
        base_molecule_smile: SMILES of the starting fragment (e.g. "C", "CCN").
        molecule_type: "cation" or "anion".
        num_molecules: Number of molecules to generate (default 10).
        temperature: Sampling diversity, 0.1–1.0 (default 0.7). Lower → more similar.

    Returns:
        Path to the sorted CSV file (columns: smiles, sa_score, similarity).

    Requires:
        ILS4GAS_MOLECULE_GEN_SCRIPT set to the pipeline script path.
        The script must accept: --ion, --start, --num, --temp, --output.
    """
    script_path = os.environ.get("ILS4GAS_MOLECULE_GEN_SCRIPT", "")
    if not script_path:
        raise ValueError(
            "ILS4GAS_MOLECULE_GEN_SCRIPT is not set. "
            "Set it to the path of the generation pipeline script."
        )

    script = Path(script_path)
    if not script.exists():
        raise FileNotFoundError(f"Script not found: {script}")

    ion_map = {"cation": "ca", "anion": "an"}
    if molecule_type not in ion_map:
        raise ValueError(
            f"molecule_type must be 'cation' or 'anion', got '{molecule_type}'"
        )

    if not (0.1 <= temperature <= 1.0):
        raise ValueError(
            f"temperature must be 0.1–1.0, got {temperature}"
        )

    work_path = generate_work_path(create=True)
    output_prefix = str(
        Path(work_path) / f"{molecule_type}_{num_molecules}_{temperature:.1f}"
    )

    cmd = (
        f"bash {script} "
        f"--ion {ion_map[molecule_type]} "
        f"--start {base_molecule_smile} "
        f"--num {num_molecules} "
        f"--temp {temperature} "
        f"--output {output_prefix}"
    )

    return_code, out, err = run_command(cmd, cwd=str(script.parent))

    if return_code != 0:
        raise RuntimeError(f"Pipeline failed (exit {return_code}):\n{err}")

    sorted_csv = f"{output_prefix}_results_sorted.csv"
    if not Path(sorted_csv).exists():
        raise FileNotFoundError(f"Output not found: {sorted_csv}")

    return sorted_csv


@mcp.tool()
def predict_binding_energy(
    complex_xyz_path: str,
    part1_xyz_path: str,
    part2_xyz_path: str,
    part1_eps_path: str,
    part2_eps_path: str,
) -> Dict[str, float]:
    """Predict binding energy between two molecular structures using a pre-trained model.

    Runs the Eb_predict.py script in single mode with the provided structure and EPS files.

    Args:
        complex_xyz_path: Path to complex structure file (xyz format).
        part1_xyz_path: Path to monomer 1 structure file (xyz format).
        part2_xyz_path: Path to monomer 2 structure file (xyz format).
        part1_eps_path: Path to monomer 1 electrostatic potential file (eps format).
        part2_eps_path: Path to monomer 2 electrostatic potential file (eps format).

    Returns:
        Dict with "binding_energy" key containing the predicted value.

    Requires:
        ILS4GAS_EB_PREDICT_SCRIPT set to the Eb_predict.py script path.
    """
    script_path = os.environ.get("ILS4GAS_EB_PREDICT_SCRIPT", "")
    if not script_path:
        raise ValueError(
            "ILS4GAS_EB_PREDICT_SCRIPT is not set. "
            "Set it to the path of the Eb_predict.py script."
        )

    script = Path(script_path)
    if not script.exists():
        raise FileNotFoundError(f"Script not found: {script}")

    # Validate all input files exist
    input_files = [
        ("complex_xyz", complex_xyz_path),
        ("part1_xyz", part1_xyz_path),
        ("part2_xyz", part2_xyz_path),
        ("part1_eps", part1_eps_path),
        ("part2_eps", part2_eps_path),
    ]
    for name, path in input_files:
        if not Path(path).exists():
            raise FileNotFoundError(f"{name} file not found: {path}")

    work_path = generate_work_path(create=True)
    output_file = Path(work_path) / "binding_energy_result.txt"

    cmd = (
        f"python {script} single "
        f"{complex_xyz_path} "
        f"{part1_xyz_path} "
        f"{part2_xyz_path} "
        f"{part1_eps_path} "
        f"{part2_eps_path} "
        f"-o {output_file}"
    )

    return_code, out, err = run_command(cmd, cwd=str(work_path))

    if return_code != 0:
        raise RuntimeError(f"Binding energy prediction failed (exit {return_code}):\n{err}")

    # Read result from output file
    if not output_file.exists():
        raise FileNotFoundError(f"Result output file not found: {output_file}")

    with open(output_file, "r") as f:
        result_line = f.read().strip()
        try:
            binding_energy = float(result_line)
        except ValueError:
            raise RuntimeError(f"Invalid result format from prediction: {result_line}")

    return {"binding_energy": binding_energy}