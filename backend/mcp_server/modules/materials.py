from pathlib import Path
from typing import Literal, Optional, Dict, Any, List, Tuple, Union
import os

from backend.mcp_server import mcp
from backend.mcp_server.util.comm import generate_work_path, run_command
from backend.mcp_server.modules.enums import (
    GasSpecies,
    CationSpecies,
    AnionSpecies,
    DFTMolecule,
    DFTFeature,
)


@mcp.tool()
def search_properties_from_database(
    gas: str,
    cation: str,
    anion: str,
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
        gas (str): The gas species to search for (e.g., 'c2h4', 'co2', 'nh3', 'so2').
        cation (str): The cation species to search for (e.g., 'amim', 'c4mim', 'bmim', 'pf6').
        anion (str): The anion species to search for (e.g., 'bf4', 'pf6', 'alcl4', 'no3').
        properties (list): List of properties to retrieve. Valid options are:
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

    if gas not in GasSpecies._value2member_map_:
        valid = ", ".join(sorted(g.value for g in GasSpecies))
        raise ValueError(
            f"Invalid gas '{gas}'. Valid gas options are: {valid}"
        )

    if cation not in CationSpecies._value2member_map_:
        valid = ", ".join(sorted(c.value for c in CationSpecies))
        raise ValueError(
            f"Invalid cation '{cation}'. Valid cation options are: {valid}"
        )

    if anion not in AnionSpecies._value2member_map_:
        valid = ", ".join(sorted(a.value for a in AnionSpecies))
        raise ValueError(
            f"Invalid anion '{anion}'. Valid anion options are: {valid}"
        )

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
    mol_name: str,
    properties: Optional[List[str]] = None,
) -> Dict[str, float]:
    """
    Search DFT features of a molecule from the database.

    Args:
        mol_name (str): The name of the molecule to search for (e.g., 'c2h4', 'co2', 'bf4', 'amim').
        properties (list): List of DFT feature names to retrieve (defaults to all features). Options include:
            - "AtomNum": The number of atoms in the molecule.
            - "Weight": The molecular weight in amu.
            - "Mol_Radius": The radius of the molecule in Angstrom.
            - "Mol_Size_Short": The short size of the molecule in Angstrom.
            - "Mol_Size_2": The second size of the molecule in Angstrom.
            - "Mol_Size_L": The long size of the molecule in Angstrom.
            - "Length_Ratio": The ratio of the long size to the short size.
            - "Len_Div_Diameter": The ratio of the long size to the diameter of the molecule.
            - "MPP": The maximum polarizability of the molecule in Angstrom^3.
            - "SDP": The static dipole polarizability of the molecule in Angstrom^3.
            - "Dipole_Moment": The dipole moment of the molecule in Debye.
            - "Volume": The volume of the molecule in Angstrom^3.
            - "Density": The density of the molecule in g/cm^3.
            - "ESPmin": The minimum electrostatic potential of the molecule in eV.
            - "ESPmax": The maximum electrostatic potential of the molecule in eV.
            - "Overall_Surface_Area": The overall surface area of the molecule in Angstrom^2.
            - "Pos_Surface_Area": The positive surface area of the molecule in Angstrom^2.
            - "Neg_Surface_Area": The negative surface area of the molecule in Angstrom^2.
            - "Overall_Average": The overall average electrostatic potential of the molecule in eV.
            - "Pos_Average": The positive average electrostatic potential of the molecule in eV.
            - "Neg_Average": The negative average electrostatic potential of the molecule in eV.
            - "Overall_Variance": The overall variance of the electrostatic potential of the molecule in eV^2.
            - "Nu": The nucleophilicity of the molecule.
            - "Pi": The electrophilicity of the molecule.
            - "MPI": The maximum polarizability index of the molecule.
            - "Nonpolar_Area": The nonpolar area of the molecule in Angstrom^2.
            - "Polar_Area": The polar area of the molecule in Angstrom^2.
            - "ALIEmin": The minimum absolute local ionization energy of the molecule in eV.
            - "ALIEmax": The maximum absolute local ionization energy of the molecule in eV.
            - "ALIE_Ave": The average absolute local ionization energy of the molecule in eV.
            - "ALIE_Var": The variance of the absolute local ionization energy of the molecule in eV^2.
            - "LEAmin": The minimum local electrophilicity of the molecule in eV.
            - "LEAmax": The maximum local electrophilicity of the molecule in eV.
            - "LEA_Ave": The average local electrophilicity of the molecule in eV.

    Returns:
        dict: A dictionary containing the requested DFT features of the specified molecule.
    """
    import pickle

    if mol_name not in DFTMolecule._value2member_map_:
        valid = ", ".join(sorted(m.value for m in DFTMolecule))
        raise ValueError(
            f"Invalid molecule '{mol_name}'. Valid molecule options are: {valid}"
        )

    if properties is not None:
        for prop in properties:
            if prop not in DFTFeature._value2member_map_:
                valid = ", ".join(sorted(f.value for f in DFTFeature))
                raise ValueError(
                    f"Invalid feature '{prop}'. Valid feature options are: {valid}"
                )

    dataset_path = os.environ.get("ILS4GAS_QM_FEATURE_PATH", None)
    if dataset_path is None:
        raise ValueError(
            "Environment variable ILS4GAS_QM_FEATURE_PATH is not set. "
            "Please set it to the path of the QM feature dataset."
        )

    dataset = pickle.load(open(dataset_path, "rb"))
    if mol_name not in dataset:
        raise ValueError(f"Data for {mol_name} not found in the QM feature dataset.")

    all_features = dataset[mol_name]

    if properties is None:
        return all_features

    return {prop: all_features.get(prop) for prop in properties}


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