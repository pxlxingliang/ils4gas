from pathlib import Path
from typing import Literal, Optional, Dict, Any, List

from backend.mcp_server import mcp
from backend.mcp_server.util.comm import generate_work_path


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
    molecule_name: str = "H2O",
    cell: Optional[List[List[float]]] = None,
    vacuum: Optional[float] = 5.0,
    output_file_format: Literal["cif", "poscar", "abacus"] = "abacus",
) -> Dict[str, Any]:
    """
    Generate molecule structure from ASE's collection of molecules or single atoms.

    Valid molecule_name can be either:
    - A chemical symbol of a single element (e.g., 'H', 'O', 'C'). Covers all elements
      in the periodic table except 'X'.
    - A molecule name in the g2 collection from ASE (e.g., 'H2O', 'CO2', 'NH3').

    Args:
        molecule_name (str): The name of the molecule or atom to generate.
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
def convert_structure_format(
    input_file: str,
    output_format: Literal["cif", "xyz", "poscar", "vasp", "pdb", "xsd", "abacus"] = "xyz",
    input_format: Optional[Literal["cif", "xyz", "poscar", "vasp", "pdb", "xsd"]] = None,
    output_filename: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Convert a structure file from one format to another using ASE.

    Args:
        input_file: Path to the input structure file.
        output_format: Output format. Options include 'cif', 'xyz', 'poscar', 'vasp', 'pdb', 'xsd', 'abacus'.
        input_format: Input format (optional). If not provided, ASE will try to auto-detect from file extension.
        output_filename: Output filename (optional). If not provided, a name will be generated automatically.

    Returns:
        A dictionary containing:
        - output_file: The absolute path to the converted structure file.
        - input_format: The format of the input file.
        - output_format: The format of the output file.
        - num_atoms: Number of atoms in the structure.
    """
    from ase.io import read, write

    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    try:
        atoms = read(str(input_path), format=input_format)
    except Exception as e:
        raise ValueError(f"Failed to read input file: {str(e)}") from e

    work_path = generate_work_path(create=True)

    if output_filename is None:
        base_name = input_path.stem
        ext = output_format
        if output_format in ["poscar", "vasp"]:
            ext = "vasp"
        elif output_format == "abacus":
            ext = "stru"
        output_path = Path(work_path) / f"{base_name}.{ext}"
    else:
        output_path = Path(work_path) / output_filename

    write_format = output_format
    if write_format == "poscar":
        write_format = "vasp"

    try:
        write(str(output_path), atoms, format=write_format)
    except Exception as e:
        raise ValueError(f"Failed to write output file: {str(e)}") from e

    detected_format = input_format or input_path.suffix.lstrip(".") or "unknown"

    return {
        "output_file": str(output_path.absolute()),
        "input_format": detected_format,
        "output_format": output_format,
        "num_atoms": len(atoms),
    }
