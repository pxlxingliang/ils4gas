import subprocess
import select
from pathlib import Path
from typing import List, Tuple, Union, Optional
import os
import time
import json
import traceback
import uuid
import glob



def run_command(
        cmd,
        shell=True
):
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=shell,
        executable='/bin/bash'
    )
    out = ""
    err = ""
    while True:
        readable, _, _ = select.select(
            [process.stdout, process.stderr], [], [])

        for fd in readable:
            if fd == process.stdout:
                line = process.stdout.readline()
                print(line.decode()[:-1])
                out += line.decode()
            elif fd == process.stderr:
                line = process.stderr.readline()
                print("STDERR:", line.decode()[:-1])
                err += line.decode()

        return_code = process.poll()
        if return_code is not None:
            break
    return return_code, out, err

def remove_comm_prefix(paths: Union[List[Path], List[str]]) -> List[str]:
    """
    Remove the common prefix from a list of paths.
    This is useful for displaying relative paths in logs.
    """
    if not paths:
        return []

    if len(paths) == 1:
        return [os.path.basename(str(paths[0]))]
    
    # Convert all paths to absolute paths
    abs_paths = [Path(p).absolute() for p in paths]
    
    # Find the common prefix
    common_prefix = os.path.commonpath(abs_paths)
    
    # Remove the common prefix from each path
    relative_paths = [str(p.relative_to(common_prefix)) for p in abs_paths]
    
    return relative_paths

def get_physical_cores():
    """
    """
    # 对于Linux系统，解析/proc/cpuinfo文件
    with open('/proc/cpuinfo', 'r') as f:
        cpuinfo = f.read()
    
    # 统计物理ID的数量和每个物理ID下的核心数
    physical_ids = set()
    cores_per_socket = {}
    
    for line in cpuinfo.split('\n'):
        if line.startswith('physical id'):
            phys_id = line.split(':')[1].strip()
            physical_ids.add(phys_id)
        elif line.startswith('cpu cores'):
            cores = int(line.split(':')[1].strip())
            if phys_id in cores_per_socket:
                cores_per_socket[phys_id] = max(cores_per_socket[phys_id], cores)
            else:
                cores_per_socket[phys_id] = cores
    
    # 计算总物理核心数
    if physical_ids and cores_per_socket:
        return sum(cores_per_socket.values())
    else:
        # 备选方法：使用lscpu命令
        output = subprocess.check_output('lscpu', shell=True).decode()
        for line in output.split('\n'):
            if line.startswith('Core(s) per socket:'):
                cores_per_socket = int(line.split(':')[1].strip())
            elif line.startswith('Socket(s):'):
                sockets = int(line.split(':')[1].strip())
        return cores_per_socket * sockets
   
def generate_work_path(create: bool = True) -> str:
    """
    Generate a unique working directory path based on call function and current time.
    
    directory = calling function name + current time + random string.
    
    Returns:
        str: The path to the working directory.
    """
    calling_function = traceback.extract_stack(limit=2)[-2].name
    current_time = time.strftime("%Y%m%d%H%M%S")
    random_string = str(uuid.uuid4())[:8]
    work_path = f"{current_time}.{calling_function}.{random_string}"
    if create:
        os.makedirs(work_path, exist_ok=True)
    
    return work_path

def xyz_to_smiles(xyz_path):
    """
    Convert an XYZ file to SMILES format using Open Babel.
    
    Args:
        xyz_path (str): The path to the XYZ file.
        
    Returns:
        str: The SMILES representation of the molecule.
    
    Example:
        >>> smiles = xyz_to_smiles("molecule.xyz") # molecule is benzene
        >>> print(smiles)
        C1=CC=CC=C1
    """
    from openbabel import pybel
    mol = next(pybel.readfile("xyz", xyz_path))
    smiles = mol.write("smi").strip().split()[0]
    
    return to_canonical_smiles(smiles)

def to_canonical_smiles(smiles):
    """ Convert a SMILES string to its canonical form using RDKit.
    Args:
        smiles (str): The SMILES string to convert.
    Returns:
        str: The canonical SMILES string, or None if conversion fails.
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return None
