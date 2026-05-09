import subprocess, select
from pathlib import Path
from typing import List, Tuple, Union, Optional
import os, time, json, traceback, uuid


def run_command(
    cmd: Union[str, List[str]],
    cwd: str = ".",
    shell: bool = True,
) -> Tuple[int, str, str]:
    process = subprocess.Popen(
        cmd,
        cwd=os.path.abspath(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=shell,
        executable="/bin/bash" if shell else None,
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
    if not paths:
        return []
    if len(paths) == 1:
        return [os.path.basename(str(paths[0]))]
    abs_paths = [Path(p).absolute() for p in paths]
    common_prefix = os.path.commonpath(abs_paths)
    relative_paths = [str(p.relative_to(common_prefix)) for p in abs_paths]
    return relative_paths


def get_physical_cores():
    with open('/proc/cpuinfo', 'r') as f:
        cpuinfo = f.read()
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
    if physical_ids and cores_per_socket:
        return sum(cores_per_socket.values())
    else:
        output = subprocess.check_output('lscpu', shell=True).decode()
        for line in output.split('\n'):
            if line.startswith('Core(s) per socket:'):
                cores_per_socket = int(line.split(':')[1].strip())
            elif line.startswith('Socket(s):'):
                sockets = int(line.split(':')[1].strip())
        return cores_per_socket * sockets


def generate_work_path(create: bool = True) -> str:
    calling_function = traceback.extract_stack(limit=2)[-2].name
    current_time = time.strftime("%Y%m%d%H%M%S")
    random_string = str(uuid.uuid4())[:8]
    work_path = f"{current_time}.{calling_function}.{random_string}"
    if create:
        os.makedirs(work_path, exist_ok=True)
    return work_path


def xyz_to_smiles(xyz_path):
    from openbabel import pybel
    mol = next(pybel.readfile("xyz", xyz_path))
    smiles = mol.write("smi").strip().split()[0]
    return to_canonical_smiles(smiles)


def to_canonical_smiles(smiles):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return None
