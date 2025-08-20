"""
Compute fingerprints for all structures in opt.traj files inside mp-* folders.
Generate one pickle file per folder.
"""

from pathlib import Path
import pickle
import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor as AAA
from pymatgen.analysis.local_env import CrystalNN
import libfp

def print_nn_dists(struct, min_dist=3.0):
    cnn = CrystalNN()
    natoms = len(struct)
    dists = []
    for i in range(natoms):
        for j in range(i):
            dist = struct.get_distance(i, j)
            dists.append(dist)
    dists = np.array(dists)/0.529177
    if np.all(dists > min_dist):
        return dists
    else:
        return None

def create_fingerprints(atom, positions, cell, natx):
    all_symbols = atom.get_chemical_symbols()
    unique_symbols = list(dict.fromkeys(all_symbols))
    type_nums = [all_symbols.count(x) for x in unique_symbols]

    types = []
    for i, num in enumerate(type_nums):
        types += [i+1]*num
    types = np.array(types, int)
    znucl = np.array(list(set(atom.get_atomic_numbers())), int)

    cutoff = 4.0
    fp1 = libfp.get_lfp((cell, positions, types, znucl), cutoff, log=False, orbital='s', natx=natx)
    return fp1

def compute_fingerprints(root_dir, natx=100):
    root = Path(root_dir)

    for folder in sorted(root.glob("mp-*")):
        traj_path = folder / "opt.traj"
        if not traj_path.exists():
            print(f"❌ {folder.name} no opt.traj, skip")
            continue

        final_atoms = []
        final_fps = []
        failcount = 0
        passcount = 0

        atoms = read(traj_path, index=":")
        for i, atom in enumerate(atoms):
            try:
                testatom = atom.copy()
                dists = print_nn_dists(AAA.get_structure(testatom))
                testatom.wrap()
                fingerprint = create_fingerprints(testatom, testatom.positions, testatom.get_cell()[:], natx=natx)
                final_fps.append(fingerprint)
                final_atoms.append(testatom)
                passcount += 1
            except:
                failcount += 1

        # save fingerprints per folder
        with open(folder / "atoms_fingerprints.pkl", "wb") as f:
            pickle.dump([final_atoms, final_fps], f)

        print(f"[{folder.name}] Passed: {passcount}, Failed: {failcount}, "
              f"Pass%: {passcount/(passcount+failcount):.2f}")

    print("✅ Finished computing fingerprints for all folders")

