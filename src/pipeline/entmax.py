"""
Generate structures through entropy maximization using MatterSim and FIRE.
"""

from pathlib import Path
import numpy as np
import torch
from ase.io import read
from ase.optimize import FIRE
from mattersim.forcefield import MatterSimCalculator
from .custom_calc_entropymaxim import EntMaxCalc, calculate_entropy

def generate_structures(root_dir, n_structs=5000, device=None):
    
    root = Path(root_dir)
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")

    for folder in sorted(root.glob("mp-*")):

        poscar_super = folder / "POSCAR_super"
        poscar_default = folder / "POSCAR"
        if poscar_super.exists():
            poscar_path = poscar_super
        elif poscar_default.exists():
            poscar_path = poscar_default
        else:
            print(f"❌ {folder.name} no POSCAR or POSCAR_super, skip")
            continue

        print(f"▶ processing {folder.name} ({poscar_path.name})")

        atom = read(poscar_path, format="vasp")
        atom.rattle(stdev=0.01)
        atom.calc = EntMaxCalc(MatterSimCalculator(device=device), 300, 5.0)

        def print_entropy(a=atom):
            entropy = calculate_entropy(a, cutoff=4.0)
            force_max = np.amax(a.get_forces())
            print(f"[{folder.name}] Entropy = {entropy:.5f}, Max Force = {force_max:.5f}")

        traj_path = folder / "opt.traj"
        dyn = FIRE(atom, trajectory=str(traj_path))
        dyn.attach(print_entropy, interval=1)
        dyn.run(fmax=0.01, steps=n_structs)  # step limit

    print("✅ finish generating entmax structures")

