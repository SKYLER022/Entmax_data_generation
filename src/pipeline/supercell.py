"""
Expand structures with small number of atoms into supercells.
"""

from pathlib import Path
from itertools import product
from pymatgen.io.vasp import Poscar

def minimal_replication(n_atoms, target=33, max_factor=4):
    
    best = None
    for a, b, c in product(range(1, max_factor + 1), repeat=3):
        total = n_atoms * a * b * c
        if total >= target:
            vol = a * b * c
            shape = max(a, b, c) - min(a, b, c)
            key = (vol, shape, a, b, c)
            if best is None or key < best[0]:
                best = (key, (a, b, c))
    return best[1] if best else (1, 1, 1)


def expand_small_structures(root_dir, min_atoms=32):
    
    root = Path(root_dir)
    for poscar_path in root.rglob("POSCAR"):
        struct = Poscar.from_file(poscar_path).structure
        natoms = len(struct)

        if natoms >= min_atoms:
            print(f"[skip] {poscar_path.parent.name}: {natoms} atoms")
            continue

        a, b, c = minimal_replication(natoms, target=min_atoms)
        super_struct = struct * (a, b, c)
        total_atoms = len(super_struct)

        out_path = poscar_path.with_name("POSCAR_super")
        Poscar(super_struct).write_file(out_path)

        print(f"[new ] {poscar_path.parent.name}: {natoms} → {total_atoms} atoms "
              f"(replicated {a}×{b}×{c})")

    print("Done – supercells written where needed.")

