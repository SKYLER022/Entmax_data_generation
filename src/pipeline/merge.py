from ase.io import read, write
from pathlib import Path

def merge_clustered_atoms(root_dir, output_file="all_clustered_atoms.extxyz"):
    root = Path(root_dir)
    all_atoms = []
    for folder in sorted(root.glob("mp-*")):
        clustered_file = folder / "clustered_atoms.extxyz"
        if clustered_file.exists():
            atoms = read(clustered_file, index=":")
            all_atoms.extend(atoms)
    if all_atoms:
        write(Path(root) / output_file, all_atoms, format="extxyz")
        print(f"✅ All clustered atoms merged into {output_file}")
    else:
        print("❌ No clustered_atoms.extxyz files found to merge.")
