"""
Download structures from Materials Project and save as POSCAR files.
"""

import os
from pathlib import Path
from tqdm import tqdm
from pymatgen.io.vasp import Poscar
from mp_api.client import MPRester

def download_structures(elements=None, num_elements=None, save_dir="data"):
    
    API_KEY = os.getenv("MP_API_KEY")
    if not API_KEY:
        raise ValueError("❌ MP_API_KEY is not set. Please export it before running.")

    OUTDIR = Path(save_dir)
    OUTDIR.mkdir(exist_ok=True)

    if elements is None:
        user_input = input("Enter elements to search (comma-separated, e.g., C,H,O): ")
        elements = [el.strip() for el in user_input.split(",") if el.strip()]

    if num_elements is None:
        num_elements_input = input("Enter number of elements to restrict (or press Enter for no restriction): ")
        num_elements = int(num_elements_input) if num_elements_input.strip() else None

    with MPRester(API_KEY) as mpr:
        atoms_docs = mpr.materials.summary.search(
            elements=elements, num_elements=num_elements, fields=["material_id"]
        )
        mpids = [d.material_id for d in atoms_docs]
        print(f"Found {len(mpids)} materials with elements {elements}.")

        # Download and save POSCAR
        for mpid in tqdm(mpids, desc="Downloading POSCARs"):
            struct = mpr.get_structure_by_material_id(mpid, conventional_unit_cell=False)
            d = OUTDIR / mpid
            d.mkdir(exist_ok=True)
            Poscar(struct).write_file(d / "POSCAR")

    print("✓ All POSCARs saved to", OUTDIR.resolve())

