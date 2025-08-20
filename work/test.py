from pipeline.download import download_structures
from pipeline.supercell import expand_small_structures
from pipeline.entmax import generate_structures
from pipeline.fingerprint import compute_fingerprints
from pipeline.cluster import select_representatives
from pipeline.merge import merge_clustered_atoms

def test():


    download_structures(save_dir="data")

    expand_small_structures("data", min_atoms=32)

    generate_structures("data", n_structs=100)

    compute_fingerprints("data")

    select_representatives("data", n_select=10)
    
    merge_clustered_atoms("data", output_file="all_clustered_atoms.extxyz")

if __name__ == "__main__":
    test()

