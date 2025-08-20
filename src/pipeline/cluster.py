import pickle
import numpy as np
from pathlib import Path
from sklearn_extra.cluster import KMedoids
from sklearn.decomposition import PCA
from ase.io import write
import matplotlib.pyplot as plt
import libfp
from tqdm import tqdm

def create_fp_dist_matrix(fps, types):
    N = len(fps)
    dist_matrix = np.zeros((N,N))
    shift = np.sqrt(fps.shape[2])
    for i in tqdm(range(N), desc="Building distance matrix"):
        for j in range(i+1):
            dist_matrix[i,j] = libfp.get_fp_dist(fps[i], fps[j], types)*shift
            dist_matrix[j,i] = dist_matrix[i,j]
    return dist_matrix

def create_and_save_fp_dist_matrix(atoms, fingerprints, filename="fp_dist_matrix.pkl"):
    all_symbols = atoms[0].get_chemical_symbols()
    unique_symbols = list(dict.fromkeys(all_symbols))
    type_nums = [all_symbols.count(x) for x in unique_symbols]

    types = []
    for i, num in enumerate(type_nums):
        types += [i+1]*num
    types = np.array(types, int)

    dist_matrix = create_fp_dist_matrix(fingerprints, types)
    with open(filename, "wb") as f:
        pickle.dump(dist_matrix, f)
    print(f"Distance matrix saved to {filename}")

def select_representatives(root_dir, n_select=10):
    root = Path(root_dir)

    for folder in sorted(root.glob("mp-*")):
        fp_file = folder / "atoms_fingerprints.pkl"
        if not fp_file.exists():
            print(f"❌ {folder.name} no atoms_fingerprints.pkl, skip")
            continue

        with open(fp_file, "rb") as f:
            atoms, fingerprints = pickle.load(f)

        if len(atoms) == 0:
            print(f"❌ {folder.name} has no atoms, skip")
            continue

        fingerprints = np.array(fingerprints)
        create_and_save_fp_dist_matrix(atoms, fingerprints, folder / "fp_dist_matrix.pkl")

        with open(folder / "fp_dist_matrix.pkl", "rb") as f:
            dist_matrix = pickle.load(f)

        # PCA for visualization
        flattened_fps = fingerprints.reshape(fingerprints.shape[0], -1)
        pca = PCA(n_components=3)
        red_fps = pca.fit_transform(flattened_fps)
        print(f"[{folder.name}] Explained variance:", pca.explained_variance_ratio_)

        clusterer = KMedoids(n_clusters=min(n_select, len(atoms)), metric="precomputed")
        clusterer.fit(dist_matrix)

        lclustdict = {}
        for i, lab in enumerate(clusterer.labels_):
            lclustdict.setdefault(lab, []).append(i)

        # plot 3D clusters
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for l in lclustdict:
            clust_nums = np.array(lclustdict[l])
            ax.scatter(red_fps[clust_nums, 0], red_fps[clust_nums, 1], red_fps[clust_nums, 2], c=np.random.rand(3))
        ax.set_title(f"Clusters in {folder.name}")
        ax.legend()
        plt.savefig(folder / "clusters.png", dpi=300)
        plt.close(fig)

        # select one representative per cluster
        selected_atoms = [atoms[np.random.choice(lclustdict[l])] for l in lclustdict]
        write(folder / "clustered_atoms.extxyz", selected_atoms, format="extxyz")
        print(f"[{folder.name}] Selected {len(selected_atoms)} representative structures saved.")

    print("✅ Finished clustering for all folders")

