#  Data-Driven Structure Generation Pipeline

This Python pipeline automates downloading crystal structures from the Materials Project, expanding small structures, generating configurations via entropy maximization, computing fingerprints, clustering structures, and merging filtered results.

## ✨ Features

- Download structures from Materials Project based on element types and number of elements.

- Expand small structures to ensure a minimum number of atoms.

- Entropy maximization: generate perturbed structures for sampling configurational space.

- Fingerprint calculation: generate structural fingerprints for clustering.

- Clustering & filtering: select representative structures using fingerprint similarity.

- Merge filtered structures into a single file: `all_clustered_atoms.extxyz`.

## 🛠️ Installation

```
git clone https://github.com/SKYLER022/Entmax_data_generation.git
```
## Dependence

```
python=3.11
pip install "mp-api>=0.45" "emmet-core>=0.45" pymatgen tqdm torch ase mattersim libfp scikit-learn-extra "numpy<2"
export MP_API_KEY="YOUR_MP_KEY"
```
## Usage

`python main.py`
For test use: `python test.py`

You will see
```
Enter elements to search (comma-separated, e.g., C,H,O):
Enter number of elements to restrict (or press Enter for no restriction):
```
For example, we can type:  Fe,F (Enter) 2(Enter)







