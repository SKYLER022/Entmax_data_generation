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
conda create -n mp_env python=3.11
conda activate mp_env
cd Entmax_data_generation/work
python main.py
```





