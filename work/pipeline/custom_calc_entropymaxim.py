"""This module defines an ASE interface to SpkEnergyOnly calculator
   that can compute forces and stresses through finite differences
"""

import re
from pathlib import Path
from subprocess import check_output
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import ase
from ase import units
from ase.constraints import FixAtoms
from ase.calculators.calculator import Calculator, all_changes
from ase.vibrations import Vibrations
import logging
from copy import deepcopy
from ase import Atoms
from ase.optimize.lbfgs import LBFGS
from ase.constraints import ExpCellFilter
from ase.io import read, write
import numpy as np
# import potential
# import fplibmod
# import rcovdata
import libfp

def from_ase(atom):
    lat = atom.get_cell()[:]
    syms = atom.get_chemical_symbols()
    types = []
    for i, s in enumerate(set(syms)):
        types += [i+1]*syms.count(s)
    types = np.array(types, int)
    rxyz = atom.positions
    znucl = np.array(list(set(atom.get_atomic_numbers())), int)
    return lat, rxyz, types, znucl

def create_fingerprints(atom, natx, cutoff):
    lat, rxyz, types, znucl = from_ase(atom)
    fp, dfp = libfp.get_dfp((lat, rxyz, types, znucl), cutoff=cutoff, natx=natx)
    return fp, dfp

def find_min_fp_dist(fp, ind):
    N = len(fp)
    minind = 0
    mindist = 100000000000
    for i in range(N):
        if i != ind:
            diff = fp[ind] - fp[i]
            dist = np.linalg.norm(diff)
            if dist < mindist:
                mindist = dist
                minind = i
    if mindist < 1e-8:
        mindist = 1e-8
    return minind, mindist

def find_min_fp_inds(fp):
    N = len(fp)
    mininds = []
    for i in range(N):
        mindist = 1e12
        minind = 0
        for j in range(N):
            if i == j:
                continue
            dist = np.linalg.norm(fp[i] - fp[j])
            if dist < mindist:
                minind = j
                mindist = dist
        mininds.append(minind)
    return mininds

def calculate_entropy(atom, cutoff):
    fingerprints, dfps = create_fingerprints(atom, natx=4*len(atom), cutoff=cutoff)
    fingerprint = fingerprints
    dfp = dfps
    N = len(fingerprint)
    entropy = 0
    for i in range(N):
        lmin, deltaqmin = find_min_fp_dist(fingerprint, i)
        entropy += np.log(N*deltaqmin)
    entropy = entropy/N
    return entropy

def calculate_entropy_derivative(atom, cutoff):
    fingerprints, dfps = create_fingerprints(atom, natx=4*len(atom), cutoff=cutoff)
    fingerprint = fingerprints
    dfp = dfps
    N = len(fingerprint)
    entropy_deriv = np.zeros((N,3))

    lvals = find_min_fp_inds(fingerprint)
    for i in range(N):
        for j in range(N):
            fpdiff = fingerprint[j] - fingerprint[lvals[j]]
            deltaq = np.linalg.norm(fpdiff)
            if deltaq == 0.0:
              continue
            for k in range(3):
                derivdiff = dfp[j, i, k] - dfp[lvals[j], i, k]
                numerator = -(1/N)*np.dot(derivdiff, fpdiff)
                denominator = (deltaq**2)
                if denominator < 1e-6:
                    factdiff = -6 - int(np.log10(denominator))
                    denominator = denominator*(10**factdiff)
                entropy_deriv[i, k] += numerator/denominator
    return entropy_deriv

class EntMaxCalc(Calculator):
    energy = "energy"
    forces = "forces"
    stress = "stress"
    implemented_properties=[energy, forces, stress]

    def __init__(
        self,
        calculator,
        kfactor,
        cutoff
        ):
        Calculator.__init__(self)
        self.kfactor = kfactor
        self.calculator = calculator
        self.fp_cutoff = cutoff

    def calculate(
        self,
        atoms: ase.Atoms = None,
        properties: list[str] = ["energy"],
        system_changes: list[str] = all_changes,
    ):
        if self.calculation_required(atoms, properties):
            newatom = atoms.copy()
            newatom.calc = self.calculator
            results = {}
            for p in properties:
                if p == "energy":
                    results[p] = newatom.get_potential_energy() - self.kfactor*calculate_entropy(atoms, self.fp_cutoff)
                elif p=="forces":
                    results[p] = (newatom.get_forces() - self.kfactor*calculate_entropy_derivative(atoms, self.fp_cutoff))/max(self.kfactor, 0)
                elif p=="stress":
                    results[p] = newatom.get_stress()
            self.results = results


#
