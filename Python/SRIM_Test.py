#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 16:13:04 2022

@author: salomonguinchard
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from srim import TRIM, Ion, Layer, Target
from srim.output import Results

# Construct a 3MeV Nickel ion
ion = Ion('Ni', energy=3.0e6)

# Construct a layer of nick 20um thick with a displacement energy of 30 eV
layer = Layer({
        'Ni': {
            'stoich': 1.0,
            'E_d': 30.0,
            'lattice': 0.0,
            'surface': 3.0
        }}, density=8.9, width=20000.0)

# Construct a target of a single layer of Nickel
target = Target([layer])

# Initialize a TRIM calculation with given target and ion for 25 ions, quick calculation
trim = TRIM(target, ion, number_ions=25, calculation=1)

# Specify the directory of SRIM.exe
# For windows users the path will include C://...
srim_executable_directory = '/tmp/srim'

# takes about 10 seconds on my laptop
results = trim.run(srim_executable_directory)
# If all went successfull you should have seen a TRIM window popup and run 25 ions!

output_directory = '/tmp/srim_outputs'
os.makedirs(output_directory, exist_ok=True)
TRIM.copy_output_files('/tmp/srim', output_directory)

srim_executable_directory = '/tmp/srim'
results = Results(srim_executable_directory)

def plot_damage_energy(folder, ax):
     results = Results(folder)
     phon = results.phonons
     dx = max(phon.depth) / 100.0 # to units of Angstroms
     energy_damage = (phon.ions + phon.recoils) * dx
     ax.plot(phon.depth, energy_damage / phon.num_ions, label='{}'.format(folder))
     return sum(energy_damage)

def plot_ionization(folder, ax):
     results = Results(folder)
     ioniz = results.ioniz
     dx = max(ioniz.depth) / 100.0 # to units of Angstroms
     ax.plot(ioniz.depth, ioniz.ions, label='Ionization from Ions')
     ax.plot(ioniz.depth, ioniz.recoils, label='Ionization from Recoils')

def plot_vacancies(folder, ax):
     results = Results(folder)
     vac = results.vacancy
     vacancy_depth = vac.knock_ons + np.sum(vac.vacancies, axis=1)
     ax.plot(vac.depth, vacancy_depth, label="Total vacancies at depth")
     return sum(vacancy_depth)

folders = ['test_files/2', 'test_files/4']
image_directory = 'examples/images'
os.makedirs(image_directory, exist_ok=True)


fig, axes = plt.subplots(1, len(folders), sharex=True, sharey=True)

for ax, folder in zip(np.ravel(axes), folders):
    energy_damage = plot_damage_energy(folder, ax)
    print("Damage energy: {} eV".format(energy_damage))
    ax.set_xlabel('Depth [Angstroms]')
    ax.set_ylabel('eV')
    ax.legend()

fig.suptitle('Damage Energy vs. Depth', fontsize=15)
fig.set_size_inches((20, 6))
fig.savefig(os.path.join(image_directory, 'damagevsdepth.png'), transparent=True)
