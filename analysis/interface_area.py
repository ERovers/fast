# Author: Evianne M. Rovers <evianne.rovers@mail.utoronto.ca>
# Contributors:
# Copywright (C) 2023, University of Toronto
# All rights reserved.
# Unauthorized copying of this file, via any medium, is strictly prohibited
# Proprietary and confidential


#######################################################################
# imports
#######################################################################


import gc
import glob
import mdtraj as md
import numpy as np
import os
from .base_analysis import base_analysis
from .. import tools


#######################################################################
# code
#######################################################################


class InterfaceWrap(base_analysis):
    """Analysis wrapper for calculating state interface area

    Parameters
    ----------
    base_struct : str or md.Trajectory,
        The base structure to compare for native contacts. This
        topology must match the structures to analyse. Can be provided
        as a pdb location or an md.Trajectory object.
    atom_indices1 : str or array,
        The atom indices for protein 1 (part1). Can be
        provided as a data file to load or an array.
    atom_indices2 : str or array,
        The atom indices for protein 2 (part2). Can be
        provided as a data file to load or an array.
    atom_indices3 : str or array,
        The atom indices for PROTAC (part3). Can be
        provided as a data file to load or an array.

    Attributes
    ----------
    output_name : str,
        The file containing rankings.
    """
    def __init__(
            self, base_struct, atom_indices1=':', atom_indices2=':', atom_indices3=':'):
        # determine base_struct
        self.base_struct = base_struct
        if type(base_struct) is md.Trajectory:
            self.base_struct_md = self.base_struct
        else:
            self.base_struct_md = md.load(base_struct)
        # determine atom indices protein 1
        self.atom_indices1 = atom_indices1
        if type(atom_indices1) is str:
            self.atom_indices_vals1 = np.loadtxt(atom_indices1, dtype=int)
        else:
            self.atom_indices_vals1 = self.atom_indices1
        # determine atom indices protein 2
        self.atom_indices2 = atom_indices2
        if type(atom_indices2) is str:
            self.atom_indices_vals2 = np.loadtxt(atom_indices2, dtype=int)
        else:
            self.atom_indices_vals2 = self.atom_indices2
        self.atom_indices3 = atom_indices3
        if type(atom_indices2) is str:
            self.atom_indices_vals3 = np.loadtxt(atom_indices3, dtype=int)
        else:
            self.atom_indices_vals3 = self.atom_indices3
    @property
    def class_name(self):
        return "InterfaceWrap"

    @property
    def config(self):
        return {
            'base_struct': self.base_struct,
        }

    @property
    def analysis_folder(self):
        return None

    @property
    def base_output_name(self):
        return "int_area_per_state"

    def run(self):
        # determine if file already exists
        if os.path.exists(self.output_name):
            pass
        else:
            # load centers
            centers_total = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md)
            centers_part1 = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md,
                atom_indices=self.atom_indices_vals1)
            centers_part2 = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md,
                atom_indices=self.atom_indices_vals2)
            centers_part3 = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md,
                atom_indices=self.atom_indices_vals3)
            # calculate SASA of individual parts
            total_sasa = md.shrake_rupley(centers_total).sum(axis=1)
            part1 = md.shrake_rupley(centers_part1).sum(axis=1)
            part2 = md.shrake_rupley(centers_part2).sum(axis=1)
            part3 = md.shrake_rupley(centers_part3).sum(axis=1)
            # calculate interface area by formula SASA(individual parts) - SASA(total)
            interface = (np.array(part1) + np.array(part2) + np.array(part3)) - np.array(total_sasa)
            np.save(self.output_name, interface)
        
