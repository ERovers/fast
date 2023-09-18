# Author: Maxwell I. Zimmerman <mizimmer@wustl.edu>
# Contributors:
# Copywright (C) 2017, Washington University in St. Louis
# All rights reserved.
# Unauthorized copying of this file, via any medium, is strictly prohibited
# Proprietary and confidential
# have reference structures (10 RMSD) maximize minimum rmsd to references

#######################################################################
# imports
#######################################################################


import gc
import glob
import mdtraj as md
import numpy as np
import os
import subprocess
from .base_analysis import base_analysis
from .. import tools


#######################################################################
# code
#######################################################################


class RMSD_PROTACWrap(base_analysis):
    """Analysis wrapper for calculating state rmsds

    Parameters
    ----------
    base_struct : str or md.Trajectory,
        The base structure to compare for native contacts. This
        topology must match the structures to analyse. Can be provided
        as a pdb location or an md.Trajectory object.
    atom_indices : str or array,
        The atom indices to use for computing native contacts. Can be
        provided as a data file to load or an array.

    Attributes
    ----------
    output_name : str,
        The file containing rankings.
    """
    def __init__(
            self, base_struct, atom_indices=None):
        # determine base_struct
        self.base_struct = base_struct
        if type(base_struct) is md.Trajectory:
            self.base_struct_md = self.base_struct
        else:
            self.base_struct_md = md.load(base_struct)
        # determine atom indices
        self.atom_indices = atom_indices
        if type(atom_indices) is str:
            self.atom_indices_vals = np.loadtxt(atom_indices, dtype=int)
        else:
            self.atom_indices_vals = self.atom_indices

    @property
    def class_name(self):
        return "RMSD_PROTACWrap"

    @property
    def config(self):
        return {
            'base_struct': self.base_struct,
            'atom_indices': self.atom_indices,
        }

    @property
    def analysis_folder(self):
        return None

    @property
    def base_output_name(self):
        return "max_rmsd_per_state"

    def run(self):
        # determine if file already exists
        if os.path.exists(self.output_name):
            pass
        else:
            if not os.path.exists("reference.xtc"):
               subprocess.call("gmx trjconv -f "+self.base_struct+" -o reference.xtc")
            # load centers
            centers = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md,
                atom_indices=self.atom_indices_vals)
            # load references
            ref_struct = md.load("reference.xtc", top=self.base_struct_md,
                atom_indices=self.atom_indices_vals)
            # get subset if necessary
            if self.atom_indices_vals is None:
                struct_sub = self.base_struct_md
            else:
                struct_sub = self.base_struct_md.atom_slice(self.atom_indices_vals)
            # calculate and save rmsds
            rmsds = []
            for frame in centers: 
                tmp = md.rmsd(centers, frame)
                tmp.sort()
                rmsds.append(rmsd[1])
            np.save(self.output_name, rmsds)
            # save the complexes that differ more than 10 A from any of the references (to be references for future runs)
            idx = np.argwhere(rmsds>10)
            tmp2 = self.base_struct_md
            if len(idx>0):
                for i in idx:
                    tmp2 = md.join([tmp2,centers[i]])
                with md.formats.XTCTrajectoryFile('reference.xtc', 'w') as f:
                    f.write(tmp2)
