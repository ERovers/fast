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


class SASAWrap(base_analysis):
    """Analysis wrapper for calculating state SASA

    Parameters
    ----------
    base_struct : str or md.Trajectory,
        The base structure to compare for native contacts. This
        topology must match the structures to analyse. Can be provided
        as a pdb location or an md.Trajectory object.

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

    @property
    def class_name(self):
        return "SASAWrap"

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
        return "SASA_per_state"

    def run(self):
        # determine if file already exists
        if os.path.exists(self.output_name):
            pass
        else:
            # load centers
            centers = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md)
            # calculate and save SASA
            SASAs = md.shrake_rupley(centers)
            total_sasa = SASAs.sum(axis=1)
            print(total_sasa)
            np.save(self.output_name, total_sasa)
        
