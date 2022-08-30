import numpy as np

from . import RootSystem

class TheAlgebra:
    '''
    TheAlgebra Class
    -----------------

    Takes a RootSystem class object. This class mainly is used to simplify
    a general computator between 2 roots of a given root system.
    '''

    def __init__(self, rootsys=None):
        if rootsys is None or not type(rootsys) == RootSystem.RootSystem:
            print('\nWARNING NO ROOT SYSTEM PROVIDED. Please provide one either\
by initializing it or at a later stage using set_root_sys method\n\n')
        self.rootsys = rootsys

    def gen_structure_constants(self):
        '''
        Generates the structure constants for the algebra.
        '''

        # Checks whether the roots are generated
        if self.rootsys.roots is None:
            self.rootsys.gen_all_roots()
        
        # Checks whether the basics commutator is generated.
        if self.rootsys.basics_commutator is None:
            self.rootsys.gen_basics_commutator()

        # Create an empty matrix for the structure constants.
        self.structure_constants = np.zeros((2*self.rootsys.n_roots, 2*self.rootsys.n_roots))