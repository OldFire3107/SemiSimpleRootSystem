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
        Generates the structure constants for the algebra, expected for the ones
        that involve the Cartan generators or have them in the result as thier 
        calculation is trivial enough.

        Antisymmetric of the structure constants can be used to get all the relationship
        between the roots. Also a positive root subtracted with a negetive root will 
        never be a root. All these are exploited here to get the structure constants.
        '''

        # Checks whether the roots are generated
        if self.rootsys.roots is None:
            self.rootsys.gen_all_roots()
        
        # Checks whether the basics commutator is generated.
        if self.rootsys.basics_commutator is None:
            self.rootsys.gen_basics_commutator()

        # Create an empty matrix for the structure constants. Not including Cartan generators as it needs the roots in vector form.
        # Also we know that Cartan Generator with the generator would just give the component of the root in the direction of the generator.
        # It is also know that positive and the negetive roots will just give the root dotted with the Cartan generators.
        self.structure_constants = np.zeros((2*self.rootsys.n_roots, 2*self.rootsys.n_roots))

        