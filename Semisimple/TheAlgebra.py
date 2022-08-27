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

    def compute_computator(self, root1, root2):
        '''
        Computes the computator for a given root.
        '''

        # Checks whether it is a valid root
        if self.rootsys.roots is None:
            self.rootsys.gen_all_roots()

        if not root1 in self.rootsys.roots and not root2 in self.rootsys.roots:
            print('\nWARNING INVALID ROOT. Please provide a valid root\n\n')
            return None
        
        # Computes the computator