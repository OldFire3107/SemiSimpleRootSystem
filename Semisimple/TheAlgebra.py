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
        self.raising_norm = self.rootsys.get_raising_norm() if self.rootsys.raising_norm is None else self.rootsys.raising_norm

    def gen_algebra_constants(self):
        '''
        Generates the algebra constants for the algebra, expected for the ones
        that involve the Cartan generators or have them in the result as thier 
        calculation is trivial enough.

        Antisymmetric of the algebra constants can be used to get all the relationship
        between the roots. Also a positive root subtracted with a negetive root will 
        never be a root. All these are exploited here to get the algebra constants.
        '''

        # Checks whether the roots are generated
        if self.rootsys.roots is None:
            self.rootsys.get_all_roots()

        # Create an empty matrix for the algebra constants. Not including Cartan generators as it needs the roots in vector form.
        # Also we know that Cartan Generator with the generator would just give the component of the root in the direction of the generator.
        # It is also know that positive and the negetive roots will just give the root dotted with the Cartan generators.
        self.algebra_constants = np.empty((self.rootsys.n_roots, self.rootsys.n_roots))
        self.algebra_constants.fill(np.nan)
        # We are storing the indices of the roots in the rootsys.roots list.

        # First we calculate the algebra constants or commuttors between positive roots as positive negetive requires them to be 
        # calculated first.
        self.__pos_struct()
        # Make all the algebra constants with nan's 0.
        self.algebra_constants = np.nan_to_num(self.algebra_constants)

        # Starting from the simple roots, we will calculate the algebra constants.
        for i in range(len(self.rootsys.PositiveRoots)):
            for j in range(len(self.rootsys.PositiveRoots[i])):
                for k in range(i+1):
                    for l in range(j+1):
                        
                        # New root in final state
                        newrootposneg = self.rootsys.PositiveRoots[i][j] - self.rootsys.PositiveRoots[k][l]

                        if i == k and l == j:
                            # For the pos neg case the end result is alpha dot H
                            # They obivously commute as they are the same terms
                            continue

                        # Computing for positive negetive first
                        if newrootposneg in self.rootsys.roots:
                            Term_C = self.rootsys.rootsyms.index(newrootposneg)
                            Term_A = self.rootsys.rootsyms.index(self.rootsys.PositiveRoots[i][j])
                            Term_Aneg = self.rootsys.rootsyms.index(-self.rootsys.PositiveRoots[i][j])
                            Term_B = self.rootsys.rootsyms.index(self.rootsys.PositiveRoots[k][l])
                            Term_Bneg = self.rootsys.rootsyms.index(-self.rootsys.PositiveRoots[k][l])
                            Term_D = self.rootsys.rootsyms.index(newrootposneg - self.rootsys.PositiveRoots[k][l])

                            # algebra constant for the positive negetive case
                            BdotC = self.rootsys.PositiveRoots[k][l] * newrootposneg
                            BdotC.expand().subs(self.rootsys.proper_repalcements)
                            struct_const = 1 / self.algebra_constants[Term_B, Term_C, Term_A] *(self.algebra_constants[Term_C, Term_Bneg, Term_D]*self.algebra_constants[Term_B, Term_D, Term_C] + BdotC)

                            # Using the antisymmetric property of the algebra constants to fill entries
                            self.algebra_constants[Term_A, Term_Bneg, Term_C] = self.algebra_constants[Term_Aneg, Term_B, Term_C] = struct_const
                            self.algebra_constants[Term_Bneg, Term_A, Term_C] = self.algebra_constants[Term_B, Term_Aneg, Term_C]= -struct_const


        return self.algebra_constants

    def __pos_struct(self):
        '''
        An internal function to calculate algebra constant or a general commutator between 2 roots, recursively.
        As of now this is not ooptimized as it is not possible to write it as tail recursion as it uses 2 calls.
        '''

        pass
