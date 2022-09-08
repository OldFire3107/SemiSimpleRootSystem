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

        # Check if norms and dot product is defined
        if self.rootsys.raising_norm is None:
            self.raising_norm = self.rootsys.get_raising_norm()

        # Create an empty matrix for the algebra constants. Not including Cartan generators as it needs the roots in vector form.
        # Also we know that Cartan Generator with the generator would just give the component of the root in the direction of the generator.
        # It is also know that positive and the negetive roots will just give the root dotted with the Cartan generators.
        self.algebra_constants = np.empty((self.rootsys.n_roots, self.rootsys.n_roots))
        self.algebra_constants.fill(np.nan)
        # We are storing the indices of the roots in the rootsys.roots list.

        # First we calculate the algebra constants or commuttors between positive roots as positive negetive requires them to be 
        # calculated first.
        
        # First find the commutators with simple roots so that it can be used for recursive implementation.
        for i in range(len(self.rootsys.PositiveRoots)-1): # We are not considering the last one as it is the highest k value root.
            for j in range(len(self.rootsys.PositiveRoots[i])):
                for k in range(len(self.rootsys.rootsyms)):
                    indexa = self.rootsys.roots.index(self.rootsys.PositiveRoots[i][j])
                    indexb = self.rootsys.roots.index(self.rootsys.PositiveRoots[0][k])
                    indexaneg = self.rootsys.roots.index(-self.rootsys.PositiveRoots[i][j])
                    indexbneg = self.rootsys.roots.index(-self.rootsys.PositiveRoots[0][k])

                    if np.isnan(self.algebra_constants[indexa][indexb]): # If the value is not calculated yet.
                        if self.rootsys.p[i][j][k] > 0:
                            j_val = (self.rootsys.q[i][j][k] + self.rootsys.p[i][j][k])/2
                            m_val = (self.rootsys.q[i][j][k] - self.rootsys.p[i][j][k])/2
                            const_term = np.sqrt((j_val + m_val + 1) * (j_val - m_val) / 2) * self.rootsys.raising_norm[k]

                            self.algebra_constants[indexa][indexb] = self.algebra_constants[indexaneg][indexbneg] = -const_term # Convention followed for basic roots
                            self.algebra_constants[indexb][indexa] = self.algebra_constants[indexbneg][indexaneg] = const_term

                        else:
                            self.algebra_constants[indexa][indexb] = self.algebra_constants[indexaneg][indexbneg] = 0
                            self.algebra_constants[indexb][indexa] = self.algebra_constants[indexbneg][indexaneg] = 0


        # Recursive implementation of the commutator calculation.
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
                            Term_C = self.rootsys.roots.index(newrootposneg)
                            Term_A = self.rootsys.roots.index(self.rootsys.PositiveRoots[i][j])
                            Term_Aneg = self.rootsys.roots.index(-self.rootsys.PositiveRoots[i][j])
                            Term_B = self.rootsys.roots.index(self.rootsys.PositiveRoots[k][l])
                            Term_Bneg = self.rootsys.roots.index(-self.rootsys.PositiveRoots[k][l])
                                
                            struct_const = 0
                            if self.algebra_constants[Term_B, Term_C] != 0:
                                try:    
                                    Term_D = self.rootsys.roots.index(newrootposneg - self.rootsys.PositiveRoots[k][l])

                                    # algebra constant for the positive negetive case
                                    BdotC = self.rootsys.PositiveRoots[k][l] * newrootposneg
                                    BdotC = BdotC.expand().subs(self.rootsys.proper_replacements)
                                    struct_const = 1 / self.algebra_constants[Term_B, Term_C] *(self.algebra_constants[Term_C, Term_Bneg]*self.algebra_constants[Term_B, Term_D] + BdotC)
                                except ValueError:
                                    BdotC = self.rootsys.PositiveRoots[k][l] * newrootposneg
                                    BdotC = BdotC.expand().subs(self.rootsys.proper_replacements)
                                    struct_const = 1 / self.algebra_constants[Term_B, Term_C] *(BdotC)

                            
                            # Using the antisymmetric property of the algebra constants to fill entries
                            self.algebra_constants[Term_A, Term_Bneg] = self.algebra_constants[Term_Aneg, Term_B] = struct_const
                            self.algebra_constants[Term_Bneg, Term_A] = self.algebra_constants[Term_B, Term_Aneg] = -struct_const
            
 
        return self.algebra_constants

    def __pos_struct(self):
        '''
        An internal function to calculate algebra constant or a general commutator between 2 roots, recursively.
        As of now this is not ooptimized as it is not possible to write it as tail recursion as it uses 2 calls.
        '''

        # Give the input for the recursive function
        for i in range(1,len(self.rootsys.PositiveRoots)): # Skipping simple roots as they are already calculated.
            for j in range(len(self.rootsys.PositiveRoots[i])):
                for k in range(1,i+1): # Skipping simple roots as they are already calculated.
                    for l in range(j+1):
                        rootA = self.rootsys.roots.index(self.rootsys.PositiveRoots[i][j])
                        rootB = self.rootsys.roots.index(self.rootsys.PositiveRoots[k][l])

                        if np.isnan(self.algebra_constants[rootA][rootB]):
                            self.__recursive_cal_pos(self.rootsys.PositiveRoots[i][j], self.rootsys.PositiveRoots[k][l])

    
    def compute_commutator(self, root1, root2):
        '''
        Computes the commutator of 2 roots. The roots are given as symbolic constants.
        '''
        if root1 == root2:
            return 0
        else:
            return self.algebra_constants[self.rootsys.roots.index(root1), self.rootsys.roots.index(root2)]


    def __recursive_cal_pos(self, root1, root2):
        '''
        The main recursive function to calculate the algebra constants between 2 positive (negetive) roots.
        THIS HAS TO BE OPTIMIZED! TAIL RECURSION
        '''

        rootA = self.rootsys.roots.index(root1)
        rootB = self.rootsys.roots.index(root2)
        rootAneg = self.rootsys.roots.index(-root1)
        rootBneg = self.rootsys.roots.index(-root2)

        i, j = self.__index_2d(root1)
        k, l = self.__index_2d(root2)

        Add_entry = False

        if np.isnan(self.algebra_constants[rootA][rootB]):
            Add_entry = True
        else:
            return self.algebra_constants[rootA][rootB]

        if (root1 + root2) not in self.rootsys.roots:
            if Add_entry:
                self.algebra_constants[rootA][rootB] = 0
                self.algebra_constants[rootB][rootA] = 0
                self.algebra_constants[rootAneg][rootBneg] = 0
                self.algebra_constants[rootBneg][rootAneg] = 0
            return 0

        mask = self.rootsys.q[i][j] != 0
        roots_to_switch = np.where(mask)[0]
        mask0 = self.rootsys.p[k][l] == 0
        roots_to_exchagne = np.where(mask0)[0]

        # The above step is to eliminate the roots which are not in the positive root system.
        # This is done by checking the q and p matrices.
        common = list(set(roots_to_switch).intersection(roots_to_exchagne))

        if len(roots_to_switch) == 0:
            if Add_entry:
                self.algebra_constants[rootA][rootB] = 0
                self.algebra_constants[rootB][rootA] = 0
                self.algebra_constants[rootAneg][rootBneg] = 0
                self.algebra_constants[rootBneg][rootAneg] = 0
            return 0

        if common:
            common = common[0]
            subcommutator_term11 = root1 - self.rootsys.rootsyms[common]
            subcommutator_term12 = subcommutator_term11 + root2
            subcommutator_term11_index = self.rootsys.roots.index(subcommutator_term11)
            try:
                subcommutator_term12_index = self.rootsys.roots.index(subcommutator_term12)
            except ValueError:
                if Add_entry:
                    self.algebra_constants[rootA][rootB] = 0
                    self.algebra_constants[rootB][rootA] = 0
                    self.algebra_constants[rootAneg][rootBneg] = 0
                    self.algebra_constants[rootBneg][rootAneg] = 0
                return 0

            in1, in2 = self.__index_2d(subcommutator_term12)
            indexsimp = common
            if self.rootsys.p[in1][in2][indexsimp] == 0:
                if Add_entry:
                    self.algebra_constants[rootA][rootB] = 0
                    self.algebra_constants[rootB][rootA] = 0
                    self.algebra_constants[rootAneg][rootBneg] = 0
                    self.algebra_constants[rootBneg][rootAneg] = 0
                return 0

            const = 1 / self.algebra_constants[indexsimp][subcommutator_term11_index]
            if np.isnan(self.algebra_constants[subcommutator_term11_index][rootB]):
                const *= self.__recursive_cal_pos(subcommutator_term11, root2)
            else:
                const *= self.algebra_constants[subcommutator_term11_index][rootB]

            const *= self.algebra_constants[indexsimp][subcommutator_term12_index]

            if Add_entry:
                self.algebra_constants[rootA][rootB] = -const
                self.algebra_constants[rootB][rootA] = const
                self.algebra_constants[rootAneg][rootBneg] = -const
                self.algebra_constants[rootBneg][rootAneg] = const

            return -const

        else:
            common = roots_to_switch[0]
            indexsimp = common

            # Seperate the roots into 2 parts
            subcommutator_term11 = root1 - self.rootsys.rootsyms[common]
            subcommutator_term12 = subcommutator_term11 + root2
            subcommutator_term11_index = self.rootsys.roots.index(subcommutator_term11)

            const = 1 / self.algebra_constants[indexsimp][subcommutator_term11_index]

            const_term_1 = 1
            const_term_2 = 1


            try:
                subcommutator_term12_index = self.rootsys.roots.index(subcommutator_term12)
                in1, in2 = self.__index_2d(subcommutator_term12)

                # Commuting the first term just like in the if part.
                if self.rootsys.p[in1][in2][indexsimp] == 0: # Only second term matters
                    const_term_1 = 0
                else:
                    if np.isnan(self.algebra_constants[subcommutator_term11_index][rootB]):
                        const_term_1 *= self.__recursive_cal_pos(subcommutator_term11, root2)
                    else:
                        const_term_1 *= self.algebra_constants[subcommutator_term11_index][root2]

                    const_term_1 *= self.algebra_constants[indexsimp][subcommutator_term12_index] # Since simple root stuff already calculated
            except ValueError:
                const_term_1 = 0
            
            # Commuting the second term
            subcommutator_term21 = root2 + self.rootsys.rootsyms[common]
            subcommutator_term21_index = self.rootsys.roots.index(subcommutator_term21)
            const_term_2 *= self.algebra_constants[indexsimp][subcommutator_term21_index] # Since simple root stuff already calculated

            if np.isnan(self.algebra_constants[subcommutator_term11_index][subcommutator_term21_index]):
                const_term_2 *= self.__recursive_cal_pos(subcommutator_term11 ,subcommutator_term21)
            else:
                const_term_2 *= self.algebra_constants[subcommutator_term11_index][subcommutator_term21_index]

            const *= const_term_1 - const_term_2

            if Add_entry:
                self.algebra_constants[rootA][rootB] = -const
                self.algebra_constants[rootB][rootA] = const
                self.algebra_constants[rootAneg][rootBneg] = -const
                self.algebra_constants[rootBneg][rootAneg] = const

            return -const



    def __index_2d(self, v):
        '''
        Finds the 2D index of a root in the positive root system.
        '''
        for i, x in enumerate(self.rootsys.PositiveRoots):
            if v in x:
                return i, x.index(v)
