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
            self.rootsys.get_all_roots()

        # Create an empty matrix for the structure constants. Not including Cartan generators as it needs the roots in vector form.
        # Also we know that Cartan Generator with the generator would just give the component of the root in the direction of the generator.
        # It is also know that positive and the negetive roots will just give the root dotted with the Cartan generators.
        self.structure_constants = np.zeros((self.rootsys.n_roots, self.rootsys.n_roots, self.rootsys.n_roots))
        # We are storing the indices of the roots in the rootsys.roots list.

        # Starting from the simple roots, we will calculate the structure constants.
        for i in range(len(self.rootsys.PositiveRoots)):
            for j in range(len(self.rootsys.PositiveRoots[i])):
                for k in range(i+1):
                    for l in range(j+1):
                        
                        # New root with both positve
                        newrootpospos = self.rootsys.PositiveRoots[i][j] + self.rootsys.PositiveRoots[k][l]
                        newrootposneg = self.rootsys.PositiveRoots[i][j] - self.rootsys.PositiveRoots[k][l]

                        if i == k and l == j:
                            # For the pos neg case the end result is alpha dot H
                            # They obivously commute as they are the same terms
                            continue

                        # Computing for positive negetive first
                        if newrootposneg in self.rootsys.roots:
                            mask = self.rootsys.q[i][j] != 0
                            roots_to_switch = np.where(mask)[0]

                        if newrootpospos in self.rootsys.roots:
                            struct_constant = 0

                            if not(self.rootsys.PositiveRoots[k][l] in self.rootsys.rootsyms):
                                mask = self.rootsys.q[i][j] != 0
                                roots_to_switch = np.where(mask)[0]
                                mask0 = self.rootsys.p[k][l] == 0
                                roots_to_exchagne = np.where(mask0)[0]

                                common = list(set(roots_to_switch).intersection(roots_to_exchagne))

                                if common:
                                    common = common[0]
                                else:
                                    common = roots_to_switch[0]
                                    print("Didn't work out!!\n\n\n\n\n") # Debugging statement 
                                
                                subcommutator_term = self.rootsys.PositiveRoots[i][j] - self.rootsys.rootsyms[common]

                                index1 = self.rootsys.roots.index(subcommutator_term)
                                index2 = self.rootsys.roots.index(self.rootsys.PositiveRoots[k][l])
                                # index3 = self.rootsys.roots.index(self.rootsys.PositiveRoots[i][j])

                                # struct_constant = 1 / self.structure_constants[common][index1][index3]

                                try:
                                    index0 = self.rootsys.roots.index(subcommutator_term + self.rootsys.PositiveRoots[k][l])

                                    struct_constant = self.structure_constants[index1][index2][index0] * struct_constant

                                    prefinal = subcommutator_term + self.rootsys.PositiveRoots[k][l]
                                    secondaryindex = self.rootsys.PositiveRoots[i+k].index(prefinal)

                                    j_val = (self.rootsys.q[i+k][secondaryindex][common] + self.rootsys.p[i+k][secondaryindex][common])/2
                                    m_val = (self.rootsys.q[i+k][secondaryindex][common] - self.rootsys.p[i+k][secondaryindex][common])/2

                                    struct_constant = struct_constant * np.sqrt((j_val + m_val + 1) * (j_val - m_val) / 2) * self.raising_norm[common]
                                except ValueError:
                                    struct_constant = 0

                            else:
                                common = self.rootsys.rootsyms.index(self.rootsys.PositiveRoots[k][l])
                                j_val = (self.rootsys.q[i][j][common] + self.rootsys.p[i][j][common])/2
                                m_val = (self.rootsys.q[i][j][common] - self.rootsys.p[i][j][common])/2
                                
                                struct_constant =  struct_constant * np.sqrt((j_val + m_val + 1) * (j_val - m_val) / 2) * self.raising_norm[common]
                            
                            a = self.rootsys.roots.index(self.rootsys.PositiveRoots[i][j])
                            b = self.rootsys.roots.index(self.rootsys.PositiveRoots[k][l])
                            c = self.rootsys.roots.index(newrootpospos)

                            print(struct_constant, a, b, c)
                           
                            self.structure_constants[a][b][c] = self.structure_constants[b][c][a] = self.structure_constants[c][a][b] = struct_constant
                            self.structure_constants[b][a][c] = self.structure_constants[a][c][b] = self.structure_constants[c][b][a] = -struct_constant

                            a = self.rootsys.roots.index(-self.rootsys.PositiveRoots[i][j])
                            b = self.rootsys.roots.index(-self.rootsys.PositiveRoots[k][l])
                            c = self.rootsys.roots.index(-newrootpospos)
                          
                            self.structure_constants[a][b][c] = self.structure_constants[b][c][a] = self.structure_constants[c][a][b] = struct_constant
                            self.structure_constants[b][a][c] = self.structure_constants[a][c][b] = self.structure_constants[c][b][a] = -struct_constant


        