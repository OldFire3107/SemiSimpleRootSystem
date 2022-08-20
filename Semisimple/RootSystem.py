from iteration_utilities import deepflatten
from sympy import symbols
import numpy as np
import itertools
from matplotlib import rcParams
import matplotlib.pyplot as plt

rcParams.update({'figure.autolayout': True})

# An important note is that while calculating p-q values using the replacements for the symbols from
# the cartan matrix initialization or when they are found while finding cartan matrix are actually
# not alpha_i * alpha_j but alpha_i * alpha_j / alpha_i ** 2 and hence should not be used for other 
# calculations as its non commutative. This is used to make things easier for finding all the roots.
class RootSystem:
    '''
    RootSystem Class
    -----------------

    The master class that does operations related to the entire system
    of roots based on provided simple roots.

    - Use fromCartanMatrix Method if needed to be initialized from Cartan Matrix
    - Use fromRootList if needed to be initialized from a root vector provided 
      as in eucleadian vector space

    In case the cartan matrix and such variables are inaccessible make sure that
    they are generated first.
    '''

    def __init__(self, rootlist=None, cart_mat=None):

        self.roots = None # Computationally intensive, only generates when function called or needed
        self.PositiveRoots = None
        self.q = None
        self.p = None
        self.basic_commutators = None

        if rootlist is not None:
            rootlist = np.array(rootlist)
            length = len(rootlist[0])
            self.dim = len(rootlist)
            self.rootsyms = symbols(f'alpha1:{self.dim+1}', commutative=False)
            # Symbols can't be replaced directly with arrays so have to store it as a 2D matrix so they can be substituted
            # in a loop multiple times to get the list later on.
            self.simpreplace = [[(self.rootsyms[i], rootlist[i][j]) for i in range(self.dim)] for j in range(self.dim)]
            for simpleroot in rootlist[1:]:
                assert len(simpleroot) == length , 'Root vector sizes inconsistent!'

            assert np.linalg.matrix_rank(rootlist) == self.dim, 'Simple Roots provided are not linearly independant!'

            self.simprootlist = rootlist
            self.gen_cartan_mat()

        else:
            self.cartan_matrix = cart_mat
            self.dim = len(cart_mat)
            self.rootsyms = symbols(f'alpha1:{self.dim+1}', commutative=False)
            self.replacements = [(self.rootsyms[i] * self.rootsyms[j], self.cartan_matrix[i][j]) for i in range(self.dim) for j in range(self.dim)]
            self.simprootlist = None # Would have to generate it
            self.simpreplace = None # Would have to generate it with the above list

    @classmethod
    def fromRootList(cls, RootList):
        '''
        Initialize the class from the Eucleadian Vector Space Roots
        '''

        assert RootList is not None, 'Provide at least one root!'
        return cls(RootList, None)

    @classmethod
    def fromCartanMatrix(cls, CartMat):
        '''
        Initialize the class from the Cartan Matrix
        '''
        assert CartMat is not None, 'Provide a non empty Cartan Matrix!'

        CartMat = np.array(CartMat, dtype=int)
        
        ## Checking it angles are valid.
        diag = np.diag(CartMat)
        diaglist = np.unique(diag)
        assert diaglist[0] == 2 and diaglist.size == 1, 'Invalid Cartan Matrix'
        diag2 = np.diag(diag)
        check_mat = CartMat - diag2
        num_list = np.unique(check_mat)
        for i in range(len(num_list)):
            assert num_list[i] in [0, -1, -2, -3], 'Invalid angle between the roots!'

        return cls(None, CartMat)

    def gen_cartan_mat(self):
        '''
        Computes and Returns the adjecency list of the system and its Cartan Matrix.
        '''

        self.adj_list = {}
        self.cartan_matrix = np.zeros((self.dim, self.dim), dtype=int)
        for i in range(len(self.simprootlist)):
            simple_root = self.simprootlist[i]
            # Calculating the column of cartan matrix
            result = 2. * np.dot(self.simprootlist, simple_root) / (np.linalg.norm(simple_root) ** 2)
            # Cartan matrix can only be integers
            resultint = np.array(result, dtype=int)
            self.cartan_matrix[:,i] =  resultint
            assert np.isclose(result, resultint).all(), 'Invalid angle between the roots!'
            # list_node = np.where(np.logical_not(np.logical_or(resultint == 0, resultint == 2)))
            # self.adj_list[i] = list_node

        ## Checking it angles are valid.
        diag2 = np.diag(np.diag(self.cartan_matrix))
        check_mat = self.cartan_matrix - diag2
        num_list = np.unique(check_mat)
        for i in range(len(num_list)):
            assert num_list[i] in [0, -1, -2, -3], 'Invalid angle between the roots!'

        ## divided by right side and multiplied by 2
        self.replacements = [(self.rootsyms[i] * self.rootsyms[j], self.cartan_matrix[i][j]) for i in range(self.dim) for j in range(self.dim)]

        return self.cartan_matrix

    def get_positive_roots_layered(self):
        '''
        Returns an array of all roots generated from the simple roots.
        '''
        PositiveRoots = []

        # Starting with the i th row of the cartan matrix as alpha_i
        # k = 0 is just the cartan matrix
        # k = 1 layer will always be there
        kzero = []
        self.q = []
        self.p = []

        for simpleroots in self.rootsyms:
            kzero.append(simpleroots)
        
        PositiveRoots.append(kzero)
        # for values k > 1
        # Now 
        q_vals = [np.zeros(self.dim, dtype=int) for i in range(self.dim)]
        for i in range(self.dim):
            q_vals[i][i] = 2
        q_vals = np.array(q_vals, dtype=int)
        k_val = 1

        base_qp = []
        # Get the q-p vals for all k=1 layer roots.
        for i in range(self.dim):
            qp_row = []
            for j in range(self.dim):
                qpval = PositiveRoots[k_val-2][i] * self.rootsyms[j]
                # We take the cartan matrix values, so the division by self.rootsyms[j]
                # square is not necessary.
                qpval = qpval.expand().subs(self.replacements)
                qp_row.append(qpval)

            base_qp.append(qp_row)
        base_qp = np.array(base_qp, dtype=int)
        qpval = base_qp

        Possible = True
        while Possible:
            Possible = False
            new_qvals = []
            new_qpvals = []
            k_val += 1
            krowlist = []

            p_vals = q_vals - qpval
            self.p.append(p_vals)
            self.q.append(q_vals)
            for i in range(len(p_vals)):
                for j in range(self.dim):
                    if p_vals[i][j] > 0:
                        Possible = True
                        newroot = PositiveRoots[k_val-2][i] + self.rootsyms[j]
                        newqp_val = qpval[i] + base_qp[j]

                        inList = False

                        for no in range(len(krowlist)):
                            if newroot == krowlist[no]:
                                inList = True
                                new_qvals[no][j] += 1

                        if not inList:
                            inList = True
                            krowlist.append(newroot)
                            new_qpvals.append(newqp_val)
                            qval = []
                            for k in range(self.dim):
                                if k == j:
                                    qval.append(q_vals[i][k] + 1)
                                else:
                                    qval.append(0)
                            
                            new_qvals.append(qval)
            
            if len(krowlist) == 0:
                break

            PositiveRoots.append(krowlist)
            q_vals = np.array(new_qvals, dtype=int)
            qpval = np.array(new_qpvals, dtype=int)
            
        self.PositiveRoots = PositiveRoots
        return PositiveRoots

    def get_all_roots(self):
        '''
        Returns all the posible roots.
        '''
        if self.PositiveRoots is None:
            self.get_positive_roots_layered()
        All_roots = list(itertools.chain(*self.Positive_roots))
        for i in range(len(All_roots)):
            nege = All_roots[i]*-1
            All_roots.append(nege.expand())

        self.roots = All_roots
        return All_roots

    def plot_2D_space(self):
        '''
        Plots the root system (only 2 dimensional system) on the cartesian 2D plane.
        Note that the first root entry in cartan matrix corresponds to the x axis.

        Returns the matplotlib figure and axis objects.
        '''
        assert self.dim == 2 , f'Not a 2 dimentional root system, it has {self.dim} simple roots instead of 2'

        # Check if the roots are defined else it is made
        # Making roots if not defined
        if self.simprootlist is None:
            root1 = (1, 0) # Fixing the first root
            root2 = None
            angle_cart = self.cartan_matrix[0][1]
            if angle_cart == -1:
                angle_cart = -1 * self.cartan_matrix[1][0]

            # Finding the 2nd root based on the 3 angles
            multiplier = None
            if abs(angle_cart) == 1: # 120 deg
                multiplier = 1
                root2 = (-0.5, np.sqrt(3)/2)
            elif abs(angle_cart) == 2: # 135 deg
                multiplier = np.sqrt(2) 
                root2 = (-np.sqrt(2), np.sqrt(2))
            elif abs(angle_cart) == 3: # 150 deg
                multiplier = np.sqrt(3) 
                root2 = (-np.sqrt(3)/2, 0.5)
            else:
                multiplier = 1 # It is 2 independant SU(2) system
                print('Here both the roots are taken to be of the same magnitude but it can be different as they are 2 independant SU(2)  systems.')
                root2 = (0, 1) # 90 deg

            if angle_cart < 0:
                root2 = [rootx / multiplier for rootx in root2]
            else:
                root2 = [rootx * multiplier for rootx in root2]

            self.simprootlist = np.array([root1, root2])
        
        if self.simpreplace is None:
            self.simpreplace = [[(self.rootsyms[i], self.simprootlist[i][j]) for i in range(self.dim)] for j in range(self.dim)]

        # Generate all the roots if not done already
        if self.roots is None:
            self.get_all_roots()

        # Substituting all the alphas in the root list
        replaced_roots = [tuple([a_root.subs(repl) for repl in self.simpreplace]) for a_root in self.roots]

        fig, ax = plt.subplots()
        ax.scatter(*zip(*replaced_roots))
        ax.scatter(0, 0)

        return fig, ax
    
    def gen_basic_commutators(self):
        '''
        Computes the basic commutators for the root system and stores them.

        Returns the basic commutators.
        '''

        self.basic_commutators = {}
        
        # Checks if the positive roots are generated
        if self.PositiveRoots is None:
            self.get_positive_roots_layered()

        for i in range(len(self.PositiveRoots)-1):
            for j in range(len(self.PositiveRoots[i])):
                for k in range(len(self.PositiveRoots[0])):
                    if self.p[i][j][k] > 0:
                        j_val = (self.q[i][j][k] + self.p[i][j][k])/2
                        m_val = (self.q[i][j][k] - self.p[i][j][k])/2
                        check_root = self.PositiveRoots[i][j] + self.PositiveRoots[0][k]

                        # Check if the root is already in the dictionary of commutators
                        if check_root in self.basic_commutators:
                            break

                        try:
                            const_term, commutator = self.basic_commutators[self.PositiveRoots[i][j]]
                        except:
                            const_term = 1
                            commutator = self.PositiveRoots[i][j]

                        const_term = const_term * 1 / np.sqrt((j_val + m_val + 1) * (j_val - m_val) / 2)
                        
                        self.basic_commutators[check_root] = (const_term, [self.PositiveRoots[0][k], commutator])

        return self.basic_commutators
        