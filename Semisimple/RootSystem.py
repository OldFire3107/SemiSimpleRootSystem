from sympy import symbols
import numpy as np

class RootSystem:
    '''
    RootSystem Class
    -----------------

    The master class that does operations related to the entire system
    of roots based on provided simple roots.
    '''

    def __init__(self, rootlist=None):
        assert rootlist is not None, 'Provide at least one root!'
        rootlist = np.array(rootlist)
        length = len(rootlist[0])
        self.dim = len(rootlist)
        self.rootsyms = symbols(f'alpha1:{self.dim+1}', commutative=False)
        for simpleroot in rootlist[1:]:
            assert len(simpleroot) == length , 'Root vector sizes inconsistent!'

        assert np.linalg.matrix_rank(rootlist) == self.dim, 'Simple Roots provided are not linearly independant!'

        self.simprootlist = rootlist
        self.gen_adj_list_and_cartan_mat()

    def gen_adj_list_and_cartan_mat(self):
        '''
        Computes and Returns the adjecency list of the system and its Cartan Matrix.
        '''

        self.adj_list = {}
        self.cartan_matrix = np.zeros((self.dim, self.dim), dtype=int)
        for i in range(len(self.simprootlist)):
            simple_root = self.simprootlist[i]
            result = 2. * np.dot(self.simprootlist, simple_root) / (np.linalg.norm(simple_root) ** 2)
            resultint = np.array(result, dtype=int)
            self.cartan_matrix[:,i] =  resultint
            assert np.isclose(result, resultint).all(), 'Invalid angle between the roots!'
            list_node = np.where(np.logical_not(np.logical_or(resultint == 0, resultint == 2)))
            self.adj_list[i] = list_node

        ## Checking it angles are valid.
        diag2 = np.diag(np.diag(self.cartan_matrix))
        check_mat = self.cartan_matrix - diag2
        num_list = np.unique(check_mat)
        for i in range(len(num_list)):
            assert num_list[i] in [0, -1, -2, -3], 'Invalid angle between the roots!'

        ## divided by right side and multiplied by 2
        self.replacements = [(self.rootsyms[i] * self.rootsyms[j], self.cartan_matrix[i][j]) for i in range(self.dim) for j in range(self.dim)]

        return self.adj_list, self.cartan_matrix

    def get_positive_roots_layered(self):
        '''
        Returns an array of all roots generated from the simple roots.
        '''
        PositiveRoots = []

        # Starting with the i th row of the cartan matrix as alpha_i
        # k = 0 is just the cartan matrix
        # k = 1 layer will always be there
        kzero = []

        for simpleroots in self.rootsyms:
            kzero.append(simpleroots)
        
        PositiveRoots.append(kzero)
        
        # Now for values k > 1
        k_val = 1
        Possible = True
        while Possible:
            Possible = False
            k_val += 1
            krowlist = []
            # Get the value of p associate with the k value
            for i in range(self.dim):
                for j in range(self.dim):
                    qpval = PositiveRoots[k_val-2][i] * self.rootsyms[j]
                    # We take the cartan matrix values, so the division by self.rootsyms[j]
                    # square is not necessary.
                    qpval = qpval.expand().subs(self.replacements)
                    print(qpval)

            

        return PositiveRoots

    def generate_simple_roots(self):
        '''
        Returns and generates simples roots from the Cartan Matrix in eucleadian geometry Vector Space.
        '''

        # ToDo
        

        pass
