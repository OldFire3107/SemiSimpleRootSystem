from sympy import symbols
import numpy as np
import itertools

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

        if rootlist is not None:
            rootlist = np.array(rootlist)
            length = len(rootlist[0])
            self.dim = len(rootlist)
            self.rootsyms = symbols(f'alpha1:{self.dim+1}', commutative=False)
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
            result = 2. * np.dot(self.simprootlist, simple_root) / (np.linalg.norm(simple_root) ** 2)
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
            

        return PositiveRoots

    def get_all_roots(self):
        Positive_roots = self.get_positive_roots_layered()
        All_roots = list(itertools.chain(*Positive_roots))
        for i in range(len(All_roots)):
            nege = All_roots[i]*-1
            All_roots.append(nege.expand())
        return All_roots

    def generate_simple_roots(self):
        '''
        Returns and generates simples roots from the Cartan Matrix in eucleadian geometry Vector Space.
        '''

        # Starts by taking a root to be (1, 0, 0, ...dims)
        rootlist = []
        root1 = np.zeros(self.dim)
        root1[0] = 1.
        rootlist.append(root1)
        
        # Now generate other roots based on this and the cartan matrix
        for i in range(1, self.dim):
            # Process all other roots
            # Rotation
            pass
