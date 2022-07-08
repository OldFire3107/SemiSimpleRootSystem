from cmath import isclose
from matplotlib.pyplot import axis
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
        dim = len(rootlist)
        for simpleroot in rootlist[1:]:
            assert len(simpleroot) == length , 'Root vector sizes inconsistent!'

        assert np.linalg.matrix_rank(rootlist) == dim, 'Simple Roots provided are not linearly independant!'

        self.simprootlist = rootlist

        # Getting the adjecency list
        self.adj_list = {}
        self.cartan_matrix = np.zeros((dim,dim), dtype=int)
        for i in range(len(self.simprootlist)):
            simple_root = self.simprootlist[i]
            result = 2. * np.dot(self.simprootlist, simple_root) / (np.linalg.norm(simple_root) ** 2)
            resultint = np.array(result, dtype=int)
            assert np.isclose(result, resultint).all(), 'Invalid angle between the roots!'
            print(resultint)
            list_node = np.where(np.logical_not(np.logical_or(resultint == 0, resultint == 2)))
            self.adj_list[i] = list_node

        print(self.adj_list)
        
        

        
    def get_simple_roots(self):
        '''
        Returns the simple roots as an array of row vectors.
        '''

        return self.simprootlist

    def get_cartan_matrix(self):
        '''
        Returns the cartan matrix of the system.
        '''

        return self.cartan_matrix

    def get_adj_list(self):
        '''
        Returns the adjecency list of the system.
        '''

        return self.adj_list