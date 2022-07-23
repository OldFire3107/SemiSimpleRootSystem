from queue import Queue
import numpy as np
from matplotlib import rcParams
import matplotlib.patches as pc
import matplotlib.pyplot as plt

from . import RootSystem

class DynkinDraw:
    '''
    DynkinDraw Class
    -----------------

    Takes a RootSystem class object. This class mainly is used to
    draw or save a dynkin diagram for a particular semi-simple Lie
    algebra to a file.
    '''

    def __init__(self, rootsys=None):
        rcParams.update({'figure.autolayout': True})
        if rootsys is None or not type(rootsys) == RootSystem.RootSystem:
            print('\nWARNING NO ROOT SYSTEM PROVIDED. Please provide one either\
by initializing it or at a later stage using set_root_sys method\n\n')
        self.rootsys = rootsys
        self.form_size_adjlist()

    def form_size_adjlist(self):
        '''
        Stores and returns the number of ode along depth and heigth or
        breadth and length of the Dynkin Diagram.
        '''
        # Form an adjcency list
        rootnums = np.linspace(0, self.rootsys.dim-1, self.rootsys.dim, dtype=int)
        # Uses the rootsys cartan matrix to make adjlist
        adjlist = {}

        breadthy_dia_num = 1

        remove_next_list = Queue()
        while not len(rootnums) == 0:
            if remove_next_list.empty():
                remove_next = 0
                row = rootnums[remove_next]
                rootnums = np.delete(rootnums, remove_next)
            else:
                remove_next = remove_next_list.get()
                row = remove_next
                rootnums = np.delete(rootnums, np.where(rootnums == remove_next))
            remove_next = 0
            layer_adj_list = {}
            count = 0
            for i in range(self.rootsys.dim):
                if i in rootnums and not i == row:
                    adjnode = self.rootsys.cartan_matrix[row][i]
                    if not adjnode == 0:
                        count += 1
                        remove_next_list.put(i)
                        if adjnode == -1 and self.rootsys.cartan_matrix[i][row] == -1:
                            layer_adj_list[i] = 1
                        elif adjnode == -1:
                            layer_adj_list[i] = -1 * self.rootsys.cartan_matrix[i][row]
                        else:
                            layer_adj_list[i] = -1 * adjnode
            if not count == 0:
                breadthy_dia_num = breadthy_dia_num * count
            adjlist[row] = layer_adj_list

        notEnd = True
        length_dia_num = 1
        listofkeys = [0]
        while notEnd:
            notEnd = False
            newlistofkeys = []
            for keystosearch in listofkeys:
                for akey in adjlist[keystosearch].keys():
                    newlistofkeys.append(akey)
                    if akey in adjlist.keys():
                        notEnd = True

            listofkeys = newlistofkeys
            if notEnd:
                length_dia_num += 1

        self.adjlist = adjlist
        self.breadth_num = breadthy_dia_num
        self.length_num = length_dia_num

    def make_fig(self, linecol = 'k', fillcol = 'w'):
        '''
        Parameters:
            - Line Colours (linecol) defaults to black
            - Fill Colours (fillcol) defaults to white
        '''

        # Find length and heigth
        maxy = 0.5 * self.breadth_num # Proportional to side ways expansion
        lowx = -0.3 # Fixed as diagram expands rightwards
        maxx = 2.7 # Based on horizontal expansion


        # make the figure
        fig, ax = plt.subplots()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_xlim(lowx, maxx)
        ax.set_ylim(-maxy, maxy)

        # set aspect ratio to make sure its not streched
        ax.set_aspect(1)

        # Loop over the adjecency list to draw dynkin diagram
        reflines = [0]
        list_of_points = []
        for i in range(len(self.adjlist)):
            newreflines = []
            minstep = None # Make stuff


        ax.add_patch(pc.PathPatch(pc.Path([(0, 0), (0.8, 0)])))
        ax.add_patch(pc.PathPatch(pc.Path([(0.8, -0.05), (1.6, -0.05)])))
        ax.add_patch(pc.PathPatch(pc.Path([(0.8, 0.05), (1.6, 0.05)])))
        ax.add_patch(pc.PathPatch(pc.Path([(1.6, 0.), (2.4, 0)])))
        ax.add_patch(pc.PathPatch(pc.Path([(1.6, -0.1), (2.4, -0.1)])))
        ax.add_patch(pc.PathPatch(pc.Path([(1.6, 0.1), (2.4, 0.1)])))

        for point in list_of_points:
            ax.add_patch(pc.Circle(point, 0.2, edgecolor=linecol, facecolor=fillcol))
        # ax.add_patch(pc.Circle((0.8, 0), 0.2, edgecolor='k', facecolor='w'))
        # ax.add_patch(pc.Circle((1.6, 0), 0.2, edgecolor='k', facecolor='w'))
        # ax.add_patch(pc.Circle((2.4, 0), 0.2, edgecolor='k', facecolor='w'))

        self.ax = ax
        self.fig = fig

        return ax, fig


    def set_root_sys(self, rootsys):
        '''
        Sets the root system for the diagram.
        Parameters:
            - RootSystem object
        '''

        assert isinstance(rootsys, RootSystem), 'Dynkin Diagram was not provided a RootSystem object. Provide a RootSystem object.'     
        self.rootsys = rootsys
        self.form_size_adjlist
