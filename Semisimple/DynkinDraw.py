from queue import Queue
import numpy as np
from matplotlib import rcParams
import matplotlib.patches as pc
import matplotlib.pyplot as plt

from . import RootSystem

rcParams.update({'figure.autolayout': True})

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
                            layer_adj_list[i] = self.rootsys.cartan_matrix[i][row]
                        else:
                            layer_adj_list[i] = -1 * adjnode
            if not count == 0:
                breadthy_dia_num = breadthy_dia_num * count
            adjlist[row] = layer_adj_list
            adjlist[row]['layer'] = 0 # Will be updated later

        notEnd = True
        length_dia_num = 1
        listofkeys = [0]
        adjkeys = list(adjlist.keys())
        while notEnd:
            notEnd = False
            newlistofkeys = []
            for keystosearch in listofkeys:
                try:
                    adjkeys.remove(keystosearch)
                except ValueError:
                    print('Warning: Root System may have a cyclic structure. It may not be fully supported at the moment.')
                    pass
                for akey in adjlist[keystosearch].keys():
                    adjlist[keystosearch]['layer'] = length_dia_num
                    if akey == 'layer':
                        continue
                    newlistofkeys.append(akey)
                    if akey in adjlist.keys():
                        notEnd = True

            # If they are not connected start antother layer
            if not notEnd and len(adjkeys) != 0:
                newlistofkeys.append(adjkeys[0])
                notEnd = True

            listofkeys = newlistofkeys
            if notEnd:
                length_dia_num += 1

        self.adjlist = adjlist
        self.breadth_num = breadthy_dia_num
        self.length_num = length_dia_num
        
        return self.adjlist, self.breadth_num, self.length_num

    def make_fig(self, linecol = 'k', fillcol = 'w'):
        '''
        Parameters:
            - Line Colours (linecol) defaults to black
            - Fill Colours (fillcol) defaults to white
        '''

        # Find length and heigth
        maxy = 0.5 * self.breadth_num # Proportional to side ways expansion
        lowx = -0.3 # Fixed as diagram expands rightwards
        maxx = 0.3 + 0.8 * (self.length_num - 1) # Based on horizontal expansion


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
        list_of_points = [(0, 0)]
        y_range = {0: (maxy, -maxy)}
        to_beadded = None
        i = 0
        for keybase, valbase in self.adjlist.items():

            if to_beadded is not None:
                list_of_points.append(to_beadded)
                to_beadded = None

            startx = 0.8 * (self.adjlist[i]['layer'] - 1)
            endx = 0.8 * self.adjlist[i]['layer']
            try:
                upper_y = y_range[keybase][0]
                lower_y = y_range[keybase][1]
            except KeyError: # If they are not connected start antother layer
                upper_y = maxy
                lower_y = -maxy
            ref = (upper_y + lower_y) / 2.
            num_nodes_next = len(self.adjlist[i]) - 1
            if num_nodes_next == 0:
                step = 1 # If next layer exists
            else:
                step = (upper_y - lower_y) / num_nodes_next
            start = upper_y - step / 2.
            j = 0
                
            for key, values in self.adjlist[i].items():
                if len(self.adjlist[i].items()) == 1:
                    to_beadded = (endx, start)
                    break
                if key == 'layer':
                    continue
                if abs(values) == 1:
                    ax.add_patch(pc.PathPatch(pc.Path([(startx, ref), (endx, start - step * j)])))
                if abs(values) == 2:
                    ax.add_patch(pc.PathPatch(pc.Path([(startx, ref+0.05), (endx, start + 0.05 - step * j)])))
                    ax.add_patch(pc.PathPatch(pc.Path([(startx, ref-0.05), (endx, start - 0.05 - step * j)])))
                if abs(values) == 3:
                    ax.add_patch(pc.PathPatch(pc.Path([(startx, ref), (endx, start - step * j)])))
                    ax.add_patch(pc.PathPatch(pc.Path([(startx, ref-0.1), (endx, start - 0.1 - step * j)])))
                    ax.add_patch(pc.PathPatch(pc.Path([(startx, ref+0.1), (endx, start + 0.1 - step * j)])))

                # Find mid point
                midx = (startx + endx) / 2
                slide_param = 0.0005 # arrow has a definite size

                # Calculates the coordinates of arrows
                midx1 = midx + slide_param
                midx2 = midx - slide_param
                dmidx1 = midx1 - midx2
                midy1 = ref + (midx1 - startx) * (start - step * j - ref) / (endx - startx)
                midy2 = ref + (midx2 - startx) * (start - step * j - ref) / (endx - startx)
                dmidy1 = midy1 - midy2

                if abs(values) > 1:
                    arrow = None
                    if values < 0:
                        midx1 = midx - 0.05 + slide_param # For arrow adjustment as the tail appears otherwise
                        midy1 = ref + (midx1 - startx) * (start - step * j - ref) / (endx - startx)
                        arrow = pc.FancyArrow(midx1, midy1, -dmidx1, -dmidy1, color=linecol, width=0,
                                head_width=0.3, head_length=0.1, length_includes_head=True, overhang=0.9)
                    else:
                        midx2 = midx + 0.05 - slide_param # For arrow adjustment as the tail appears otherwise
                        midy2 = ref + (midx2 - startx) * (start - step * j - ref) / (endx - startx)
                        arrow = pc.FancyArrow(midx2, midy2, dmidx1, dmidy1, color=linecol, width=0,
                                head_width=0.3, head_length=0.1, length_includes_head=True, overhang=0.9)

                    ax.add_patch(arrow) # plots the arrow

                y_range[key] = (upper_y, upper_y - step * num_nodes_next)
                list_of_points.append((endx, start - step * j))
                j += 1
            
            i += 1

        for point in list_of_points:
            ax.add_patch(pc.Circle(point, 0.2, edgecolor=linecol, facecolor=fillcol))

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
