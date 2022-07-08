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
        if rootsys is None:
            print('\nWARNING NO ROOT SYSTEM PROVIDED. Please provide one either\
                 by initializing it or at a later stage using set_root_sys method\n\n')
        self.rootsys = rootsys

    def make_fig(self, linecol = 'k', fillcol = 'w'):
        '''
        Parameters:
            - Line Colours (linecol) defaults to black
            - Fill Colours (fillcol) defaults to white
        '''
        # Form an adjcency list

        # Find lenght and heigth
        maxy = 0.5 # Proportional to side ways expansion
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
        ax.add_patch(pc.PathPatch(pc.Path([(0, 0), (0.8, 0)])))
        ax.add_patch(pc.PathPatch(pc.Path([(0.8, -0.05), (1.6, -0.05)])))
        ax.add_patch(pc.PathPatch(pc.Path([(0.8, 0.05), (1.6, 0.05)])))
        ax.add_patch(pc.PathPatch(pc.Path([(1.6, 0.), (2.4, 0)])))
        ax.add_patch(pc.PathPatch(pc.Path([(1.6, -0.1), (2.4, -0.1)])))
        ax.add_patch(pc.PathPatch(pc.Path([(1.6, 0.1), (2.4, 0.1)])))
        ax.add_patch(pc.Circle((0, 0), 0.2, edgecolor='k', facecolor='w'))
        ax.add_patch(pc.Circle((0.8, 0), 0.2, edgecolor='k', facecolor='w'))
        ax.add_patch(pc.Circle((1.6, 0), 0.2, edgecolor='k', facecolor='w'))
        ax.add_patch(pc.Circle((2.4, 0), 0.2, edgecolor='k', facecolor='w'))

        return ax, fig

    def set_root_sys(self, rootsys):
        '''
        Sets the root system for the diagram.
        Parameters:
            - RootSystem object
        '''

        assert(isinstance(rootsys, RootSystem), 'Dynkin Diagram was not provided a RootSystem object. Provide a RootSystem object.')        
        self.rootsys = rootsys


