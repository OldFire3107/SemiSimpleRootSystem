from Semisimple import RootSystem, DynkinDraw
from time import sleep
import numpy as np

if __name__ == '__main__':
    rootsys = RootSystem.RootSystem.fromRootList([(0,1), (np.sqrt(3)/2,-3/2)])
    diagram = DynkinDraw.DynkinDraw(rootsys)
    diagram.make_fig()
    diagram.fig.savefig('dia.png')

    rootsys2 = RootSystem.RootSystem.fromCartanMatrix(np.array([[2, -1, 0], [-1, 2 , -1], [0, -2, 2]]))
    diagram2 = DynkinDraw.DynkinDraw(rootsys2)
    diagram2.make_fig()
    diagram2.fig.savefig('dia2.png')