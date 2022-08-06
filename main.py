from Semisimple import RootSystem, DynkinDraw
import numpy as np

if __name__ == '__main__':
    rootsys = RootSystem.RootSystem.fromRootList([(np.sqrt(3)/2,-3/2), (0,1)])
    print(rootsys.get_all_roots())
    diagram = DynkinDraw.DynkinDraw(rootsys)
    diagram.make_fig()
    diagram.fig.savefig('dia.png')

    rootsys2 = RootSystem.RootSystem.fromCartanMatrix(np.array([[2, -1, 0], [-1, 2 , -1], [0, -2, 2]]))
    print(rootsys2.get_all_roots())
    diagram2 = DynkinDraw.DynkinDraw(rootsys2)
    diagram2.make_fig()
    diagram2.fig.savefig('dia2.png')

    rootsys3 = RootSystem.RootSystem.fromCartanMatrix(np.array([[2, -1], [-1, 2]]))
    fig, ax = rootsys3.plot_2D_space()
    fig.savefig('SU(3).png')