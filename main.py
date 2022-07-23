from os import wait
from time import sleep
from Semisimple import RootSystem, DynkinDraw
from time import sleep
import numpy as np

if __name__ == '__main__':
    rootsys = RootSystem.RootSystem.fromCartanMatrix(np.array([[2, -1, 0], [-1, 2 , -1], [0, -2, 2]]))
    diagram = DynkinDraw.DynkinDraw(rootsys)
    fig, ax = diagram.make_fig()
    ax.savefig('dia.png')