from os import wait
from time import sleep
from Semisimple import RootSystem, DynkinDraw
from time import sleep

if __name__ == '__main__':
    rootsys = RootSystem.RootSystem([(2,1)])
    diagram = DynkinDraw.DynkinDraw(rootsys)
    fig, ax = diagram.make_fig()
    ax.savefig('dia.png')