import argparse
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from random import randint


# parámetros del graficador
rcParams['toolbar'] = 'None'
fig, ax = plt.subplots(num='FRUCHTERMAN', figsize = (8, 8))

# función graficadora
def plot(nodesPosDisp, listNodes, listEdges, fa, ps):
    flen = int(np.sqrt(fa))
    textOffset = 0.2
    ax.axis([-flen, flen, -flen, flen])

    for edge in listEdges:
        x1, y1 = nodesPosDisp[edge[0]][0]
        x2, y2 = nodesPosDisp[edge[1]][0]
        ax.plot([x1, x2], [y1, y2], color = 'blue', marker = 'o',
                 markeredgecolor = 'red', markerfacecolor = 'red')

    for node in listNodes:
        x, y = nodesPosDisp[node][0]
        ax.text(x + textOffset , y + textOffset, node)

    plt.pause(ps)
    ax.clear()

# inicializa los nodos en una posición aleatoria
def randomize_positions(listNodes, fa):
    nodesPosDisp = {}
    flen = int(np.sqrt(fa))
    listX = [randint(-flen, flen) for node in listNodes]
    listY = [randint(-flen, flen) for node in listNodes]

    for i, node in enumerate(listNodes):
        nodesPosDisp[node] = [np.array((listX[i], listY[i])), np.array((0, 0))]

    return nodesPosDisp

# calcula fuerza de atracción
def fatr(x, k): return x**2 / k

# calcula fuerza de repulsión
def frep(x, k): return k**2 / x

# retorna posición aleatoria
def random_position(pos, flen):
    pos[0] += randint(-flen / 2, flen / 2)
    pos[1] += randint(-flen / 2, flen / 2)

    return np.array((pos[0], pos[1]))

# función principal
def layout(G, verb, it, fa, t, ce, mfg, ps):
    listNodes, listEdges = G
    flen = int(np.sqrt(fa))
    k = np.sqrt(fa / len(listNodes))
    eps = 0.1

    """nodesPosDisp es la estructura de datos principal que usa el agoritmo.
    Las claves son los nodos del grafo y un valor asociado a una clave (nodo)
    es una lista con dos elementos tipo numpy.array cada uno, el primero guarda
    la posición actual del nodo y el segundo sirve de acumulador"""
    nodesPosDisp = {} 
    nodesPosDisp = randomize_positions(listNodes, fa)

    difPos = x = y = 0
    n1Pos = n2Pos = np.array((0, 0))
    n1Disp = n2Disp = np.array((0, 0))
    norm = 0

    for i in range(it):
        if verb == 1: print("\n---------------iteracion:",i, "---------------\n")
        for n1 in listNodes:
			# reiniciar acumuladores
            n1Disp = np.array((0, 0))
            n1Pos = nodesPosDisp[n1][0]

			# calcular fuerzas de repulsión
            for n2 in listNodes:
                if n2 != n1:
                    n2Pos = nodesPosDisp[n2][0]
                    difPos = n1Pos - n2Pos
                    norm = np.linalg.norm(difPos)
                    if(norm < eps): # prevenir nodos cercanos
                        n1Pos = random_position(nodesPosDisp[n1][0], flen)
                        nodesPosDisp[n1][0] = n1Pos
                        difPos = n1Pos - n2Pos
                        norm = np.linalg.norm(difPos)
                    n1Disp = n1Disp + (difPos / norm) * frep(norm, k)
                    if verb == 1: print("-frep(", n1,", ", n2,") = ",frep(norm,k), sep='')
            nodesPosDisp[n1][1] = n1Disp

        # calcular fuerzas de atracción
        for edge in listEdges:
            n1Pos = nodesPosDisp[edge[0]][0]
            n2Pos = nodesPosDisp[edge[1]][0]
            n1Disp = nodesPosDisp[edge[0]][1]
            n2Disp = nodesPosDisp[edge[1]][1]
            difPos = n1Pos - n2Pos
            norm = np.linalg.norm(difPos)
            nodesPosDisp[edge[0]][1] = n1Disp - (difPos / norm) * fatr(norm, k)
            nodesPosDisp[edge[1]][1] = n2Disp + (difPos / norm) * fatr(norm, k)
            if verb == 1: print("-fatr(", edge,") = ", fatr(norm,k), sep='')

        # limitar despazamiento máximo con t
        # prevenir que los nodos salgan del area
        # añadir gravedad
        # actualizar posiciones
        for n1 in listNodes:
            n1Pos = nodesPosDisp[n1][0]
            n1Disp = nodesPosDisp[n1][1]

            norm = np.linalg.norm(n1Disp)
            nodesPosDisp[n1][0] = n1Pos + (n1Disp / norm) * \
            min(norm, t) - (n1Pos / np.linalg.norm(n1Pos)) * \
            min(np.linalg.norm(n1Pos) * mfg, t)
            x = min(flen, max(-flen, nodesPosDisp[n1][0][0]))
            y = min(flen, max(-flen, nodesPosDisp[n1][0][1]))
            n1Pos = np.array((x, y))
            nodesPosDisp[n1][0] = n1Pos

        plot(nodesPosDisp, listNodes, listEdges, fa, ps)
        t = ce * t
        if verb == 1: print("-temp:",t)

# lee archivo de archivo y lo retorna en forma de lista
def lee_grafo_archivo(file_path):
    count, cantVer = 0, 0
    listaNodos, listaAristas =  [], []

    with open(file_path, 'r') as f:
        first_line = f.readline()
        cantVer = int(first_line)

        for line in f:
            node = ""
            edge = ["", ""]
            count += 1
            line = line.strip('\n')
            if count <= cantVer:
                listaNodos.append(line)
            else:
                    space = line.find(' ')
                    edge[0] = line[:space]
                    edge[1] = line[space+1:]
                    listaAristas.append(tuple(edge))
    return (listaNodos, listaAristas)


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("file", help="Ruta del archivo donde se encuntra el grafo")
    parser.add_argument("-v", type=int, choices=[0,1], default=0, help="Modo ""verbose")
    parser.add_argument("-it", nargs='?', type=int, default=40, help="Asigna cantidad de iteraciones, valor defecto: 40")
    parser.add_argument("-fa", nargs='?', type=int, default=100, help="Asigna area del marco o frame, valor defecto: 100")
    parser.add_argument("-t", nargs='?', type=int, default=30, help="Asigna temperatura inicial, valor defecto: 30")
    parser.add_argument("-ce", nargs='?', type=int, default=0.9, help="Asigna coeficiente de enfriamiento, valor defecto: 0.9")
    parser.add_argument("-mfg", nargs='?', type=int, default=0.1, help="Asigna modulo de la fuerza de gravedad, valor defecto: 0.1")
    parser.add_argument("-ps", nargs='?', type=int, default=0.01, help="Asigna tiempo límite en que se refrezca el dibujo: 0.01")
    return parser.parse_args()


def main():
    args = parser()
    G = lee_grafo_archivo(args.file)
    layout(G, args.v, args.it, args.fa, args.t, args.ce, args.mfg, args.ps)
    input("<presione alguna tecla para salir>")

if __name__ == "__main__":
    main()
