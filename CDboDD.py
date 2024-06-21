from math import sin
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

distMatrix = pd.read_csv("distMatrix.csv", encoding="utf-8")
distMatrix = distMatrix.values.tolist()

adjMatrix = pd.read_csv("adjMatrix.csv", encoding="utf-8")
adjMatrix = adjMatrix.values.tolist()

def deg(u, adjMatrix):
    return sum(adjMatrix[u])


def neighbor(u, adjMatrix):
    neighborList_no_self = []
    for i in range(len(adjMatrix[u])):
        if adjMatrix[u][i] == 1:
            neighborList_no_self.append(i)
    neighborList_with_self = neighborList_no_self
    neighborList_with_self.append(u)
    return set(neighborList_no_self), set(neighborList_with_self)


def distance(u, v, adjMatrix):
    up = down = 0 
    for x in list(neighbor(u, adjMatrix)[1].union(neighbor(v, adjMatrix)[1])):
        for y in list(neighbor(u, adjMatrix)[1].union(neighbor(v, adjMatrix)[1])):
            if x != y:
                down += adjMatrix[x][y]
    for x in list(neighbor(u, adjMatrix)[1].intersection(neighbor(v, adjMatrix)[1])):
        up += (adjMatrix[u][x] + adjMatrix[v][x])

    dist = 1 - up/down
    return dist


def direct_index(u, v, adjMatrix):
    return -(sin(1 - distance(u, v, adjMatrix))/deg(u, adjMatrix) + sin(1 - distance(u, v, adjMatrix))/deg(v, adjMatrix))


def common_index(u, v, adjMatrix):
    CI = 0
    for x in list(neighbor(u, adjMatrix)[0].intersection(neighbor(v, adjMatrix)[0])):
        CI -= ((1/deg(u, adjMatrix)) * sin(1 - distance(x, u, adjMatrix)) * (1 - distance(x, v, adjMatrix)) + 1/deg(v, adjMatrix) * sin(1 - distance(x, v, adjMatrix)) * (1 - distance(x, u, adjMatrix)))
    return CI


def exclusive_index(u, v, adjMatrix, lamda=0.5):
    EI = 0
    for x in list(neighbor(u, adjMatrix)[1].difference(neighbor(v, adjMatrix)[1])):
        if (1 - distance(x, v, adjMatrix)) >= lamda:
            rou = 1 - distance(x, v, adjMatrix)
        else:
            rou = 1 - distance(x, v, adjMatrix) - lamda
        EI -= (1/deg(u, adjMatrix)) * sin(1 - distance(x, u, adjMatrix)) * rou
    for y in list(neighbor(v, adjMatrix)[1].difference(neighbor(u, adjMatrix)[1])):
        if (1 - distance(y, u, adjMatrix)) >= lamda:
            rou = 1 - distance(y, u, adjMatrix)
        else:
            rou = 1 - distance(y, u, adjMatrix) - lamda
        EI -= (1/deg(v, adjMatrix)) * sin(1 - distance(y, v, adjMatrix)) * rou
    return EI


# Initialization of distances
for i in range(len(adjMatrix)):
    for j in range(len(adjMatrix[i])):
        if adjMatrix[i][j] != 1:
            continue
        else:
            distMatrix[i][j] = distance(i, j, adjMatrix)


# Dynamic Interaction
Flag = True
while Flag:
    Flag = False
    for i in range(len(adjMatrix)):
        for j in range(len(adjMatrix[i])):
            if adjMatrix[i][j] != 1:
                continue
            else:
                # print(i, j)
                DI = direct_index(i, j, adjMatrix)
                CI = common_index(i, j, adjMatrix)
                EI = exclusive_index(i, j, adjMatrix)
                delta = DI + CI + EI
                if delta != 0:
                    distMatrix[i][j] += delta
                    if distMatrix[i][j] > 1:
                        distMatrix[i][j] = 1
                    if distMatrix[i][j] < 0:
                        distMatrix[i][j] = 0
                    FLag = True


for i in range(len(distMatrix)):
    print(distMatrix[i])

# Find community
for i in range(len(adjMatrix)):
    for j in range(len(adjMatrix[i])): 
        if adjMatrix[i][j] != 1:
            continue
        else:
            if distMatrix[i][j] == 1:
                adjMatrix[i][j] = 0

G = nx.from_numpy_array(np.array(adjMatrix))
nx.draw(G, node_size=500, with_labels=True)
plt.show()

