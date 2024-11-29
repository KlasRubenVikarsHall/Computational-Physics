import copy

import numpy as np
from numpy.random import rand, randn
from matplotlib import pyplot as plt


def initialize_system(M=100): #Creates a town of M houses on an identity square.
    city = [None for j in range(M)]
    for i in range(M):
        city[i] = [rand(), rand()]
    for i in range(M):
        plt.scatter(city[i][0],city[i][1])
    plt.show()
    return city


def theoretical_solution(city): #Solve by combinatoric
    pass


def total_length(route):
    length = 0
    for i in range(len(route) - 1):
        length += metric(route[i], route[i + 1])
    length += metric(route[-1], route[0])
    return length


def metric(a, b):
    distance = np.sqrt((a[0]-b[0]) ** 2 + (a[1]-b[1]) ** 2)
    return distance


def exchange(route, i, j):
    a = copy.deepcopy(route[i])
    b = copy.deepcopy(route[j])
    route[i] = copy.deepcopy(b)
    route[j] = copy.deepcopy(a)
    return route


if __name__ == "__main__":
    initialize_system(5)