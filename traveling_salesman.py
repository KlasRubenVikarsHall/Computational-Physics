import copy

import numpy as np
from numpy.random import rand, randn
from matplotlib import pyplot as plt
import pickle


def initialize_system(M=100): #Creates a town of M houses on an identity square.
    city = [None for j in range(M)]
    for i in range(M):
        city[i] = [rand(), rand()]
    # for i in range(M):
    #     plt.scatter(city[i][0],city[i][1])
    # plt.show()
    distances = np.zeros((M,M))
    for i in range(M):
        for j in range(M):
            if i == j:
                distances[i][j] = 10
            else:
                distances[i][j] = metric(city[i], city[j])
    return city, distances


def optimal_solution(city):
    pass


def metric(a, b):
    distance = np.sqrt(np.min((np.abs(a[0]-b[0]), 1 - np.abs(a[0]-b[0]))) ** 2 + np.min((np.abs(a[1]-b[1]), 1 - np.abs(a[1]-b[1]))) ** 2) 
    return distance


def initialize_route(M):
    route = [i for i in range(M)]
    return route


def nnsolution(M, distances):
    distances = np.array(distances)
    route = np.zeros(M,dtype=np.int8)
    counter = 0
    current_index = 0
    while counter < M - 1:
        counter += 1
        route[counter] = int(np.argmin(distances[current_index]))
        distances[:, current_index] = 10
        current_index = route[counter]
    return route


def exchange(M, prev_route, a, b): # a: cut between a and a+1. b cut between b and b+1
    new_route = np.zeros(M,dtype=np.int8)
    new_route[0:a+1] = prev_route[0:a+1]
    new_route[a+1:b+1] = np.flip(prev_route[a+1:b+1])
    new_route[b+1:] = prev_route[b+1:]
    return new_route


def length(M, distances, route):
    total_length = 0
    for i in range(M-1):
        total_length += distances[route[i],route[i+1]]
    total_length += distances[0,route[-1]]
    return total_length


def MC():
    pass


def delta_E(D, a, b):
    dE = D[a,b] + D[a+1,b+1] 
    dE -= D[a,a+1] + D[b,b+1] 
    return dE


def annealing():
    pass


if __name__ == "__main__":
    # M = 200
    # Create and save new city
    # city, dist = initialize_system(M)
    # with open('city_200_1.pkl', 'wb') as file:
    #     pickle.dump([city, dist], file)
    with open('city_200_1.pkl', 'rb') as file:
        [city, dist] = pickle.load(file)
    M = len(city)
    route_1 = initialize_route(M)
    print(length(M, dist, route_1))
    # route_nn = nnsolution(M, dist)
    # print(length(M, dist, route_nn))
    # print(exchange(M, initialize_route(M), 0, 5))