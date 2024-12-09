import numpy as np
from numpy.random import rand
from random import randrange
from matplotlib import pyplot as plt
import pickle
from copy import deepcopy
from itertools import permutations


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


def optimal_solution(M, D):
    routes = [i for i in range(M)]
    routes = list(permutations(routes))
    distances = [length(M, D, route) for route in routes]
    optimal = np.min(distances)
    return optimal



def metric(a, b):
    distance = np.sqrt(np.min((np.abs(a[0]-b[0]), 1 - np.abs(a[0]-b[0]))) ** 2 + np.min((np.abs(a[1]-b[1]), 1 - np.abs(a[1]-b[1]))) ** 2) 
    return distance


def initialize_route(M):
    route = [i for i in range(M)]
    return route


def nnsolution(M, distances):
    distances = np.array(distances)
    route = np.zeros(M,dtype=np.int16)
    counter = 0
    current_index = 0
    while counter < M - 1:
        counter += 1
        route[counter] = int(np.argmin(distances[current_index]))
        distances[:, current_index] = 10
        current_index = route[counter]
    return route


def exchange(M, prev_route, a, b): # a: cut between a and a+1. b cut between b and b+1
    new_route = np.zeros(M,dtype=np.int16)
    if a < b:
        new_route[0:a+1] = prev_route[0:a+1]
        new_route[a+1:b+1] = np.flip(prev_route[a+1:b+1])
        new_route[b+1:] = prev_route[b+1:]
    else: 
        new_route[0:b+1] = prev_route[0:b+1]
        new_route[b+1:a+1] = np.flip(prev_route[b+1:a+1])
        new_route[a+1:] = prev_route[a+1:]
    return new_route


def length(M, distances, route):
    total_length = 0
    for i in range(M - 1):
        total_length += distances[route[i], route[i+1]]
    total_length += distances[route[0], route[-1]]
    return total_length


def MC(N, M, T, distances, route):
    distance = [length(M, distances, route)]
    for i in range(N-1):
        a = randrange(M)
        b = randrange(M)
        while a == b:
            b = randrange(M)
        ksi = rand()
        d_E = delta_E(M, distances, route, a, b)
        d_e = np.exp(-d_E / T)
        if ksi < np.min((d_e, 1)):
            route = exchange(M, route, a, b)
            distance.append(distance[-1] + d_E)
        else:
            distance.append(distance[-1])
    return route, distance


def delta_E(M, D, route, a, b):
    adist = route[a]
    if a == M - 1:
        adist2 = route[0]
    else:
        adist2 = route[a+1]

    bdist = route[b]
    if b == M - 1:
        bdist2 = route[0]
    else:
        bdist2 = route[b+1]

    dE = D[adist,bdist] + D[adist2,bdist2] # new
    dE -= D[adist,adist2] + D[bdist,bdist2] # old -> longer -> positive
    return dE


def annealing(M, D, route, T_tot, T_step, T_0, r):
    # M: #of cities
    # D: Matrix of distances
    # route
    # T_tot: Total amount of MCMC steps
    # T_step: Amount of steps of MC before lowering T
    # T_0: Initial T
    # r: Cooling rate, r < 1
    dist_list = []
    T_protocol = [T_0 * r ** i for i in range(int(T_tot / T_step))]
    for T in T_protocol:
        route, tot_dist = MC(T_step, M, T, D, route)
        dist_list.extend(tot_dist)
    return route, dist_list


def tempering(M, D, route, N_tot, N_step, T_max, T_min, T_n):
    # N_tot: Total number of MC steps
    # N_step: number of steps before exchanging temps
    # T_max: highest T
    # T_min: lowest temp
    # T_n: number of simulations at different T
    T_step = T_max / T_n
    T_list = [T_min + T_step * i for i in range(T_n)]
    route_list = [deepcopy(route) for i in range(T_n)]
    dist_matrix = [[length(M, D, route)] for route in route_list]
    for j in range(int(N_tot / N_step)):
        for i in range(len(route_list)):
            route_list[i], tot_dist = MC(N_step, M, T_list[i], D, route_list[i])
            dist_matrix[i].extend(tot_dist)
        # Change places 
        T_rand = randrange(T_n - 2) # Randomly pick a Temperature between T_0 to T_max - 1
        ksi = rand()
        beta = + 1 / T_list[T_rand] - 1 / T_list[T_rand + 1]
        dE =  + length(M, D, route_list[T_rand]) - length(M, D, route_list[T_rand + 1])
        check = np.exp(beta * dE)
        if ksi < np.min((1, check)):
            new_route1 = deepcopy(route_list[T_rand])
            new_route2 = deepcopy(route_list[T_rand + 1])
            route_list[T_rand] = new_route2
            route_list[T_rand + 1] = new_route1
    len_list = [length(M, D, route_list[i]) for i in range(T_n)]
    return route_list, dist_matrix
    

if __name__ == "__main__":
    M = 200
    #N_tot = 10000
    city, dist= initialize_system(M)
    # opt = optimal_solution(M, D)
    # random = initialize_route(M)
    # nn = nnsolution(M, D)
    # print(length(M, D, random), length(M, D, nn), opt)
    # Create and save new city
    with open('city_200_4.pkl', 'wb') as file:
        pickle.dump([city, dist], file)
    # with open('city_200_1.pkl', 'rb') as file:
    #     [city, dist] = pickle.load(file)
    # M = len(city)
    # route_1 = initialize_route(M)
    # route_2 = nnsolution(M, dist)
    # route_annealing, tot_dist = annealing(M, dist, route_1, N_tot, 1000, 0.1, 0.5)
    # route_3, tot_dist = MC(1000, M, 0.1, dist, route_1)
    # plt.plot(range(N_tot), tot_dist)
    # plt.show()
    # temp_route, dist_matrix = tempering(M, dist, route_2, 1000000, 50000, 0.01, 0.0005, 8)
    # print(length(M, dist, route_1))
    # print(length(M, dist, route_2))
    #print(length(M, dist, route_3))
    #print(length(M, dist, route_annealing))
    # print(temp_route)
    # route_nn = nnsolution(M, dist)
    # print(length(M, dist, route_nn))
    
    # for i in range(8):
    #     plt.plot(range(1000001), dist_matrix[i][:], label=f"v_{i+1}")

    # # Add labels and legend
    # plt.xlabel("Iteration step")
    # plt.ylabel("Total length")
    # plt.title("Parallel tempering")
    # plt.legend()
    # plt.show()