import numpy as np
from numpy.random import rand
from random import randrange
import matplotlib
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
    print(routes)
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


def MC(N, M, T, distances, route, k=1):
    distance = [length(M, distances, route)]
    for i in range(N-1):
        a = randrange(M)
        b = randrange(M)
        while a == b:
            b = randrange(M)
        ksi = rand()
        d_E = delta_E(M, distances, route, a, b)
        d_e = np.exp(-d_E / (k * T))
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


def tempering(M, D, route, N_tot, N_step, T_max, T_min, L, k=1):
    # N_tot: Total number of MC steps
    # N_step: number of steps before exchanging temps
    # T_max: highest T
    # T_min: lowest temp
    # L: number of simulations at different T
    T_step = (T_max - T_min) / (L - 1)
    T_list = [T_min + T_step * i for i in range(L)]
    route_list = [deepcopy(route) for i in range(L)]
    dist_matrix = [[length(M, D, route)] for route in route_list]
    for j in range(int(N_tot / N_step)):
        for i in range(len(route_list)):
            route_list[i], tot_dist = MC(N_step, M, T_list[i], D, route_list[i], k)
            dist_matrix[i].extend(tot_dist)
        # Change places 
        T_rand = randrange(0, L - 1) # Randomly pick a Temperature between T_0 to T_max - 1
        # while dist_matrix[T_rand,-1] == np.max(dist_matrix[:,-1]):
        #     T_rand = randrange(0, L - 1)
        diff = np.zeros(L)
        for i in range(L):
            diff[i] = dist_matrix[T_rand][-1] - dist_matrix[i][-1]
            if diff[i] == 0:
                diff[i] = 10 ** 6
        T_close = np.argmin([np.abs(element) for element in diff])
        ksi = rand()
        beta = + 1 /( T_list[T_rand] * k) - 1 /( T_list[T_close] * k)
        dE =  length(M, D, route_list[T_rand]) - length(M, D, route_list[T_close])
        check = np.exp(beta * dE)
        print(beta, dE)
        print(check, T_rand, T_close, j)
        if ksi < np.min((1, check)):
            print(j)
            new_route1 = deepcopy(route_list[T_rand])
            new_route2 = deepcopy(route_list[T_close])
            route_list[T_rand] = new_route2
            route_list[T_close] = new_route1
    len_list = [length(M, D, route_list[i]) for i in range(L)]
    return route_list, dist_matrix
    

if __name__ == "__main__":
    # M = 10
    # N_tot = 1000000
    # N_r = 40000
    # r = 0.95
    # T_0 = 0.1
    # city, dist= initialize_system(M)
    # random = initialize_route(M)
    # nn = nnsolution(M, D)
    # print(length(M, D, random), length(M, D, nn), opt)
    # Create and save new city
    # with open('city_10_1.pkl', 'wb') as file:
    #     pickle.dump([city, dist], file)


    # with open('city_200_1.pkl', 'rb') as file:
    #     [city, dist] = pickle.load(file)
    # M = len(city)
    # route_2 = nnsolution(M, dist)

    # OPTIMAL VS ANNEALING
    # route_1 = initialize_route(M)
    # route_annealing, tot_dist = annealing(M, dist, route_1, N_tot, N_r, T_0, r)

    # plt.plot(range(N_tot),tot_dist, label=f"Simulated annealing, d = {np.round(np.min(tot_dist), 4)}")
    # plt.plot(range(N_tot), [opt for i in range(N_tot)], '--', label=f"Optimal solution, d = {np.round(opt, 4)}")
    # plt.legend()
    # plt.xlabel("MC Iteration step")
    # plt.ylabel("Total length")
    # plt.show()
    # plt.savefig(f"Optimal vs annealing M=20, d = {np.round(np.min(tot_dist), 4)}.pdf")

    # OPTIMAL VS TEMPERING
    # route_1 = initialize_route(M)
    # opt = 2.2737415617330017
    # T_max = 0.1
    # L = 3
    # T_min = 0.001
    # temp_route, dist_matrix = tempering(M, dist, route_1, N_tot, N_r, T_max, T_min, L)
    # T_step = (T_max - T_min) / (L - 1)
    # T_list = [T_min + T_step * i for i in range(L)]
    # for i in range(L):
    #     plt.plot(range(N_tot +1), dist_matrix[i][:], label=f"R(T = {T_list[i]}), d = {np.round(np.min(dist_matrix[i][-1]), 4)}")
    # plt.plot(range(N_tot +1 ), [opt for i in range(N_tot + 1)], '--', label=f"Optimal solution, d = {np.round(opt, 4)}")
    # plt.legend()
    # plt.xlabel("MC Iteration step")
    # plt.ylabel("Total length")
    # plt.show()
    # plt.savefig(f"Optimal vs tempering M=20, N_tot={N_tot}, N_r={N_r}, (T_max,L,T_min) = {T_max, L, T_min}).pdf")

    # ANNEALING 200
    # with open('city_200_1.pkl', 'rb') as file:
    #     [city, dist] = pickle.load(file)
    # M = len(city)

    # N_tot = 1000000
    # N_r = 10000
    # r = 0.98

    # basic_route = nnsolution(M, dist)
    # T_0 = np.sqrt(length(M, dist, basic_route)) / M
    # route, tot_dist = annealing(M, dist, basic_route, N_tot, N_r, T_0, r)

    # plt.plot(range(N_tot),tot_dist, label=f"Simulated annealing, d = {np.round(np.min(tot_dist), 4)}")
    # plt.legend()
    # plt.xlabel("MC Iteration step")
    # plt.ylabel("Total length")
    # # plt.savefig(f"Annealing_M=200_r={r}_city2_d={np.round(np.min(tot_dist), 4)}.pdf", format="pdf")
    # plt.show()

    # GRID SEARCH ANNEALING
    # with open('city_200_1.pkl', 'rb') as file:
    #     [city, dist] = pickle.load(file)
    # M = len(city)
    # basic_route = nnsolution(M, dist)
    # T_0 = np.sqrt(length(M, dist, basic_route)) / M

    # N_sim = 1
    # N_tot = 1000000
    # r_list = [0.5, 0.8, 0.9, 0.95, 0.99]
    # N_r_list = [int(N_tot/5), int(N_tot/10), int(N_tot/20), int(N_tot/50), int(N_tot/100)]
    # means = np.zeros((len(r_list), len(N_r_list)))
    # stds = np.zeros((len(r_list), len(N_r_list)))
    # for rid,r in enumerate(r_list):
    #     for N_rid, N_r in enumerate(N_r_list):
    #         lngth = []
    #         for i in range(N_sim):
    #             print(rid, N_rid, i)
    #             route, tot_dist = annealing(M, dist, basic_route, N_tot, N_r, T_0, r)
    #             lngth.append(np.min(tot_dist))
    #         means[rid,N_rid] = np.mean(lngth)
    #         stds[rid,N_rid] = np.std(lngth)
    # print(means)
    # print(stds)
    # r_list = np.array(r_list)
    # N_r_list = np.array(N_r_list)
    # # Define x and y coordinates for the grid
    # r_listedge = np.zeros(len(r_list) + 1)  # One more edge than center points
    # r_listedge[1:-1] = (r_list[:-1] + r_list[1:]) / 2  # Midpoints between consecutive centers
    # r_listedge[0] = r_list[0] - (r_list[1] - r_list[0]) / 2  # Extrapolate for the first edge
    # r_listedge[-1] = r_list[-1] + (r_list[-1] - r_list[-2]) / 2  # Extrapolate for the last edge
    # N_r_listedge = np.zeros(len(N_r_list) + 1)  # One more edge than center points
    # N_r_listedge[1:-1] = (N_r_list[:-1] + N_r_list[1:]) / 2  # Midpoints between consecutive centers
    # N_r_listedge[0] = N_r_list[0] - (N_r_list[1] - N_r_list[0]) / 2  # Extrapolate for the first edge
    # N_r_listedge[-1] = N_r_list[-1] + (N_r_list[-1] - N_r_list[-2]) / 2  # Extrapolate for the last edge
    # X, Y = np.meshgrid(r_listedge, N_r_listedge)

    # # Create a contour plot
    # plt.figure(figsize=(6, 5))
    # contour = plt.pcolormesh(X, Y, means, cmap="viridis_r", shading='flat')  # Filled contour plot

    # # Add a color bar to indicate data values
    # plt.colorbar(contour, label="Value")

    # # Add labels and a title
    # plt.xlabel("r")
    # plt.ylabel("N_r")
    # plt.title(f"Contour Plot of distances")
    # plt.show()

    # # Define x and y coordinates for the grid
    # r_listedge = np.zeros(len(r_list) + 1)  # One more edge than center points
    # r_listedge[1:-1] = (r_list[:-1] + r_list[1:]) / 2  # Midpoints between consecutive centers
    # r_listedge[0] = r_list[0] - (r_list[1] - r_list[0]) / 2  # Extrapolate for the first edge
    # r_listedge[-1] = r_list[-1] + (r_list[-1] - r_list[-2]) / 2  # Extrapolate for the last edge
    # N_r_listedge = np.zeros(len(N_r_list) + 1)  # One more edge than center points
    # N_r_listedge[1:-1] = (N_r_list[:-1] + N_r_list[1:]) / 2  # Midpoints between consecutive centers
    # N_r_listedge[0] = N_r_list[0] - (N_r_list[1] - N_r_list[0]) / 2  # Extrapolate for the first edge
    # N_r_listedge[-1] = N_r_list[-1] + (N_r_list[-1] - N_r_list[-2]) / 2  # Extrapolate for the last edge
    # X, Y = np.meshgrid(r_listedge, N_r_listedge)

    # # Create a contour plot
    # plt.figure(figsize=(6, 5))
    # contour = plt.pcolormesh(X, Y, stds, cmap="viridis_r", shading='flat')  # Filled contour plot

    # # Add a color bar to indicate data values
    # plt.colorbar(contour, label="Value")

    # # Add labels and a title
    # plt.xlabel("r")
    # plt.ylabel("N_r")
    # plt.title(f"Contour Plot of distance std for {N_sim} simulations")

    # # Display the plot
    # plt.show()


    # TEMPERING 200
    with open('city_200_1.pkl', 'rb') as file:
        [city, dist] = pickle.load(file)
    M = len(city)
    N_tot = 500000
    N_r = 5000
    T_max = 0.01
    L = 10
    T_min = 0.0001

    # route = initialize_route(M)
    route = nnsolution(M, dist)
    temp_route, dist_matrix = tempering(M, dist, route, N_tot, N_r, T_max, T_min, L)

    T_step = (T_max - T_min) / (L - 1)
    T_list = [T_min + T_step * i for i in range(L)]
    for i in range(L):
        plt.plot(range(N_tot +1), dist_matrix[i][:], label=f"R(T = {np.round(T_list[i],5)})")
        # plt.plot(range(N_tot +1), dist_matrix[i][:], label=f"R(T = {np.round(T_list[i],4)}), d = {np.round(np.min(dist_matrix[i][-1]), 2)}")
    plt.legend()
    plt.xlabel("MC Iteration step")
    plt.ylabel("Total length")
    # plt.savefig(f"Tempering_M={M}_N_tot={N_tot}_N_r={N_r}_(T_max,L,T_min)={T_max, L, T_min}).pdf")
    plt.show()



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
    
    ##for i in range(tot_dist):
    #    plt.plot(range(1000001), dist_matrix[i][:], label=f"v_{i+1}")

    # Add labels and legend
    # plt.plot(range(N_tot), tot_dist, label="Simulated annealing")
    # #plt.plot(range(N_tot+1), [opt for i in range(N_tot+1)], label="Exact ground state")
    # plt.xlabel("MC Iteration step")
    # plt.ylabel("Total length")
    # # plt.title("Simulated annealing vs optimal sol")
    # plt.legend()
    # plt.show()
    # plt.savefig(f"Simulated_annealingN={N_tot},r={r},T_0={T_0}.pdf", format="pdf")