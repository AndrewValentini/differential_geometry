from sympy import *

def christoffel_symbols(metric_tensor, inverse_metric_tensor, coordinates):
    num_coords = len(coordinates)
    # Creating a rank 3 tensor(3D matrix) where all of the Christoffel symbols will be stored
    christoffel = [[[0 for _ in range(num_coords)] for _ in range(num_coords)] for _ in range(num_coords)]

    for i in range(num_coords):
        for j in range(num_coords):
            for k in range(num_coords):
                sum_term = 0.0
                for l in range(num_coords):
                    term1 = 0.5 * inverse_metric_tensor[i, l] * (
                        diff(metric_tensor[l, j], coordinates[k]) +
                        diff(metric_tensor[l, k], coordinates[j]) -
                        diff(metric_tensor[j, k], coordinates[l])
                    )
                    sum_term += term1
                christoffel[i][j][k] = sum_term

    return christoffel


def Riemann(metric_tensor, inverse_metric_tensor, coordinates):
    num_coordinates = len(coordinates)
    # Creating a 4D matrix to store the Riemann curvature tensor values in
    riemann = [[[[0 for _ in range(num_coordinates)] for _ in range(num_coordinates)] for _ in range(num_coordinates)] for _ in range(num_coordinates)]

    for i in range(num_coordinates):
        for j in range(num_coordinates):
            for k in range(num_coordinates):
                for l in range(num_coordinates):
                    sum_term = 0
                    for m in range(num_coordinates):
                        term1 = .5* inverse_metric_tensor[i,m] * (
                            diff(diff(metric_tensor[m,l], coordinates[j]), coordinates[k]) -
                            diff(diff(metric_tensor[m,k], coordinates[j]), coordinates[l]) +
                            diff(diff(metric_tensor[j,k], coordinates[m]), coordinates[l]) -
                            diff(diff(metric_tensor[j,l], coordinates[m]), coordinates[k])
                        ) 
                        sum_term += term1
                    riemann[i][j][k][l] = sum_term
    return riemann


def Ricci_tensor(metric_tensor, inverse_metric_tensor, coordinates):
    num_coordinates = len(coordinates)
    ricci = [[0 for _ in range(num_coordinates)] for _ in range(num_coordinates)]

    for i in range(num_coordinates):
        for j in range(num_coordinates):
            sum_term = 0
            for k in range(num_coordinates):
                for l in range(num_coordinates):
                    for m in range(num_coordinates):
                        term1 = .5 * inverse_metric_tensor[i, m] * (
                            diff(diff(metric_tensor[m, l], coordinates[j]), coordinates[k]) -
                            diff(diff(metric_tensor[m, k], coordinates[j]), coordinates[l]) +
                            diff(diff(metric_tensor[j, k], coordinates[m]), coordinates[l]) -
                            diff(diff(metric_tensor[j, l], coordinates[m]), coordinates[k])
                        )
                        sum_term += term1
            ricci[i][j] = sum_term

    return ricci



    