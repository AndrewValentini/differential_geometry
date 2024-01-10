from sympy import *

def christoffel_symbols(metric_tensor, inverse_metric_tensor, coordinates):
    num_coords = len(coordinates)
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
    riemann = [[[[0 for _ in range(num_coordinates)] for _ in range(num_coordinates)] for _ in range(num_coordinates)] for _ in range(num_coordinates)]

    for i in range(num_coordinates):
        for j in range(num_coordinates):
            for k in range(num_coordinates):
                for l in range(num_coordinates):
                    sum_term = 0
                    for m in range(num_coordinates):
                        term1 = .5 * inverse_metric_tensor[i, m] * (
                            diff(diff(metric_tensor[m, l], coordinates[j]), coordinates[k]) -
                            diff(diff(metric_tensor[m, k], coordinates[j]), coordinates[l]) +
                            diff(diff(metric_tensor[j, k], coordinates[m]), coordinates[l]) -
                            diff(diff(metric_tensor[j, l], coordinates[m]), coordinates[k])
                        ) 
                        sum_term += term1
                    riemann[i][j][k][l] = sum_term
    return riemann
    
#Here, the Riemann tensor is only defined in terms of the Christoffel symbols since they seem to be computed correctly:

def Riemann_from_christoffel(christoffel, coordinates):
    num_coordinates = len(coordinates)
    riemann = [[[[0 for _ in range(num_coordinates)] for _ in range(num_coordinates)] for _ in range(num_coordinates)] for _ in range(num_coordinates)]

    for i in range(num_coordinates):
        for j in range(num_coordinates):
            for k in range(num_coordinates):
                for l in range(num_coordinates):
                    term1 = diff(christoffel[i][l][k], coordinates[j]) - diff(christoffel[i][k][j], coordinates[l])
                    term2 = sum(christoffel[m][k][j] * christoffel[i][m][l] for m in range(num_coordinates))
                    term3 = sum(christoffel[m][l][k] * christoffel[i][m][j] for m in range(num_coordinates))
                    riemann[i][j][k][l] = term1 + term2 - term3

    return riemann
    

def Ricci_tensor(metric_tensor, inverse_metric_tensor, coordinates):
    num_coordinates = len(coordinates)
    ricci_tensor = [[0 for _ in range(num_coordinates)] for _ in range(num_coordinates)]

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
            ricci_tensor[i,j] = sum_term

    return ricci_tensor

def Ricci_scalar(metric_tensor, inverse_metric_tensor, coordinates):
    num_coordinates = len(coordinates)
    ricci_scalar = 0
    
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
            ricci_scalar += sum_term

    return ricci_scalar

def Einstein_tensor(metric_tensor, inverse_metric_tensor, coordinates):
    num_coordinates = len(coordinates)
    
    # Compute Ricci tensor
    ricci_tensor = Ricci_tensor(metric_tensor, inverse_metric_tensor, coordinates)
    
    # Compute Ricci scalar
    ricci_scalar = Ricci_scalar(metric_tensor, inverse_metric_tensor, coordinates)

    # Compute Einstein tensor using the formula: G_{mu nu} = R_{mu nu} - (1/2) R g_{mu nu}
    einstein_tensor = [[0 for _ in range(num_coordinates)] for _ in range(num_coordinates)]
    for mu in range(num_coordinates):
        for nu in range(num_coordinates):
            einstein_tensor[mu,nu] = ricci_tensor[mu,nu] - (1/2) * ricci_scalar * metric_tensor[mu, nu]

    return einstein_tensor











    



    