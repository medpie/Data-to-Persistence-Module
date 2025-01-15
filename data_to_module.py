import gudhi
import numpy as np
from sklearn.metrics import pairwise_distances
print("#####################################################################")

def homology_generators(data, r):
     generators = []
     
     def loop_checker(array):
        x, y = zip(*array)
        list_of_endpoints = list(x) + list(y)
        set_of_vertices = set(list_of_endpoints)
        array_sum = sum(list_of_endpoints)
        set_sum = sum(set(list_of_endpoints))
        if array_sum == 2*set_sum:
             return True
        else:
             False

     rips = gudhi.RipsComplex(points=data, max_edge_length=r)
     simplex_tree = rips.create_simplex_tree(max_dimension=2)
     persistence = simplex_tree.persistence(homology_coeff_field=2)
     one_skeleton_total = list(simplex_tree.get_skeleton(1))
     one_skeleton = [pair for pair in one_skeleton_total if len(pair[0]) == 2]
     filtration = list(set(pair[1] for pair in one_skeleton))
     simplices_dic = {}
     for parameter in filtration:
          simplices_par = []
          for pair in one_skeleton:
               if pair[1] == parameter:
                    simplices_par.append(pair[0])
                    simplices_dic[parameter] = simplices_par
     for parameter in filtration:
          if loop_checker(simplices_dic[parameter]) == True:
                   list_of_simplices = []
                   for simplex in simplices_dic[parameter]:
                       list_of_simplices.append([data[simplex[0]], data[simplex[1]]])
                   simplices_dic[parameter] = list_of_simplices
                   generators.append(simplices_dic[parameter])

     return generators


def biRips_a(data, a, p):
    def density_function(data, p):
      gamma_dictionary = {}
      distance_matrix = pairwise_distances(data)
      for i, point in enumerate(data):
        point_row = distance_matrix[i]
        gamma_dictionary[tuple(point)] = np.count_nonzero(point_row <= p)
      return gamma_dictionary
    
    observed_datapoints = set()
    gamma_function = density_function(data,p)
    for v in gamma_function.values():
        if a <= v:
          for data_point in gamma_function.keys():
              if gamma_function[data_point] == v and data_point not in observed_datapoints:
                observed_datapoints.add(data_point)

    return list(observed_datapoints)


def biRips(data, p):
  N = len(data) + 1
  distance_matrix = pairwise_distances(data)
  r_max = int(np.max(distance_matrix)) + 1
  biRips_dictionary = {}

  for a in range(N):
    Rips = biRips_a(data, a, p)  

    if len(Rips) > 0:  
      for r in range(r_max):
        index_a_r = homology_generators(Rips, r)
        biRips_dictionary.update({(a, r): index_a_r}) 

  return biRips_dictionary

#Genrate a linear transformation H_1(Rips(a, r)) -> H_1(Rips(a, r+1))
def vertical_homology_mapping(data, a, p, r):
    points_1 = biRips(data, p)[(a, r)]
    points_2 = biRips(data, p)[(a, r + 1)]
    if len(points_2) == 0 and len(points_1) != 0:
       mapping = np.zeros((1, len(points_1)))
    elif len(points_1) == 0 and len(points_2) != 0:
       mapping = np.zeros((len(points_2), 1))
    elif len(points_1) == 0 and len(points_2) == 0:
      mapping = np.zeros((1, 1))
    else:
      mapping = np.zeros((len(points_2), len(points_1)))
      for j, sub_list in enumerate(points_1):
         if sub_list in points_2:
            i = points_2.index(sub_list)
            mapping[i, j] = 1

    return mapping

#Genrate a linear transformation H_1(Rips(a, r)) -> H_1(Rips(a+1, r))
def horizontal_homology_mapping(data, a, p, r):
    points_1 = biRips(data, p)[(a, r)]
    points_2 = biRips(data, p)[(a + 1, r)]

    if len(points_2) == 0 and len(points_1) != 0:
      mapping = np.zeros((1, len(points_1)))
    elif len(points_1) == 0 and len(points_2) !=0:
      mapping = np.zeros((len(points_2), 1))
    elif len(points_1) == 0 and len(points_2) ==0:
       mapping = np.zeros((1,1))
    else:          
      mapping = np.zeros((len(points_2), len(points_1)))
      for j, sub_list in enumerate(points_1):
         if sub_list in points_2:
            i = points_2.index(sub_list)
            mapping[i, j] = 1

    return mapping


def data_to_pModule(data, p):
    N = len(data)
    distance_matrix = pairwise_distances(data)
    r_max = int(np.max(distance_matrix))
    dictionary = {}

    biRips_dict = biRips(data, p)  # Call once and reuse

    for a in range(N):
        for r in range(r_max):
            if (a, r) in biRips_dict:
                dim_a = len(biRips_dict[(a, r)])
            else:
                dim_a = 0

            map_1 = horizontal_homology_mapping(data, a, p, r)
            map_2 = vertical_homology_mapping(data, a, p, r)

            if dim_a != 0:
                identity_a_r = np.eye(dim_a)
                dictionary.update({(a, r): [identity_a_r, [map_1, map_2]]})
            else:
                dictionary.update({(a, r): [np.zeros((1, 1)), [map_1, map_2]]})

    return dictionary

data = [[1, 0], [0.7, 0.7], [0, 1], [-0.7, 0.7], [-0.7, -0.7], [0.7, -0.7], [5, 1], [4.7, 4.7], [4, 5], [3.3, 4.7], [3.3, 3.3], [4.7, 3.3]]

#print(data_to_pModule(data, 6))
