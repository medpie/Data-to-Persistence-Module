import gudhi
import numpy as np
from sklearn.metrics import pairwise_distances
print("#####################################################################")

""" List of functions:
    - homology_generators(data, p)
    - biRips_a(data, a, p)
    - biRips(data, p)
    - vertical_homology_mapping(data, a, p, r)
    - horizontal_homology_mapping(data, a, p, r)
    - data_to_pModule(data, p)
"""

def homology_generators(data, r):
     generators = []
     def loop_finder(edges):
      from collections import defaultdict

      graph = defaultdict(list)
      for a, b in edges:
          graph[a].append(b)
          graph[b].append(a)

      visited = set()
      loops = []

      def dfs(vertex, parent, path):
          visited.add(vertex)
          path.append(vertex)

          for neighbor in graph[vertex]:
              if neighbor == parent:  
                  continue
              if neighbor in path:  
                  loop_start_index = path.index(neighbor)
                  loop_edges = [(path[i], path[i + 1]) for i in range(loop_start_index, len(path) - 1)]
                  loop_edges.append((path[-1], neighbor))  
                  loops.append(loop_edges)
              elif neighbor not in visited:
                  dfs(neighbor, vertex, path)

          path.pop()  

      for vertex in graph:
          if vertex not in visited:
              dfs(vertex, None, [])

      unique_loops = []
      seen_loops = set()
      for loop in loops:
          loop_tuple = tuple(sorted(loop))
          if loop_tuple not in seen_loops and len(loop) > 2:
              seen_loops.add(loop_tuple)
              unique_loops.append(loop)

      return unique_loops

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
          loops = loop_finder(simplices_dic[parameter])
          if len(loops) != 0:
             for loop in loops:
                   generators.append([[data[simplex[0]], data[simplex[1]]] for simplex in loop])

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
        if v <= a + 1:
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
    else:
        for r in range(r_max):
           biRips_dictionary.update({(a, r): []})

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


#======================================================================================
#======================================================================================
#======================================================================================
# Example1

#data = [[1, 0], [0.7, 0.7], [0, 1], [-0.7, 0.7], [-0.7, -0.7], [0.7, -0.7], [5, 1], [4.7, 4.7], [4, 5], [3.3, 4.7], [3.3, 3.3], [4.7, 3.3]]

#print(homology_generators(data, 21))
#print(data_to_PModule(data, 7))

#======================================================================================
# Exmaple 2

import matplotlib.pyplot as plt

# This function creates a loop of a given radius and center
def create_loop(center, radius, num_points=100):
    angles = np.linspace(0, 2 * np.pi, num_points)
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)
    return np.column_stack((x, y))


small_loops = []
for i in range(3):  
    center = (np.random.uniform(-5, 5), np.random.uniform(-5, 5))
    small_loops.append(create_loop(center, radius=1))


big_loops = []
for i in range(2):  
    center = (np.random.uniform(-5, 5), np.random.uniform(-5, 5))
    big_loops.append(create_loop(center, radius=3))

data = np.vstack(small_loops + big_loops)

# Optionally, plot the loops
plt.figure(figsize=(8, 8))
plt.scatter(data[:, 0], data[:, 1], s=5)
plt.title("Point Data with Small and Big Loops")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.axis('equal')
plt.grid()
#plt.show() 

# Finding the first radius `r` for which the set of 1-homology generators is non-empty
r = 0
hom = homology_generators(data, r)
while homology_generators(data, r) == []:
   r += 1
   hom = homology_generators(data, r)
   print(f'r = {r}: pending ..')
   if hom != []:
      print(f'The 1-homology generators for r = {r}:')
      print(hom)
      break
