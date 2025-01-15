# Data-to-Persistence-Module
## Copyright 

Copyright (c) 2024 Mehdi Nategh  
This code is licensed under the MIT License. See the LICENSE file for details.


## Overview

This repository contains a Python module for performing persistent homology computations using the GUDHI library. The code implements several functions to analyze the topological features of point cloud data through the construction of Rips complexes and biRips complexes. It also provides mappings between homology groups of different dimensions.

## Requirements

To run this code, you need to have the following Python libraries installed:

- `gudhi`
- `numpy`
- `scikit-learn`

You can install these dependencies using pip:

```bash
pip install gudhi numpy matplotlib scikit-learn
```

## Functions

### `homology_generators(data, r)`

Calculates the generators of the 1-homology group for a given Rips complex defined by the input data and maximum edge length `r`.

**Parameters:**

- `data`: A list of points (2D coordinates).
- `r`: The maximum edge length for the Rips complex.

**Returns:**

- A list of homology generators. Every generator is a list of at least 3 edges.

### `biRips_a(data, a, p)`

Identifies points in the dataset that meet certain density conditions defined by parameters `a` and `p`.

**Parameters:**

- `data`: A list of points.
- `a`: A lower threshold.
- `p`: An upper threshold.

**Returns:**

- A list of 2D points. The belongs to the list if it has at least a and at most p datapoints in its neighborhood. 

### `biRips(data, p)`

Constructs a biRips complex and computes the homology generators for various parameters.

**Parameters:**

- `data`: A list of points.
- `p`: An upper distance threshold.

**Returns:**

- A dictionary with keys as tuples `(a, r)` and values as corresponding homology generators of the dataset `biRips(data, a, p)`.

### `vertical_homology_mapping(data, a, p, r)`

Generates a mapping between the 1-dimensional homology groups of two consecutive biRips complexes along the vertical direction (i., `r` direction)

**Parameters:**

- `data`: A list of points.
- `a`: The parameter for the current Rips complex.
- `p`: An upper distance threshold.
- `r`: The dimension of the Rips complex.

**Returns:**

- A matrix representing a linear transformation from `biRips(data, p)[(a, r)]` to `biRips(data, p)[(a, r+ 1)]`

### `horizontal_homology_mapping(data, a, p, r)`

Generates a mapping between the 1-dimensional homology groups of two consecutive biRips complexes along the horizontal direction (i.e., `a` direction).

**Parameters:**

- `data`: A list of points.
- `a`: The parameter for the first Rips complex.
- `p`: An upper distance threshold.
- `r`: The dimension of the Rips complex.

**Returns:**

- A matrix representing from `biRips(data, p)[(a, r)]` to `biRips(data, p)[(a + 1, r)]`.

### `data_to_pModule(data, p)`

Creates a data structure that represents the persistent module for the input data.

**Parameters:**

- `data`: A list of points.
- `p`: An upper distance threshold.

**Returns:**

- A dictionary representing the persistent module `pModule`. 
- `pModule` is a dictionary where keys are d-dimensional tuples representing indices (grades) 
  `z = (z_1, z_2, ..., z_d)`. 
- Values are lists containing:
  - The identity matrix of dimension `dim(M_z)`.
  - A list of linear transformations:
    - `M_z -> M_{z + (1, 0, 0, ..., 0)}`
    - `M_z -> M_{z + (0, 1, 0, ..., 0)}`
    - ...
    - `M_z -> M_{z + (0, 0, ..., 0, 1)}`


## Example Usage

Here’s a simple example of how to use the provided functions:

```python
import numpy as np

# Sample data
data = np.array([[1, 0], [0.7, 0.7], [0, 1], [-0.7, 0.7], 
                 [-0.7, -0.7], [0.7, -0.7], [5, 1], 
                 [4.7, 4.7], [4, 5], [3.3, 4.7], 
                 [3.3, 3.3], [4.7, 3.3]])

# Compute homology generators
generators = homology_generators(data, r=1)

`

## Acknowledgments

This project utilizes the GUDHI library for topological data analysis. For more information, please visit https://gudhi.inria.fr/.
