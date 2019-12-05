"""
This is the library to store existed (blocked)MDS matrix
"""

from sympy import symbols, Matrix
import string
from sympy.matrices import randMatrix
import random
# import scratch_2
import time
from itertools import combinations
z = symbols('z')


# original 8*4 k = 2, provided by zyh
# G1 = Matrix(([[z, 1, 1, z],[z, z, 1, z+1],[z, 0, 1, z],[z, z+1, 0, 1],[z, z+1, 0, 0],[1, z+1, z, z+1], [z+1, 1, 1, 0], [0, 0, z+1, 1]]))
G1_zyh = Matrix(([[z+1, 1, 1, z],[1, 1, z+1, 0],[z, z, z+1, z+1],[z+1, 0, z, 1],[z,z,z+1,1],[0,1,z,0],[1,1,0,z],[0,0,0,z]]))

# new 8*4 k = 2 renewed!
G8_4 = Matrix(([z + 1, 1, 1, z], [z + 1, 0, 1, 0], [z, z + 1, 0, 1], [z + 1, 1, z + 1, 1], [1, 1, z + 1, z + 1],
             [0, z + 1, z + 1, z], [1, z + 1, 0, 0], [0, 1, 0, z + 1]))

# 10*6 k = 3 block_mds
G10_6_block = Matrix([[0, z, 1, z, 1, z], [1, 0, 1, 0, 1, 0], [0, z, 1, z, 1, 0], [z + 1, z, 1, z, z, z], [z, z, z + 1, 1, 1, z], [z + 1, 0, z + 1, z + 1, z, z], [z, z + 1, 0, 1, z + 1, z], [0, z, 0, z + 1, 1, 0], [1, z, 1, z, z, 1], [z, 1, 0, z, z + 1, 1]])

# 12*8 k = 4 block_mds
G12_8_block = Matrix([[1, z, 0, z, 1, z + 1, z, z + 1], [0, z + 1, 1, z, z, z, 1, 0], [1, z, 0, z, z, 0, 1, 1], [1, 1, z + 1, 1, z + 1, z, 0, z], [0, 0, 0, z, z + 1, z, z + 1, z + 1], [z + 1, 0, z, z + 1, 0, z + 1, z + 1, 0], [0, 0, 0, 1, 1, z, 0, z], [z, 1, z + 1, z + 1, 0, z + 1, 0, 1], [0, 1, 0, z, z + 1, 1, 1, z + 1], [z + 1, 0, z, z + 1, 0, z + 1, 1, 0], [z + 1, z, z, 0, z, z, 0, z], [0, 1, 1, 0, 0, z, z, 0]])


# 14*8 k = 4 block_mds
G14_8_block = Matrix(([1,0,1,z,1,z,1,z+1],[1,1,z,z,1,0,z+1,z+1],[z+1,z+1,1,1,z+1,0,z,0],[z+1,0,0,0,0,1,0,1],[0,z,0,z+1,z+1,0,z,z],[z,1,z,z,0,z+1,z+1,z+1],[1,0,0,0,0,z,1,1],[z,z+1,z,z+1,z,0,z,0],[z,0,1,0,z,z+1,z,z],[z+1,1,0,z+1,z,z,z,1],[0,1,z,1,z+1,z,z+1,z],[z,0,z+1,z+1,z,0,z+1,0],[1,1,z+1,z,z,1,z,1],[0,0,z,1,z,1,1,1]))


# 14*10 k = 5 block_mds  //
G14_10_block = Matrix([[1, 0, 1, z, 1, z, z + 1, z, 0, 0], [0, 1, z, z, z, z + 1, z + 1, 1, 0, z + 1], [z, 0, 1, z, z, 0, 1, 0, z + 1, 1], [0, z, z + 1, 0, z + 1, z, z + 1, 1, 0, z + 1], [z, z + 1, 0, z, z, 0, 0, 0, z + 1, 1], [z, z + 1, z + 1, 1, z + 1, z, 1, z + 1, 0, 1], [z, z + 1, z + 1, z + 1, z + 1, 0, 1, z + 1, z + 1, z + 1], [z + 1, z, 1, z + 1, z + 1, z, z + 1, 1, 0, 1], [1, 1, 1, z + 1, 1, z + 1, z + 1, 1, 0, 1], [z + 1, z + 1, z + 1, 0, z + 1, 0, z, z + 1, z + 1, z + 1], [0, 0, 1, 0, z + 1, z + 1, z, 0, z, 1], [z, z, z + 1, 1, z + 1, 1, z, z + 1, 0, 0], [z, 0, z, z + 1, 0, 0, 1, 1, z, z], [z + 1, 1, z + 1, z, z + 1, 1, z + 1, z + 1, 0, 1]])


# 16 * 12 k = 6 block_mds  //
G16_12_block = Matrix([[1, 0, z + 1, 0, 1, 0, z, z + 1, 1, z + 1, 0, z + 1], [z + 1, z + 1, 1, 1, z, 1, 1, z, z, 0, z, 1], [z + 1, 1, z + 1, 1, 0, 1, z, 1, 0, z, 0, 1], [z + 1, z + 1, 1, 0, 1, 1, 0, z, 0, 1, 0, 0], [1, z + 1, z + 1, z + 1, z, 1, 0, 1, z, 1, z + 1, 0], [1, 0, z + 1, z + 1, z, z, z + 1, z + 1, 0, z, 1, 1], [1, 0, z, 1, 1, 0, z + 1, z + 1, 1, z + 1, 1, 1], [z, 1, z, 0, 1, z, z, z + 1, 0, z, z + 1, z], [z, 1, 1, 0, z + 1, z, z + 1, z, z, z, 1, 1], [z, z + 1, z + 1, z, 0, 0, z, 0, z, z + 1, z + 1, z], [z + 1, z + 1, z, 1, z, 0, z + 1, 0, z, z + 1, 0, 0], [1, 1, 0, z, 1, 0, z, 0, 0, z + 1, z, 0], [z, 0, z, z + 1, z + 1, z, z + 1, 1, 1, 1, 1, 1], [z, z, z, z, z + 1, 0, 0, z + 1, z + 1, z, 0, z + 1], [z + 1, z, z, z, 0, z, 0, z, z, 1, z, 1], [1, 0, 0, z, z + 1, 1, 1, 1, 0, z + 1, 1, z]])

# 18 * 14 k = 7 block_mds  //
G18_14_block = Matrix([[z, z, 0, 0, 0, z, z + 1, 1, 0, z, z, z, z + 1, 0], [z + 1, z, z, z + 1, z, z + 1, z, z + 1, z + 1, z + 1, 1, z + 1, z, 1], [z + 1, z + 1, 1, 0, 0, 1, 0, 0, z + 1, 1, 0, 1, 1, 1], [0, 1, z, z, 1, z, 0, 0, 0, 1, 0, 1, z + 1, 1], [0, z, 1, z, z + 1, 0, z + 1, 1, z, z, 1, z + 1, z, z + 1], [0, z, 0, z + 1, 1, 0, z + 1, 0, 1, z, 1, z + 1, z, z], [z, z + 1, 0, 1, z + 1, z, z, z + 1, 1, 0, 0, z + 1, 0, z + 1], [z + 1, 0, z, 0, 1, z, z, 0, 1, 0, 1, 1, 1, 0], [z, 0, 0, z, 1, 0, z + 1, z, 0, 0, z, 0, 1, z], [z, 0, z, z + 1, z, 0, 0, 1, z, 1, z, z + 1, z + 1, 1], [z, z, z, z, 0, z + 1, 1, z + 1, z, 0, 0, z, z, z], [z + 1, z + 1, 0, z, z + 1, 0, z, z + 1, z, z, 0, z, z, z], [1, 0, 0, 1, z, 0, z, z, 0, z + 1, z, z + 1, z, z], [z + 1, 0, 1, z, z + 1, z + 1, 1, z, z, 1, 1, z + 1, z + 1, z], [z + 1, 1, z + 1, z, z + 1, z + 1, z, z + 1, 0, 1, z, 0, z + 1, 1], [1, 0, 0, 0, 1, 0, z + 1, 1, 1, z, z + 1, 0, 0, z], [z, z, z, z, 0, z, z, 0, z, 1, 1, z + 1, z, z], [z, z, z, z + 1, 0, z, 1, 1, 1, z, 0, 1, 0, z]])



# 20 * 16 k = 8 block_mds
G20_16_block = Matrix([[1, z + 1, 1, 1, z, 1, z, 1, 1, z, z + 1, 0, 0, 1, 1, z], [z, 0, 0, 0, 0, 0, 1, 0, z + 1, 1, 1, z + 1, z, 0, z + 1, 1], [z, z + 1, z + 1, 0, 1, 0, 1, z, z + 1, 0, z + 1, 0, z, z + 1, z + 1, z + 1], [1, z, z + 1, 0, 0, 0, z, 0, z, z, 1, 0, 0, 0, z, 0], [z + 1, z, z, 0, z + 1, z + 1, z + 1, 0, z, z + 1, 0, z, z + 1, z, z, 1], [0, z + 1, 0, 1, 0, 0, z + 1, z, z, 1, z, 0, 1, z, z + 1, z + 1], [0, z + 1, z + 1, z, z + 1, 0, z, z + 1, z, z, 0, 1, 0, 1, z, 0], [z, 1, 1, 1, z, 1, 1, z + 1, z, 1, z + 1, 0, 0, z + 1, z, z], [0, 1, 1, z + 1, 1, z + 1, z, z, 1, 1, z + 1, 0, z, z, 1, z], [0, z, z, 0, z, z + 1, z + 1, 1, z + 1, z + 1, z, 0, z, 0, z, z + 1], [1, z + 1, 1, 0, 0, z, 0, z + 1, 1, z, z + 1, 1, 0, 1, z + 1, z], [z, 1, z + 1, 1, 1, z, 1, z + 1, z, 1, 1, 1, z + 1, 0, z, z], [0, z, 1, z, 0, z, z + 1, z + 1, 1, z + 1, 1, 0, z, z, z, 0], [z + 1, 1, 1, 0, z, 1, 0, 1, 1, 0, z + 1, z, 1, z, z, z + 1], [1, z, z, 0, z + 1, z, 1, z, 0, 1, z, z + 1, 1, 0, z, z + 1], [1, z, 1, 0, 0, z + 1, z + 1, z, z + 1, 1, z, 1, z + 1, 1, z + 1, z + 1], [z, 1, 1, 0, 0, z, 0, z + 1, z + 1, 0, 1, z, z + 1, z + 1, z, z + 1], [0, 0, 1, z, z + 1, z + 1, z + 1, z, z, 1, z + 1, z + 1, z, 1, z + 1, z + 1], [z, z + 1, z, z + 1, 1, 0, 1, 0, z, 1, z + 1, z + 1, 0, z, 0, z + 1], [1, 0, 1, z + 1, z + 1, z + 1, z + 1, z, z, z + 1, z, 1, z, z, z + 1, 0]])



# 15 * 6 k = 2 block_mds
G15_6_block = Matrix([[0, 1, 1, z + 1, z, 0], [1, 1, z, 1, z, z + 1], [0, z, 0, z + 1, z + 1, 1], [z + 1, 0, 1, 0, z + 1, z], [z, z + 1, 1, z + 1, z + 1, z], [1, 1, 0, z + 1, 1, 0], [z + 1, 1, z + 1, z + 1, z + 1, z + 1], [z + 1, 1, 1, z, 0, z], [z, 1, z + 1, z, z, z], [z + 1, z, 0, z, 0, 0], [0, z, 1, 1, z, z], [1, 1, z + 1, z, 1, z + 1], [z, 0, 0, z, z, 1], [z, z, z, z + 1, z + 1, z + 1], [0, z + 1, z, z + 1, 0, z + 1]])


# 18 * 9 k = 3 block_mds
G18_9_block = Matrix([[0, z, z + 1, z, 1, z, 1, z + 1, 0], [z + 1, 1, 0, z + 1, 1, z + 1, z, 1, z], [z, 1, 0, 0, z, z + 1, 1, 0, z + 1], [z, 1, 1, 0, 1, z, 0, 0, 0], [1, z, 1, 0, z + 1, 0, 0, 0, z + 1], [1, 1, 0, 1, 0, 0, z, 0, 1], [1, z, 0, z, 0, z + 1, z, z + 1, z + 1], [1, z + 1, 0, z, 1, z, 0, 0, 1], [z, z + 1, z + 1, 1, 0, z + 1, z, z + 1, z + 1], [z, 1, 1, 0, z, z, z + 1, z, z], [0, 1, z + 1, z + 1, 1, 0, z, 0, z], [z + 1, 1, 0, 0, 0, z + 1, z + 1, z + 1, z], [0, 0, z + 1, z, 1, z, 0, 1, 0], [0, z + 1, 1, z, 1, z, 1, 1, z], [0, 1, z + 1, z + 1, z + 1, z, 1, z + 1, 0], [z, 0, 0, z, z, 1, z + 1, z, z], [1, 1, z, 1, z, z + 1, z, z, 0], [z + 1, 0, z + 1, z + 1, z, z + 1, z + 1, z + 1, z + 1]])

# 21 * 12 k = 4 block_mds
G21_12_block = Matrix([[0, z, z + 1, z + 1, 1, z + 1, z + 1, 0, 1, z, 1, z + 1], [z + 1, z + 1, 0, 1, z + 1, 0, z, 1, 1, z, z + 1, z], [0, z, z, 0, 0, z, 0, 1, z + 1, 0, 0, z], [z, 0, 0, 0, 1, z, z, z, z + 1, z + 1, z, 1], [z, 0, z, z + 1, z, z + 1, z, 0, z + 1, 0, z + 1, z], [z + 1, z + 1, z, 0, 0, z, z, z, z + 1, 1, z + 1, z + 1], [z + 1, z, z + 1, 1, z, z + 1, z, z + 1, z, z, 1, 1], [0, z + 1, z, z, 1, 0, 0, z + 1, z, z, z, 0], [z + 1, 0, 0, z + 1, 0, z, 1, z + 1, z + 1, z + 1, 0, z], [0, z, z, z, z + 1, 1, z + 1, 0, z + 1, 0, z, z], [0, z + 1, z + 1, z + 1, 0, 0, z + 1, 1, z + 1, 1, z + 1, 1], [z, z + 1, 0, 0, 0, z + 1, z, 1, z + 1, z + 1, z, 1], [z, z, 0, 0, 1, 1, z + 1, z, z + 1, 0, 0, 0], [z, z, z + 1, 1, 0, z + 1, 1, 1, z, z + 1, 0, z], [z + 1, z + 1, z, 0, 1, 0, 0, z, z + 1, 0, z, z + 1], [z + 1, 0, z, z + 1, z, 0, 1, 1, z + 1, 0, 1, 0], [z + 1, z, z + 1, z, z + 1, z, 0, 0, z + 1, z + 1, z, z], [0, z, 1, 0, 1, z + 1, z + 1, 0, 0, 0, 1, 0], [z + 1, z, z, z, z, z + 1, 0, 1, z + 1, 1, z, 1], [1, 1, 1, z, z, z + 1, 0, z, z, 1, 1, 0], [z + 1, z, 1, z + 1, 1, 0, 1, 0, 0, 0, 0, 0]])

# 24 * 15 k = 5 block_mds  //
G24_15_block = Matrix([[1, 0, 0, 0, z, z + 1, 1, z + 1, 0, z, 0, 1, 0, 0, 1], [z, 0, 0, z, z, 0, 0, 0, 1, z + 1, z + 1, 0, 0, 1, 1], [1, 0, 1, 0, z + 1, 1, z, z, 1, z, z + 1, z + 1, 0, 1, 1], [z + 1, z + 1, 0, z, 1, 0, z, 1, z, z + 1, 1, 1, 0, 0, 0], [0, z + 1, z, z, 1, z, 0, z, 1, 1, z + 1, 1, 1, z, 1], [0, 1, z, z + 1, 0, 1, z + 1, 0, 0, z, 0, 1, 1, z + 1, z], [1, 0, z + 1, z, z, z, z, z, 0, 0, 0, 0, 0, z + 1, 0], [z, 0, 0, z, z + 1, 1, z, 1, z, z + 1, z, z, z, z + 1, z + 1], [1, 0, z + 1, z, z + 1, 0, z, 1, z, z, 0, 0, 1, z + 1, 0], [z + 1, 1, 1, z + 1, z, 0, 1, 0, z + 1, 0, 1, z, 0, z + 1, z], [1, 0, z + 1, 0, z, 1, z, z + 1, 1, 0, 0, 0, 1, z, 0], [z, z, z + 1, z, 1, z, 1, z, z, z, z + 1, 0, 0, 1, z], [z, z + 1, z, 1, z + 1, 1, z + 1, z + 1, z, z + 1, 0, z + 1, z, z, 0], [0, z + 1, 0, 1, 1, z, z + 1, z, z + 1, z, 0, z, z + 1, z + 1, z], [z + 1, z + 1, 0, 0, 1, z + 1, 1, 0, z + 1, 1, 0, 1, 0, z, z + 1], [0, 0, 1, z + 1, z, 1, 1, 1, z + 1, z, 1, z, z + 1, 1, z + 1], [1, 0, z + 1, z, z, 0, z + 1, 0, 1, z + 1, 1, z + 1, 0, 0, 1], [z + 1, 1, 1, 0, 1, 1, 1, z + 1, 0, z, 0, 0, 0, 0, 0], [z + 1, z + 1, z, z + 1, z, 0, z, z, z, z, 0, z, 0, z + 1, z], [1, 1, z + 1, 0, 1, z, 0, z, z, 1, 1, z + 1, 1, 0, z + 1], [1, z, z + 1, 0, z + 1, 0, 0, z + 1, 1, z, z, 0, z, z + 1, z], [z, 0, z + 1, 1, z + 1, 0, 1, z, z, 0, 0, 1, z, 0, z], [z, z, 0, 0, z + 1, 0, z + 1, z + 1, z, 0, 0, z + 1, 1, 0, 1], [1, 1, z, z, z + 1, z, z + 1, 1, z + 1, 0, 0, 1, 1, z, 0]])


# 32 * 16 k = 4 block_mds
G32_16_block = Matrix([[0, 0, z + 1, 0, 1, z + 1, 1, z + 1, z + 1, 0, 1, z + 1, z, 1, z, 1], [z + 1, 1, z + 1, 0, 0, z + 1, 1, 1, z, 1, 0, 0, z, 0, 0, 1], [z, z, z + 1, 0, z, z + 1, z, z + 1, 1, 0, z + 1, 0, z, z + 1, z, z], [0, 1, z, 1, z + 1, z, z, 1, z + 1, 0, 0, 0, z, 0, 1, z + 1], [z, 1, 1, z + 1, z, z + 1, z + 1, z + 1, 1, 0, z + 1, z, z + 1, 0, z, z + 1], [0, z + 1, 0, z + 1, z, z, 0, 0, 1, z + 1, 0, z + 1, z, z, 1, 1], [z, z + 1, 0, 0, z + 1, z, z + 1, 1, 0, z + 1, z, 1, 0, 1, 0, 0], [z, 0, z + 1, z, z, 0, z, z, 1, z, 1, z + 1, 0, z, z + 1, 0], [1, 0, 1, 1, 0, 1, 1, z, 1, 0, 1, 1, z + 1, z + 1, z, z + 1], [z, 1, 0, 0, z, z, z + 1, 1, 0, 0, 0, 0, z, 0, z + 1, 1], [z + 1, 1, 1, 1, z, 0, z, 0, z + 1, z + 1, 1, 1, z, 0, z, z], [z + 1, z + 1, 1, z, z + 1, z + 1, z + 1, 1, z + 1, 1, 1, 1, 0, z + 1, z, z + 1], [z + 1, 1, 0, z + 1, z + 1, 1, 1, 1, 1, 1, 0, z, 0, 0, z, 1], [0, z + 1, 1, z + 1, 1, 1, 1, 1, z, z + 1, z + 1, 0, 0, z, z + 1, 1], [z, z, 1, z + 1, 0, 0, 1, 1, 0, z, z + 1, z, 1, z, z, 0], [1, z + 1, 1, 1, z + 1, 1, z + 1, z, z + 1, 1, z, 1, 1, 1, z + 1, z], [0, z, 1, z + 1, z + 1, z + 1, z, z + 1, 1, z + 1, z + 1, 0, z, 1, z + 1, 0], [z + 1, 0, z + 1, z, 0, z + 1, 0, 1, 1, 1, 1, 0, 1, 1, z, 0], [z, 0, 1, 0, 1, z + 1, 1, z, 1, 0, z + 1, z + 1, z + 1, 1, 1, 0], [z + 1, 0, 0, z, 1, 1, 1, z, z, z, z + 1, 0, z, z + 1, z + 1, z + 1], [1, z, 0, z, 0, z + 1, 1, 0, 0, 1, 1, z, 0, z, 0, 0], [z + 1, 0, 0, 0, 1, z + 1, 0, 1, z, 0, 1, 0, 0, z + 1, 1, z], [z + 1, z, z, z + 1, z + 1, z, z + 1, 0, z, 1, z, z + 1, 0, 1, 1, 0], [0, z, z + 1, 1, z + 1, z, 0, 0, z, z, z + 1, 0, z + 1, z, 0, z + 1], [z, 0, 0, 1, z, z, 0, z, z, 1, z, z + 1, z, 0, z, z + 1], [z + 1, 0, 0, 1, z, z + 1, 1, z, 0, z, 1, 0, 0, z + 1, 1, z + 1], [z, 1, z, 1, 0, z, 0, 1, 0, z + 1, z + 1, 1, 1, z, 1, 1], [z + 1, 0, 1, z + 1, 0, 1, z + 1, 1, 0, 1, 0, 0, 1, 1, z + 1, 1], [0, 1, z + 1, z, 0, 1, z, z + 1, z, z + 1, 1, 0, 1, z, 0, 1], [0, 0, z + 1, z, 1, z + 1, 0, 1, 0, z, z + 1, z, 1, z, 0, z], [1, z, 1, z + 1, 1, z + 1, 1, 1, z + 1, 1, z + 1, 0, z + 1, 0, z + 1, z], [z, 0, z, 0, z, z, 0, z, z, z, 1, z, z, 1, 0, z + 1]])

# 24 * 8 k = 2 block_mds
G24_8_block = Matrix([[z, 1, z, z + 1, 1, 0, 1, 0], [1, 1, 1, z, 0, z, z + 1, z + 1], [0, z, z, 1, z + 1, z + 1, 0, z + 1], [0, 0, 1, z, z, 1, z + 1, 0], [z, z + 1, 0, 0, 1, z, z + 1, 0], [0, z + 1, z + 1, 0, 0, 0, 1, z + 1], [0, z + 1, z, z + 1, z + 1, z + 1, 0, 1], [1, z, 0, 0, z, 0, 1, 0], [z + 1, 1, z + 1, z + 1, 0, 0, 0, z], [0, 1, 0, z + 1, 1, z, 1, 0], [z, 1, 0, z + 1, z + 1, z + 1, z, 0], [1, z + 1, 1, 1, z + 1, 0, z + 1, z], [1, 1, 0, z, z, 1, 1, 1], [0, z + 1, 0, z, z, z + 1, z + 1, z + 1], [1, 1, z + 1, 0, 1, 1, z, 0], [z + 1, 0, z + 1, 0, z, z + 1, 1, z], [1, 1, 0, z + 1, 1, 1, z, z + 1], [z + 1, 0, z + 1, 1, z + 1, 1, z + 1, z + 1], [z, z, z, 0, 0, z, z, z], [z, z + 1, 1, z, z, 1, 1, z], [z, 1, z + 1, z + 1, z + 1, z + 1, 1, 0], [1, z + 1, z + 1, z + 1, 0, 0, 0, 0], [z, 1, 1, 1, z + 1, z, 0, z + 1], [0, 1, 0, z + 1, z + 1, 1, 1, z]])



# 28 * 12 k = 3 block_mds
G28_12_block = Matrix([[z + 1, 0, z + 1, 0, 0, 0, z, z + 1, 1, 1, z + 1, z], [z + 1, z + 1, 0, z, 0, z, 1, z + 1, 0, 1, z + 1, 1], [z + 1, z, 1, z, 1, z + 1, 0, z, 0, z + 1, z + 1, 1], [z, z, 0, 0, z + 1, 0, 1, 0, 1, 1, z, z], [z, z, 0, 1, 1, z, 0, 0, z + 1, z, z + 1, 1], [z + 1, z + 1, z + 1, 1, 0, z, 0, 1, z, z + 1, 0, z + 1], [1, z + 1, z, z + 1, 0, 1, 0, 0, z, z, 0, z], [z + 1, 0, z, 1, z, 0, 1, z + 1, 0, z + 1, 0, z + 1], [1, 0, z + 1, z, 1, z + 1, z + 1, 1, z + 1, z + 1, z, 1], [1, z + 1, z, 1, 1, z, z + 1, z, 1, 1, z, z + 1], [1, 1, z + 1, z, 1, z + 1, z + 1, 0, 0, z, z, 0], [1, z, z, 1, z + 1, 0, z + 1, z + 1, z + 1, z + 1, z, 0], [z + 1, 0, 0, 1, 0, 1, 1, 1, 0, z, z, 0], [z + 1, z, z, 1, 0, z + 1, z, z + 1, 0, z, z + 1, z], [0, z + 1, z, z + 1, 1, 0, 1, 1, z + 1, 0, z, z + 1], [z, 0, z + 1, z + 1, z, z, 1, 0, z, 1, z, z], [z, z, 1, z + 1, z, z + 1, 1, z, z + 1, z, z, 0], [z + 1, 0, 0, z, 0, 1, 0, z, 0, 0, 0, 1], [z + 1, 0, z + 1, z + 1, 0, z + 1, 1, z, z, 1, z, 1], [z, z + 1, z, 1, z, 0, 1, z + 1, z, 0, z + 1, 0], [0, 0, z + 1, 1, 1, z, z + 1, z + 1, 1, 1, z, 1], [z + 1, 0, 0, 1, z, 1, 0, z + 1, z + 1, 0, 0, 1], [1, 0, z + 1, z, z + 1, 0, z + 1, z, 1, z, z, 1], [z, z + 1, z, 0, 1, z + 1, 1, 0, z, 0, 1, 0], [1, 1, 1, 1, z + 1, z + 1, z + 1, z + 1, 0, z + 1, z + 1, z + 1], [z, z + 1, 1, z, 0, 0, 1, z, 1, z + 1, z, z + 1], [1, 1, z, z + 1, 0, z + 1, 0, 1, 1, 1, 1, 0], [0, 1, 1, z, z, 0, 0, 0, 1, 1, 0, 0]])



# 36 * 20 k = 5 block_mds
G36_20_block = Matrix([[0, z, 1, 0, 1, 0, 0, z + 1, 1, 0, z, 1, z, z + 1, z, z + 1, z + 1, z + 1, 0, 1], [1, z + 1, 0, 1, 0, z + 1, 0, z + 1, 1, z, 0, 1, 1, z, z, z, 1, 1, z + 1, z], [1, z, 0, z + 1, 0, 1, 1, 1, 0, z, 0, 1, 1, 0, z + 1, 1, z + 1, z, 0, 1], [z + 1, 0, z, z + 1, 1, z, 1, 1, z + 1, z, 1, z + 1, 1, z + 1, 1, 0, 1, z + 1, 1, z], [z + 1, 0, 1, 0, z, 0, z, 0, 1, 0, z + 1, 1, 0, z, 1, 0, z, z + 1, 0, 1], [0, 1, 1, z, 0, 0, 0, z, z + 1, z, z + 1, 1, 0, 1, 0, 0, z + 1, 0, z, 0], [z + 1, z + 1, z + 1, z, z, 0, z + 1, 1, 1, z, z, z, z, 1, 0, 0, z + 1, 1, z, 1], [z, 0, z + 1, 0, 0, z, 0, 1, 0, 1, 1, 1, 0, 0, 1, z + 1, z + 1, z, 1, z], [z, 0, z + 1, z + 1, z, z + 1, 0, 1, z + 1, 1, z + 1, 1, 0, 0, 0, z + 1, z + 1, z + 1, 0, z + 1], [z + 1, 0, 1, z, 1, 1, 1, z + 1, 0, 0, z + 1, z + 1, 1, z, z + 1, z, 0, z + 1, 1, z], [1, z + 1, z, 0, 1, z, z, z, z, 0, z, z + 1, z + 1, z, z, 1, 0, 1, 0, 1], [z + 1, z, z + 1, z, z, 0, 1, 0, z, z + 1, 0, z, z, 0, z, 0, z, z, 1, z], [0, 0, z, z + 1, z + 1, 0, 1, z, z, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, z], [z + 1, z + 1, 1, 1, 1, z, 1, 1, z + 1, 1, z + 1, z + 1, z, z, z + 1, 1, 1, z + 1, z + 1, 1], [z, z + 1, 0, 1, 1, z, 0, 0, 1, 0, 1, 1, 0, z, z + 1, 1, 0, 0, 1, z + 1], [z, z, z + 1, z + 1, z, z, 1, z + 1, z, z + 1, 0, 1, z, 1, 1, z, 0, z + 1, 0, 1], [z + 1, 1, 1, z, 0, z, z + 1, z, z + 1, 1, 0, 1, z + 1, 0, 0, z + 1, z, 1, z, z], [z, 1, z, z, z, 1, 1, 1, 0, z + 1, 1, z + 1, 0, 1, z + 1, z, 1, 1, z, 0], [0, 1, z, 0, z, 1, z + 1, 0, z, z, z, z + 1, 0, 0, z + 1, z + 1, z, z, z + 1, 0], [1, z + 1, 0, 0, 0, z, 0, z, z, 0, z, z + 1, 0, z, z + 1, 1, 1, 0, 0, z], [1, 0, 1, 1, z, z + 1, 1, 1, 0, 0, z, 0, 0, z + 1, 1, 1, z, z, z + 1, 0], [1, 0, z, 1, 1, 1, 0, 1, 0, z, z, 0, 1, 0, 0, 0, 0, z, z + 1, z + 1], [0, z + 1, 0, 1, 0, 1, 1, 0, z + 1, 0, z + 1, 1, 0, 0, 0, 1, 0, z, z, 1], [1, 1, 0, z + 1, 0, 0, 0, 1, 0, 1, z, z, z, z + 1, z, z + 1, z + 1, 0, z + 1, z + 1], [0, z, 0, 0, z + 1, z, z, z + 1, z, z, 1, z + 1, 0, z, 0, 0, 0, z, 1, z], [z + 1, z, z + 1, z + 1, 1, 1, z, z + 1, 1, z + 1, 1, 0, 0, 1, 0, 0, 1, 0, z, z], [1, z, z, z + 1, z, z, z, 1, 1, z + 1, z, z + 1, 0, z, 1, 0, z, 1, z, z + 1], [0, z, 0, 0, 1, 0, 0, 0, z + 1, 1, 1, z + 1, 0, 1, z + 1, 1, 0, z, z + 1, z + 1], [1, 0, 0, z, 0, z, z, 0, 0, 0, z + 1, z + 1, z, z + 1, 1, z, z + 1, 1, 1, 0], [0, 0, 1, z, z + 1, z, 0, 0, 1, z + 1, 0, 1, 0, z + 1, z, z, z, 0, z + 1, 0], [1, 1, 0, 0, z, 0, z, z, 1, z + 1, 1, 1, 1, z + 1, 1, 0, 1, z + 1, 0, 1], [z, 1, 1, z, 0, z + 1, 0, 1, 0, 0, z, z, z, z + 1, z + 1, z, 0, 0, z + 1, z], [0, 1, z + 1, z, 0, 1, z + 1, 0, z + 1, 0, z + 1, z, z + 1, z, z, z + 1, z + 1, 1, 0, 1], [0, z + 1, z + 1, z + 1, 0, 1, z, z + 1, z + 1, z, z, 0, z + 1, 0, 1, 1, 1, 1, 1, z + 1], [0, 0, z + 1, 1, z, 0, 0, z, 1, z, z, z, 1, z, z + 1, z + 1, 1, z, 0, z + 1], [0, 1, z + 1, z, 0, z, 0, 0, z + 1, z, z + 1, z, 1, z + 1, 1, 1, 1, 1, 0, 0]])


