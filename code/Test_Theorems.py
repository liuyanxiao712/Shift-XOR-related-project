"""
This is a set of tests for some useful Theorems.
Normally if tests are passed, we continue to focusing on theoretical proofs.
"""

from sympy import symbols, Matrix
from sympy.matrices import randMatrix, zeros, eye
from itertools import combinations
import random
import numpy as np
import Functional_regenerating_code_test
import MDS_matrix_library

z = symbols('z')


def generate_strong_Z(alpha, d):
    """
    This Z satisfies that there are no less than k element "1"  in its each row
    :param alpha:
    :param d:
    :return:
    """
    k = d+1-alpha
    Z = zeros(alpha, d)
    # print("Z",Z)
    for i in range(alpha):
        N = random.randint(k, d)   # N is the number of element "1" in this row
        randNs = random.sample(range(d), N)
        # print(randNs)
        for j in range(N):
            t = randNs[j]
            Z[i, t] = 1
    return Z


def repair_1_block_strong(G, k):
    """ Totally random b, Z, d and fail_node, try to see if success
    :param G: matrix to be repaired, failed node is randomly chosen in this marix
    :param k: originally we have k nodes
    :return: 0 to fail and 1 to success
    """
    # print("G", G)
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    fail_node = random.randint(1,n)
    Z = generate_strong_Z(alpha, d)
    while Z.rank() != alpha:
        Z =generate_strong_Z(alpha, d)
    # print(row, col, alpha, n, d, fail_node, Z)
    print("Z", Z)

    access_nodes = set()
    while len(access_nodes) < d:
        access_nodes.add(random.randint(1,n))
        if fail_node in access_nodes:
            access_nodes.remove(fail_node)
    access_nodes = list(access_nodes)
    # print("access_nodes", access_nodes)

    access_matrix = Functional_regenerating_code_test.generate_b_row_vector(alpha) * G[(access_nodes[0]-1) * alpha: access_nodes[0] * alpha, :]

    for i in range(1, d):
        access_matrix = Matrix([access_matrix, Functional_regenerating_code_test.generate_b_row_vector(alpha) * G[(access_nodes[i]-1)*alpha: access_nodes[i]*alpha, :]])
    # print("access_matrix", access_matrix)
    newcomer = Z * access_matrix
    newcomer = Functional_regenerating_code_test.tiny_binary_operation(newcomer)
    print("newcomer", newcomer)


    fail = 1
    tempG = Matrix([G[0:(fail_node - 1)*alpha, :], G[fail_node*alpha:, :]])
    # print("tempG", tempG)
    for p in combinations(list(range(1, n)), k-1):
        temp_matrix = tempG[(p[0]-1)*alpha:p[0]*alpha, :]
        for i in range(1, k-1):
            temp_matrix = Matrix([temp_matrix, tempG[(p[i]-1)*alpha:alpha*p[i], :]])
        temp_matrix = Matrix([temp_matrix, newcomer])

        det = temp_matrix.det()
        print("det", det)

        if det == 0:
            fail = 0
            break
    # print(fail)

    if fail == 0:
        return 0
    else:
        print("success")
        return 1

def test_repair_1_block(num, G, k):
    print("This time we try ", num, "times.")
    count = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if repair_1_block_strong(G, k) == 1:
            count = count + 1
    return count/num


def Theorm2(G, k):
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    fail_node = random.randint(1, n)
    Z = Functional_regenerating_code_test.generate_Z(alpha, d)



    access_nodes = set()
    while len(access_nodes) < d:
        access_nodes.add(random.randint(1, n))
        if fail_node in access_nodes:
            access_nodes.remove(fail_node)
    access_nodes = list(access_nodes)

    b = []
    b1 = Functional_regenerating_code_test.generate_b_row_vector(alpha)
    # access_matrix =  b1* G[access_nodes[0] * alpha: (access_nodes[0] + 1) * alpha, :]
    b.append(b1)

    for i in range(1, d):
        bi =Functional_regenerating_code_test.generate_b_row_vector(alpha)
        b.append(bi)
        # access_matrix = Matrix([access_matrix,
        #                         bi * G[(access_nodes[i] - 1) * alpha: access_nodes[i] * alpha,
        #                                                        :]])
    original_H = eye(row)

    H_sub = b[0] * original_H[access_nodes[0] * alpha: (access_nodes[0] + 1) * alpha, :]
    for i in range(1, d):
        H_sub = Matrix([H_sub, b[i] * original_H[(access_nodes[i] - 1) * alpha: access_nodes[i] * alpha, :]])

    H = Z * H_sub

    tempH = Matrix([original_H[0:(fail_node - 1) * alpha, :], original_H[fail_node * alpha:, :]])

    H = Matrix([tempH, H])

    # print("H", H)
    fail = 1
    for p in combinations(list(range(n)), k):
        # print(p)
        temp_matrix = H[(p[0]-1)*alpha:p[0]*alpha, :]
        for i in range(1, k):
            temp_matrix = Matrix([temp_matrix, H[(p[i] - 1) * alpha:alpha * p[i], :]])
        # print("temp matrix", temp_matrix, temp_matrix.rows, temp_matrix.cols)

        rank = temp_matrix.rank()

        if rank != k*alpha:
            break
            fail = 0
    if fail == 0:
        return 0
    else:
        return 1

def test_Theorm2(num, G, k):
    print("This time we try ", num, "times.")
    count = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if Theorm2(G, k) == 1:
            count = count + 1
    return count/num




if __name__ == "__main__":
    print("The success probability for G8_4 is: ", test_Theorm2(50, MDS_matrix_library.G8_4, 2))


