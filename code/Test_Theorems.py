"""
This is a set of tests for some useful Theorems.
Normally if tests are passed, we continue to focusing on theoretical proofs.
Globally, 1 means success and 0 means fail
"""


from sympy import symbols, Matrix
from sympy.matrices import randMatrix, zeros, eye
from itertools import combinations
import random
from random import  choice
import numpy as np
import Functional_regenerating_code_test
import MDS_matrix_library

z = symbols('z')


def generate_MDS_H(n,alpha):
    H = eye(n*alpha)


def generate_strong_Z(alpha, d):
    """
    Z satisfies there are at least one "1" in each column
    """
    Z = zeros(alpha, d)
    while Z.rank() != alpha:
        Z = zeros(alpha, d)
        for i in range(d):
            N = random.randint(1, alpha)  # N means randomly choose how many "1" in this column
            randNs = random.sample(range(alpha), N)
            for j in range(N):
                t = randNs[j]
                Z[t, i] = 1
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

# def test_repair_1_block(num, G, k):
#     print("This time we try ", num, "times.")
#     count = 0
#     row = G.rows
#     col = G.cols
#     for i in range(num):
#         print("This is turn ", i, "in calculating matrix with size ", row, col)
#         if repair_1_block_strong(G, k) == 1:
#             count = count + 1
#     return count/num


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
        bi = Functional_regenerating_code_test.generate_b_row_vector(alpha)
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


# Test for multi-repair
def multi_repair_1(G, k, N):
    """
    This is to test multi-round repairs, corresponding to repair_1
    Each round random once, to check the success rate after N rounds
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1

    fail = 1  # fail = 1 means success and fail = 0 means fail

    repaired_matrix = G

    for t in range(N):
        fail_node = random.randint(1, n)  # randomly fail a node in each round

        access_nodes = set()
        while len(access_nodes) < d:
            access_nodes.add(random.randint(1, n))
            if fail_node in access_nodes:
                access_nodes.remove(fail_node)
        access_nodes = list(access_nodes)  # nodes we access when we repair

        access_matrix = Functional_regenerating_code_test.generate_b_full_rank(alpha) * repaired_matrix[(access_nodes[0] - 1) * alpha: access_nodes[0] * alpha, :]
        for i in range(1, d):
            access_matrix = Matrix([access_matrix, Functional_regenerating_code_test.generate_b_full_rank(alpha) * repaired_matrix[(access_nodes[i] - 1) * alpha:
                                                                                    access_nodes[i] * alpha, :]])

        # Z = Functional_regenerating_code_test.generate_Z(alpha, d)
        Z = generate_strong_Z(alpha, d)
        newcomer = Z * access_matrix
        newcomer = Functional_regenerating_code_test.tiny_binary_operation(newcomer)

        tempG = Matrix([repaired_matrix[0:(fail_node - 1) * alpha, :], repaired_matrix[fail_node * alpha:, :]])
        # tempG = Matrix([tempG, newcomer])

        # TODO: NEED DOUBLE CHECK
        for p in combinations(list(range(1, (n-1)*alpha)), (k - 1)*alpha):
            temp_matrix = tempG[(p[0] - 1):p[0], :]
            for i in range(1, (k - 1)*alpha):
                temp_matrix = Matrix([temp_matrix, tempG[(p[i] - 1):p[i], :]])
            temp_matrix = Matrix([temp_matrix, newcomer])
            det = temp_matrix.det()
            if det == 0:
                fail = 0
                break
        if fail == 0:
            break
        else:
            repaired_matrix = Matrix([tempG, newcomer])

    if fail == 0:
        return 0
    elif fail == 1:
        return 1


def test_Multi_repair_1(num, G, k, N):
    print("This time we try ", num, "times.")
    count = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if multi_repair_1(G, k, N) == 1:
            count = count + 1
    return count / num


def multi_repair_block(G, k, N):
    """
    *************************************************************************
    Note: this is the final version I use to test.
    On a certain matrix G, we fail its blocks once by once,
    to see if we can always find a way to repair after N rounds.
    *************************************************************************
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    trytime = 100
    Z_set = [Matrix([[1, 1, 0], [0, 1, 1]]), Matrix([[1, 0, 1], [1, 1, 0]]), Matrix([[1, 0, 1], [0, 1, 1]])]

    # fail = 0
    fail2 = 1   # if fail2 = 0, it means we try many times but cannot repair


    repaired_matrix = G
    for t in range(N):
        print("This is round", t)
        fail = 0  # If fail=1 it means we success in this round
        # fail_node = random.randint(1, n)
        fail_node = 3
        print("The failed node is node ", fail_node)

        access_nodes = set()
        while len(access_nodes) < d:
            access_nodes.add(random.randint(1, n))
            if fail_node in access_nodes:
                access_nodes.remove(fail_node)
        access_nodes = list(access_nodes)

        access_matrix = Functional_regenerating_code_test.generate_b_full_rank(alpha) * repaired_matrix[(access_nodes[0]- 1) * alpha: access_nodes[0] * alpha, :]
        for i in range(1, d):
            access_matrix = Matrix([access_matrix, Functional_regenerating_code_test.generate_b_full_rank(alpha) * repaired_matrix[(access_nodes[i] - 1) * alpha:
                                                                                       access_nodes[i] * alpha, :]])
        access_matrix = Functional_regenerating_code_test.tiny_binary_operation(access_matrix)

        print("access_matrix is", access_matrix)
        tempG = Matrix([repaired_matrix[0:(fail_node - 1) * alpha, :], repaired_matrix[fail_node * alpha:, :]])
        print("tempG we temporarily use after failure is", tempG)

        r = 0  # 有两种情况你会跳出while，要么你repair成功了，要么试了很多次发现再也没法repair了
        while fail == 0 and r <= trytime:  # r is round number
            print("This is trytime ", r, "in round", t)
            # Z = generate_strong_Z(alpha, d)
            Z = choice(Z_set)
            print("This time Z is", Z)

            newcomer = Z * access_matrix
            newcomer = Functional_regenerating_code_test.tiny_binary_operation(newcomer)
            print("This time newcomer is", newcomer)

            fail1 = 1  # to denote if repair success. If it's 0, means we fail
            for p in combinations(list(range(1, n)), k - 1):
                print("This is p", p)

                temp_matrix = tempG[(p[0] - 1) * alpha:p[0] * alpha, :]
                for i in range(1, k - 1):
                    temp_matrix = Matrix([temp_matrix, tempG[(p[i] - 1) * alpha:alpha * p[i], :]])

                temp_matrix = Matrix([temp_matrix, newcomer])
                print("After adding the newcomer, the temp_matrix is ", temp_matrix)
                det = temp_matrix.det()
                print("Before GF2, the temp_matrix's det is", det)
                det = Functional_regenerating_code_test.det_GF2(det)
                print("After GF2, the temp_matrix's det is", det)
                if det == 0:  # det = 0 means fail this time
                    fail1 = 0
                    print("No, this newcomer cannot fit")
                    break

            if fail1 == 0:  # means we fail this time
                r = r + 1
                if r == trytime - 2:
                    fail2 = 0

            if fail1 == 1:
                repaired_matrix = Matrix([tempG, newcomer])
                print("In round", t, "you success, let's torture this matrix again!")
                print("Now the repaired matrix is ", repaired_matrix)
                fail = 1

        if fail2 == 0:
            print("This stupid system is permanently down")
            break

    if fail2 == 1:
        print("We survive after", N, "rounds !")
        print("The final matrix we obtain is", repaired_matrix)



def multi_repair(G, k, N):
    """
    Note: this is the final version I use to test FOR MDS INSTEAD OF MDS_BLOCK
    On a certain matrix G, we fail its blocks once by once,
    to see if we can always find a way to repair after N rounds.
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    trytime = 100

    # fail = 0
    fail2 = 1   # if fail2 = 0, it means we try many times but cannot repair

    repaired_matrix = G
    for t in range(N):
        print("This is round", t)
        fail = 0  # If fail=1 it means we success in this round
        fail_node = random.randint(1, n)
        print("THe failed node is node", fail_node)

        access_nodes = set()
        while len(access_nodes) < d:
            access_nodes.add(random.randint(1, n))
            if fail_node in access_nodes:
                access_nodes.remove(fail_node)
        access_nodes = list(access_nodes)
        access_matrix = Functional_regenerating_code_test.generate_b_full_rank(alpha) * repaired_matrix[(access_nodes[0] - 1) * alpha: access_nodes[0] * alpha, :]
        for i in range(1, d):
            access_matrix = Matrix([access_matrix, Functional_regenerating_code_test.generate_b_full_rank(alpha) * repaired_matrix[(access_nodes[i] - 1) * alpha:
                                                                                       access_nodes[i] * alpha, :]])

        tempG = Matrix([repaired_matrix[0:(fail_node - 1) * alpha, :], repaired_matrix[fail_node * alpha:, :]])

        r = 0 # 有两种情况你会跳出while，要么你repair成功了，要么试了很多次发现再也没法repair了
        while fail == 0 and r <= trytime:  # r is round number
            print("This is trytime ", r, "in round", t)
            Z = generate_strong_Z(alpha, d)
            print("This time Z is", Z)
            newcomer = Z * access_matrix
            newcomer = Functional_regenerating_code_test.tiny_binary_operation(newcomer)
            print("This time newcomer is", newcomer)

            fail1 = 1  # to denote if repair success. If it's 0, means we fail
            if Functional_regenerating_code_test.test_mds(Matrix([newcomer, tempG]), k) == 0:
                fail1 = 0
                print("No, this newcomer cannot fit")

            if fail1 == 0:  # means we fail this time
                r = r + 1
                if r == trytime - 2:
                    fail2 = 0

            if fail1 == 1:
                repaired_matrix = Matrix([tempG, newcomer])
                print("In round", t, "you success, let's torture this matrix again!")
                print("Now the repaired matrix is ", repaired_matrix)
                fail = 1

        if fail2 == 0:
            print("This stupid system is permanently down")
            break

    if fail2 == 1:
        print("We survive after", N, "rounds !")
        print("The final matrix we obtain is", repaired_matrix)




def multi_repair_block_1(G, k, N):
    """
    This is to test multi-turn repairs
    Note: The test of MDS is still on block level
    :param G:
    :param k:
    :param N: number of repair turns
    :return:
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1

    fail = 1

    repaired_matrix = G

    for t in range(N):
        fail_node = random.randint(1, n)

        access_nodes = set()
        while len(access_nodes) < d:
            access_nodes.add(random.randint(1, n))
            if fail_node in access_nodes:
                access_nodes.remove(fail_node)
        access_nodes = list(access_nodes)

        access_matrix = Functional_regenerating_code_test.generate_b_row_vector(alpha) * repaired_matrix[(access_nodes[0] - 1) * alpha: access_nodes[0] * alpha, :]
        for i in range(1, d):
            access_matrix = Matrix([access_matrix, Functional_regenerating_code_test.generate_b_row_vector(alpha) * repaired_matrix[(access_nodes[i] - 1) * alpha:
                                                                                    access_nodes[i] * alpha, :]])
        Z = Functional_regenerating_code_test.generate_Z(alpha, d)
        newcomer = Z * access_matrix
        newcomer = Functional_regenerating_code_test.tiny_binary_operation(newcomer)

        tempG = Matrix([repaired_matrix[0:(fail_node - 1) * alpha, :], repaired_matrix[fail_node * alpha:, :]])
        for p in combinations(list(range(1, n)), k - 1):
            temp_matrix = tempG[(p[0] - 1) * alpha:p[0] * alpha, :]
            for i in range(1, k - 1):
                temp_matrix = Matrix([temp_matrix, tempG[(p[i] - 1) * alpha:alpha * p[i], :]])
            temp_matrix = Matrix([temp_matrix, newcomer])

            det = temp_matrix.det()
            if det == 0:
                fail = 0
                break
        if fail == 0:
            break
        else:
            repaired_matrix = Matrix([tempG, newcomer])
    if fail == 0:
        return 0
    else:
        return 1


def coverall_hamming_w(H, n):
    """
    This is used to randomly generate a new H
    whose newcomer's all Hamming weight is covered,
    after randomly one node fails
    :param H: Original H, H should be a square matrix
    :param n: Number of all nodes
    :return: New H
    """
    row = H.rows
    alpha = row // n
    fail_node = random.randint(1, n)
    remained_H = Matrix([H[0:(fail_node - 1) * alpha, :], H[fail_node * alpha:, :]])

    H_new = zeros(alpha, n * alpha)

    for i in range(n*alpha):
        temp_set = []
        while any(temp_set) is False:
            temp_set = np.random.randint(2, size=alpha)
            temp_set = list(temp_set)
        for j in range(alpha):
            H_new[i+j*row] = temp_set[j]

    H = Matrix([H_new, remained_H])
    return H


def test_coverall_hamming_w(G, k, N):
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    try_time = 40

    H = eye(row)
    count = 1
    for t in range(N):
        fail_node = random.randint(1, n)
        print("failed node", fail_node)

        access_nodes = set()
        while len(access_nodes) < d:
            access_nodes.add(random.randint(1, n))
            if fail_node in access_nodes:
                access_nodes.remove(fail_node)
        access_nodes = list(access_nodes)
        print("access_nodes", access_nodes)

        b = []
        b1 = Functional_regenerating_code_test.generate_b_row_vector(alpha)
        b.append(b1)
        for i in range(1, d):
            bi = Functional_regenerating_code_test.generate_b_row_vector(alpha)
            b.append(bi)
        print("bs", b)

        Z = Functional_regenerating_code_test.generate_Z(alpha, d)
        print("Z", Z)

        remained_H = Matrix([H[0:(fail_node - 1)*alpha, :], H[fail_node*alpha: , :]])
        print("remained_H", remained_H)

        # alpha * n(alpha) matrix denoting "helpers"
        H_helpers = b[0] * H[(access_nodes[0]-1) * alpha: access_nodes[0] * alpha, :]
        for i in range(1, d):
            H_helpers = Matrix([H_helpers, b[i] * H[(access_nodes[i] - 1) * alpha: access_nodes[i] * alpha, :]])
        print("H_helpers", H_helpers)

        #   Now use Z*H_helpers to calculate H_new denoting the "newcomer"
        H_new = Z * H_helpers

        #   Update H until we success(or until we reach largest try time)
        H = coverall_hamming_w(H, n)
        count1 = 0
        while Functional_regenerating_code_test.test_mds_block(H, n) == 0 and count1 <= try_time:
            H = coverall_hamming_w(H, n)
            count1 = count1 + 1

        print("Now we are using this H", H)

        # test H is MDS
        # fail = 1
        # for p in combinations(list(range(1, (n - 1))), (k - 1)):
        #     print(p)
        #     # print("remained_H", remained_H)
        #     temp_matrix = remained_H[(p[0]-1)*alpha:p[0]*alpha, :]
        #     for i in range(1, (k - 1) ):
        #         temp_matrix = Matrix([temp_matrix, remained_H[(p[i] - 1)*alpha:p[i]*alpha, :]])
        #     temp_matrix = Matrix([temp_matrix, H_new])
        #     rank = temp_matrix.rank()
        #     if rank != k * alpha:
        #         fail = 0
        #         print("This H is not MDS, its rank is", H.rank())
        #         break

        # test if H is MDS directly
        if Functional_regenerating_code_test.test_mds_block(H, n) == 0:
            print("This H is not MDS, its rank is", H.rank())
            print("We fail in round ", t + 1)
            break
        else:
            print("We success in round ", t + 1, "Let's go to next round")
            count = count + 1

    if count >= N:
        print("We survive! Cheers!")


def multi_repair_3_block(G, k, N):
    """
    This is to test multi-round repairs, corresponding to repair_3
    Each round random once, to check the success rate after N rounds
    :param G: Original encoding matrix
    :param k: Number of nodes needed to recover
    :param N: number of repair rounds
    :return: 0 to fail and 1 to success
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    count = 0

    try_time = 40

    # 在每一轮中，对于每一个node out of n，是不是ANY d nodes can repair，尝试足够多的random repair。如果是，再把这个node替换掉进入下一轮。
    # 确认是否有可能对于某个node死掉，在下一轮存在一个node修不回来了

    current_matrix = G
    for t in range(N):
        print("This is round ", t+1)
        # print("Now the current matrix we are testing is ", current_matrix)
        print("Row of current_matrix is ", current_matrix.rows, "Col of current_matrix is ", current_matrix.cols)
        for i in range(1, n + 1):
            print("This is to test node ", i)
            for p in combinations(list(range(1, n)), d):
                count2 = 0
                fail = 0
                print("p", p)
                while fail == 0 and count2 <= try_time:
                    print("In combination ", p, "we are trying time ", count2+1)
                    Z = Functional_regenerating_code_test.generate_Z(alpha, d)

                    temp_G = Matrix([current_matrix[0:(i - 1) * alpha, :], current_matrix[i * alpha:, :]])
                    access_matrix = Functional_regenerating_code_test.generate_b_row_vector(alpha) * temp_G[(p[0] - 1) * alpha: p[0] * alpha, :]
                    for j in range(1, d):
                        access_matrix = Matrix(
                            [access_matrix, Functional_regenerating_code_test.generate_b_row_vector(alpha) * temp_G[(p[j] - 1) * alpha: p[j] * alpha, :]])

                    newcomer = Functional_regenerating_code_test.tiny_binary_operation(Z * access_matrix)
                    # print("The newcomer is ", newcomer)

                    fail2 = 1
                    for q in combinations(list(range(1, (n-1)*alpha)), (k - 1)*alpha):
                        temp_matrix = temp_G[(q[0] - 1):q[0], :]
                        for j in range(1, (k - 1) * alpha):
                            temp_matrix = Matrix([temp_matrix, temp_G[(q[j] - 1):q[j], :]])
                        temp_matrix = Matrix([temp_matrix, newcomer])
                        det = temp_matrix.det()
                        if det == 0:
                            fail2 = 0
                            break
                    if fail2 == 0:   # repair fails
                        fail = 0
                        count2 = count2 + 1
                    elif fail2 == 1:
                        # print("MDS is protected")
                        fail = 1

                if fail == 0:
                    print("Now we cannot fix this node")
                    break
                else:
                    print("In turn p: ", p, "we succeed")
            if fail == 1:
                count = count + 1

        if count == n:
            print("In round ", t+1, "MDS is protected, now we can going to next round")
            fail_node = random.randint(1, n)
            access_nodes = set()
            while len(access_nodes) < d:
                access_nodes.add(random.randint(1, n))
                if fail_node in access_nodes:
                    access_nodes.remove(fail_node)
            access_nodes = list(access_nodes)

            access_matrix = Functional_regenerating_code_test.generate_b_row_vector(alpha) * current_matrix[(access_nodes[0] - 1) * alpha: access_nodes[0] * alpha, :]
            for i in range(1, d):
                access_matrix = Matrix([access_matrix, Functional_regenerating_code_test.generate_b_row_vector(alpha) * current_matrix[(access_nodes[i] - 1) * alpha:
                                                                                      access_nodes[i] * alpha, :]])

            temp_G = Matrix([current_matrix[0:(fail_node - 1) * alpha, :], current_matrix[fail_node * alpha:, :]])

            failure = 0  # to check if this repair success or not
            count3 = 0
            while failure == 0 and count3 <= try_time:
                Z = Functional_regenerating_code_test.generate_Z(alpha, d)
                newcomer = Functional_regenerating_code_test.tiny_binary_operation(Z * access_matrix)
                print("To come to next round, newcomer is", newcomer)

                fail3 = 1
                for p in combinations(list(range(1, (n-1)*alpha)), (k - 1)*alpha):
                    temp_matrix = temp_G[(p[0] - 1):p[0], :]
                    for i in range(1, (k - 1)*alpha):
                        temp_matrix = Matrix([temp_matrix, temp_G[(p[i] - 1):p[i], :]])
                    temp_matrix = Matrix([temp_matrix, newcomer])
                    # print(temp_matrix.rows, temp_matrix.cols)

                    det = temp_matrix.det()
                    if det == 0:
                        fail3 = 0
                        break

                if fail3 == 1:
                    failure = 1
                    current_matrix = Matrix([temp_G, newcomer])
                else:
                    failure = 0
                    print("Oh this newcomer is bad, let's try one more time")
                    count3 = count3 + 1

        if count3 >= try_time:
            print("In round", t + 1, "we fail")
            break

        if t == N-1:
            print("After ", N, "round we survive! Cheers!")


# def test_cover_hamming_weight_on_matrix_H(num, n, k, alpha):
#     """
#     This is to test d=n-1 case regenerating code on representing matrix H,
#     based on theorem that H MDS rank iff G MDS
#     :param num:
#     :param n:
#     :param k:
#     :param alpha:
#     :return:
#     """
#     H = eye(n*alpha)
#     d = alpha + k - 1
#
#     fail_block = random.randint(1, n)
#
#     access_nodes = set()
#     for i in range(1, 1 + n):
#         access_nodes.add(i)
#     access_nodes.remove(fail_block)
#
#
# access_matrix = generate_b_full_rank(alpha) * G[(access_nodes[0]-1) * alpha: access_nodes[0] * alpha, :]


def test_Multi_repair_block_1(num, G, k, N):
    print("This time we try ", num, "times.")
    count = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if multi_repair_block_1(G, k, N) == 1:
            count = count + 1
    return count/num


if __name__ == "__main__":

    # print("The success probability for G8_4 is: ", test_Multi_repair_1(50, MDS_matrix_library.G8_4, 2, 3))
    # print("The success probability for G8_4 is: ", test_Multi_repair_1(500, MDS_matrix_library.G8_4, 2, 2))
    # print("The success probability for G8_4 is: ", test_Multi_repair_1(500, MDS_matrix_library.G8_4, 2, 3))


    # print("The success probability fgenerate_b_full_rankor G10_6(repair for 3 turn) is: ", test_Multi_repair_1(50, MDS_matrix_library.G10_6_block, 3, 3))


    # test Multi_3
    # print("For a matrix G8_4", multi_repair_3_block(MDS_matrix_library.G10_6_block, 3, 200))
    # test_coverall_hamming_w(MDS_matrix_library.G10_6_block, 3, 100)
