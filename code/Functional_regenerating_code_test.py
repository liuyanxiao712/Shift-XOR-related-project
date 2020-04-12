"""
This is the library for our use in manipulating MDS matrices

THIS IS THE MAIN PROGRAM I USE NOW

Please mainly use repair_1_block, repair_2_block, repair_3_block
"""

from sympy import symbols, Matrix
from sympy.matrices import randMatrix, zeros, ones

from itertools import combinations
import random
import MDS_matrix_library

z = symbols('z')

def output_para(row, col, k):
    alpha = col // k
    n = row//alpha
    d = alpha + k - 1
    print(alpha, n, d)

def generate_b_col_vector(alpha):
    b = randMatrix(alpha, 1, 0, 1)
    while b.rank() == 0:
        b = randMatrix(alpha, 1, 0, 1)
    return b

def binary_operation(G):
    G = G.subs(z ** 18, z ** 51)
    G = G.subs(z ** 16, z ** 53)
    G = G.subs(z ** 14, z ** 55)
    G = G.subs(z ** 12, z ** 57)
    G = G.subs(z ** 10, z ** 59)
    G = G.subs(z ** 8, z ** 61)
    G = G.subs(z ** 6, z ** 63)
    G = G.subs(z ** 4, z ** 65)
    G = G.subs(z ** 2, z ** 67)
    # G = G.subs(lambda x: x % 2 == 0, 0)
    return G

def tiny_binary_operation(G):
    G = G.subs(2, 0)
    G = G.subs(4, 0)
    G = G.subs(6, 0)
    G = G.subs(8, 0)
    G = G.subs(3, 1)
    G = G.subs(5, 1)
    G = G.subs(7, 1)
    G = G.subs(9, 1)
    return G

def GF2(G):
    """
    Make the target operate on GF(2)
    """
    G = binary_operation(G)
    G = tiny_binary_operation(G)
    return G


def generate_b_full_rank(alpha):
    """
    Generate a row vector with all entries 1 and length alpha
    """
    b = ones(1, alpha)
    return b


def generate_b_row_vector(alpha):
    b = randMatrix(1, alpha, 0, 1)
    while b.rank() == 0:
        b = randMatrix(1, alpha, 0, 1)
    return b


def generate_Z(alpha, d):
    Z = randMatrix(alpha, d, 0, 1)
    while Z.rank() != alpha:
        Z = randMatrix(alpha, d, 0, 1)
    return Z


def test_mds(G):
    row = G.rows
    col = G.cols
    items = list(range(1, row+1))
    failure = 1
    for q in combinations(items, col):
        temp_matrix = G[q[0]-1:q[0], :]
        for i in range(1, col):
            temp_matrix = Matrix([temp_matrix, G[q[i]-1:q[i], :]])
        det = temp_matrix.det()
        if det == 0:
            failure = 0   # fail
            break
            return 0
    if failure == 1:
        return 1   # success

def generate_strong_Z(alpha, d):
    """
    This Z satisfies that there are no less than k element "1"  in its each row
    :param alpha:
    :param d:
    :return:
    """
    k = d+1-alpha
    Z = zeros(alpha, d)
    for i in range(alpha):
        N = random.randint(k, d)   # N is the number of element "1" in this row
        randNs = random.sample(range(d), N)
        for j in range(N):
            t = randNs[j]
            Z[i, t] = 1
    while Z.rank() != alpha:
        generate_strong_Z(alpha, d)
    return Z




# Z84 = Matrix(2, 3, [1, 1, 1, 1, 1, 0])
# Z106 = Matrix(2, 4, [1, 1, 1, 1, 1, 1, 1, 0])
# Z128 = Matrix(2, 5, [1, 1, 1, 1, 1, 1, 1, 0,1,1])
# Z1410 =Matrix(2, 6, [1, 1, 1, 1, 1, 1,1,1, 1, 0,1,1])
Z156= Matrix(3,4, [0, 1, 1, 1, 1, 0, 1, 1, 1,1,0,1])
# Z189 =Matrix(3,5, [0, 1, 1, 1, 1, 1, 1, 1, 1,1,0,1,1,1,0])
# Z248=Matrix([[1, 0, 0, 1, 1], [0, 0, 1, 1, 1], [0, 1, 1, 1, 0], [0, 0, 1, 1, 0]])

def test_mds_block(G, k):
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    # d = alpha + k - 1
    fail = 1
    for q in combinations(list(range(1, n+1)), k):
        # print(q)
        new_matrix = G[(q[0]-1)*alpha:q[0]*alpha, :]
        for i in range(1, k):
            new_matrix = Matrix([new_matrix, G[(q[i]-1)*alpha:q[i]*alpha, :]])
        det = new_matrix.det()
        # print(det)
        if binary_operation(det) == 0:
            fail = 0
    if fail == 0:
        return 0 # fail
    else:
        return 1 # success

def generate_mds_f(i,j):
    if i == j:
        return random.choice([0, 1, z, 1 + z])
    else:
        return random.choice([0, 1, z, 1 + z])


def generate_mds(row, col):
    a = 0
    while a == 0:
        MDS = Matrix(row, col, generate_mds_f)
        if test_mds(MDS) == 0:
            a = 0;
        else:
            a = 1
    # print("for matrix", row, "*", col, ", the MDS matrix is ", MDS)
    return MDS

def generate_mds_block(row, col, k):
    alpha = col // k
    a = 0
    while a == 0:
        MDS = Matrix(row, col, generate_mds_f)
        if test_mds_block(MDS, k) == 0:
            a = 0
        else:
            a = 1
    return MDS


def repair_1(G, k):
    """ Totally random b, Z, d and fail_node, try to see if success
    :param G: matrix to be repaired, failed node is randomly chosen in this marix
    :param k: originally we have k nodes
    :return: 0 to fail and 1 to success
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    fail_node = random.randint(1,n)
    Z = generate_Z(alpha, d)
    # print(row, col, alpha, n, d, fail_node, Z)

    access_nodes = set()
    while len(access_nodes) < d:
        access_nodes.add(random.randint(1,n))
        if fail_node in access_nodes:
            access_nodes.remove(fail_node)
    access_nodes = list(access_nodes)
    # print("access_nodes", access_nodes)

    access_matrix = generate_b_row_vector(alpha) * G[access_nodes[0]*alpha : (access_nodes[0] + 1) * alpha, :]

    for i in range(1,d):
        access_matrix = Matrix([access_matrix, generate_b_row_vector(alpha) * G[(access_nodes[i]-1)*alpha : access_nodes[i]*alpha, :]])
    # print("access_matrix", access_matrix)
    newcomer = Z * access_matrix

    fail = 1
    tempG = Matrix([G[0:(fail_node - 1)*alpha, :], G[fail_node*alpha: , :]])
    # for i in range(alpha):
    #     tempG.row_del((fail_node - 1)*alpha)
    for p in combinations(list(range(1, tempG.rows+1)), col - newcomer.rows):
        # print(p)
        temp_matrix = tempG[p[0]-1:p[0], :]
        for i in range(1, col - newcomer.rows):
            temp_matrix = Matrix([temp_matrix, G[p[i]-1:p[i], :]])
        temp_matrix = Matrix([temp_matrix, newcomer])
        temp_matrix = tiny_binary_operation(temp_matrix)
        if binary_operation(temp_matrix.det()) == 0:
            fail = 0
    if fail == 0:
        print("This time fails, Z is", Z, "newcomer is: ", newcomer)
        return 0
    else:

        return 1

def repair_1_block(G, k):
    """ Totally random b, Z, d and fail_node, try to see if success
    :param G: matrix to be repaired, failed node is randomly chosen in this marix
    :param k: originally we have k nodes
    :return: 0 is fail and 1 is success
    """
    # print("G", G)
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    fail_node = random.randint(1, n)
    # Z = generate_strong_Z(alpha, d)
    Z = Z156
    # print(row, col, alpha, n, d, fail_node, Z)

    access_nodes = set()
    while len(access_nodes) < d:
        access_nodes.add(random.randint(1,n))
        if fail_node in access_nodes:
            access_nodes.remove(fail_node)
    access_nodes = list(access_nodes)
    # print("access_nodes", access_nodes)

    access_matrix = generate_b_full_rank(alpha) * G[(access_nodes[0]-1) * alpha: access_nodes[0] * alpha, :]

    for i in range(1, d):
        access_matrix = Matrix([access_matrix, generate_b_full_rank(alpha) * G[(access_nodes[i]-1)*alpha: access_nodes[i]*alpha, :]])
    # print("access_matrix", access_matrix)
    newcomer = Z * access_matrix
    newcomer = tiny_binary_operation(newcomer)
    # print("newcomer", newcomer)


    fail = 1
    tempG = Matrix([G[0:(fail_node - 1)*alpha, :], G[fail_node*alpha:, :]])
    # print("tempG", tempG)
    for p in combinations(list(range(1, n)), k-1):
        # print(p)
        temp_matrix = tempG[(p[0]-1)*alpha:p[0]*alpha, :]
        for i in range(1, k-1):
            temp_matrix = Matrix([temp_matrix, tempG[(p[i]-1)*alpha:alpha*p[i], :]])
        temp_matrix = Matrix([temp_matrix, newcomer])
        # print("temp_matrix before", temp_matrix)
        # temp_matrix = tiny_binary_operation(temp_matrix)
        # print("temp_matrix after", temp_matrix)

        det = temp_matrix.det()
        # print("det", det)


        # This is to determine if all determinant components are even
        # if (det/2) == (det/2).subs(z**1/2, 0):
        #     det = 0
        #     print("oh")

        if det == 0:
            fail = 0
            break

    if fail == 0:
        return 0
    else:
        return 1


def repair_2(G, k):
    """
    b, d, fail_node are set beforehand and we try Z to see how many times we need to reapir successfully
    :param G:G: Original encoding matrix
    :param k:Number of nodes needed to recover
    :return: number of success cases
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    fail_node = random.randint(1, n)
    b = generate_b_row_vector(alpha)

    access_nodes = set()
    while len(access_nodes) < d:
        access_nodes.add(random.randint(1, n))
        if fail_node in access_nodes:
            access_nodes.remove(fail_node)
    access_nodes = list(access_nodes)

    access_matrix = b * G[access_nodes[0]*alpha : (access_nodes[0] + 1) * alpha, :]
    for i in range(1, d):
        access_matrix = Matrix([access_matrix, b * G[(access_nodes[i] - 1) * alpha: access_nodes[i] * alpha, :]])
    num1 = 1

    tempG = Matrix([G[0:(fail_node - 1) * alpha, :], G[fail_node * alpha:, :]])
    fail = False

    Z = generate_Z(alpha, d)
    newcomer = binary_operation(Z * access_matrix)

    for p in combinations(list(range(1, tempG.rows+1)), col - alpha):
        temp_matrix = tempG[p[0] - 1:p[0], :]
        for i in range(1, col - alpha):
            temp_matrix = Matrix([temp_matrix, G[p[i] - 1:p[i], :]])
        temp_matrix = Matrix([temp_matrix, newcomer])
        if binary_operation(temp_matrix) == 0:
            fail = True
    while fail is True:
        Z = generate_Z(alpha, d)
        newcomer = binary_operation(Z * access_matrix)
        for p in combinations(list(range(1, tempG.rows + 1)), col - alpha):
            temp_matrix = tempG[p[0] - 1:p[0], :]
            for i in range(1, col - alpha):
                temp_matrix = Matrix([temp_matrix, G[p[i] - 1:p[i], :]])
            temp_matrix = Matrix([temp_matrix, newcomer])
            if binary_operation(temp_matrix) == 0:
                fail = True
            else:
                fail = False
        num1 = num1 + 1
    return num1


def repair_2_block(G, k):
    """
    b vectors are set, try ALL possibilities of b(each is set beforehand), random Z to see if we can always fix random failures
    :param G: Original encoding matrix
    :param: Number of nodes needed to recover
    :return: 0 to fail and 1 to success
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    b_set = list(range(n))
    b_set[0] = generate_b_row_vector(alpha)
    tempG = b_set[0] * G[0:alpha, :]
    for i in range(1, n):
        b_set[i] = generate_b_row_vector(alpha)
        tempG = Matrix([tempG, b_set[i]*G[alpha*i:(alpha)*(i+1), :]])

    fail = 0
    count = 0

    while fail == 0 and count <= 40:
        fail_node = random.randint(1, n)   # random failure
        Z = generate_Z(alpha, d)   # random Z, transformation matrix

        access_nodes = set() # random access nodes to repair
        while len(access_nodes) < d:
            access_nodes.add(random.randint(1, n))
            if fail_node in access_nodes:
                access_nodes.remove(fail_node)
        access_nodes = list(access_nodes)

        # accessing matrix we use to repair
        access_matrix = tempG[(access_nodes[0] - 1):access_nodes[0], :]
        for i in range(1,d):
            access_matrix = Matrix([access_matrix, tempG[(access_nodes[i] - 1):access_nodes[i], :]])

        newcomer = tiny_binary_operation(Z * access_matrix)

        new_G = Matrix([G[0:(fail_node - 1) * alpha, :], G[fail_node * alpha:, :]])
        # new_G = Matrix([new_G, newcomer])

        fail2 = 1
        for p in combinations(list(range(1, n)), k - 1):
            temp_matrix = new_G[(p[0] - 1) * alpha:p[0] * alpha, :]
            for i in range(1, k - 1):
                temp_matrix = Matrix([temp_matrix, new_G[(p[i] - 1) * alpha:alpha * p[i], :]])
            temp_matrix = tiny_binary_operation(Matrix([temp_matrix, newcomer]))
            det = temp_matrix.det()
            if det == 0:
                fail2 = 0
                break
        if fail2 == 0:
            fail = 0
        else:
            fail = 1
        count = count + 1
    if fail == 0:
        print("For matrix", G.rows, "*", G.cols, "we cannot fix it")
        print("The set of b vectors here is", b_set)
        return 0
    if fail == 1:
        print("In", count, "attempt we success to repair", G.rows, "*", G.cols, "matrix")
        return 1








def repair_3_block(G, k):
    """
    Make sure for a matrix, whichever a node fails, we can random b and Z to repair it.
    Iterate all failure possibilities to see if we can always successfully repair.
    Detail: For any 1 node fails, we access ANY d out of remaining nodes, can repair it in several random rounds.
    :param G: Original encoding matrix
    :param k: Number of nodes needed to recover
    :return: 0 to fail and 1 to success
    """
    row = G.rows
    col = G.cols
    alpha = col // k
    n = row // alpha
    d = alpha + k - 1
    count = 0
    fail = 0
    for i in range(1, k+1):
        temp_G = Matrix([G[0:(i-1)*alpha, :], G[i*alpha:, :]])
        for p in combinations(list(range(1, n)), d):
            count2 = 0
            while fail == 0 and count2 <= 50:
                print("count2 in while loop", count2, "fail", fail)
                Z = generate_Z(alpha, d)
                access_matrix = generate_b_row_vector(alpha) * temp_G[(p[0]-1)*alpha : p[0]*alpha, :]
                for j in range(1, d):
                    access_matrix = Matrix([access_matrix, generate_b_row_vector(alpha)*temp_G[(p[j]-1)*alpha : p[j]*alpha, :]])
                newcomer = Z * access_matrix
                temp_matrix = tiny_binary_operation(Matrix([temp_G, newcomer]))
                fail = test_mds_block(temp_matrix, k)
                count2  = count2 + 1
            if fail == 1:
                break
        if fail == 0:
            break
        elif fail == 1:
            count = count + 1
    if count == k:
        print("---------------------For matrix ", G.rows, "*", G.cols, "we can successfully repair above matrix, no matter which node fails----------------------")
    else:
        print("---------------------For matrix ", G.rows, "*", G.cols, "sometimes we cannot fix above matrix-----------------------")


def repair_4_block(G, k):
    """ Only difference with repair_1_block is: this is repair-by-transfer, vector b is to select only one packet in a node
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
    Z = generate_Z(alpha, d)
    # print(row, col, alpha, n, d, fail_node, Z)

    access_nodes = set()
    while len(access_nodes) < d:
        access_nodes.add(random.randint(1,n))
        if fail_node in access_nodes:
            access_nodes.remove(fail_node)
    access_nodes = list(access_nodes)
    # print("access_nodes", access_nodes)

    access_matrix = generate_b_1_vector(alpha) * G[(access_nodes[0]-1) * alpha: access_nodes[0] * alpha, :]

    for i in range(1,d):
        access_matrix = Matrix([access_matrix, generate_b_1_vector(alpha) * G[(access_nodes[i]-1)*alpha : access_nodes[i]*alpha, :]])
    # print("access_matrix", access_matrix)
    newcomer = Z * access_matrix
    # print("newcomer", newcomer)

    fail = 1
    tempG = Matrix([G[0:(fail_node - 1)*alpha, :], G[fail_node*alpha: , :]])
    # print("tempG", tempG)
    for p in combinations(list(range(1, n)), k-1):
        # print(p)
        temp_matrix = tempG[(p[0]-1)*alpha:p[0]*alpha, :]
        for i in range(1, k-1):
            temp_matrix = Matrix([temp_matrix, tempG[(p[i]-1)*alpha:alpha*p[i], :]])
        temp_matrix = Matrix([temp_matrix, newcomer])
        temp_matrix = tiny_binary_operation(temp_matrix)
        # print("temp_matrix", temp_matrix)
        det = temp_matrix.det()
        # print("det", det)
        if det == 0:
            fail = 0
            # print("before binary operation", det)
            # print("after binary operation", binary_operation(det))
            break
    # print(fail)
    if fail == 0:
        # print("This time fails")
        # print("This time fails, Z is", Z, "newcomer is: ", newcomer)
        return 0
    else:
        return 1



def test_repair_1(num, G, k):
    print("This time we try ", num, "times.")
    count = 0
    for i in range(num):
        print("This is turn ", i)
        if repair_1(G, k) == 1:
            count = count + 1
    return count/num

def test_repair_1_block(num, G, k):
    print("This time we try ", num, "times.")
    count = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if repair_1_block(G, k) == 1:
            count = count + 1
    return count/num

def test_repair_2_block(num, G, k):
    print("This time we try ", num, "times.")
    count1 = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if repair_2_block(G, k) == 1:
            count1 = count1 + 1
    return count1/num


def test_repair_2(num, G, k):
    print("This time we try ", num, "times.")
    poss_set = []
    for i in range(num):
        print("This is turn ", i)
        num1 = repair_2(G, k)
        poss_set.append(num1)
    result = {}
    for key in poss_set:
        result[key] = result.get(key, 0) + 1
    print("There are B times we need random Z for A times to repair, A:B", result)

def test_repair_4_block(num, G, k):
    # print("This time we try ", num, "times.")
    count = 0
    row = G.rows
    col = G.cols
    for i in range(num):
        print("This is turn ", i, "in calculating matrix with size ", row, col)
        if repair_4_block(G, k) == 1:
            count = count + 1
    return count/num





if __name__ == "__main__":
    print(Matrix([[1, 0, 0, 0], [0, 1, 1, 1]]).rank())
    # test 1

    # print("The success probability for G8_4 is: ", test_repair_1_block(1000, MDS_matrix_library.G8_4, 2))   # 0.11
    # print("The success probability for G10_6_block is: ", test_repair_1_block(1000, MDS_matrix_library.G10_6_block, 3))   # 0.525 in 200 tries
    # print("The success probability for G12_8_block is: ", test_repair_1_block(1000, MDS_matrix_library.G12_8_block, 4))   # 0.675 in 200 tries
    # print("The success probability for G14_10_block is: ", test_repair_1_block(1000, MDS_matrix_library.G14_10_block, 5))   # 0.88 in 200 tries
    # print("The success probability for G16_12_block is: ", test_repair_1_block(100, MDS_matrix_library.G16_12_block, 6))   # 0.925 in 200 tries
    #
    # print("The success probability for G14_8_block is: ", test_repair_1_block(100, MDS_matrix_library.G14_8_block, 4))   # 0.655 in 200 tries
    # print("The success probability for G24_8_block is: ", test_repair_1_block(50, MDS_matrix_library.G24_8_block, 2))  # 0.79
    # print("The success probability for G24_15_block is: ", test_repair_1_block(50, MDS_matrix_library.G24_15_block, 5))   # 0.94

    # print("The success probability for G15_6 is: ", test_repair_1_block(1000, MDS_matrix_library.G15_6_block, 2))  #   0.485 in 200 tries
    # print("The success probability for G18_9 is: ", test_repair_1_block(100, MDS_matrix_library.G18_9_block, 3))  #   0.815 in 200 tries
    # print("The success probability for G20_16_block is: ", test_repair_1_block(50, MDS_matrix_library.G20_16_block, 8))
    # print("The success probability for G21_12_block is: ", test_repair_1_block(100, MDS_matrix_library.G21_12_block, 4))   # 0.91
    # print("The success probability for G18_14_block is: ", test_repair_1_block(50, MDS_matrix_library.G18_14_block, 7))   # 0.95
    # print("The success probability for G32_16_block is: ", test_repair_1_block(1, MDS_matrix_library.G32_16_block, 4))   # 0.96
    # print("The success probability for G28_12_block is: ", test_repair_1_block(100, MDS_matrix_library.G28_12_block, 3))      # 0.84
    # print("The success probability for G36_20_block is: ", test_repair_1_block(100, MDS_matrix_library.G36_20_block, 5))





    # test 2

    # print("The success probability for G8_4 is: ", test_repair_2_block(100, MDS_matrix_library.G8_4, 2))# mostly in 1-3 attempts and a few 10-20 attempts we success
    # print("The success probability for G10_6_block is: ", test_repair_2_block(200, MDS_matrix_library.G10_6_block, 3))# mostly in 1-5 attempts we success
    # print("The success probability for G12_8_block is: ", test_repair_2_block(200, MDS_matrix_library.G12_8_block, 4))   # mostly in 1-3 attempts we success
    # print("The success probability for G14_8_block is: ", test_repair_2_block(50, MDS_matrix_library.G14_8_block, 4))# mostly in 1-2 attempts we success
    # print("The success probability for G14_10_block is: ", test_repair_2_block(30, MDS_matrix_library.G14_10_block, 5))# mostly 1 and a few 2 attempts we success
    # print("The success probability for G16_12_block is: ", test_repair_2_block(30, MDS_matrix_library.G16_12_block, 6))# mostly in 1 attempt we success
    # print("The success probability for G15_6 is: ", test_repair_2_block(30, MDS_matrix_library.G15_6_block, 2))# mostly in 1-4 attempts we success
    # print("The success probability for G18_9 is: ", test_repair_2_block(30, MDS_matrix_library.G18_9_block, 3))# mostly in 1-2 attempts we success
    # print("The success probability for G21_12_block is: ", test_repair_2_block(30, MDS_matrix_library.G21_12_block, 4))# mostly 1 and a few 2 attempt we success
    # print("The success probability for G24_15_block is: ", test_repair_2_block(30, MDS_matrix_library.G24_15_block, 5))# mostly in 1 attempt we success
    # print("The success probability for G28_12_block is: ", test_repair_2_block(30, MDS_matrix_library.G28_12_block, 3))
    # print("The success probability for G32_16_block is: ", test_repair_2_block(30, MDS_matrix_library.G32_16_block, 4))# mostly 1 and a few 2 attempt we success
    # print("The success probability for G36_20_block is: ", test_repair_2_block(30, MDS_matrix_library.G36_20_block, 5))
    # print("The success probability for G20_16_block is: ", test_repair_2_block(30, MDS_matrix_library.G20_16_block, 8))






    # test 3

    # print("For a matrix G8_4", repair_3_block(MDS_matrix_library.G8_4, 2))   # success
    # print("For a matrix G10_6", repair_3_block(MDS_matrix_library.G10_6_block, 3))    # success
    # print("For a matrix G12_8", repair_3_block(MDS_matrix_library.G12_8_block, 4))    # success
    # print("For a matrix G14_8", repair_3_block(MDS_matrix_library.G14_8_block, 4))    # success
    # print("For a matrix G14_10", repair_3_block(MDS_matrix_library.G14_10_block, 5))    # success
    # print("For a matrix G16_12", repair_3_block(MDS_matrix_library.G16_12_block, 6))    # success
    # print("For a matrix G15_6", repair_3_block(MDS_matrix_library.G15_6_block, 2))    # success
    # print("For a matrix G18_9", repair_3_block(MDS_matrix_library.G18_9_block, 3))    # success
    # print("For a matrix G20_16", repair_3_block(MDS_matrix_library.G20_16_block, 8))    # success
    # print("For a matrix G21_12", repair_3_block(MDS_matrix_library.G21_12_block, 4))    # success
    # print("For a matrix G32_16", repair_3_block(MDS_matrix_library.G32_16_block, 4))    # success
    # print("For a matrix G28_12", repair_3_block(MDS_matrix_library.G28_12_block, 3))    # success
    # print("For a matrix G36_20", repair_3_block(MDS_matrix_library.G36_20_block, 5))


    # test 4

    # print("The success probability for G8_4 is: ", test_repair_4_block(200, MDS_matrix_library.G8_4, 2))   #
    # print("The success probability for G10_6_block is: ", test_repair_4_block(200, MDS_matrix_library.G10_6_block, 3))   #
    # print("The success probability for G12_8_block is: ", test_repair_4_block(200, MDS_matrix_library.G12_8_block, 4))   #
    # print("The success probability for G14_10_block is: ", test_repair_4_block(200, MDS_matrix_library.G14_10_block, 5))   #
    # print("The success probability for G16_12_block is: ", test_repair_4_block(100, MDS_matrix_library.G16_12_block, 6))   #
    #
    # print("The success probability for G24_8_block is: ", test_repair_4_block(100, MDS_matrix_library.G24_8_block, 2))  #
    # print("The success probability for G24_15_block is: ", test_repair_4_block(100, MDS_matrix_library.G24_15_block, 5))  #
    # print("The success probability for G18_14_block is: ", test_repair_4_block(100, MDS_matrix_library.G18_14_block, 7))  #


    # print("The success probability for G15_6 is: ", test_repair_4_block(100, MDS_matrix_library.G15_6_block, 2))  #
    # print("The success probability for G18_9 is: ", test_repair_4_block(100, MDS_matrix_library.G18_9_block, 3))  #
    # print("The success probability for G20_16_block is: ", test_repair_4_block(100, MDS_matrix_library.G20_16_block, 8))
    # print("The success probability for G21_12_block is: ", test_repair_4_block(100, MDS_matrix_library.G21_12_block, 4))   #
    # print("The success probability for G28_12_block is: ", test_repair_4_block(40, MDS_matrix_library.G28_12_block, 3))      #
    # print("The success probability for G24_15_block is: ", test_repair_4_block(40, MDS_matrix_library.G24_15_block, 5))  #
    # print("The success probability for G32_16_block is: ", test_repair_4_block(30, MDS_matrix_library.G32_16_block, 4))   #
    #
    # print("The success probability for G36_20_block is: ", test_repair_4_block(30, MDS_matrix_library.G36_20_block, 5))