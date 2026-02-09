"""
Module: GaloisTools
"""
import galois
import numpy as np
from mip import Model, BINARY, minimize, xsum


##  CHAPTER: ====Solve Linear Equations====
def solve(matrix, vector):
    """
    input.matrix：galois.GF(2) 或 01矩阵，形状 (m, n)
    input.vector：galois.GF(2) 或 01向量，形状 (n,)
    output：galois.GF(2) 向量解；若无解返回 None
    """

    # PART: Details
    assert matrix.shape[1] == vector.shape[0]
    A = matrix.T
    b = vector.reshape(-1, 1)

    # PART: Details
    aug = np.concatenate((A, b), axis=1)
    n, m_plus_1 = aug.shape
    m = m_plus_1 - 1
    rank = 0

    for col in range(m):

        pivot_row = None
        for i in range(rank, n):
            if aug[i, col] == 1:
                pivot_row = i
                break
        if pivot_row is None:
            continue  # 该列没有主元，跳过

        if pivot_row != rank:
            aug[[rank, pivot_row], :] = aug[[pivot_row, rank], :]

        for i in range(n):
            if i != rank and aug[i, col] == 1:
                aug[i, :] ^= aug[rank, :]
        rank += 1

    for i in range(rank, n):
        if aug[i, -1] == 1:
            return None  # 无解的情况

    GF=galois.GF(2)
    solution = GF.Zeros(m)
    leading_cols = []

    for i in range(rank):
        for j in range(m):
            if aug[i, j] == 1:
                leading_cols.append(j)
                break

    for i in range(rank):
        col = leading_cols[i]
        solution[col] = aug[i, -1]

        for j in range(col + 1, m):
            if aug[i, j] == 1:
                solution[col] ^= solution[j]

    # PART: Details
    return solution


##  CHAPTER: ====Subspace Intersection====
def cap(matrix1, matrix2):
    """""
    input.matrix1：GF(2)矩阵，行张成子空间1
    input.matrix2：GF(2)矩阵，行张成子空间2
    output：GF(2)矩阵，表示交集的行空间基
    """""

    ##  PART: Details
    m = matrix1.shape[0]

    ## basis1^T * x + basis2^T * y = 0 (), [basis1^T | basis2^T]
    aug_matrix = galois.GF(2)(np.concatenate((matrix1.T, matrix2.T), axis=1))

    nullspace = aug_matrix.null_space()

    ab_space = nullspace[:, :m]

    ## basis1
    if len(ab_space) == 0:
        return galois.GF(2).Zeros((0, matrix1.shape[1]))

    intersection_vectors = ab_space @ matrix1

    rref_intersection = intersection_vectors.row_reduce()

    ##  0
    nz_mask = np.any(rref_intersection != 0, axis=1)
    rref_basis = rref_intersection[nz_mask]

    # PART: Details
    return rref_basis


##  CHAPTER: ====Subspace Difference====
def minus(matrix1, matrix2):
    """""
    input.matrix1：GF(2)矩阵，行张成子空间1
    input.matrix2：GF(2)矩阵，行张成子空间2
    output：列表[ndarray]，为子空间差的基向量（01行向量）
    """""

    # PART: Details
    if matrix1 is None:
        return None
    if matrix2 is None or len(matrix2) == 0:
        return galois.GF(2)(matrix1)
    assert matrix1.shape[1] == matrix2.shape[1]

    intersect = cap(matrix1, matrix2)
    if len(intersect) == 0:
        return galois.GF(2)(matrix1)

    result = []
    for i in range(len(matrix1)):
        rank = np.linalg.matrix_rank(intersect)
        intersect = np.vstack((intersect, matrix1[i]))
        if np.linalg.matrix_rank(intersect) > rank:
            result.append(matrix1[i])

    # PART: Details
    return galois.GF(2)(np.array(result,dtype=int))


##  CHAPTER: ====Subspace Direct Sum====
def direct_sum(matrix1, matrix2):
    """
    input.matrix1：GF(2)矩阵
    input.matrix2：GF(2)矩阵
    output：GF(2)矩阵，直和后的独立基
    """

    # PART: Details
    assert isinstance(matrix1, galois.GF(2))
    assert isinstance(matrix2, galois.GF(2))

    # PART: Details
    result = matrix1[0].copy()

    for i in range(1, len(matrix1)):
        if len(result.shape)>1:
            rank = np.linalg.matrix_rank(result)
        else:
            rank=1
        temp = np.vstack((result, matrix1[i]))
        if np.linalg.matrix_rank(temp) > rank:
            result = temp
    for i in range(len(matrix2)):
        if len(result.shape)>1:
            rank=np.linalg.matrix_rank(result)
        else:
            rank=1
        temp = np.vstack((result, matrix2[i]))
        if np.linalg.matrix_rank(temp) > rank:
            result = temp

    # PART: Details
    return galois.GF(2)(result)


##  CHAPTER: ====Orthogonalization====
def orthogonalize(matrix):
    """""
    input.matrix：可迭代的GF(2)行向量列表/矩阵
    output：列表[GF(2)行向量]，两两正交的基
    """""

    B_i=[v.copy() for v in matrix]
    ortho_basis=[]  # 存储正交基

    while True:

 ## B_i
        length=len(B_i)
        for i in range(length):
            if np.mod(np.count_nonzero(B_i[i]), 2)==1:
                for j in range(length):
                    if j!=i and np.mod(np.count_nonzero(B_i[j]), 2)==0:
                        B_i[j]=B_i[j]+B_i[i]
                break
        flag=0
        b1=B_i[0]
        o_i=b1
        ortho_basis.append(o_i)
        next_B=[]
        for j in range(len(B_i)):
            if j!=flag:
                b=B_i[j]
                coef=np.dot(b, o_i)
                b_new=b+coef*o_i
                next_B.append(b_new)
        B_i=next_B
        if len(B_i)==0:
            break

    return ortho_basis


##  CHAPTER: ====Minimum Hamming Distance====
def mip_distance_caculator(H, logicOp, code_type='pauli')->int:
    """""
    input.H：01矩阵或GF(2)矩阵
    input.logicOp：01矩阵或GF(2)矩阵，逻辑算子候选
    input.code_type：'pauli' 或 'majorana'
    output：int，最小权重
    """""

    ##  PART: ----Return code distance----
    H = np.array(H, dtype=int)  # 转换为整数类型的numpy数组
    logicOp = np.array(logicOp, dtype=int)  # 转换为整数类型的numpy数组
    d = H.shape[1]  # 初始化距离为量子比特数量（最大可能距离）

    ##  PART: Details
    for i in range(logicOp.shape[0]):
        logicOp_i = logicOp[i, :]
        n = H.shape[1]  # 量子比特数量（稳定子矩阵的列数）
        m = H.shape[0]  # 稳定子数量（稳定子矩阵的行数）

        num_anc_stab = int(np.ceil(np.log2(2 * n + 1)))
        num_anc_logical = int(np.ceil(np.log2(2 * n + 1)))
        num_var = n + m * num_anc_stab + num_anc_logical  # 总变量数量

        model = Model()
        model.verbose = 0  # 关闭详细输出
        x = [model.add_var(var_type=BINARY) for i in range(num_var)]  # 创建二进制变量数组
        model.objective = minimize(xsum(x[i] for i in range(n)))  # 目标函数：最小化Hamming权重

        for row in range(m):
            weight = [0] * num_var
            supp = np.nonzero(H[row, :])[0]
            w_stab = len(supp)

            if code_type == 'majorana':
                w_s_mod2 = w_stab % 2
                for q in range(n):
                    coeff = 0
                    if q in supp:
                        coeff += 1
                    if w_s_mod2:
                        coeff += 1
                    weight[q] = coeff
            else:
                for q in supp:
                    weight[q] = 1

            cnt = 1
            for q in range(num_anc_stab):
                weight[n + row * num_anc_stab + q] = -(1 << cnt)
                cnt += 1
            model += xsum(weight[i] * x[i] for i in range(num_var)) == 0
        supp = np.nonzero(logicOp_i)[0]
        w_log = len(supp)
        weight = [0] * num_var

        if code_type == 'majorana':
            w_l_mod2 = w_log % 2
            for q in range(n):
                coeff = 0
                if q in supp:
                    coeff += 1
                if w_l_mod2:
                    coeff += 1
                weight[q] = coeff
        else:
            for q in supp:
                weight[q] = 1

        cnt = 1
        for q in range(num_anc_logical):
            weight[n + m * num_anc_stab + q] = -(1 << cnt)
            cnt += 1

        model += xsum(weight[i] * x[i] for i in range(num_var)) == 1

        model.optimize()
        
        if model.num_solutions > 0:
            opt_val = np.sum([x[i].x for i in range(n)])
            d = min(d, int(opt_val))

    ##  PART: Return code distance
    return d

def random_distance_caculator_css(H_X,H_Z,number):
    F = galois.GF(2)
    w_X = F(H_Z.null_space())  # 计算Z稳定子的零空间（X逻辑算子空间）
    w_Z = F(H_X.null_space())  # 计算X稳定子的零空间（Z逻辑算子空间）
    logical_matrix_X=minus(w_X,H_X)
    logical_matrix_Z=minus(w_Z,H_Z)
    rows_wx, cols_wx = w_X.shape  # 获取零空间矩阵的维度信息
    rows_wz, cols_wz = w_Z.shape  # 获取零空间矩阵的维度信息
    dist_bound = cols_wx + 1  # 初始化距离上界为最大可能值（列数+1）
    vec_count = 0  # 计数器：记录找到当前最小权重的向量数量

    for i in range(number):
        per = np.random.permutation(cols_wx)  # 生成随机排列，用于随机化搜索顺序
        wx1 = w_X[:, per]  # 对零空间矩阵的列进行随机排列
        wx2 = wx1.row_reduce()  # 对排列后的矩阵进行行约简（高斯消元）
        wx2 = wx2[:, np.argsort(per)]  # 将列顺序恢复为原始顺序

        for j in range(rows_wz):
            temp_vec = wx2[j, :]  # 获取当前行向量
            temp_weight = np.count_nonzero(temp_vec)  # 计算向量的Hamming权重（非零元素个数）

            if 0 < temp_weight <= dist_bound:
                if np.any(logical_matrix_X @ temp_vec):

                    if temp_weight < dist_bound:
                        dist_bound = temp_weight
                        vec_count = 1

                    elif temp_weight == dist_bound:
                        vec_count += 1

            if dist_bound <= 2:
                return 2

    for i in range(number):
        per = np.random.permutation(cols_wz)  # 生成随机排列，用于随机化搜索顺序
        wz1 = w_Z[:, per]  # 对零空间矩阵的列进行随机排列
        wz2 = wz1.row_reduce()  # 对排列后的矩阵进行行约简（高斯消元）
        wz2 = wz2[:, np.argsort(per)]  # 将列顺序恢复为原始顺序

        for j in range(rows_wz):
            temp_vec = wz2[j, :]  # 获取当前行向量
            temp_weight = np.count_nonzero(temp_vec)  # 计算向量的Hamming权重（非零元素个数）

            if 0 < temp_weight <= dist_bound:
                if np.any(logical_matrix_Z @ temp_vec):

                    if temp_weight < dist_bound:
                        dist_bound = temp_weight
                        vec_count = 1

                    elif temp_weight == dist_bound:
                        vec_count += 1

            if dist_bound <= 2:
                return 2

    return dist_bound

##  CHAPTER: ====Random Heuristic====
def random_distance_caculator(stabilizers_matrix, gauge_matrix, number)->int:
    """""
    input.stabilizers_matrix：np.ndarray of GF(2), stabilizer check matrix
    input.gauge_matrix：np.ndarray of GF(2), gauge check matrix
    input.number：int, number of random search
    output：int, estimated code distance
    """""

    F = galois.GF(2)
    w = F(stabilizers_matrix.null_space())  # 计算Z稳定子的零空间（X逻辑算子空间）
    logical_matrix=minus(w, gauge_matrix)
    logical_matrix=minus(logical_matrix,stabilizers_matrix)
    rows_wz, cols_wz = w.shape  # 获取零空间矩阵的维度信息
    dist_bound = cols_wz + 1  # 初始化距离上界为最大可能值（列数+1）
    vec_count = 0  # 计数器：记录找到当前最小权重的向量数量

    for i in range(number):
        per = np.random.permutation(cols_wz)  # 生成随机排列，用于随机化搜索顺序
        wz1 = w[:, per]  # 对零空间矩阵的列进行随机排列
        wz2 = wz1.row_reduce()  # 对排列后的矩阵进行行约简（高斯消元）
        wz2 = wz2[:, np.argsort(per)]  # 将列顺序恢复为原始顺序

        for j in range(rows_wz):
            temp_vec = wz2[j, :]  # 获取当前行向量
            temp_weight = np.count_nonzero(temp_vec)  # 计算向量的Hamming权重（非零元素个数）

            if 0 < temp_weight <= dist_bound:
                if np.any(logical_matrix @ temp_vec):

                    if temp_weight < dist_bound:
                        dist_bound = temp_weight
                        vec_count = 1

                    elif temp_weight == dist_bound:
                        vec_count += 1

            if dist_bound <= 2:
                return 2

    return dist_bound

