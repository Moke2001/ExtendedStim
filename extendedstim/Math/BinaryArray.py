import copy
import galois
import numpy as np
from mip import Model, xsum, minimize, BINARY


class BinaryArray:
    GF = galois.GF(2)

    #%%  USER：构造方法
    def __init__(self, array):
        if isinstance(array, BinaryArray):
            self._array = copy.deepcopy(array._array)  # 拷贝生成二进制数组
        else:
            self._array = self.GF(np.array(array, dtype=int))  # 转换为二进制数组
        self.shape = self._array.shape  # 数组形状

    #%%  USER：重载运算符
    ##  USER：返回数组的长度
    def __len__(self):
        return self.shape[0]

    ##  USER：判断两个数组是否相等
    def __eq__(self, other):
        if isinstance(other, BinaryArray):
            return self._array == other._array
        else:
            return self._array == other

    ##  USER：获取数组的元素
    def __getitem__(self, item):
        return self._array[item]

    ##  USER：设置数组的元素
    def __setitem__(self, item, value):
        self._array[item] = value

    ##  USER：返回数组的字符串表示
    def __str__(self):
        return str(self._array)

    ##  USER：数组相加
    def __add__(self, other):
        assert (isinstance(other, BinaryArray))
        return BinaryArray(self._array + other._array)

    ##  USER：数组相减
    def __sub__(self, other):
        assert (isinstance(other, BinaryArray))
        return BinaryArray(self._array - other._array)

    ##  USER：数组相乘
    def __mul__(self, other):
        assert (isinstance(other, BinaryArray)) or (isinstance(other, int))
        if isinstance(other, int):
            return BinaryArray(self._array * other)
        else:
            return BinaryArray(self._array * other._array)

    ##  USER：数组右乘
    def __rmul__(self, other):
        assert (isinstance(other, BinaryArray)) or (isinstance(other, int))
        if isinstance(other, int):
            return BinaryArray(self._array * other)
        else:
            return BinaryArray(self._array * other._array)

    ##  USER：数组矩阵乘
    def __matmul__(self, other):
        assert (isinstance(other, BinaryArray))
        return BinaryArray(self._array @ other._array)

    ##  USER：数组右矩阵乘
    def __rmatmul__(self, other):
        assert (isinstance(other, BinaryArray))
        return BinaryArray(other._array @ self._array)

    ##  USER：数组幂运算
    def __pow__(self, power):
        return BinaryArray(np.linalg.matrix_power(self._array, power))

    #%%  USER：属性方法
    ##  USER：返回数组的转置
    @property
    def T(self):
        return BinaryArray(self._array.T)

    ##  USER：求占据的位置
    @property
    def occupy(self):
        return np.where(self._array != 0)[0]

    ##  USER：返回数组的权重
    @property
    def weight(self):
        return np.count_nonzero(self._array)

    ##  USER：返回数组的零空间
    @property
    def null_space(self):
        return BinaryArray(self._array.null_space())

    ##  USER：返回数组的秩
    @property
    def rank(self):
        return np.linalg.matrix_rank(self._array)

    #%%  USER：对象方法
    ##  USER：复制函数
    def copy(self):
        return copy.deepcopy(self)

    #%%  USER：静态方法
    ##  USER：矩阵求和
    @staticmethod
    def sum(array_list):
        assert (isinstance(array_list, list))
        assert (isinstance(array_list[0], BinaryArray))
        return BinaryArray(np.sum([array._array for array in array_list], axis=0))

    ##  USER：求数组的水平堆叠
    @staticmethod
    def hstack(array, *args):
        assert len(args) > 0
        assert isinstance(array, BinaryArray)
        for temp in args:
            assert isinstance(temp, BinaryArray)
        result = array
        for temp in args:
            result = np.hstack((result, temp._array))
        return BinaryArray(result)

    ##  USER：求数组的垂直堆叠
    @staticmethod
    def vstack(array, *args):
        assert len(args) > 0
        assert isinstance(array, BinaryArray)
        for temp in args:
            assert isinstance(temp, BinaryArray)
        result = array
        for temp in args:
            result = np.vstack((result, temp._array))
        return BinaryArray(result)

    ##  USER：求解线性系统
    @staticmethod
    def solve(matrix, vector):
        assert isinstance(matrix, BinaryArray)
        assert isinstance(vector, BinaryArray)
        assert matrix.shape[1] == vector.shape[0]

        ##  使用GF(2)数组
        A = matrix._array.T
        b = vector._array.reshape(-1, 1)

        # 构造增广矩阵
        aug = np.concatenate((A, b), axis=1)
        n, m_plus_1 = aug.shape
        m = m_plus_1 - 1
        rank = 0

        ##  前向消元
        for col in range(m):
            ##  寻找主元
            pivot_row = None
            for i in range(rank, n):
                if aug[i, col] == 1:
                    pivot_row = i
                    break

            if pivot_row is None:
                ##  该列没有主元，跳过
                continue

            ##  交换行
            if pivot_row != rank:
                aug[[rank, pivot_row], :] = aug[[pivot_row, rank], :]

            ##  消去当前列下方和上方的元素
            for i in range(n):
                if i != rank and aug[i, col] == 1:
                    aug[i, :] ^= aug[rank, :]
            rank += 1

        ##  检查是否有解
        for i in range(rank, n):
            if aug[i, -1] == 1:
                ##  矛盾方程，无解
                return None

        ##  构造解向量
        solution = BinaryArray.GF.Zeros(m)
        leading_cols = []

        ##  找出主元列
        for i in range(rank):
            for j in range(m):
                if aug[i, j] == 1:
                    leading_cols.append(j)
                    break

        ##  回代求解
        for i in range(rank):
            col = leading_cols[i]
            solution[col] = aug[i, -1]

            ##  消去当前行中主元列右侧的元素对解的影响
            for j in range(col + 1, m):
                if aug[i, j] == 1:
                    solution[col] ^= solution[j]

        ##  返回解向量
        return BinaryArray(solution)

    ##  USER：返回两个线性空间的交集
    @staticmethod
    def cap(matrix1, matrix2):
        assert isinstance(matrix1, BinaryArray)
        assert isinstance(matrix2, BinaryArray)
        GF2 = galois.GF(2)

        ##  转化为GF2数组
        matrix1 = matrix1._array
        matrix2 = matrix2._array
        m = matrix1.shape[0]
        k = matrix2.shape[0]

        ##  解方程 basis1^T * x + basis2^T * y = 0 (表示交集向量), 构造矩阵 [basis1^T | basis2^T]
        aug_matrix = GF2(np.concatenate((matrix1.T, matrix2.T), axis=1))

        ##  计算零空间（解空间）
        nullspace = aug_matrix.null_space()

        ##  从零空间中取对应于每个解向量的系数（对应于方程组中的解向量）
        ab_space = nullspace[:, :m]

        ##  将系数乘以 basis1 得到具体解（交集）
        if len(ab_space) == 0:
            return GF2.Zeros((0, matrix1.shape[1]))

        intersection_vectors = ab_space @ matrix1

        ##  产生基向量的线性无关组合，即每种向量的不同的高斯行组合样本
        rref_intersection = intersection_vectors.row_reduce()

        ##  除去全0行的行（即为没有有用基底捕获的时候）
        nz_mask = np.any(rref_intersection != 0, axis=1)
        rref_basis = rref_intersection[nz_mask]

        ##  返回补集
        return BinaryArray(rref_basis)

    ##  USER：返回两个线性空间的差集
    @staticmethod
    def minus(matrix1, matrix2):
        assert isinstance(matrix1, BinaryArray)
        assert isinstance(matrix2, BinaryArray)

        ##  转化为GF2数组
        matrix1_origin = matrix1
        matrix2_origin = matrix2
        matrix1 = matrix1._array
        matrix2 = matrix2._array

        intersect = BinaryArray.cap(matrix1_origin, matrix2_origin)
        if len(intersect) == 0:
            return matrix1

        result = []
        for i in range(len(matrix1)):
            rank = np.linalg.matrix_rank(intersect)
            intersect = np.vstack((intersect, matrix1[i]))
            if np.linalg.matrix_rank(intersect) > rank:
                result.append(matrix1[i])

        ##  返回差集
        return BinaryArray(result)

    ##  USER：返回两个线性空间的直积
    @staticmethod
    def direct_sum(matrix1, matrix2):
        assert isinstance(matrix1, BinaryArray)
        assert isinstance(matrix2, BinaryArray)

        ##  转化为GF2数组
        matrix1 = matrix1._array
        matrix2 = matrix2._array
        result = matrix1[0].copy()

        ##  直接拼接
        for i in range(1, len(matrix1)):
            rank = np.linalg.matrix_rank(result)
            temp = np.vstack((result, matrix1[i]))
            if np.linalg.matrix_rank(temp) > rank:
                result = temp
        for i in range(len(matrix2)):
            rank = np.linalg.matrix_rank(result)
            temp = np.vstack((result, matrix2[i]))
            if np.linalg.matrix_rank(temp) > rank:
                result = temp

        ##  返回结果
        return BinaryArray(result)

    ##  USER：正交化基矢组
    @staticmethod
    def orthogonalize(matrix):
        def gf2_gram_schmidt(vectors, bilinear_form):
            k = len(vectors)
            ##  复制向量以避免修改原始数据
            B_i = [v.copy() for v in vectors]
            O = []  # 存储正交基
            for i in range(k):
                if np.mod(np.count_nonzero(B_i[i]), 2) == 1:
                    for j in range(k):
                        if j != i and np.mod(np.count_nonzero(B_i[j]), 2) == 0:
                            B_i[j] = B_i[j] + B_i[i]
                    break
                if i == k - 1:
                    print('There is no odd-weight vector')
                    return False

            while True:
                length = 10000
                b1 = None
                flag = 0
                b1 = B_i[0]
                # for i in range(len(B_i)):
                #     if length > np.count_nonzero(B_i[i]):
                #         b1 = B_i[i]
                #         flag = i
                #         length = np.count_nonzero(B_i[i])

                o_i = b1

                O.append(o_i)

                next_B = []
                for j in range(len(B_i)):
                    if j != flag:
                        b = B_i[j]
                        coef = bilinear_form(b, o_i)
                        b_new = b + coef * o_i
                        next_B.append(b_new)

                B_i = next_B
                if len(B_i) == 0:
                    return O
                length = len(B_i)
                for i in range(length):
                    if np.mod(np.count_nonzero(B_i[i]), 2) == 1:
                        for j in range(length):
                            if j != i and np.mod(np.count_nonzero(B_i[j]), 2) == 0:
                                B_i[j] = B_i[j] + B_i[i]
                        break
                    if i == len(B_i):
                        return O

        ##  双线性形式矩阵
        bilinear_form1 = lambda u, v: np.dot(u, v)
        ortho_basis1 = gf2_gram_schmidt(matrix, bilinear_form1)

        ##  输出结果
        return BinaryArray(ortho_basis1)

    ##  USER：计算张量积
    @staticmethod
    def kron(matrix1, matrix2):
        assert isinstance(matrix1, BinaryArray)
        assert isinstance(matrix2, BinaryArray)
        return BinaryArray(np.kron(matrix1._array, matrix2._array))

    ##  USER：计算code distance
    @staticmethod
    def distance(H, method):
        assert isinstance(H, BinaryArray)
        if method == 'mip':
            logic_op = BinaryArray.minus(H.null_space, H)
            return mip_distance_caculator(H._array, logic_op._array)
        elif method == 'random':
            return random_distance_caculator(H._array, H._array, 500, mindist=1, debug=1, field=2)

    ##  USER：计算子系统距离
    @staticmethod
    def subsystem_distance(stabilizers, gauges, method):
        pass

    ##  USER：根据占用位置生成二进制数组
    @staticmethod
    def FromOccupy(occupy, *args):
        if len(args) == 0:
            assert isinstance(occupy, list) or isinstance(occupy, np.ndarray)
            temp = np.zeros(np.max(occupy) + 1)
            temp[occupy] = 1
            return BinaryArray(temp)
        elif len(args) == 1:
            assert isinstance(occupy, list) or isinstance(occupy, np.ndarray)
            assert isinstance(args[0], int)
            temp = np.zeros(args[0], dtype=int)
            temp[occupy] = 1
            return BinaryArray(temp)
        else:
            raise ValueError("FromOccupy() takes 1 or 2 positional arguments but {} were given".format(len(args) + 1))

    ##  USER：根据数组生成二进制数组
    @staticmethod
    def FromArray(array):
        return BinaryArray(array)

    ##  USER：生成全零数组
    @staticmethod
    def zeros(number):
        return BinaryArray(np.zeros(number, dtype=int))

    ##  USER：生成单位矩阵
    @staticmethod
    def eye(number):
        return BinaryArray(np.eye(number, dtype=int))

    ##  USER：生成循环左移矩阵
    @staticmethod
    def shift(number, shift):
        S = np.zeros((number, number), dtype=int)
        for i in range(number):
            S[i, (i + shift) % number] = 1
        return BinaryArray(S)

    ##  USER：生成全一数组
    @staticmethod
    def ones(number):
        return BinaryArray(np.ones(number, dtype=int))


#%%  KEY：计算code distance
def mip_distance_caculator(H, logicOp):
    """"
    计算量子纠错码的代码距离

    该函数通过计算每个逻辑算子的最小权重，来确定量子纠错码的代码距离。
    代码距离是衡量量子码纠错能力的重要指标。

    Args:
        H: 稳定子生成器矩阵（二进制矩阵），表示量子码的稳定子群
        logicOp: 逻辑算子矩阵（二进制矩阵），每一行代表一个逻辑算子

    Returns:
        d: 量子码的代码距离，即所有逻辑算子中最小的非平凡权重
    """""
    ##  格式化输入
    H = np.array(H, dtype=int)  # 转换为整数类型的numpy数组
    logicOp = np.array(logicOp, dtype=int)  # 转换为整数类型的numpy数组
    d = H.shape[1]  # 初始化距离为量子比特数量（最大可能距离）

    ##  遍历每个逻辑算子
    for i in range(logicOp.shape[0]):
        logicOp_i = logicOp[i, :]
        n = H.shape[1]  # 量子比特数量（稳定子矩阵的列数）
        m = H.shape[0]  # 稳定子数量（稳定子矩阵的行数）
        wstab = np.max([np.sum(H[i, :]) for i in range(m)])  # 计算最大稳定子权重（单个稳定子中非零元素的最大数量）
        wlog = np.count_nonzero(logicOp_i)  # 计算逻辑算子的权重
        num_anc_stab = int(np.ceil(np.log2(wstab)))  # 计算稳定子约束所需的辅助变量数量（基于最大稳定子权重的对数）
        num_anc_logical = int(np.ceil(np.log2(wlog)))  # 计算逻辑算子约束所需的辅助变量数量（基于逻辑算子权重的对数）
        num_var = n + m * num_anc_stab + num_anc_logical  # 总变量数量 = 量子比特变量 + 稳定子辅助变量 + 逻辑算子辅助变量

        ##  创建混合整数规划模型
        model = Model()
        model.verbose = 0  # 关闭详细输出
        x = [model.add_var(var_type=BINARY) for i in range(num_var)]  # 创建二进制变量数组
        model.objective = minimize(xsum(x[i] for i in range(n)))  # 目标函数：最小化前n个变量（量子比特变量）的和（即最小化Hamming权重）

        # 为每个稳定子添加正交性约束（模2）
        for row in range(m):
            weight = [0] * num_var  # 初始化权重向量
            supp = np.nonzero(H[row, :])[0]  # 获取当前稳定子的支持集（非零元素的位置）

            ##  设置qubit变量的权重为1
            for q in supp:
                weight[q] = 1

            ##  添加辅助变量来处理模2约束
            cnt = 1
            for q in range(num_anc_stab):
                weight[n + row * num_anc_stab + q] = -(1 << cnt)  # 设置辅助变量的权重为负的2的幂次方
                cnt += 1
            model += xsum(weight[i] * x[i] for i in range(num_var)) == 0  # 添加约束：权重向量与变量向量的点积等于0

        ##  添加逻辑算子的奇数重叠约束
        supp = np.nonzero(logicOp_i)[0]  # 获取逻辑算子的支持集
        weight = [0] * num_var  # 初始化权重向量

        ##  设置qubit变量的权重为1
        for q in supp:
            weight[q] = 1

        ##  添加辅助变量来处理模2约束
        cnt = 1
        for q in range(num_anc_logical):
            # 设置辅助变量的权重为负的2的幂次方
            weight[n + m * num_anc_stab + q] = -(1 << cnt)
            cnt += 1

        ##  添加约束：权重向量与变量向量的点积等于1（奇数重叠）
        model += xsum(weight[i] * x[i] for i in range(num_var)) == 1

        ##  求解优化问题，计算最优解中前n个变量的和，即最小Hamming权重
        model.optimize()
        # noinspection PyTypeChecker
        opt_val = sum([x[i].x for i in range(n)])
        d = min(d, int(opt_val))

    ##  返回最小距离
    return d


#%%  KEY：计算随机距离
def random_distance_caculator(gx, gz, num, mindist=2, debug=0, field=None):
    """"
    使用随机算法计算量子纠错码的代码距离
    该函数通过随机排列和行约简技术来寻找满足特定约束的最小权重向量，
    从而估计量子码的代码距离。适用于CSS码（Calderbank-Shor-Steane codes）。

    Args:
        gx: X稳定子生成器矩阵（二进制矩阵）
        gz: Z稳定子生成器矩阵（二进制矩阵）
        num: 随机迭代次数，控制算法的精度
        mindist: 最小距离阈值，如果找到的距离小于等于此值则提前终止
        debug: 调试标志位，控制调试信息的输出
        field: 有限域，默认为GF(2)（二进制域）
        maxav: 是否记录最大权重向量的详细信息

    Returns:
        int: 估计的代码距离
    """""
    ##  设置默认有限域为GF(2)（二进制域）
    if field is None:
        F = galois.GF(2)
    else:
        F = galois.GF(field)

    ##  调试模式检查：验证X和Z稳定子是否正交（CSS码的基本要求）
    if debug & 2:
        if np.any(gx @ gz.T):
            raise ValueError("Input matrices must have orthogonal rows")

    wz = F(gx.null_space())  # 计算X稳定子的零空间（Z逻辑算子空间）
    wx = F(gz.null_space())  # 计算Z稳定子的零空间（X逻辑算子空间）
    rows_wz, cols_wz = wz.shape  # 获取零空间矩阵的维度信息
    dist_bound = cols_wz + 1  # 初始化距离上界为最大可能值（列数+1）
    vec_count = 0  # 计数器：记录找到当前最小权重的向量数量

    ##  主循环：进行num次随机迭代
    for i in range(num):
        per = np.random.permutation(cols_wz)  # 生成随机排列，用于随机化搜索顺序
        wz1 = wz[:, per]  # 对零空间矩阵的列进行随机排列
        wz2 = wz1.row_reduce()  # 对排列后的矩阵进行行约简（高斯消元）
        wz2 = wz2[:, np.argsort(per)]  # 将列顺序恢复为原始顺序

        ##  遍历行约简后的每一行
        for j in range(rows_wz):
            temp_vec = wz2[j, :]  # 获取当前行向量
            temp_weight = np.count_nonzero(temp_vec)  # 计算向量的Hamming权重（非零元素个数）

            ##  检查权重是否在有效范围内且小于等于当前最小距离
            if 0 < temp_weight <= dist_bound:

                ##  检查向量是否与X逻辑算子空间有非零重叠（即是否为非平凡逻辑算子）
                if np.any(wx @ temp_vec):

                    ##  如果找到更小的权重，更新最小距离
                    if temp_weight < dist_bound:
                        dist_bound = temp_weight
                        vec_count = 1

                    ##  如果权重等于当前最小距离，增加计数器
                    elif temp_weight == dist_bound:
                        vec_count += 1

            ##  检查是否达到最小距离阈值，如果达到则提前终止
            if dist_bound <= mindist:
                return mindist

    ##  返回找到的最小距离
    return dist_bound
