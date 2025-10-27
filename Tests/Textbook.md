# ExtendedStim：A Python Library for Addressing Both Fermionic and Bosonic Quantum Error-Correction

---

## 1. Clifford method

### 1.1 Gottesman-Knill theorem

A quantum platform's operations are limited to: resets, Clifford gates, string-operator measurements and classical controls based on outputs of quantum platforms can be simulated with a classical computers efficiently under complexity about $O(n^3m)$. Here, $n$ is the number of qubits, and $m$ is the number of operations.

### 1.2 Pauli-Clifford algebra

用Pauli group可以表示一个Hilbert space中的一些特殊的态，我们将这些态称为stabilizer state，也就是在一组Pauli group群元中与Hilbert space中一个state建立了一个一一映射。Stabilizer state记为$|S\rangle$，$S$是group的群元形成的集合，要求如下。

1. Stabilizers满秩，$\text{rank}(S) = N$;
2. Stabilizers对易，$\forall\hat S_i,\hat S_j \in S:[\hat S_i,\hat S_j]=0$;
3. Stabilizers是Hermitian operators，$\forall S \in S,\hat S= \hat S^\dagger$;
4. Stabilizers位于Pauli group中，$S=\{\hat S_0,\cdots,\hat S_{N-1}\},\forall\hat S\in S:\hat S\in\mathcal{P}^N$.

Quantum state与stabilizer state的关系是stabilizer state是它所有stabilizers本征值为$1$的本征空间的交形成的一维子空间，因此代表了一个态。

$$
\begin{equation}
|S\rangle=\{|S\rangle\mid\forall\hat S\in S,\hat S|S\rangle=1|\S\rangle\}
\end{equation}
$$

定义一个quantum platform，含$N$个qubits，其允许的操作集合如下。

1. Initialization: $\forall|\psi_0\rangle\in\text{QP}:|\psi_0\rangle\in\{|S\rangle\}$，这里$\{|S\rangle\}$是所有stabilizer state的集合，出于方便我们一般假设初始化为$|0\cdots0\rangle$;
2. Gate set: $X, Y, Z, S, H, CX, CY, CZ, \text{SWAP}$;
3. Measurements: $\{\hat M\mid\hat M\in \mathcal{P}^N\land\hat M=\hat M^\dagger\}$.

Gate作用效果可用stabilizer mapping实现，Clifford gate将一个stabilizer state总是映射到另一个stabilizer state，因此计算复杂度不是指数增长的。

$$
\begin{equation}
\forall\hat U\in\text{Clifford}(\mathcal{P}^N):\hat U|S\rangle=|\{\hat U\hat S_0\hat U^\dagger,\cdots,\hat U\hat S_{N-1}\hat U^\dagger\}\rangle
\end{equation}
$$

我们只需要计算Clifford gate不对易的Pauli group子群的生成元在Clifford gate下如何变化即可，剩余的由于满足对易关系不会受Clifford gate的影响，在计算时只需要用对易关系跳过即可。

Clifford gate mapping relationship如下。

| Gate               | \hat X_0               | \hat Z_0    | \hat X_1           | \hat Z_1           |
|--------------------|------------------------|-------------|--------------------|--------------------|
| $\text{X}_0$       | $\hat X_0$             | $-\hat Z_0$ | $\hat X_1$         | $-Z_1$             |
| $\text{Y}_0$       | $-\hat X_0$            | $-\hat Z_0$ | $\hat X_1$         | $\hat Z_1$         |
| $\text{Z}_0$       | $-\hat X_0$            | $\hat Z_0$  | $\hat X_1$         | $\hat Z_1$         |
| $\text{H}_0$       | $\hat Z_0$             | $\hat X_0$  | $\hat X_1$         | $\hat Z_1$         |
| $\text{S}_0$       | $i\hat X_0\hat Z_0$    | $\hat Z_0$  | $\hat X_1$         | $\hat Z_1$         |
| $\text{CX}(0,1)$   | $\hat X_0\hat X_1$     | $\hat Z_0$  | $\hat X_1$         | $\hat Z_0\hat Z_1$ |
| $\text{CY}(0,1)$   | $i\hat X_0\hat X_1Z_1$ | $\hat Z_0$  | $\hat Z_0\hat X_1$ | $\hat Z_0\hat Z_1$ |
| $\text{CZ}(0,1)$   | $\hat X_0\hat Z_1$     | $\hat Z_0"  | $\hat Z_0\hat X_1$ | $\hat Z_1$         |
| $\text{SWAP}(0,1)$ | $\hat X_1$             | $\hat Z_1$  | $\hat X_0$         | $\hat Z_0$         |

举个例子$|\{\hat X_0\hat Z_1\hat X_2,\hat X_0\hat X_1\hat Z_2,\hat X_0\}\rangle$，作用Hadamard gate计算过程如下。

$$
\begin{equation}
\begin{aligned}
&\text{H}_0|\{\hat X_0\hat Z_1\hat X_2,\hat X_0\hat X_1\hat Z_2,\hat X_0\}\rangle=\\
&|\{\text{H}_0\hat X_0\hat Z_1\hat X_2\text{H}_0^\dagger,\text{H}_0\hat X_0\hat X_1\hat Z_2\text{H}_0^\dagger,\text{H}_0\hat X_0\text{H}_0^\dagger\}\rangle=\\
&|\{\text{H}_0\hat X_0\text{H}_0^\dagger\hat Z_1\hat X_2,\text{H}_0\hat X_0\text{H}_0^\dagger\hat X_1\hat Z_2,\text{H}_0\hat X_0\text{H}_0^\dagger\}\rangle=\\
&|\{\hat Z_0\hat Z_1\hat X_2,\hat Z_0\hat X_1\hat Z_2,\hat Z_0\}\rangle
\end{aligned}
\end{equation}
$$

作用$\text{CNOT}(0,2)$效果如下。

$$
\begin{equation}
\begin{aligned}
&\text{CNOT}(0,2)|\{\hat X_0\hat Z_1\hat X_2,\hat X_0\hat X_1\hat Z_2,\hat X_0\}\rangle=\\
&|\{\text{CNOT}(0,2)\hat X_0\hat Z_1\hat X_2\text{CNOT}(0,2)^\dagger,\text{CNOT}(0,2)\hat X_0\hat X_1\hat Z_2\text{CNOT}(0,2)^\dagger,\text{CNOT}(0,2)\hat X_0\text{CNOT}(0,2)^\dagger\}\rangle=\\
&|\{\text{CNOT}(0,2)\hat X_0\text{CNOT}(0,2)^\dagger\hat Z_1\text{CNOT}(0,2)\hat X_2\text{CNOT}(0,2)^\dagger,\text{CNOT}(0,2)\hat X_0\text{CNOT}(0,2)^\dagger\hat X_1\text{CNOT}(0,2)\hat Z_2\text{CNOT}(0,2)^\dagger,\text{CNOT}(0,2)\hat X_0\text{CNOT}(0,2)^\dagger\}\rangle=\\
&|\{\hat X_0\hat X_2\hat Z_1\hat X_2,\hat X_0\hat X_2\hat X_1\hat Z_0\hat Z_2,\hat X_0\hat X_2\}\rangle=\\
&|\{\hat X_0\hat Z_1,\hat X_0\hat Z_0\hat X_2\hat X_1\hat X_2\hat Z_2,\hat X_0\hat X_2\}\rangle
\end{aligned}
\end{equation}
$$

Measurements作用效果是对stabilizers的遍历与处理。

第一种情况是stabilizer state处于observable的本征态上，$\forall\hat S \in S:[\hat S,\hat M]=0$，那么有$\hat M=\pm\prod_{i=0}^{N-1}\hat S_i^{\eta_i}$，$\pm1$决定了测量结果。

第一种情况是stabilizer state不处于observable的本征态上，也就是$\exists\hat S\in S:[\hat S,\hat M]_+=0$，我们将这些stabilizers记为$\hat S_{j_0},\cdots,\hat S_{j_{K-1}}$，$|S\rangle$等价改写为$|\{\cdots,\hat S_{j_0},\hat S_{j_0}\hat S_{j_1},\cdots,\hat S_{j_0}\hat S_{j_{K-1}}\}\rangle$，测量结果为$|\{\cdots,\pm\hat M,\hat S_{j_0}\hat S_{j_1},\cdots,\hat S_{j_0}\hat S_{j_{K-1}}\}\rangle$，$\pm1$各有$50\%$几率生成，决定了测量结果与坍缩态。

重置与classical control是以上两种操作的classical combination，因此也可以用经典计算机高效模拟。

### 1.3 Majorana-Clifford algebra

基于Majorana group $\mathcal{M}^N$中stabilizers可以表示含$N$个fermionic sites的平台，与Pauli group情况类似，区别是群元的对易关系不同。

Majorana-Clifford gates有$\text{U},\text{V}, \text{N}, \text{FSWAP}, \text{BRAID}, \text{CN},\text{FS}, \text{TUNNEL}$。

| Gate                 | $\hat\gamma_0$                           | $\hat\gamma_0'$                          | $\hat\gamma_1$                           | $\hat\gamma_1'$                           |
|----------------------|------------------------------------------|------------------------------------------|------------------------------------------|-------------------------------------------|
| $\text{U}_0$         | $\hat\gamma_0$                           | $-\hat\gamma_0$                          | -$\hat\gamma_1$                          | $-\hat\gamma_1'$                          |
| $\text{V}_0$         | $-\hat\gamma_0$                          | $\hat\gamma_0$                           | $-\hat\gamma_1$                          | $-\hat\gamma_1'$                          |
| $\text{N}_0$         | $-\hat\gamma_0$                          | $-\hat\gamma_0$                          | $\hat\gamma_1$                           | $\hat\gamma_1'$                           |
| $\text{FS}_0$        | $\hat\gamma_0'$                          | $-\hat\gamma_0$                          | $\hat\gamma_0'$                          | $\hat\gamma_1'$                           |
| $\text{FSWAP}(0,1)$  | $\hat\gamma_1$                           | $\hat\gamma_1'$                          | $\hat\gamma_0$                           | $\hat\gamma_0'$                           |
| $\text{BRAID}(0,1)$  | $\hat\gamma_0$                           | $-\hat\gamma_1'$                         | $\hat\gamma_0'$                          | $\hat\gamma_1'$                           |
| $\text{CN}(0,1)$     | $i\hat\gamma_0\hat\gamma_1\hat\gamma_1'$ | $i\hat\gamma_0\hat\gamma_1\hat\gamma_1'$ | $i\hat\gamma_0\hat\gamma_0'\hat\gamma_1$ | $i\hat\gamma_0\hat\gamma_0'\hat\gamma_1'$ |
| $\text{TUNNEL}(0,1)$ | $-\gamma'_1$                             | $\gamma_1$                               | $-\gamma_0'$                             | $\gamma_0$                                |

Fermionic sites 上的测量与 qubits 情况完全相同。

### 1.4 Mixed algebra

一个既含fermionic sites又含qubits的quantum platform可以用$\mathcal{P}^{N_Q}\otimes\mathcal{M}^{N_F}$表示态。

Mixed platforms的额外的Clifford gate建立在fermionic sites and qubits之间，有$\text{CDX},\text{CVX},\text{CNX}\text{CNZ}$等。

| Gate       | $\hat\gamma$     | $\hat\gamma'$     | $\hat X$                     | $\hat Z$                     |
|------------|------------------|-------------------|------------------------------|------------------------------|
| \text{CUX} | \hat\gamma       | \hat\gamma'\hat X | \hat X                       | \hat\gamma\hat Z             |
| \text{CVX} | \hat\gamma\hat X | \hat\gamma'       | \hat X                       | \hat\gamma'\hat Z            |
| \text{CNX} | \hat\gamm\hat X  | \hat\gamma'\hat X | \hat X                       | i\hat\gamma\hat\gamma'\hat Z |
| \text{CNZ} | \hat\gamma\hat Z | \hat\gamma'\hat Z | i\hat\gamma\hat\gamma'\hat X | \hat Z                       |

---

## 2. Monte-Carlo method

### 2.1 Quantum error sampling

Mixed states等价于量子系综，因此如果非纯态过程只将stabilizer states投射到若干stabilizer states组成的系综上，我们还是可以用Clifford algebra的方式高效模拟这样的随机过程。

当输入mixed state是一个ensemble of stabilizer states时，Clifford unitary noise将它投射到另一组ensemble of stabilizer states上，这一过程等价于我们按概率将每个可能的态分别投射，形成另一组可能的态。一个基本的Clifford unitary noise是Pauli noise和Majorana noise，包括depolarization, dephasing, fermionic depolarization and fermionic dephasing等。

Measurement error的效果一般是影响feedback，因此也可以归纳到unitary noise中。

考虑最简单的一种模拟，输入态为pure stabilizer state，经过了一系列错误过程结果为$\hat\rho$，这个态可以用Monte-Carlo方法实现，也就是在计算过程中我们每一步unitary noise都以对应概率$p$作用对应的unitary，这样就完全是pure stabilizer state's evolution，最终的结果是一系列pure stabilizer state，而它们形成的系综就是最后的$\hat\rho$，若要计算保真度等可以直接用这种方式sampling得到。

### 2.2 Independent error model

考虑解码问题时，我们需要更强的error model，也就是independent error model，它把所有错误源视为独立发生的。一般而言Markovian noise都完全等价地映射到这样地模型上，但non-Markovian noise的情况则不清楚。

举例，depolarization是$\{\sqrt\frac{p}{3}\hat X, \sqrt\frac{p}{3}\hat Y, \sqrt\frac{p}{3}\hat Z\}$，这三种噪声并不是独立发生的，而是互斥的，但是我们可以将它拆成$X, Y, Z$三种error独立发生，每个error发生的概率为$\frac{9}{26}p$，这样error model与depolarization是等价的。

### 2.3 Tanner graph

Independent error model下的quantum circuit可以表示为一个Tanner graph。现在我们考虑一个理想线路$\text{QC}$和一个含某种噪声的线路$\text{QC}'$，$\text{QC}$具有纠错机制，也就是具有一系列的measurement output $\{m_i\}$；$\text{QC}'$则是$\text{QC}$加上了某种噪声，对应输出为$\{m'_i\}$。二者最终态对应的测量量记为$\{l_i\}$和$\{l'_i\}$

定义纠错机制中的比较项detectors为$\{d_i=m_{i_0}==m_{i_1}\}$，在无噪声线路中$\{d_i\}$应该是始终稳定的且都是$\text{True}$，没有随机性。噪声引入后每个sample的$\{d'_i\}$都可能出现错误，也就是一个alert，代表了某种error的发生。在独立错误模型下，错误是独立的，因此我们将每个错误描述为$\{e_i\}$，它们每个单独发生可能会产生的alert记为$\{(e_i,d_{i_0},\cdots,d_{i_{W-1}})\}$，我们还需要知道每个error是如何影响logical information的，与理想态相比是否翻转，也就是$\{(e_i,l_{i_0},\cdots,l_{i_{K-1}})\}$，这样依赖，我们得到了以error,logical information and detectors作为节点的无向图，每个error节点都连接logical information节点和detectors节点，logical information节点和detectors节点之间没有相连，这就是一个Tanner graph.

`Stim`遍历quantum circuit中的错误生成Tanner graph，即`stim.DetectorErrorModel`。算法流程是：取一个路径，假设只有其中一个error发生，观测detectors and logical information是否变化，生成对应edge信息。

Tanner graph在decoders中用来解码syndrome信息，如果解码出来的logical information与模拟出来的值一样，那么就说明这个error是被成功检测到的，记为error-free events，反之记为error，我们可以求出logical error rate出来。目前`stimbposd`可以根据Tanner graph用BP-OSD algorithm解码syndrome，不需要对Tanner graph做进一步的适配。

这样的Monte-Carlo方法能够模拟的过程是有限的，但是根据纠错条件定理和错误离散化定理，我们这样计算出来的logical error rate和真实值应该是同一量级甚至接近一样的。

## 3. Structure

### 3.1 Framework

ExtendedStim自底向上有四个levels.

1. Math level: `BinaryArray`，用于$\mathbb{F}_2$上的线性空间运算；
2. Physics level: `Operator`,`MajoranaOperator`,`PauliOperator`，表示物理上Majorana group和Pauli group及算符之间的关系；
3. Code level: `LinearCode`,`QuantumCode`,`QuantumCSSCode`,`MajoranaCode`,`MajoranaCSSCode`,`PauliCode`,`PauliCSSCode`，表示纠错码，用stabilizer表示；
4. Platform level: `Platform`，表示预设的平台结构，包括对平台state的刻画采用stabilizer state，相应的操作在平台内用stabilizer map或者measurement map实现，是合成纯数学的error model的物理基础；
5. Circuit level: `Circuit`，用于表示quantum circuit，可以生成`stim.DetectorErrorModel`，与`stim`的用法基本一致。

### 3.2 BinaryArray

#### 3.2.1 构造方法

`BinaryArray(self, array)`
> **输入**：`array->np.ndarray or array->list`
>
> **输出**：`BinaryArray`
>
> 构造方法，输入一个未格式化的二进制数组，输出一个格式化的对象，后续调用该对象可以不需要考虑格式问题。

#### 3.2.2 属性

1. `_array`：`np.ndarray of GF(2)`对象，用于存储二进制数组；
2. `shape`：`tuple`对象，数组形状。

#### 3.2.3 重载运算符

重载运算符改写运算符对对象的作用效果。

`__len__(self)`
> **输出**：`int`
>
> 返回数组的长度，即数组中元素的个数。

`__eq__(self, other)`
> **输入**：`other->BinaryArray`
>
> **输出**：`bool`
>
> 判断两个数组是否相等，返回一个布尔值。

`__getitem__(self, item)`
> **输入**：`item->int or item->slice`
>
> **输出**：`int`
>
> 返回数组中指定位置的元素。

`__setitem__(self, item, value)`
> **输入**：`item->int or item->slice`, `value->int`
>
> 将数组中指定位置的元素设置为给定值。

`__str__(self)`
> **输出**：`str`
>
> 返回数组的字符串表示，方便打印。

`__add__(self, other)`
> **输入**：`other->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 实现数组的$\mod 2$加法。

`__sub__(self, other)`
> **输入**：`other->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 实现数组的$\mod 2$减法。

`__mul__(self, other)`
> **输入**：`other->BinaryArray or other->int`
>
> **输出**：`BinaryArray`
>
> 实现数组的$Z_2$乘法。

`__rmul__(self, other)`
> **输入**：`other->BinaryArray or other->int`
>
> **输出**：`BinaryArray`
>
> 实现数组的右$Z_2$乘法。

`__matmul__(self, other)`
> **输入**：`other->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 实现数组的矩阵乘法。


`__rmatmul__(self, other)`
> **输入**：`other->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 实现数组的右矩阵乘法。

`__pow__(self, power)`
> **输入**：`int`
>
> **输出**：`BinaryArray`
>
> 实现数组的幂运算，主要是数组的矩阵幂运算。

#### 3.2.4 属性方法

属性方法即`@property`修饰的方法，用于获取属性值。

`T(self)`
> **输出**：`BinaryArray`
>
> 求数组的转置$A^T$.

`occupy(self)`
> **输出**：`list of int or list of tuple`
>
> 求数组中$1$的位置。

`weight(self)`
> **输出**：`int`
>
> 求数组中$1$的个数。

`null_space(self)`
> **输出**：`BinaryArray`
>
> 求数组的nullspace $\text{ker}(A)$，我们总是以行向量作为向量组中的向量。

`rank(self)`
> **输出**：`int`
>
> 求数组的秩。

#### 3.2.5 对象方法

对象方法是一般的成员方法，被对象调用，这里只有一个。

`copy(self)`
> **输出**：`BinaryArray`
>
> 复制数组。

#### 3.2.6 静态方法

静态方法是类所属的方法，不属于单个的对象，一般用于处理多个对象。

`sum(array_list)`
> **输入**：`array_list->list of BinaryArray`
>
> **输出**：`BinaryArray`
>
> 对多个数组进行数组间的$\mod 2$加法。

`hstack(array, *args)`
> **输入**：`array->BinaryArray`, `args->tuple of BinaryArray`
>
> **输出**：`BinaryArray`
>
> 水平拼接数组。

`vstack(array, *args)`
> **输入**：`array->BinaryArray`, `args->tuple of BinaryArray`
>
> **输出**：`BinaryArray`
>
> 垂直拼接数组。

`solve(matrix, vector)`
> **输入**：`matrix->BinaryArray`, `vector->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 解线性方程组$(x_0,...,x_{n-1})A=(b_0,...,b_{n-1})$，返回一个新的`BinaryArray`对象$x$，相当于对矩阵的每一行求组合方式得到`b`;

`cap(matrix1, matrix2)`
> **输入**：`matrix1->BinaryArray`, `matrix2->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 把矩阵的每一行作为基矢，求两个基矢组的交空间。

`minus(matrix1, matrix2)`
> **输入**：`matrix1->BinaryArray`, `matrix2->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 把矩阵的每一行作为基矢，求两个基矢组的差空间。

`direct_sum(matrix1, matrix2)`
> **输入**：`matrix1->BinaryArray`, `matrix2->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 把矩阵的每一行作为基矢，求两个基矢组的直和。

`orthogonalize(matrix)`
> **输入**：`matrix->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 把矩阵的每一行作为基矢，对矩阵进行正交化，返回一个新的`BinaryArray`对象，其中$B[i]$都是对应的基矢；

`kron(matrix1, matrix2)`
> **输入**：`matrix1->BinaryArray`, `matrix2->BinaryArray`
>
> **输出**：`BinaryArray`
>
> 把矩阵的每一行作为向，量求矩阵的张量积，返回一个新的`BinaryArray`对象，其中$B[i]$都是对应的向量；

`distance(H, method)`
> **输入**：`H->BinaryArray`, `method->str`
>
> **输出**：`int`
>
> 把矩阵视为code的check matrix，求$\text{ker}(H)-\text{im}(H)$中weight最小的'l的weight。计算方法可选值为`'randmo'`和`'mip'`，分别对应随机求解和最小正数规划精确求解。

`subsystem_distance(stabilizers, gauges, method)`
> **输入**：`stabilizers->BinaryArray`, `gauges->BinaryArray`, `method->str`
>
> **输出**：`int`
>
> 把矩阵视为subsystem code的check matrix，求$\text{ker}(H)-\text{im}(H)-\text{G}$中weight最小的'l的weight。计算方法可选值为`'randmo`和`'mip'`，分别对应随机求解和最小正数规划精确求解。

`FromOccupy(occupy, *args)`
> **输入**：`occupy->list of int`, `args->tuple of int`
>
> **输出**：`BinaryArray`
>
> 根据$1$的位置构造数组，`occupy`是`list of int`对象，`*args`是其他`int`对象，返回一个新的`BinaryArray`对象；

`FromArray(array)`
> **输入**：`array->np.ndarray or array->list`
>
> **输出**：`BinaryArray`
>
> 根据数组构造`BinaryArray`格式化对象。

`zeros(shape)`
> **输入**：`shape->tuple of int`
>
> **输出**：`BinaryArray`
>
> 根据形状构造一个全$0$数组。

`eye(number)`
> **输入**：`number->int`
>
> **输出**：`BinaryArray`
>
> 构造一个单位矩阵。

`shift(number, shift)`
> **输入**：`number->int`, `shift->int`
>
> **输出**：`BinaryArray`
>
> 构造一个循环左移矩阵。

`ones(shape)`
> **输入**：`shape->tuple of int`
>
> **输出**：`BinaryArray`
>
> 根据形状构造一个全$1$数组。
> 

#### 3.3 Graph

`Graph`模块提供图相关的工具。常用函数应包括图的连通性、独立环查找等，用于逻辑测量关系的图论分析。

### 3.4 Operator

#### 3.4.1 构造方法

`Operator(self, occupy_x, occupy_z, coff)`
> 输入：`occupy_x->np.ndarray or list`，`occupy_z->np.ndarray or list`，`coff->complex or int`
>
> 输出：`Operator`（抽象类实例通常由子类构造）
>
> 保存占据索引并做类型检查、排序。

#### 3.4.2 属性

1. `occupy_x`：`np.ndarray of int`对象，算符在X上的占据；
2. `occupy_z`：`np.ndarray of int`对象，算符在Z上的占据；
3. `coff`：`1,-1,1j or -1j`，系数

#### 3.4.3 重载运算符

`__matmul__(self, other)`
> 输入：`other->Operator`（子类实例）
>
> 输出：`Operator`（乘积）
>
> 说明：子类需实现具体代数规则（Pauli/Majorana）。

`__rmatmul__(self, other)`
> 输入：`other->Operator`
>
> 输出：`Operator`

`__mul__(self, other)` / `__rmul__(self, other)`
> 输入：`other->int or complex`
>
> 输出：`Operator`


`__str__(self)` / `__eq__(self, other)` / `__neg__(self)`
> 说明：字符串表示、相等判断、取负运算由子类实现。

####  3.4.4属性方法

`weight(self) -> int`
> 输出：`int`对象。
>
> 输出string算符中含有算符的数目。

`type(self)`（抽象）
> 输出：`str`，表示算符类型（例如 'MajoranaOperator' 或 'PauliOperator'）。

`is_hermitian(self) -> bool`（抽象）
> 输出：布尔值，判断算符是否是厄米算符。

`dual(self)`（抽象）
> 输出：`Operator`，算符的对偶（一般互换 x/z 占据）。

#### 3.4.5 对象方法

`index_map(self, index)`
> 输入：`index->np.ndarray or list`
>
> 输出：`tuple (x, z, coff)` 或子类实例（子类可能重载以返回子类对象）
>
> 说明：将占据索引按给定映射重新编码，返回映射后的索引与系数。

`get_vector(self, number)`
> 输入：`number->int`（物理位数）
>
> 输出：`BinaryArray`，长度为 `2*number` 的二进制向量，偶数位对应 X，奇数位对应 Z。

`is_exist_occupy_x(self, index)` / `is_exist_occupy_z(self, index)`
> 输入：`index->int`
>
> 输出：`int or None`（若存在返回索引位置，否则 None）。

`copy(self)`（抽象）
> 输出：对象深拷贝（子类实现）。

`split(self, index)`
> 输入：`index->int`
>
> 输出：`left_x, left_z, middle_x, middle_z, right_x, right_z`（基础拆分；子类可返回子类实例）

#### 3.4.6 静态方法

`get_matrix(ops, number)`
> 输入：`ops->list of Operator`，`number->int`
>
> 输出：`BinaryArray` （每行对应一个算符的向量表示）。

`HermitianOperatorFromOccupy(occupy_x, occupy_z)` / `HermitianOperatorFromVector(vector)`（抽象）
> 输出：构造的厄米算符对象（子类实现）。

`commute(A, B)`（抽象）
> 输入：`A,B->Operator`
>
> 输出：`bool`，判断对易性。

### 3.5 MajoranaOperator

#### 3.5.1 构造方法

`MajoranaOperator(self, occupy_x, occupy_z, coff)`
> 输入：`occupy_x->list or np.ndarray`，`occupy_z->list or np.ndarray`，`coff->complex or int`
>
> 输出：`MajoranaOperator`
>
> 说明：将占据向量排序，并保持 coff。用于表示 Majorana 群的乘积项。

#### 3.5.2 重载运算符

`__matmul__(self, other)`
> 输入：`other->MajoranaOperator`
>
> 输出：`MajoranaOperator`
>
> 说明：采用 setxor 计算占据结果，并统计按照扩展索引（x->2*i, z->2*i+1）排列时与另一算符的交换次数来确定附加 (-1) 相位。

`__rmatmul__(self, other)` / `__mul__` / `__rmul__` / `__str__` / `__eq__` / `__neg__`
> 与 PauliOperator 类似，由类直接实现，遵循 Majorana 代数规则。

#### 3.5.3 属性方法

`type` -> "MajoranaOperator"

`is_hermitian(self)`
> 输出：`bool`，基于算符的 weight 与组合的三角数奇偶性判断 coff 是否应为 1/-1 或 ±i。

`dual(self)`
> 输出：MajoranaOperator（交换 occupy_x 与 occupy_z）。

#### 3.5.3 对象方法

`split(self, index)`
> 输出：三个部分（left, middle, right），子类实现返回 MajoranaOperator 的三段。

`index_map(self, index)`
> 输出：MajoranaOperator（映射后）。

`copy(self)`
> 输出：MajoranaOperator（深拷贝）。

#### 3.5.4 静态方法

`HermitianOperatorFromOccupy(occupy_x, occupy_z)` / `HermitianOperatorFromVector(vector)`
> 输出：Majornana 算符（厄米形式），自动确定 coff 的相位为 1 或 ±i。

`commute(A, B)`
> 输出：`bool`，判断是否对易（基于占据重叠与权重计算）。

### 3.6 PauliOperator

#### 3.6.1 构造方法

`PauliOperator(self, occupy_x, occupy_z, coff)`
> 输入：`occupy_x->list or np.ndarray`（X 占据），`occupy_z->list or np.ndarray`（Z 占据），`coff->complex/int`
>
> 输出：`PauliOperator`

#### 3.6.2 重载运算符
`__matmul__(self, other)`
> 输入：`other->PauliOperator`
>
> 输出：`PauliOperator`
>
> 说明：占据使用对称差，额外相位由 A.z 与 B.x 的交叉次数产生 (-1)；返回新的 PauliOperator。

`__rmatmul__`, `__mul__`, `__rmul__`, `__str__`, `__eq__`, `__neg__` 等由类实现，满足 Pauli 代数。

#### 3.6.3 属性方法
`type` -> "PauliOperator"

`is_hermitian(self)`
> 输出：`bool`，若在同位上同时出现 X 和 Z 且个数为奇数，则 coff 应为 ±i，否则 ±1。

`dual(self)` 返回 PauliOperator（x/z 互换）。

#### 3.6.3 静态方法
`HermitianOperatorFromOccupy` / `HermitianOperatorFromVector`
> 用于从占据或向量生成厄米 Pauli 算符。

`commute(A,B)`
> 输出：`bool`，基于 X/Z 的交叉计算对易性。

### 3.7 LinearCode

#### 3.7.1 构造方法
`LinearCode(self, check_matrix)`
> 输入：`check_matrix->2D array-like`（校验矩阵）
>
> 输出：`LinearCode` 对象（内部将校验矩阵转为 `BinaryArray`）

#### 3.7.2 属性方法
`rank` -> 校验矩阵秩（使用 `BinaryArray.rank`）

`distance` -> 当前为占位（返回 1），将来可由具体算法实现。

`logical_number` -> 信息位数 = `number_bit - rank`

`codewords` -> `BinaryArray.null_space`（码字集合）

`dual` -> 返回 `LinearCode(self.check_matrix.null_space)`

`is_dual_containing` -> 布尔值（检查 H @ H.T 是否全 0）

### 3.8 BicycleCode (extendedstim.Code.LinearCode.BicycleCode)

#### 3.8.1 构造方法
`BicycleCode(self, N, k, M, seed)`
> 输入：`N,k,M,seed->int`
>
> 输出：`BicycleCode`（继承自 `LinearCode`）

> 说明：生成循环矩阵并平衡列权以返回 `LinearCode`；构造过程中会检查 `H @ H.T == 0`。

### 3.9 QuantumCode

#### 3.9.1 构造方法

`QuantumCode(self, generators, physical_number)`
> 输入：`generators->list or np.ndarray of Operator`，`physical_number->int`
>
> 输出：`QuantumCode`（抽象基类实例由子类生成）

#### 3.9.2 属性方法

`check_matrix` -> 通过 `Operator.get_matrix(self.generators, self.physical_number)` 得到 `BinaryArray` 校验矩阵。

`rank` / `logical_number` 等同前述定义。

抽象属性：`distance`、`logical_operators`

#### 3.9.3 对象方法

`index_map(self, index_map)` -> 映射所有生成器的物理索引

`copy(self)` -> 返回深拷贝

### 3.10 QuantumCSSCode (extendedstim.Code.QuantumCode.QuantumCSSCode)

#### 3.10.1 构造方法

`QuantumCSSCode(self, generators_x, generators_z, physical_number)`
> 输入：`generators_x->list of Operator`, `generators_z->list of Operator`, `physical_number->int`
>
> 输出：`QuantumCSSCode`

#### 3.10.2 属性/抽象属性

`generators_x`, `generators_z`, `check_matrix_x`, `check_matrix_z`, `rank_x`, `rank_z`

抽象：`logical_operators_x`, `logical_operators_z`, `distance_x`, `distance_z`

### 3.11 PauliCode

#### 3.11.1 构造方法
`PauliCode(self, generators, physical_number)`
> 输入：`generators->list of MajoranaOperator-like objects`, `physical_number->int`
>
> 输出：`PauliCode`

#### 3.11.2 静态方法
`FromCheckMatrix(check_matrix)`
> 输入：`check_matrix->BinaryArray-like or 2D array`（按 Pauli 的 2 列为一位编码）
>
> 输出：`PauliCode`（从行构造 MajoranaOperator 并返回实例）

### 3.12 PauliCSSCode

#### 3.12.1 构造方法
`PauliCSSCode(self, generators_x, generators_z, physical_number)`
> 输入：分别为 X/Z 方向的生成器列表；输出为 `PauliCSSCode` 对象。

#### 3.12.2 属性/方法
继承自 `PauliCode` 与 `QuantumCSSCode`，但许多具体方法（distance, logical_operators）为 TODO。

### 3.13 MajoranaCode (extendedstim.Code.QuantumCode.MajoranaCode)

#### 3.13.1 构造方法
`MajoranaCode(self, generators, physical_number)`
> 输入：`generators->list of MajoranaOperator`，`physical_number->int`
>
> 输出：`MajoranaCode`

#### 3.13.2 属性
`distance` -> 调用 `BinaryArray.distance(self.check_matrix, 'mip')`（MIP 求解，需要 `mip` 依赖）

`logical_operators` -> 通过校验矩阵的零空间提取线性无关基并映射为 `MajoranaOperator` 列表返回

`even_or_odd` -> 依据 H * ones 是否有解判断码的偶/奇性

`FromCheckMatrix` -> 从二进制校验矩阵构造 MajoranaCode

### 3.14 MajoranaCSSCode

#### 3.14.1 构造方法

`MajoranaCSSCode(self, generators_x, generators_z, physical_number)`
> 输入：X/Z 方向生成器列表及物理位数
>
> 输出：`MajoranaCSSCode`

#### 3.14.2 属性方法

`distance_x`, `distance_z` -> 分别调用 `BinaryArray.distance(...,'mip')`

`logical_operators_x`, `logical_operators_z` -> 从对应方向 null_space 提取并映射为 MajoranaOperator 列表

### 3.15 Platform (detailed API)

下面严格按照 BinaryArray 部分的格式，逐条列出 `extendedstim.Platform.Platform` 的每个属性与方法（输入/输出/影响/物理含义），基于仓内 `Platform.py` 源码逐项展开，确保没有遗漏。

#### 属性（构造后可读）
`pauli_number`
> **类型**：int
>
> **含义**：当前平台中 classical qubit（Pauli qubit）的数目。
>
> **初始值**：构造函数 `__init__` 中默认 0，调用 `initialize` 后由用户设置或自动填充。

`majorana_number`
> **类型**：int
>
> **含义**：当前平台中 fermionic sites（Majorana pairs 计作 1 site）的数目。
>
> **初始值**：构造函数中默认 0，调用 `initialize` 后设置。

`stabilizers_pauli`
> **类型**：list of `PauliOperator`
>
> **含义**：当前平台维护的一组 Pauli 稳定子（表示平台当前的 stabilizer group 中与 qubit 相关的生成元）；用于测量、重置与门操作的稳定子追踪和坍缩处理。

`stabilizers_majorana`
> **类型**：list of `MajoranaOperator`
>
> **含义**：当前平台维护的一组 Majorana 稳定子（表示 fermionic 站点相关的生成元）；用于 fermionic 测量、门与重置的稳定子追踪。

#### 构造方法
`Platform(self)`
> **输入**：无
>
> **输出**：`Platform` 对象
>
> **效果**：创建一个空平台对象并初始化内部计数与稳定子列表（`pauli_number=0`,`majorana_number=0`,`stabilizers_pauli=[]`,`stabilizers_majorana=[]`）。

#### 对象方法（按源码逐项详述）
`initialize(self, majorana_number, pauli_number, *args)`
> **输入**：
> - `majorana_number->int`：期望的 fermionic 站点数目；
> - `pauli_number->int`：期望的 qubit 数量；
> - `*args`（可选）：若提供两项，视为用户给定的稳定子列表 `(stabilizers_majorana_list, stabilizers_pauli_list)`。
>
> **输出**：None
>
> **效果**：
> - 当 `*args` 为空：按照默认策略为每个 Majorana 站点加入一个局部的 Majorana 稳定子（`MajoranaOperator([i],[i],1j)`）——代表站点上的 i*gamma*gamma'；为每个 qubit 加入对应的 Pauli Z 稳定子（`PauliOperator([], [i], 1)`），同时在 `stabilizers_majorana` / `stabilizers_pauli` 中补齐空白；并更新 `majorana_number` 与 `pauli_number`；
> - 当 `*args` 提供两个列表：直接使用用户提供的稳定子列表替代默认的初始化方式。
>
> **物理含义**：定义平台的物理 Hilbert 空间初始稳定子群；选择不同的初始稳定子相当于指定了平台的编码或占据基态。

`measure(self, op)`
> **输入**：`op->PauliOperator or MajoranaOperator`（要求 `op.is_hermitian` 为 True）
>
> **输出**：`int`，测量结果 `+1` 或 `-1`
>
> **效果（算法实现细节）**：
> 1. 根据 `op` 的类型选取 `stabilizers_now`（Pauli 类型用 `stabilizers_pauli`，Majorana 类型用 `stabilizers_majorana`）。
> 2. 遍历当前稳定子，检查与 `op` 的对易性：若有不对易的稳定子，使用高斯消去风格的替换策略（在源码中：找到首个不对易的稳定子作为 `first_pauli` 并记录索引；若存在第二个不对易稳定子则用 `stabilizers_now[i] = stabilizers_now[i] @ first_pauli` 消去）——这相当于把不对易的稳定子转换为一组互相对易的基并处理坍缩方向；
> 3. 若最终没有不对易稳定子（`first_index == -1`）：表示被测算符在当前稳定子生成群可由稳定子组合表示，代码通过把稳定子矩阵拼接为 `BinaryArray` 矩阵并解线性方程 `matrix * x = vector(op)` 来找出用于表示 `op` 的稳定子组合；根据组合算符的系数与 `op.coff` 比较，返回测量结果 `+1` 或 `-1`（对应确定性的测量结果）；
> 4. 若存在不对易的稳定子：表示被测算符不在稳定子本征态上，源码实现中以 50% 概率把第一个不对易稳定子替换为 `op` 或 `-op`（对应测量结果分别为 `+1` 或 `-1`），并返回相应结果；
>
> **物理含义**：实现了对 stabilizer 模型的测量：若被测算符可由稳定子组合表示则测量决定性；若不对易则产生随机坍缩，并通过替换/重构稳定子集合反映测量的后态。

`X(self, qubit_index: int)`
> **输入**：`qubit_index->int`
>
> **输出**：None
>
> **效果**：遍历 `stabilizers_pauli`，若某稳定子在 `occupy_z` 上占据了该 `qubit_index`（即该稳定子包含 Z 在该位），则取反该稳定子的系数 `coff = -coff`。
>
> **物理含义**：在稳定子表示中，施加 X 门会把 Z 的符号翻转（相位翻转），对应经典稳定子系数的改变；这是对 stabilizer group 的等价变换表现。

`Y(self, qubit_index: int)`
> **输入**：`qubit_index->int`
>
> **输出**：None
>
> **效果**：遍历 `stabilizers_pauli`，若稳定子在 `occupy_x` 或 `occupy_z` 包含该位，则对其系数取反（在源码中分别做了两次取反以对应 Y 的 X 与 Z 同时作用导致的相位）。
>
> **物理含义**：Y = i X Z 的复合作用在稳定子系数上表现为相位改变（可能产生 ±1 或 ±i 的变化），源码以符号翻转近似实现了该代数在稳定子表示下的映射。

`Z(self, qubit_index: int)`
> **输入**：`qubit_index->int`
>
> **输出**：None
>
> **效果**：遍历 `stabilizers_pauli`，若某稳定子在 `occupy_x` 上包含该位，则对其系数取反。
>
> **物理含义**：Z 门将 X 项的相位翻转，稳定子表示中体现为相应生成元系数的反号。

`H(self, qubit_index: int)`
> **输入**：`qubit_index->int`
>
> **输出**：None
>
> **效果（源码流程）**：
> - 对 `stabilizers_pauli` 中每个稳定子 `stabilizer` 调用 `split(qubit_index)` 将其分成 left/middle/right 部分；
> - 根据 `middle` 部分的占据情况（X,Z 在该位是否存在）执行映射：
>   - 如果 middle 为仅 Z 占据（len(middle.occupy_x)==0 && len(middle.occupy_z)==1）: 将其转换为 X 占据（对应 H: Z->X）；
>   - 若 middle 为仅 X 占据：将其转换为 Z 占据（X->Z）；
>   - 若 middle 同时包含 X 和 Z 或均不包含：源码处理为跳过或抛错（具体逻辑见实现）。
>
> **物理含义**：Hadamard 在 Pauli 基之间交换 X 与 Z 分量，稳定子在局部重写为等价生成元体现了这种映射。

`U(self, majorana_index: int)`
> **输入**：`majorana_index->int`
>
> **输出**：None
>
> **效果**：构造 `op = MajoranaOperator([majorana_index], [], 1)`（对应 γ_i），遍历 `stabilizers_majorana`，对与 `op` 不对易的稳定子取反系数。
>
> **物理含义**：实现 Majorana 平台上的单站点γ 操作，其对稳定子群的影响等价于相应生成元相位翻转。

`V(self, majorana_index: int)`
> **输入**：`majorana_index->int`
>
> **输出**：None
>
> **效果**：构造 `op = MajoranaOperator([], [majorana_index], 1)`（对应 γ'_i），对不对易稳定子取反系数。
>
> **物理含义**：实现 Majorana 上的另一类单体操作（γ'），对稳定子相位进行对应翻转。

`N(self, majorana_index: int)`
> **输入**：`majorana_index->int`
>
> **输出**：None
>
> **效果**：构造 `op = MajoranaOperator([majorana_index], [majorana_index], 1j)`（相当于 i γ_i γ'_i），对不对易稳定子取反系数。
>
> **物理含义**：实现费米子站点上的占据数或局域相位操作（iγγ'），这在费米子代数中是一个常见的局域对角算符。

`P(self, majorana_index: int)`
> **输入**：`majorana_index->int`
>
> **输出**：None
>
> **效果**：源码中该方法为 `pass`（未实现）。
>
> **物理含义**（预期）：通常用作 fermionic phase gate 的占位；需要进一步实现以反映特定费米子相位变换。

`S(self, pauli_index: int)`
> **输入**：`pauli_index->int`
>
> **输出**：None
>
> **效果**：遍历 `stabilizers_pauli`，将每个稳定子在 `pauli_index` 位的 middle 部分进行以下映射：
> - 若 middle 没有 X 占据：跳过；
> - 若 middle 只有 X：将其扩展为 X 和 Z 占据，并将中间系数设置为 `1j`；
> - 若 middle 同时含 X 与 Z：移除 Z 并将系数置为 `1j`（跟随实现细节）。
>
> **物理含义**：S 门在 Pauli 群中将 X -> Y（含相位 i），该方法在稳定子表示中通过修改占据与系数实现该映射。

`CX(self, control_index, target_index)`
> **输入**：`control_index->int`, `target_index->int`（两者均为 qubit 索引）
>
> **输出**：None
>
> **效果（源码流程）**：对每个 Pauli 稳定子：
> - 使用 `split(control_index)`、`split(target_index)` 获取左右/中间段；
> - 若控制位 middle 含 X，则 middle 替换为 X on control 与 X on target（即控制 X 扩展为目标 X）；
> - 同时处理目标位处 Z 占据对控制位的影响（通过在左右片段插入对应的 Z 组合）；
> - 最终将稳定子更新为 left @ middle @ ... 的组合。
>
> **物理含义**：在稳定子表示上实现 CNOT（CX）对生成元的代数变换；这是 Clifford 门族的核心映射之一。

`CNX(self, control_index, target_index)`
> **输入**：`control_index->int`（Majorana 站点），`target_index->int`（Pauli qubit）
>
> **输出**：None
>
> **效果**：遍历索引 i 的稳定子对（Pauli 与 Majorana 对应位置），把控制 Majorana 的 X/Z 占据映射为目标位的 Pauli 操作（并同时对 Majorana 稳定子加入由目标位产生的 Majorana 因子）；实现细节见源码。
>
> **物理含义**：实现一个混合控制门（fermionic control -> qubit target），用于 fermion-qubit 交互的 Clifford 操作。

`CUX(self, control_index, target_index)`
> **输入**：`control_index->int`（Majorana 站点），`target_index->int`（Pauli qubit）
>
> **输出**：None
>
> **效果**：类似 `CNX` 但映射规则不同（依据控制位的 z 成分进行目标 X 的生成等），参见源码。

`CVX(self, control_index, target_index)`
> **输入**：`control_index->int`（Majorana 站点），`target_index->int`（Pauli qubit）
>
> **输出**：None
>
> **效果**：另一个 Majorana->Pauli 的混合控制门映射，根据控制位的 x 成分映射到目标 Pauli，并在 Majorana 稳定子上加入相应的 Majorana 因子。

`x_error(self, pauli_index, p)` / `y_error(self, pauli_index, p)` / `z_error(self, pauli_index, p)`
> **输入**：`pauli_index->int`, `p->float`（概率，0<=p<=1）
>
> **输出**：None
>
> **效果**：以概率 `p` 调用对应的门（`X`/`Y`/`Z`）来模拟随机 Pauli 错误。
>
> **物理含义**：实现独立错误模型里 qubit 层次的随机错位采样。

`u_error(self, majorana_index, p)` / `v_error(self, majorana_index, p)` / `n_error(self, majorana_index, p)`
> **输入**：`majorana_index->int`, `p->float`
>
> **输出**：None
>
> **效果**：以概率 `p` 分别调用 `U`/`V`/`N` 来模拟 Majorana 层次的随机错误。

`clear(self, op)`
> **输入**：`op->PauliOperator or MajoranaOperator`（厄米）
>
> **输出**：None
>
> **效果（源码流程）**：与 `measure` 类似，但用于把系统重置到 `op` 的 +1 本征子空间：
> - 检查对易性以决定是否需要做稳定子替换或线性求解；
> - 若 `op` 可由现有稳定子线性组合表示，则找到参与该线性组合的一个稳定子索引并对其系数取反以实现重置（保持其余稳定子不变）；
> - 若存在不对易稳定子，则直接用 `op` 替换首个不对易稳定子（对应坍缩到该本征子之一）。
>
> **物理含义**：实现测量并将系统投射/重置到指定的算符本征态（常用于初始化或测量+反馈回路）。

`reset(self, pauli_index)`
> **输入**：`pauli_index->int`
>
> **输出**：None
>
> **效果**：构造 `PauliOperator([], [pauli_index], 1)`（Z_i）并调用 `clear(op)` 将该 qubit 重置为 |0> 对应的稳定子（+1 本征子）。
>
> **物理含义**：在稳定子表述下将指定 qubit 重置为 |0> 态。

`fermionic_reset(self, majorana_index)`
> **输入**：`majorana_index->int`
>
> **输出**：None
>
> **效果**：构造 Majorana 的局域算符 `MajoranaOperator([i],[i],1j)` 并调用 `clear(op)`，从而把该 fermionic 站点投射到空/本征态（+1/−1 由实现约定）。

#### 例外与边界情况
- 所有操作均假设输入索引在当前 `pauli_number` / `majorana_number` 的范围内；否则可能抛出 IndexError 或导致不一致的稳定子映射；
- `measure` 在存在不对易稳定子的情况下使用伪随机（源码为 `np.random.rand() < 0.5`）决定坍缩方向；如果需要确定性复现，应设置随机种子或改为确定性策略；
- `P` 方法为占位，需要实现以完整支持 fermionic phase gate；

### 3.16 Circuit (完整接口规范 — 明确标注为仓库中未实现，下面为精确接口规范/实现建议)

注意：仓库中当前未包含 `extendedstim.Circuit.Circuit` 的实现文件（我已搜索仓内源码）；因此下面的文档为严格、逐项的接口规范（API 设计），以完全满足你要求“每个属性、每个方法都要分开具体地写、没有遗漏”。如果你希望我在仓库中实现该接口，我可以接着基于下面的规范实现代码并将其加入仓库。

下面的接口设计参照 `Stim`/`ExtendedStim` 的用途（生成 Tanner graph、构建测量序列、支持混合 fermion/qubit 门、输出 `stim.DetectorErrorModel` 等），并以非常详细的方式列出每一个属性与方法的输入/输出/副作用/物理含义。

#### 假设与总体说明
- 假设 `Circuit` 需要表达一条量子电路（含 qubits 与 fermionic sites），记录一系列时序化操作（gates，measurements，resets，error insertions，classical control），并能导出以 `stim.DetectorErrorModel` 形式的检测器/错误模型；
- 假设 `Circuit` 与 `Platform` 使用相兼容的数据结构（`PauliOperator` / `MajoranaOperator` / `BinaryArray`），并能接收 `Platform` 的索引映射以生成对应的稳定子变化；

#### 属性（存储电路描述与派生表征）
`number_qubits`
> **类型**：int
>
> **含义**：电路中 qubit 的数量（与 `Platform.pauli_number` 对齐）。

`number_majorana`
> **类型**：int
>
> **含义**：电路中 Majorana 站点的数量（与 `Platform.majorana_number` 对齐）。

`instructions`
> **类型**：list of dict/object
>
> **含义**：按时间顺序存放电路中的所有操作（每条操作为一个结构体，包含字段例如 `type`（'X','H','CX','MEASURE','RESET','CNX' 等）、`targets`（索引或索引对）、`params`（概率、标号等）、`condition`（可选的经典条件）等）。
>
> **格式建议**：每个 instruction 至少包含 `{'type': str, 'targets': list[int], 'ancilla': optional, 'prob': optional}`。

`detectors`
> **类型**：list
>
> **含义**：从电路中派生出的 detector（检测器）列表，每个 detector 描述一组 measurement outputs 的比较关系（例如 m_i == m_j）；用于构造 `stim.DetectorErrorModel`。

`logical_observables`
> **类型**：list
>
> **含义**：电路中被跟踪的逻辑产出（logical measurement），例如用于计算 logical error 的测量算符或比特位映射。

`errorspecs`
> **类型**：list
>
> **含义**：定义独立错误模型中各操作对应的错误节点（每个 entry 包含错误类型、作用的位置、发生概率以及其产生的 detector/ logical 翻转影响）。

#### 构造方法
`Circuit(self, number_qubits=0, number_majorana=0)`
> **输入**：`number_qubits->int`（可选），`number_majorana->int`（可选）
>
> **输出**：`Circuit` 对象
>
> **效果**：初始化内部结构（`instructions=[]`,`detectors=[]`,`logical_observables=[]`,`errorspecs=[]`）并设置 qubit / majorana 的初始计数。

#### 核心方法（逐个列出）
`add_qubit(self, count=1)`
> **输入**：`count->int`
>
> **输出**：None（更新 `number_qubits`）
>
> **效果**：为电路增加 `count` 个 qubit（更新 `number_qubits` 并返回新索引范围以便后续引用）。

`add_majorana(self, count=1)`
> **输入**：`count->int`
>
> **输出**：None（更新 `number_majorana`）
>
> **效果**：为电路增加 `count` 个 Majorana 站点（更新计数并返回新索引范围）。

`append_instruction(self, instr)`
> **输入**：`instr->dict or object`（详见 `instructions` 结构说明）
>
> **输出**：None（把该 instruction 添加到 `instructions` 列表末尾）
>
> **效果**：以顺序方式构建电路时间轴，`instr` 可表示门/测量/重置/错误注入/条件控制等。

`add_gate(self, gate_type, targets, params=None)`
> **输入**：
> - `gate_type->str`（如 'X','Y','Z','H','S','CX','U','V','N','CNX','CUX','CVX' 等），
> - `targets->list[int]`（门作用的索引，Majorana/qubit 索引需与 `gate_type` 对应），
> - `params->dict`（可选，包含 gate 特有参数）。
>
> **输出**：None（在 `instructions` 中追加对应的 gate instruction）
>
> **效果**：添加一个门操作到电路，并在需要时为该门生成对应的错误节点（放入 `errorspecs`，例如随机 Pauli 错误模型）。

`add_measure(self, op, targets, key=None)`
> **输入**：
> - `op->str or Operator`（若为字符串则为 'Z'/'X' 等基础测量，若为 Operator 则为任意可测算符），
> - `targets->list[int]`（作用索引），
> - `key->str or int`（可选，测量结果的标签，用于构造 detector）。
>
> **输出**：`measurement_id`（用于后续引用）
>
> **效果**：在 `instructions` 中加入测量指令并在 `detectors` 或 `logical_observables` 中注册相应条目（若测量用于构建 detector 比较则把其 output 与其他测量 output 关联）。

`add_reset(self, target, basis='Z')`
> **输入**：`target->int`，`basis->str`（'Z' 或 'fermionic' 等）
>
> **输出**：None
>
> **效果**：添加 reset 指令到 `instructions`，并在语义上把该目标的后续测量值设定为已知的初始值（同时可生成与 reset 相关的检测器条目）。

`add_conditional(self, condition_key, true_instr, false_instr=None)`
> **输入**：`condition_key->str/int`（来自测量的 key），`true_instr->instr`, `false_instr->instr`（可选）
>
> **输出**：None
>
> **效果**：在电路中插入条件控制逻辑；当测量结果满足 `condition_key` 时执行 `true_instr`，否则执行 `false_instr`（若提供）。

`index_map(self, mapping)`
> **输入**：`mapping->list or np.ndarray`（长度应至少覆盖 `number_qubits`/`number_majorana` 的索引范围）
>
> **输出**：None
>
> **效果**：对电路中所有 instruction 的索引进行重新映射（用于把电路搬迁到不同物理布局或将子电路合并到更大的电路），并更新 `instructions` 与 `detectors` 中存储的索引信息。

`generate_detectors(self)`
> **输入**：None
>
> **输出**：`detectors->list`（派生的 detector 描述）
>
> **效果**：遍历 `instructions` 中的测量与比较逻辑，按定义生成 `detectors` 列表（每个 detector 描述一个比较关系，例如 m_i == m_j），这一步是构造 `stim.DetectorErrorModel` 的前置步骤。

`generate_error_specs(self, independent_error_model)`
> **输入**：`independent_error_model->dict or object`（定义每种门/操作对应的错误类型与概率）
>
> **输出**：`errorspecs->list`
>
> **效果**：基于给定的独立错误模型，将每个电路位置的可能错误映射为 `errorspecs` 条目（包含错误事件会引起哪些 detector 翻转和哪些 logical 翻转）。

`to_stim(self)`
> **输入**：None
>
> **输出**：`stim.DetectorErrorModel`（或可序列化到该格式的结构）
>
> **效果**：使用 `detectors` 与 `errorspecs` 生成符合 `stim` 格式的错误模型，包含 detector 定义与 error 行（每条 error 行包含触发该 error 的概率与它影响的 detector/ logical 列表）。
>
> **物理含义**：把电路级别的门/测量/错误语义转换为解码器可用的 Tanner graph 表示。

`sample(self, shots, noise_model=None)`
> **输入**：`shots->int`, `noise_model->object`（可选）
>
> **输出**：模拟结果（包括每次 shot 的 detector outcome、logical outcome）
>
> **效果**：对电路在指定噪声模型下进行蒙特卡洛采样；对每个 shot 按时间顺序遍历 `instructions` 并在需要时用 `Platform` 的方法更新稳定子，从而生成测量输出；返回一个包含统计的结果对象。

`run_decoder(self, decoder, shots=1, **kwargs)`
> **输入**：`decoder->callable or object`（接收 `detectors` / `errorspecs` 的解码器接口），`shots->int`
>
> **输出**：解码结果统计（logical error rate 等）
>
> **效果**：基于 `to_stim()` 生成的 `DetectorErrorModel` 调用指定解码器进行误差恢复并统计 logical error rate。

`copy(self)`
> **输入**：None
>
> **输出**：`Circuit` 的深拷贝（包含 `instructions`、`detectors`、`errorspecs` 的独立副本）
>
> **效果**：便于在不修改原电路的情况下生成修改后的版本（例如 index_map 后或合并子电路后）。

#### 额外建议的支持方法（便于和 Platform/Code 层配合）
- `from_platform(self, platform: Platform)`：根据 `Platform` 的稳定子与编号初始化 `Circuit` 中的索引与默认准备/测量指令；
- `merge(self, other: Circuit, index_offset_map)`：把另一个电路拼接到当前电路，自动调整索引并合并 detectors/errorspecs；
- `visualize(self, mode='text'|'graphviz')`：生成电路时序图或检测器 Tanner 图，便于检视。

#### 错误模型与精确性注释
- `to_stim()`/`generate_error_specs()` 需要精确定义“错误事件如何映射到 detector / logical 翻转”的规则；这是构造正确 Tanner graph 的关键；
- 若需要精确复现 Stim 风格的行为，应保证 `instructions` 中测量结果的编号与 `stim` 中 detector 编号策略一致（包括历史测量/时间层的索引分配）；

## 4 Workflow

### 4.1 Code theory

在处理code theory相关问题时，基本的工作流程如下。

1. 在QuantumCode类或LinearCode类中基于某种方式构造code;
2. 基于对象方法计算code parameters.

### 4.2 Circuit simulation

在处理quantum circuits相关问题时，主要的工作流程如下。

1. 在`Circuit`类中构造线路
2. 将线路生成错误模型
3. Monte-Carlo模拟，计算logical error rate，并由此得到pseudo-threshold或threshold等物理量。