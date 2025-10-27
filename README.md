# Moke Quantum Computing Research Framework

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen)](docs/)

本程序基于Python 3.12+开发，主要用于量子纠错码和量子线路的构造与测试。

## 🔨 1 项目依赖

- [QuTiP](https://qutip.org/) - 量子工具包
- [Stim](https://github.com/quantumlib/Stim) - 量子纠错模拟器
- [Stimbposd](https://github.com/quantumlib/Stim/blob/main/docs/bposd.md) - 基于Stim的BPOSD译码器
- [Galois](https://galois.readthedocs.io/) - 提供$$\mathbb{F}_2$$上的代数计算
- [NumPy](https://numpy.org/) - 数值计算库
- [SciPy](https://scipy.org/) - 科学计算库
- [Matplotlib](https://matplotlib.org/) - 绘图库
- [Qiskit](https://qiskit.org/) - 提供量子线路图的绘制
- [Mip](https://www.mipengine.org/) - 整数规划求解器，用于code distance的计算

## 📁 2 项目结构

```
├── extendedstim/               # 核心代码目录
│   ├── Circuit/                # 量子线路相关模块
│   │   ├── Circuit.py          # 量子线路实现
│   │   └── __init__.py         # 模块初始化文件
│   ├── Code/                   # 量子码和线性码模块
│   │   ├── LinearCode/         # 线性码实现
│   │   │   ├── BicycleCode.py  # 自行车码实现
│   │   │   ├── LinearCode.py   # 线性码基类
│   │   │   └── __init__.py     # 模块初始化文件
│   │   ├── QuantumCode/        # 量子码实现
│   │   │   ├── LatticeSurgery.py      # 格点手术
│   │   │   ├── MajoranaCSSCode.py     # Majorana CSS码
│   │   │   ├── MajoranaCode.py        # Majorana码
│   │   │   ├── PauliCSSCode.py        # Pauli CSS码
│   │   │   ├── PauliCode.py           # Pauli码
│   │   │   ├── QuantumCSSCode.py      # 量子CSS码
│   │   │   ├── QuantumCode.py         # 量子码基类
│   │   │   └── __init__.py            # 模块初始化文件
│   │   └── __init__.py         # 模块初始化文件
│   ├── Math/                   # 数学工具模块
│   │   ├── BinaryArray.py      # 二进制数组操作
│   │   ├── Graph.py            # 图论工具
│   │   └── __init__.py         # 模块初始化文件
│   └── Physics/                # 物理操作模块
│       ├── MajoranaOperator.py # Majorana算符
│       ├── Operator.py         # 算符基类
│       ├── PauliOperator.py    # Pauli算符
│       └── __init__.py         # 模块初始化文件
├── .git/                       # Git版本控制
├── .idea/                      # IDE配置文件
└── README.md                   # 项目说明文档
```

## 📖 3 基本工作流

### 3.1 计算code parameters

1. 构造量子纠错码
2. 计算量子纠错码的code parameters

### 3.2 计算logical error rate

1. 构造量子线路
2. 执行Monte-Carlo模拟，对比预测正确与否得到logical error rate

## 📄 4 许可证

本项目采用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 📞 5 联系方式

- **作者**: Moke
- **邮箱**: Moke2001@whu.edu.cn
- **地址**: 北京市海淀区清华大学蒙民伟科技楼S219
- **电话**: +86 130-3373-6868