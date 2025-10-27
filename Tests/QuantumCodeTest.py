# QuantumCode 简单示例（使用 MajoranaCode/PauliCode 子类更有意义）
try:
    from extendedstim.Code.QuantumCode.MajoranaCode import MajoranaCode
    # 构造时通常使用 FromCheckMatrix 或 传入生成器列表
    print('QuantumCode 示例提示：请使用 MajoranaCode.FromCheckMatrix 或 PauliCode.FromCheckMatrix')
except Exception as e:
    print('QuantumCode 示例不可用，原因：', e)

# QuantumCSSCode 示例（提示）
print('QuantumCSSCode 是抽象基类；请查看 MajoranaCSSCode / PauliCSSCode 的实现示例。')

# PauliCode 示例
try:
    import numpy as np
    from extendedstim.Code.QuantumCode.PauliCode import PauliCode
    H = np.array([[1,0,0,1]])  # 占位：实际按 2 列为一位编码
    pc = PauliCode.FromCheckMatrix(H)
    print('PauliCode 构造成功，physical_number=', pc.physical_number)
except Exception as e:
    print('PauliCode 示例未运行，原因：', e)

print('PauliCSSCode: 请参考源码以补全 distance / logical_operators 的实现。')

# MajoranaCode 示例
try:
    import numpy as np
    from extendedstim.Code.QuantumCode.MajoranaCode import MajoranaCode
    H = np.array([[1,1,0,0]])
    mc = MajoranaCode.FromCheckMatrix(H)
    print('MajoranaCode physical_number=', mc.physical_number)
except Exception as e:
    print('MajoranaCode 示例未运行，原因：', e)