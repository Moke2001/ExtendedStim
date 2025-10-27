# LinearCode 测试示例
try:
    import numpy as np
    from extendedstim.Code.LinearCode.LinearCode import LinearCode
    H = np.array([[1,0,1],[0,1,1]])
    code = LinearCode(H)
    print('rank=', code.rank)
    print('logical_number=', code.logical_number)
    print('is_dual_containing=', code.is_dual_containing)
except Exception as e:
    print('LinearCode 测试未运行，原因：', e)

# BicycleCode 示例（构造可能耗时且依赖 galois）
try:
    from extendedstim.Code.LinearCode.BicycleCode import BicycleCode
    bc = BicycleCode(6, 2, 2, seed=0)
    print('BicycleCode 构造成功')
except Exception as e:
    print('BicycleCode 示例未运行，原因：', e)