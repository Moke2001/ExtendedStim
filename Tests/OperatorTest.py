# 示例 / 测试（Operator 是抽象类，示例使用 PauliOperator 作为实现）
try:
    from extendedstim.Physics.PauliOperator import PauliOperator
    # 构造两个 Pauli 算符
    a = PauliOperator([0], [1], 1)
    b = PauliOperator([], [0], 1)
    print('a:', a)
    print('b:', b)
    prod = a @ b
    print('a @ b =', prod)
    print('a.weight =', a.weight)
except Exception as e:
    print('Operator示例未运行，原因：', e)
# MajoranaOperator 测试示例
try:
    from extendedstim.Physics.MajoranaOperator import MajoranaOperator
    a = MajoranaOperator([0], [], 1)
    b = MajoranaOperator([], [0], 1)
    print('Majorana a:', a)
    print('Majorana b:', b)
    print('commute(a,b)=', MajoranaOperator.commute(a,b))
    print('a @ b =', a @ b)
except Exception as e:
    print('MajoranaOperator 测试未运行，原因：', e)
# PauliOperator 测试示例
try:
    from extendedstim.Physics.PauliOperator import PauliOperator
    a = PauliOperator([0], [], 1)
    b = PauliOperator([], [0], 1)
    print('Pauli a:', a)
    print('Pauli b:', b)
    print('commute(a,b)=', PauliOperator.commute(a,b))
    print('a.dual =', a.dual)
    print('a @ b =', a @ b)
except Exception as e:
    print('PauliOperator 测试未运行，原因：', e)