import copy

import galois
import networkx as nx
import numpy as np
from .MajoranaCSSCode import MajoranaCSSCode
from .. import MajoranaCode
from ...Physics import MajoranaOperator


# %%  USER：实现fermionic lattice surgery
def MajoranaLatticeSurgery(code_A, code_B, index_A, index_B):
    code_A = code_A.copy()
    code_B = code_B.copy()
    assert isinstance(code_A, MajoranaCSSCode)
    assert isinstance(code_B, MajoranaCSSCode)
    assert isinstance(index_A, int)
    assert isinstance(index_B, int)

    # %%  SECTION：数据预处理
    majorana_number_A = code_A.physical_number
    majorana_number_B = code_B.physical_number
    code_B.index_map(np.arange(majorana_number_A, majorana_number_B+majorana_number_A))
    code_A.index_map(np.arange(majorana_number_A))
    logical_operator_0 = code_A.logical_operators_x[index_A]
    logical_operator_1 = code_B.logical_operators_x[index_B]
    support_index_vector_A = np.array(logical_operator_0.x_vector, dtype=int)
    support_index_vector_B = np.array(logical_operator_1.x_vector, dtype=int)

    ##  目标初始化
    check_origin_list = []
    check_ancilla_list = []
    check_modify_list = []
    check_stable_list = []
    check_gauge_list = []
    check_z_list = []
    check_measure_list = []

    vertex_qubit_list = []
    vertex_check_list = []

    edge_list = []
    for i in range(len(code_A.generators)):
        if len(code_A.generators[i].x_vector) > 0:
            if len(set(code_A.generators[i].x_vector) & set(support_index_vector_A)) > 0:
                check_origin_list.append(code_A.generators[i].copy())
                check_modify_list.append(code_A.generators[i].copy())
            else:
                check_stable_list.append(code_A.generators[i].copy())
        else:
            check_z_list.append(code_A.generators[i].copy())

    for i in range(len(code_B.generators)):
        if len(code_B.generators[i].x_vector) > 0:
            if len(set(code_B.generators[i].x_vector) & set(support_index_vector_B)) > 0:
                check_origin_list.append(code_B.generators[i].copy())
                check_modify_list.append(code_B.generators[i].copy())
            else:
                check_stable_list.append(code_B.generators[i].copy())
        else:
            check_z_list.append(code_B.generators[i].copy())

    ##  为待修改的stabilizers增加ancilla及其索引
    code = MajoranaCode()
    code.define_qubit(majorana_number_A + majorana_number_B)
    ancilla_list_list = []
    for i in range(len(check_modify_list)):
        x_vector_temp = check_modify_list[i].x_vector.tolist()
        z_vector_temp = []
        number_ancilla_temp = len(
            set(check_modify_list[i].x_vector) & set(np.append(support_index_vector_A, support_index_vector_B))) // 2
        temp = []
        for j in range(number_ancilla_temp):
            if check_modify_list[i].x_vector[0] < majorana_number_A:
                vertex_check_list.append((str(len(vertex_check_list)) + 'A'))
            else:
                vertex_check_list.append((str(len(vertex_check_list)) + 'B'))
            code.push_qubit(1)
            temp.append((code.qubit_list[-1], 'x'))
            temp.append((code.qubit_list[-1], 'z'))
            x_vector_temp.append(code.qubit_list[-1])
            z_vector_temp.append(code.qubit_list[-1])
            edge_list.append((str(code.qubit_list[-1]) + 'x', vertex_check_list[-1]))
            edge_list.append((str(code.qubit_list[-1]) + 'z', vertex_check_list[-1]))
        check_modify_list[i] = MajoranaOperator(x_vector_temp, z_vector_temp, 1)
        ancilla_list_list.append(temp)
    for i in range(majorana_number_A + majorana_number_B, code.number_qubit):
        vertex_qubit_list.append(str(i) + 'x')
        vertex_qubit_list.append(str(i) + 'z')

    # %%  SECTION：加入测量稳定子
    single_point = None
    single_qubit_list = []
    ##  右边logical更长的情况
    if len(support_index_vector_A) >= len(support_index_vector_B):
        ##  先将两边对齐的部分连起来
        for i in range(len(support_index_vector_B)):
            vertex_check_list.append((str(len(vertex_check_list)) + 'M'))
            x_vector_temp = [support_index_vector_A[i], support_index_vector_B[i]]
            z_vector_temp = []
            for j in range(len(check_modify_list)):
                if support_index_vector_A[i] in check_modify_list[j].x_vector or support_index_vector_B[i] in \
                        check_modify_list[j].x_vector:
                    temp = ancilla_list_list[j][-1]
                    if temp[1] == 'x':
                        x_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'x', vertex_check_list[-1]))
                    elif temp[1] == 'z':
                        z_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'z', vertex_check_list[-1]))
                    ancilla_list_list[j].pop(-1)

            if np.mod(len(x_vector_temp) + len(z_vector_temp), 2) != 0:
                if single_point is None:
                    code.push_qubit(1)
                    single_point = code.qubit_list[-1]
                    vertex_qubit_list.append(str(single_point) + 'x')
                    vertex_qubit_list.append(str(single_point) + 'z')
                    edge_list.append((str(single_point) + 'x', vertex_check_list[-1]))
                    single_qubit_list.append(code.qubit_list[-1])
                    x_vector_temp.append(code.qubit_list[-1])
                else:
                    z_vector_temp.append(single_point)
                    edge_list.append((str(single_point) + 'z', vertex_check_list[-1]))
                    single_point = None

            ##  引入新的measurement stabilizer
            check_measure_list.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))

        ##  将右边剩余的部分连起来
        length_B = len(support_index_vector_B)
        length_A = len(support_index_vector_A)
        for i in range((length_A - length_B) // 2):
            index_0 = length_B + 2 * i
            index_1 = index_0 + 1
            vertex_check_list.append((str(len(vertex_check_list)) + 'M'))
            x_vector_temp = [support_index_vector_A[index_0], support_index_vector_A[index_1]]
            z_vector_temp = []
            for j in range(len(check_modify_list)):
                if support_index_vector_A[index_0] in check_modify_list[j].x_vector:
                    temp = ancilla_list_list[j][-1]
                    if temp[1] == 'x':
                        x_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'x', vertex_check_list[-1]))
                    elif temp[1] == 'z':
                        z_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'z', vertex_check_list[-1]))
                    ancilla_list_list[j].pop(-1)
                if support_index_vector_A[index_1] in check_modify_list[j].x_vector:
                    temp = ancilla_list_list[j][-1]
                    if temp[1] == 'x':
                        x_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'x', vertex_check_list[-1]))
                    elif temp[1] == 'z':
                        z_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'z', vertex_check_list[-1]))
                    ancilla_list_list[j].pop(-1)

            if np.mod(len(x_vector_temp) + len(z_vector_temp), 2) != 0:
                if single_point is None:
                    code.push_qubit(1)
                    single_point = code.qubit_list[-1]
                    vertex_qubit_list.append(str(single_point) + 'x')
                    vertex_qubit_list.append(str(single_point) + 'z')
                    edge_list.append((str(single_point) + 'x', vertex_check_list[-1]))
                    single_qubit_list.append(code.qubit_list[-1])
                    x_vector_temp.append(code.qubit_list[-1])
                else:
                    z_vector_temp.append(single_point)
                    edge_list.append((str(single_point) + 'z', vertex_check_list[-1]))
                    single_point = None

            ##  引入新的measurement stabilizer
            check_measure_list.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))
    else:
        ##  先将两边对齐的部分连起来
        for i in range(len(support_index_vector_A)):
            vertex_check_list.append((str(len(vertex_check_list)) + 'M'))
            x_vector_temp = [support_index_vector_A[i], support_index_vector_B[i]]
            z_vector_temp = []
            for j in range(len(check_modify_list)):
                if support_index_vector_A[i] in check_modify_list[j].x_vector or support_index_vector_B[i] in \
                        check_modify_list[j].x_vector:
                    temp = ancilla_list_list[j][-1]
                    if temp[1] == 'x':
                        x_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'x', vertex_check_list[-1]))
                    elif temp[1] == 'z':
                        z_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'z', vertex_check_list[-1]))
                    ancilla_list_list[j].pop(-1)

            if np.mod(len(x_vector_temp) + len(z_vector_temp), 2) != 0:
                if single_point is None:
                    code.push_qubit(1)
                    single_point = code.qubit_list[-1]
                    vertex_qubit_list.append(str(single_point) + 'x')
                    vertex_qubit_list.append(str(single_point) + 'z')
                    edge_list.append((str(single_point) + 'x', vertex_check_list[-1]))
                    single_qubit_list.append(code.qubit_list[-1])
                    x_vector_temp.append(code.qubit_list[-1])
                else:
                    z_vector_temp.append(single_point)
                    edge_list.append((str(single_point) + 'z', vertex_check_list[-1]))
                    single_point = None

            ##  引入新的measurement stabilizer
            check_measure_list.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))

        ##  将右边剩余的部分连起来
        length_B = len(support_index_vector_B)
        length_A = len(support_index_vector_A)
        for i in range((length_B - length_A) // 2):
            index_0 = length_A + 2 * i
            index_1 = index_0 + 1
            vertex_check_list.append((str(len(vertex_check_list)) + 'M'))
            x_vector_temp = [support_index_vector_B[index_0], support_index_vector_B[index_1]]
            z_vector_temp = []
            for j in range(len(check_modify_list)):
                if support_index_vector_B[index_0] in check_modify_list[j].x_vector or support_index_vector_B[
                    index_1] in check_modify_list[j].x_vector:
                    temp = ancilla_list_list[j][-1]
                    if temp[1] == 'x':
                        x_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'x', vertex_check_list[-1]))
                    elif temp[1] == 'z':
                        z_vector_temp.append(temp[0])
                        edge_list.append((str(temp[0]) + 'z', vertex_check_list[-1]))
                    ancilla_list_list[j].pop(-1)

            if np.mod(len(x_vector_temp) + len(z_vector_temp), 2) != 0:
                if single_point is None:
                    code.push_qubit(1)
                    single_point = code.qubit_list[-1]
                    vertex_qubit_list.append(str(single_point) + 'x')
                    vertex_qubit_list.append(str(single_point) + 'z')
                    edge_list.append((str(single_point) + 'x', vertex_check_list[-1]))
                    single_qubit_list.append(code.qubit_list[-1])
                    x_vector_temp.append(code.qubit_list[-1])
                else:
                    z_vector_temp.append(single_point)
                    edge_list.append((str(single_point) + 'z', vertex_check_list[-1]))
                    single_point = None

            ##  引入新的measurement stabilizer
            check_measure_list.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))

    # %%  SECTION：图论计算规范稳定子
    for i in range(code.number_qubit - majorana_number_A - majorana_number_B):
        temp = MajoranaOperator([i + majorana_number_B + majorana_number_A], [i + majorana_number_B + majorana_number_A], 1)
        check_ancilla_list.append(temp)

    ##  获取关键参数

    for i in range(len(single_qubit_list)):
        vertex_check_list.append(str(i) + 'D')
        edge_list.append((str(single_qubit_list[i]) + 'x', str(i) + 'D'))
        edge_list.append((str(single_qubit_list[i]) + 'z', str(i) + 'D'))

    graph = nx.Graph()
    graph.add_nodes_from(vertex_check_list)
    graph.add_nodes_from(vertex_qubit_list)
    graph.add_edges_from(edge_list)
    check_gauge_list = IndependentCyclesFinder(graph)

    target = logical_operator_0.mul(logical_operator_1, code.number_qubit)

    # %%  SECTION：返回结果

    code.check_list = check_stable_list + check_modify_list+ check_gauge_list + check_z_list+check_measure_list
    code.number_checker = len(code.check_list)
    code.target=target
    assert code.commute_judge()  # 检查对易性
    return code


def ZechuanLatticeSurgery(code_ldpc, index)->MajoranaCode:

    #%%  SECTION：数据标准化
    code_ldpc=code_ldpc.copy()
    assert isinstance(code_ldpc, MajoranaCSSCode)
    assert isinstance(index, int)


    #%%  SECTION：数据预处理
    number_qubit_ldpc=code_ldpc.number_qubit
    logical_operator_ldpc=code_ldpc.logical_operator_list_x[index]
    support_index_vector_ldpc = np.array(logical_operator_ldpc.x_vector, dtype=int)
    code_color=MajoranaCode.ColorCode(len(support_index_vector_ldpc))
    number_qubit_color=code_color.number_qubit
    logical_operator_color = code_color.logical_operator_list_x[0]
    code_color.index_map(number_qubit_ldpc+number_qubit_color,np.arange(number_qubit_color)+number_qubit_ldpc)
    code_ldpc.index_map(number_qubit_ldpc+number_qubit_color,np.arange(number_qubit_ldpc))
    support_index_vector_color=np.array(logical_operator_color.x_vector,dtype=int)
    support_index_vector_color=support_index_vector_color.tolist()
    support_index_vector_ldpc=support_index_vector_ldpc.tolist()

    ##  提取与这些费米子相关联的校验子
    check_list_x_fix_ldpc=[]
    check_list_x_unfix_ldpc=[]
    check_list_z_ldpc=[]
    check_list_x_fixed_ldpc=[]
    check_list_x_fix_color=[]
    check_list_x_unfix_color=[]
    check_list_z_color=[]
    check_list_x_fixed_color=[]

    for i in range(len(code_ldpc.check_list)):
        if len(code_ldpc.check_list[i].x_vector)>0:
            if len(set(code_ldpc.check_list[i].x_vector)&set(support_index_vector_ldpc))>0:
                check_list_x_fix_ldpc.append(code_ldpc.check_list[i].copy())
            else:
                check_list_x_unfix_ldpc.append(code_ldpc.check_list[i].copy())
        else:
            check_list_z_ldpc.append(code_ldpc.check_list[i].copy())

    for i in range(len(code_color.check_list)):
        if len(code_color.check_list[i].x_vector)>0:
            if len(set(code_color.check_list[i].x_vector)&set(support_index_vector_color))>0:
                check_list_x_fix_color.append(code_color.check_list[i].copy())
            else:
                check_list_x_unfix_color.append(code_color.check_list[i].copy())
        else:
            check_list_z_color.append(code_color.check_list[i].copy())

    code=MajoranaCode()
    code.define_qubit(number_qubit_ldpc+number_qubit_color)
    code.check_origin_list=check_list_x_fix_color+check_list_x_fix_ldpc
    ##  记录与support关联的ancilla的索引
    ancilla_dict_ldpc={}
    ancilla_dict_color={}
    for i in range(len(support_index_vector_ldpc)):
        ancilla_dict_ldpc[support_index_vector_ldpc[i]]=[]
        ancilla_dict_color[support_index_vector_color[i]]=[]

    ##  对color code索引
    check_index_color=code_color.logical_plaqutte=[1,0,2,3]
    single_plaqutte=[0,1,2,3]
    def find_plaqutte(index_0,index_1):
        op_0=check_list_x_fix_color[0]
        op_1=check_list_x_fix_color[1]
        op_2=check_list_x_fix_color[2]
        op_3=check_list_x_fix_color[3]
        if index_0==0 and index_1==1:
            return op_1,[1],None
        elif index_0==0 and index_1==2:
            return op_1.mul(op_0,code.number_qubit),[0,1],None
        elif index_0==0 and index_1==3:
            return op_1.mul(op_2,code.number_qubit),[1,2],[support_index_vector_color[1],support_index_vector_color[2]]
        elif index_0==0 and index_1==4:
            return op_1.mul(op_3,code.number_qubit),[1,3],[support_index_vector_color[1],support_index_vector_color[3]]
        elif index_0==1 and index_1==2:
            return op_0,[0],None
        elif index_0==1 and index_1==3:
            return op_0.mul(op_2,code.number_qubit),[0,2],None
        elif index_0==1 and index_1==4:
            return op_0.mul(op_3,code.number_qubit),[0,3],[support_index_vector_color[2],support_index_vector_color[3]]
        elif index_0==2 and index_1==3:
            return op_2,[2],None
        elif index_0==2 and index_1==4:
            return op_2.mul(op_3,code.number_qubit),[2,3],None
        elif index_0==3 and index_1==4:
            return op_3,[3],None

    check_list_gauge=[]
    ##  对support组队
    couple_list=[]
    for i in range(len(check_list_x_fix_ldpc)):
        temp=set(check_list_x_fix_ldpc[i].x_vector) & set(support_index_vector_ldpc)
        temp=list(temp)
        couple_list.append([])
        for j in range(len(temp)//2):
            couple_list[i].append((temp[2*j],temp[2*j+1]))

    ##  修改LDPC code的check
    for i in range(len(check_list_x_fix_ldpc)):
        for j in range(len(couple_list[i])):
            x_vector_temp=check_list_x_fix_ldpc[i].x_vector.copy()
            z_vector_temp=check_list_x_fix_ldpc[i].z_vector.copy()
            code.push_qubit(1)
            ancilla_dict_ldpc[couple_list[i][j][0]].append((code.qubit_list[-1],'x'))
            ancilla_dict_ldpc[couple_list[i][j][1]].append((code.qubit_list[-1],'z'))
            x_vector_temp=np.append(x_vector_temp,code.qubit_list[-1])
            z_vector_temp=np.append(z_vector_temp,code.qubit_list[-1])
            check_list_x_fix_ldpc[i]=MajoranaOperator(x_vector_temp, z_vector_temp,1)
            index_0=couple_list[i][j][0]
            index_1=couple_list[i][j][1]
            index_min=support_index_vector_ldpc.index(index_0)
            index_max=support_index_vector_ldpc.index(index_1)
            if index_min>index_max:
                temp=index_min
                index_min=index_max
                index_max=temp
            operator_color_temp,index_list,associate_list=find_plaqutte(index_min,index_max)
            x_vector_temp = operator_color_temp.x_vector.copy()
            z_vector_temp = operator_color_temp.z_vector.copy()
            for i in range(len(index_list)):
                single_plaqutte[index_list[i]]=None
            if associate_list is not None:
                index_temp_0=associate_list[0]
                index_temp_1=associate_list[1]
                code.push_qubit(2)

                ancilla_dict_color[index_temp_0].append((code.qubit_list[-1],'x'))
                ancilla_dict_color[index_temp_0].append((code.qubit_list[-2],'x'))
                ancilla_dict_color[index_temp_1].append((code.qubit_list[-1],'z'))
                ancilla_dict_color[index_temp_1].append((code.qubit_list[-2],'z'))
                x_vector_temp_gauge=[code.qubit_list[-2],code.qubit_list[-1]]
                z_vector_temp_gauge=[code.qubit_list[-2],code.qubit_list[-1]]
                check_list_gauge.append(MajoranaOperator(x_vector_temp_gauge, z_vector_temp_gauge, 1))
                x_vector_temp=np.append(x_vector_temp,code.qubit_list[-2])
                z_vector_temp=np.append(z_vector_temp,code.qubit_list[-2])
            index_0_color = support_index_vector_color[support_index_vector_ldpc.index(index_0)]
            index_1_color = support_index_vector_color[support_index_vector_ldpc.index(index_1)]
            code.push_qubit(1)
            x_vector_temp=np.append(x_vector_temp,code.qubit_list[-1])
            z_vector_temp=np.append(z_vector_temp,code.qubit_list[-1])
            ancilla_dict_color[index_0_color].append((code.qubit_list[-1],'x'))
            ancilla_dict_color[index_1_color].append((code.qubit_list[-1],'z'))
            check_list_x_fixed_color.append(MajoranaOperator(x_vector_temp, z_vector_temp,1))
            if associate_list is not None:
                x_vector_temp=[code.qubit_list[-4],code.qubit_list[-1]]
                z_vector_temp=[code.qubit_list[-4],code.qubit_list[-1]]
            else:
                x_vector_temp=[code.qubit_list[-2],code.qubit_list[-1]]
                z_vector_temp=[code.qubit_list[-2],code.qubit_list[-1]]
            check_list_gauge.append(MajoranaOperator(x_vector_temp, z_vector_temp,1))

    for i in range(len(single_plaqutte)):
        if single_plaqutte[i] is not None:
            operator_color_temp=check_list_x_fix_color[single_plaqutte[i]]
            overlap=list(set(operator_color_temp.x_vector)&set(support_index_vector_color))
            x_vector_temp = operator_color_temp.x_vector.copy()
            z_vector_temp = operator_color_temp.z_vector.copy()
            code.push_qubit(2)
            ancilla_dict_color[overlap[0]].append((code.qubit_list[-1], 'x'))
            ancilla_dict_color[overlap[0]].append((code.qubit_list[-2], 'x'))
            ancilla_dict_color[overlap[1]].append((code.qubit_list[-1], 'z'))
            ancilla_dict_color[overlap[1]].append((code.qubit_list[-2], 'z'))
            x_vector_temp=np.append(x_vector_temp,code.qubit_list[-2])
            z_vector_temp=np.append(z_vector_temp,code.qubit_list[-2])
            check_list_x_fixed_color.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))
            x_vector_temp = [code.qubit_list[-2], code.qubit_list[-1]]
            z_vector_temp = [code.qubit_list[-2], code.qubit_list[-1]]
            check_list_gauge.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))
    check_list_measurement = []
    for i in range(len(support_index_vector_ldpc)):
        x_vector_temp=[support_index_vector_ldpc[i],support_index_vector_color[i]]
        z_vector_temp=[]
        for j in range(len(ancilla_dict_ldpc[support_index_vector_ldpc[i]])):
            if ancilla_dict_ldpc[support_index_vector_ldpc[i]][j][1]=='x':
                x_vector_temp=np.append(x_vector_temp,ancilla_dict_ldpc[support_index_vector_ldpc[i]][j][0])
            else:
                z_vector_temp = np.append(z_vector_temp, ancilla_dict_ldpc[support_index_vector_ldpc[i]][j][0])
        for j in range(len(ancilla_dict_color[support_index_vector_color[i]])):
            if ancilla_dict_color[support_index_vector_color[i]][j][1]=='x':
                x_vector_temp=np.append(x_vector_temp,ancilla_dict_color[support_index_vector_color[i]][j][0])
            else:
                z_vector_temp = np.append(z_vector_temp, ancilla_dict_color[support_index_vector_color[i]][j][0])
        check_list_measurement.append(MajoranaOperator(x_vector_temp, z_vector_temp, 1))
    check_list_x_fixed_ldpc=check_list_x_fix_ldpc
    code.check_list=check_list_x_unfix_ldpc+check_list_x_unfix_color
    code.check_stable_list=check_list_x_unfix_ldpc+check_list_x_unfix_color+check_list_z_ldpc+check_list_z_color
    code.check_list+=check_list_z_ldpc+check_list_z_color
    code.check_list+=check_list_x_fixed_color+check_list_x_fixed_ldpc
    code.check_fixed_list=check_list_x_fixed_color+check_list_x_fixed_ldpc
    code.check_gauge_list=check_list_gauge
    code.check_list+=check_list_gauge
    code.ancilla_list=[]
    code.stop=len(code.check_list)
    code.check_list+=check_list_measurement
    code.check_measure_list=check_list_measurement
    for i in range(number_qubit_ldpc+number_qubit_color,code.number_qubit):
        temp=MajoranaOperator([i],[i],1)
        code.ancilla_list.append(temp)
    code.number_checker=len(code.check_list)

    #%%  SECTION：返回结果

    assert code.commute_judge()  # 检查对易性
    return code


def IndependentCyclesFinder(G):
    qubit_list=[]
    check_list=[]
    for index,value in enumerate(G.nodes()):
        if value[-1]=='A' or value[-1]=='B' or value[-1]=='D' or value[-1]=='M':
            check_list.append(value)
        else:
            qubit_list.append(value)
    matrix=np.zeros((len(check_list),len(qubit_list)),dtype=int)
    for index0,check in enumerate(check_list):
        for index1,qubit in enumerate(qubit_list):
            if (qubit,check) in G.edges():
                matrix[index0,index1]=1
    GF=galois.GF(2**1)
    matrix=GF(matrix)
    result=matrix.null_space()
    check_gauge_list=[]
    for i in range(len(result)):
        index_list=np.where(result[i]!=0)[0]
        temp_x=[]
        temp_z=[]
        for j in range(len(index_list)):
            temp=qubit_list[index_list[j]]
            if temp[-1]=='x':
                temp_x.append(int(temp[:-1]))
            elif temp[-1]=='z':
                temp_z.append(int(temp[:-1]))
            else:
                raise ValueError
        check_gauge_list.append(MajoranaOperator(temp_x,temp_z,1))
    return check_gauge_list