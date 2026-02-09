import numpy as np
from galois import GF2
from extendedstim.Code.PrimitiveCode.MajoranaCSSCode import MajoranaCSSCode
from extendedstim.Code.PrimitiveCode.PauliCSSCode import PauliCSSCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator


##  CHAPTER: ====Check matrix of color code====
class ColorCode:
    __slots__ = ['d','H']

    ##  SECTION: ----Constructor----
    def __init__(self, d: int) -> None:
        """
        input.d：odd number representing the distance of the code
        output：None
        """

        ##  Prepare
        self.d=d  # code distance
        x_0=0
        y_0=0
        cell_vector=[np.array([0, 0], dtype=float), np.array([2, 0], dtype=float), np.array([0.5, np.sqrt(3)/2], dtype=float),
                     np.array([1.5, np.sqrt(3)/2], dtype=float), np.array([0.5, -np.sqrt(3)/2], dtype=float),
                     np.array([1.5, -np.sqrt(3)/2], dtype=float)]
        cells=[]
        x_now=0
        y_now=0

        ##  Computing cells
        for i in range(d//2):
            for j in range(3):
                for k in range(d//2-i):
                    x_now_temp=x_now+k*3
                    cells.append([np.array([x_now_temp, y_now], dtype=float)+cell_vector[k] for k in range(6)])
                if j==0:
                    y_now+=np.sqrt(3)/2
                    x_now+=1.5
                elif j==1:
                    y_now+=np.sqrt(3)/2
                    x_now-=1.5
                elif j==2:
                    y_now+=np.sqrt(3)/2
                    x_now+=1.5

        ##  computing real cells
        cells_in_triangle=[]
        for i in range(len(cells)):
            cells_in_triangle.append([])
            for j in range(6):
                if judge_in_triangle(cells[i][j][0], cells[i][j][1],x_0,y_0,d):
                    cells_in_triangle[i].append(cells[i][j])

        ##  Indexing bottom sites
        position_list=[]
        stabilizer_list=[]
        for i in range(len(cells_in_triangle)):
            for j in range(len(cells_in_triangle[i])):
                if abs(cells_in_triangle[i][j][1])<0.01:
                    flag=True
                    for k in range(len(position_list)):
                        if equal(cells_in_triangle[i][j], position_list[k]):
                            flag=False
                            break
                    if flag:
                        position_list.append(cells_in_triangle[i][j])

        ## Generating stabilizers and other sites
        for i in range(len(cells_in_triangle)):
            stabilizer_temp=[]
            for j in range(len(cells_in_triangle[i])):
                flag=True
                index=0
                for k in range(len(position_list)):
                    if equal(cells_in_triangle[i][j], position_list[k]):
                        flag=False
                        index=k
                        break
                if flag:
                    stabilizer_temp.append(len(position_list))
                    position_list.append(cells_in_triangle[i][j])
                else:
                    stabilizer_temp.append(index)
            stabilizer_list.append(stabilizer_temp)

        ##  Generating check matrix
        self.H=GF2(np.zeros((len(stabilizer_list),len(position_list)),dtype=int))
        for i in range(len(stabilizer_list)):
            self.H[i][stabilizer_list[i]]=1


##  CHAPTER: ====Majorana color codes====
class MajoranaColorCode(MajoranaCSSCode,ColorCode):

    ##  SECTION: ----Constructor----
    def __init__(self,d:int):
        ColorCode.__init__(self,d)
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H)):
            generators_x.append(MajoranaOperator.HermitianOperatorFromOccupy(np.where(self.H[i]!=0)[0],[]))
            generators_z.append(MajoranaOperator.HermitianOperatorFromOccupy([],np.where(self.H[i]!=0)[0]))
        MajoranaCSSCode.__init__(self,generators_x,generators_z,self.H.shape[1])
        self._logical_operators_x=[MajoranaOperator.HermitianOperatorFromOccupy([i for i in range(d)], [])]
        self._logical_operators_z=[MajoranaOperator.HermitianOperatorFromOccupy([], [i for i in range(d)])]


##  CHAPTER: ====Pauli color codes====
class PauliColorCode(PauliCSSCode,ColorCode):

    ##  SECTION: ----Constructor----
    def __init__(self,d:int):
        ColorCode.__init__(self,d)
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H)):
            generators_x.append(PauliOperator.HermitianOperatorFromOccupy(np.where(self.H[i]!=0)[0],[]))
            generators_z.append(PauliOperator.HermitianOperatorFromOccupy([],np.where(self.H[i]!=0)[0]))
        PauliCSSCode.__init__(self,generators_x,generators_z,self.H.shape[1])
        self._logical_operators_x=[PauliOperator.HermitianOperatorFromOccupy([i for i in range(d)], [])]
        self._logical_operators_z=[PauliOperator.HermitianOperatorFromOccupy([], [i for i in range(d)])]


##  CHAPTER: Check if a point is in a triangle
def judge_in_triangle(x, y,x_0,y_0,d) -> bool:
    """
    input.x：float，x position
    input.y：float，y position
    output：bool, whether the point is in the triangle
    """

    x_1=x_0+d//2+d-1
    x_2=x_1/2
    y_2=np.sqrt(3)*x_2

    k_0=0
    b_0=0
    k_1=(y_2-y_0)/(x_2-x_0)
    b_1=0
    k_2=-k_1
    b_2=2*y_2
    if k_0*x+b_0<y+0.01 and k_1*x+b_1>y-0.01 and k_2*x+b_2>y-0.01:
        return True
    else:
        return False


##  CHAPTER: Check if two points are equal
def equal(a, b) -> bool:
    """
    input.a：float，a position
    input.b：float，b position
    output：bool, whether the two points are equal
    """
    if abs(a[0]-b[0])+abs(a[1]-b[1])<0.01:
        return True
    else:
        return False
