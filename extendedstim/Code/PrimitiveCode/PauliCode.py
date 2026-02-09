"""
Module: PauliCode
A class for Pauli code, each of which has stabilizer generators and physical number of modes as properties.
"""
import galois
import numpy as np
from extendedstim.Code.PrimitiveCode.QuantumCode import QuantumCode
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.GaloisTools import mip_distance_caculator


class PauliCode(QuantumCode):
    """
    An example:
    $S_0=i\\hat X_0\\hat X_1\\hat X_2\\hat X_3$
    $S_1=i\\hat X_1\\hat X_2\\hat X_5\\hat X_6$
    $S_2=i\\hat X_0\\hat X_2\\hat X_4\\hat X_6$
    $S_3=i\\hat Z_0\\hat Z_1\\hat Z_2\\hat Z_3$
    $S_4=i\\hat Z_1\\hat Z_2\\hat Z_5\\hat Z_6$
    $S_5=i\\hat Z_0\\hat Z_2\\hat Z_4\\hat Z_6$
    n=7
    """

    ##  SECTION: ----Constructor----
    def __init__(self, generators:list[PauliOperator], physical_number:int):
        """
        input.generators：list of PauliOperator, Pauli stabilizer generators
        input.physical_number：number of physical qubits
        """

        ##  Check input
        try:
            generators=[generator if isinstance(generator, PauliOperator) else PauliOperator(generator) for generator in generators]
        except Exception as e:
            raise ValueError("Invalid generators") from e
        
        try:
            physical_number=int(physical_number)
            assert physical_number>=0, "physical_number must be non-negative"            
        except Exception as e:
            raise ValueError("Invalid physical_number") from e

        ##  Constructor for super class
        super().__init__(generators, physical_number)


    ##  SECTION: ----get code distance----
    @property
    def distance(self) -> int:
        """
        output：int, code distance by mip algorithm
        """
        return mip_distance_caculator(self.check_matrix, PauliOperator.get_matrix(self.logical_operators, self.physical_number))


    ##  SECTION: ----get logical operators----
    @property
    def logical_operators(self) -> list[PauliOperator]:
        """
        output：list of [(p0,q0),(p1,q1),...], a pair of logical Pauli operators corresponding to one logical qubit
        """

        ##  GF(2) numpy 0/1
        Hgf=self.check_matrix
        H_np=np.array(Hgf, dtype=int)
        m, n2=H_np.shape
        n=n2//2

        ##  J NullSpace(A) A = H * J
        J=symplectic_form(n)
        A=(H_np@J)%2
        A_gf=galois.GF2(A)
        null=A_gf.null_space()  # 返回 r x 2n 的矩阵（GF(2)）
        if null.size==0:
            return np.zeros((0, n2), dtype=int), []

        null_np=np.array(null, dtype=int)
        span=np.array(H_np, dtype=int) if H_np.shape[0]>0 else np.zeros((0, n2), dtype=int)
        span_rank=np.linalg.matrix_rank(span) if span.size else 0
        reps=[]
        for i in range(null_np.shape[0]):
            v=null_np[i:i+1, :]
            new_span=np.vstack((span, v)) if span.size else v.copy()
            new_rank=np.linalg.matrix_rank(new_span)
            if new_rank>span_rank:
                reps.append(v.flatten())
                span=new_span
                span_rank=new_rank

        if len(reps)==0:
            return np.zeros((0, n2), dtype=int), []

        logical_basis=np.vstack(reps)

        basis=[logical_basis[i].copy() for i in range(logical_basis.shape[0])]
        pairs=[]
        i=0
        while i<len(basis):
            v=basis[i]
            partner_idx=None
            for j in range(i+1, len(basis)):
                if symplectic_product(v, basis[j])==1:
                    partner_idx=j
                    break
            if partner_idx is None:
            # None
                pairs.append((i, None))
                i+=1
                continue
            w=basis[partner_idx]
            for t in range(len(basis)):
                if t==i or t==partner_idx:
                    continue
                if symplectic_product(basis[t], v)==1:
                    # basis[t] <- basis[t] + w
                    basis[t]=(basis[t]+w)%2
                if symplectic_product(basis[t], w)==1:
                    # basis[t] <- basis[t] + v
                    basis[t]=(basis[t]+v)%2
            pairs.append((i, partner_idx))
            i+=1

        ##  basis pairs basis
        logical_basis_final=np.vstack(basis)
        results=[]
        for i in range(len(logical_basis_final)):
            results.append(PauliOperator.HermitianOperatorFromVector(logical_basis_final[i]))
        return results

    ##  SECTION: ----Constructed from check matrix----
    @staticmethod
    def FromCheckMatrix(check_matrix):
        """
        input.H：np.ndarray of GF(2), check matrix of (m,2n)
        output：PauliCode, a Pauli CSS code instance
        """
        ##  Check input
        try:
            check_matrix=galois.GF2(np.array(check_matrix,dtype=int))
        except Exception as e:
            raise ValueError("Invalid check_matrix") from e

        ##  Get Pauli stabilizer generators from checker in check matrix
        generators = []
        for temp in range(check_matrix.shape[0]):
            occupy_x = np.where(check_matrix[temp, 0::2] == 1)[0]
            occupy_z = np.where(check_matrix[temp, 1::2] == 1)[0]
            generators.append(PauliOperator.HermitianOperatorFromOccupy(occupy_x, occupy_z))
        physical_number = check_matrix.shape[1] // 2

        ##  Return PauliCode instance
        return PauliCode(generators, physical_number)


##  CHAPTER: ====Symplectic form====
def symplectic_form(n):
    """
    input.n：number of physical qubits
    output: 2n*2n np.ndarray, symplectic form J
    """

    I=np.eye(n, dtype=int)  # noqa: E741
    Z=np.zeros((n, n), dtype=int)
    top=np.hstack((Z, I))
    bot=np.hstack((I, Z))
    return np.vstack((top, bot))


##  CHAPTER: ====Symplectic product====
def symplectic_product(a, b) -> int:
    """
    input.a,b: length 2n numpy.ndarray, 0/1 vectors
    output：0 or 1 as result of symplectic product
    """

    n2=a.shape[0]
    n=n2//2
    x_a=a[:n]
    z_a=a[n:]
    x_b=b[:n]
    z_b=b[n:]
    return int((x_a@z_b+z_a@x_b)%2)