"""
Module: QuantumCSSCode
An abstract class for quantum CSS code, each of which has X and Z stabilizer generators and physical number of modes as properties.
"""
from abc import abstractmethod
import numpy as np
from extendedstim.Code.PrimitiveCode.QuantumCode import QuantumCode
from extendedstim.Physics.Operator import Operator


class QuantumCSSCode(QuantumCode):
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
    def __init__(self, generators_x:list, generators_z:list, physical_number:int):
        """
        input.generators_x：list of Operator, X stabilizers
        input.generators_z：list of Operator, Z stabilizers
        input.physical_number：number of physical modes
        """

        self.generators_x=[temp for temp in generators_x]
        self.generators_z=[temp for temp in generators_z]
        self.checker_number_x=len(generators_x)
        self.checker_number_z = len(generators_z)
        self._logical_operators_x=None
        self._logical_operators_z=None
        super().__init__(generators_x+generators_z, physical_number)


    ##  SECTION: ----Get X logical operators----
    @property
    @abstractmethod
    def logical_operators_x(self)->list[Operator]:
        """
        output：list of Operator, logical X operators
        """
        
        pass


    ##  SECTION: ----Get Z logical operators----
    @property
    @abstractmethod
    def logical_operators_z(self)->list[Operator]:
        """
        output：list of Operator, logical Z operators
        """

        pass


    ##  SECTION: ----Get rank of X generators----
    @property
    def rank_x(self)->int:
        """
        output：int, rank(Hx)
        """

        ##  Get rank of X generators
        return np.linalg.matrix_rank(self.check_matrix_x)


    ##  SECTION: ----Get rank of Z generators----
    @property
    def rank_z(self)->int:
        """
        output：int, rank(Hz)
        """

        ##  Get rank of Z generators
        return np.linalg.matrix_rank(self.check_matrix_z)


    ##  SECTION: ----Get check matrix of X generators----
    @property
    def check_matrix_x(self)->np.ndarray:
        """
        output：np.ndarray of GF(2), check matrix of X generators
        """

        ##  Compute check matrix of X generators
        matrix=Operator.get_matrix(self.generators_x, self.physical_number)
        matrix=matrix[:,0::2]

        ##  Return check matrix of X generators
        return matrix


    ##  SECTION: ----Get check matrix of Z generators----
    @property
    def check_matrix_z(self)->np.ndarray:
        """
        output：np.ndarray of GF(2), check matrix of Z generators
        """

        ##  Compute check matrix of Z generators
        matrix=Operator.get_matrix(self.generators_z, self.physical_number)
        matrix = matrix[:, 1::2]

        ##  Return check matrix of Z generators
        return matrix


    ##  SECTION: ----Get code distance of X generators----
    @property
    @abstractmethod
    def distance_x(self)->int:
        """
        output：int, code distance of X generators
        """

        pass


    ##  SECTION: ----Get code distance of Z generators----
    @property
    @abstractmethod
    def distance_z(self)->int:
        """
        output：int, code distance of Z generators
        """

        pass
