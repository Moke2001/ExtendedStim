"""
Module: QuantumCode
An abstract class for quantum code, each of which has stabilizer generators and physical number of modes as properties.
"""
import copy
from abc import ABC, abstractmethod
import numpy as np
from extendedstim.Physics.Operator import Operator


class QuantumCode(ABC):
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
    def __init__(self,generators:list[Operator],physical_number:int):
        """
        input.generators：list of Operator, stabilizer generators
        input.physical_number：int, number of physical modes
        """

        ##  Assignment
        self.generators=generators  # stabilizer generators
        self.physical_number=physical_number  # number of physical modes
        self.checker_number=len(self.generators)  # number of stabilizer generators


    ##  SECTION: ----Get check matrix----
    @property
    def check_matrix(self):
        """
        output：np.ndarray of GF(2), shape=(checker_number,physical_number)
        """

        ##  Get check matrix
        return Operator.get_matrix(self.generators, self.physical_number)


    ##  SECTION: ----Get rank----
    @property
    def rank(self)->int:
        """
        output：int，rank(H)
        """

        ##  Get rank
        return np.linalg.matrix_rank(self.check_matrix)


    ##  SECTION: ----Get number of logical modes----
    @property
    def logical_number(self)->int:
        """
        output：int = physical_number - rank
        """

        ##  Get number of logical modes
        return self.physical_number-self.rank


    ##  SECTION: ----Get code distance----
    @property
    @abstractmethod
    def distance(self)->int:
        """
        output：int，minimal weight of non-trivial logical operators
        """

        pass

    ##  SECTION: ----Get logical operators----
    @property
    @abstractmethod
    def logical_operators(self)->list[Operator]:
        """
        output：list of Operator, non-trivial logical operators
        """

        pass


    ##  SECTION: ----Mapping index to new one----
    def index_map(self, index_map,physical_number):
        """
        input.index_map：list or array, shape=(physical_number,), indexing new one to each old
        output：None (in-place modify generator indices)
        """
        for generator in self.generators:
            generator.index_map(index_map)
        self.physical_number=physical_number


    ##  SECTION: ----Deepcopy----
    def copy(self):
        """
        output：deepcopy of new instance
        """
        return copy.deepcopy(self)
