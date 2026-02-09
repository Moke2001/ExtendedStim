"""
Module: TypingTools
"""
import numpy as np


##  CHAPTER: Judge if a number is integer
def isinteger(num) -> bool:
    """""
    input.num：Any, any object
    output：bool, True if num is an integer type (including built-in int/float and numpy integer types)
    influence：Used to check if a parameter is integer
    """""
    return isinstance(num, (int, float,np.int8,np.int16,np.int32, np.int64))


##  CHAPTER: Judge if a number is a list, tuple, numpy.ndarray, or range object
def islist(num) -> bool:
    """""
    input.num：Any, any object
    output：bool, True if num is a list, tuple, numpy.ndarray, or range object
    influence：Used to distinguish between single values and sequence inputs
    """""
    return isinstance(num, (list, tuple, np.ndarray,range))
