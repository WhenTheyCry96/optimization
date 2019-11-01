# -*- coding: utf-8 -*- 
import sys
import math
import numpy as np
from inspect import signature


class Ops():
    def __init__(self, func, init, MaxIter, tolerance):
        sig = signature(func)
        params = sig.parameters
        self.dimension = len(params)
        self.cur_point = init
        if len(self.cur_point) is not self.dimension:
            print("DIFFERENT DIMENSION BETWEEN INITIAL POINT and FUNCTION")
            sys.exit(-1)
    
class SteepestDescent(Ops):
    def __init__(self, func, init, delta=1e-10, MaxIter=1e6, tolerance=1e-6):
        super().__init__(func, init, MaxIter, tolerance)
        self.delta = delta
        self.maxiter = MaxIter
        self.tolerance = tolerance
        pass
    

##########
def funct(x,y):
    return math.pow(1.5, 3) + 3
sd1 = SteepestDescent(funct, 1)