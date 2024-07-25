# -*- coding: utf-8 -*-

from __future__ import annotations
from toolkit.pvector import PVector

# The class to represent a line segment in 3D. A line segment is defined by its
# end points A and B. A point in the line segment is P(t) = (1 - t)A + tB, 
# where t is a scalar in range [0, 1].
class LineSegment:
    
    # Constructor of the class.
    def __init__(self, A=PVector(), B=PVector()):
        assert A is not None, "Missing A"
        assert isinstance(A, PVector), "A is not a PVector"
        assert B is not None, "Missing B"
        assert isinstance(B, PVector), "B is not a PVector"
        self.A = A.copy()
        self.B = B.copy()
    
    
    def __getattr__(self, name:str):
        return self.__dict__[f"{name}"]
    
    
    def __repr__(self):
        return f"LineSegment:(A:{self.A}, B:{self.B})"
    
    
    def __setattr__(self, name:str, value):
        assert isinstance(value, PVector), "value not a PVector"
        self.__dict__[f"{name}"] = value
    
    
    def __str__(self):
        return f"(A:{self.A}, B:{self.B})"
    
    
    def at(self, t) -> PVector:
        assert isinstance(t, (int, float)), "t not a real number"
        return ((1.0 - t) * self.A) + (t * self.B)
    
    
    def clone(self) -> LineSegment:
        return self.copy()
    
    
    def copy(self) -> LineSegment:
        return LineSegment(self.A, self.B)
    
    
    def direction(self) -> PVector:
        return self.B - self.A
    
    
    
    
    
    def ispoint(self, P:PVector, e:int|float=1e-8) -> bool:
        assert P is not None, "Missing P"
        assert isinstance(P, PVector), "P not a PVector"
        assert isinstance(e, (int, float)), "e not a real number"
        
        # Get the vector from the start vertex of the line segment and point P
        AP = P - self.A
        
        # Get the direction vector of the line segment
        AB = self.direction()
        
        # Get the cross product between AP and the direction vector of the line
        # segment
        test = AP.cross(AB).fixzeros(e)
        
        # If the test result is not the zero vector then P does not lie along 
        # the line defined by the end points of the line segment
        if not (test.x == 0.0 and test.y == 0.0 and test.z == 0.0):
            #return False, None
            return False
        
        # Get the dot product between AP and AB, and AB with itself
        dotAP = AP.dot(AB)
        dotAB = AB.dot(AB)
        if dotAP <= e:
            dotAP = 0
        if dotAB <= e:
            dotAB = 0
        
        # If dotAP is equal to 0 then P matches with the start vertex of the 
        # line segment. If it is equal to dotAB then P matches with the end 
        # vertex of the line segment. If it is greater than zero but less than 
        # dotAB then P lies in the line segment between the end points of the 
        # line segment. In such cases return true. Otherwise, return false
        return dotAP == 0.0 or dotAP == dotAB or (dotAP > 0.0 and dotAP < dotAB)
    
    
    def fixzeros(self, e:int|float=1e-8) -> LineSegment:
        assert isinstance(e, (int, float)), "e not a real number"
        self.A.fixzeros(e)
        self.B.fixzeros(e)
        return self
    
    
    def midpoint(self) -> PVector:
        return self.A.copy().add(self.B).scale(0.5)
    
    
    def set(self, ax:int|float, ay:int|float, az:int|float, bx:int|float, by:int|float, bz:int|float) -> LineSegment:
        assert isinstance(ax, (int, float)), "ax not a real number"
        assert isinstance(ay, (int, float)), "ay not a real number"
        assert isinstance(az, (int, float)), "az not a real number"
        assert isinstance(bx, (int, float)), "bx not a real number"
        assert isinstance(by, (int, float)), "by not a real number"
        assert isinstance(bz, (int, float)), "bz not a real number"
        self.setA(ax, ay, az)
        self.setB(bx, by, bz)
        return self
    
    
    def setA(self, x:int|float, y:int|float, z:int|float) -> LineSegment:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.A.set(x, y, z)
        return self
    
    
    def setB(self, x:int|float, y:int|float, z:int|float) -> LineSegment:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.B.set(x, y, z)
        return self
    