# -*- coding: utf-8 -*-

from __future__ import annotations
from toolkit.pvector import PVector

# The class to represent a plane in 3D. A plane is defined by a normal vector N
# and a point P, the implicit definition of the plane is N * (Q - P) = 0, where
# Q is any point in the plane and * is the dot product.
class Plane:
    
    # Constructor of the class.
    def __init__(self, N:PVector=PVector(), P:PVector=PVector()):
        assert N is not None, "Missing N"
        assert isinstance(N, PVector), "N not a PVector"
        assert P is not None, "Missing P"
        assert isinstance(P, PVector), "P not a PVector"
        self.normal = N.copy()
        self.point = P.copy()
    
    
    def __getattr__(self, name:str):
        return self.__dict__[f"{name}"]
    
    
    def __repr__(self):
        return f"Plane:(normal:{self.normal}, point:{self.point})"
    
    
    def __setattr__(self, name:str, value):
        assert isinstance(value, PVector), "value not a PVector"
        self.__dict__[f"{name}"] = value
    
    
    def __str__(self):
        return f"(normal:{self.normal}, point:{self.point})"
    
    
    def clone(self) -> Plane:
        return self.copy()
    
    
    def copy(self) -> Plane:
        return Plane(self.normal, self.point)
    
    
    def fixzeros(self, e:int|float=1e-8) -> Plane:
        assert isinstance(e, (int, float)), "e not a real number"
        self.normal.fixzeros(e)
        self.point.fixzeros(e)
        return self
    
    
    def ispoint(self, Q:PVector, e:int|float=1e-8, r:int=10) -> bool:
        assert Q is not None, "Missing Q"
        assert isinstance(Q, PVector), "Q not a PVector"
        assert isinstance(e, (int, float)), "e not a real number"
        assert isinstance(r, int), "r not an integer number"
        return self.pointlocation(Q, e, r) == 0
    
    
    # Checks the location of a point Q with respect of the plane. A 1 value 
    # is in the positive half space of the plane (in the direction of the 
    # normal vector), a 0 value indicates the point is at the plane, a -1 value
    # indicates the point is at the negative half space of the plane (in the 
    # opposite direction of the normal vector). The indicator comes from 
    # Dot(N, Q - P).
    def pointlocation(self, Q:PVector, e:int|float=1e-8, r:int=10) -> int:
        assert Q is not None, "Missing Q"
        assert isinstance(Q, PVector), "Q not a PVector"
        assert isinstance(e, (int, float)), "e not a real number"
        assert isinstance(r, int), "r not an integer number"
        dot = round(self.normal.dot(Q - self.point), r)
        return 0 if abs(dot) <= e else (1 if dot > 0 else -1)
    
    
    # Returns the values a, b, c, d representing the coefficients of the 
    # equation ax + by + cz = d.
    def equationcoefficients(self) -> tuple:
        return self.normal.x, self.normal.y, self.normal.z, self.normal.dot(self.point)
    
    
    def set(self, nx:int|float, ny:int|float, nz:int|float, px:int|float, py:int|float, pz:int|float) -> Plane:
        assert isinstance(nx, (int, float)), "nx not a real number"
        assert isinstance(ny, (int, float)), "ny not a real number"
        assert isinstance(nz, (int, float)), "nz not a real number"
        assert isinstance(px, (int, float)), "px not a real number"
        assert isinstance(py, (int, float)), "py not a real number"
        assert isinstance(pz, (int, float)), "pz not a real number"
        self.normal.set(nx, ny, nz)
        self.point.set(px, py, pz)
        return self
    
    
    def setnormal(self, x:int|float, y:int|float, z:int|float) -> Plane:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.normal.set(x, y, z)
        return self
    
    
    def setpoint(self, x:int|float, y:int|float, z:int|float) -> Plane:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.point.set(x, y, z)
        return self
    