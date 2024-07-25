# -*- coding: utf-8 -*-

from __future__ import annotations
from toolkit.pvector import PVector

# The class representing a ray in 3D. A ray is defined by a point P and a 
# direction vector D.
class Ray:
    
    # Constructor of the class.
    def __init__(self, P=PVector(), D=PVector()):
        assert P is not None, "Missing P"
        assert isinstance(P, PVector), "P is not a PVector"
        assert D is not None, "Missing D"
        assert isinstance(D, PVector), "D is not a PVector"
        self.point = P.copy()
        self.direction = D.copy()
    
    
    def __getattr__(self, name:str):
        return self.__dict__[f"{name}"]
    
    
    def __hash__(self):
        return hash((self.point.x, self.point.y, self.point.z, 
                     self.direction.x, self.direction.y, self.direction.z))
    
    
    def __repr__(self):
        return f"Ray:(point:{self.point}, direction:{self.direction})"
    
    
    def __setattr__(self, name:str, value):
        assert isinstance(value, PVector), "value not a PVector"
        self.__dict__[f"{name}"] = value
    
    
    def __str__(self):
        return f"(point:{self.point}, direction:{self.direction})"
    
    
    # Returns the point at parameter t along the ray.
    def at(self, t:int|float) -> PVector:
        assert isinstance(t, (int, float)), "t is not a real number"
        return self.point.copy().add(self.direction * t)
    
    
    def clone(self) -> Ray:
        return self.copy()
    
    
    def copy(self) -> Ray:
        return Ray(self.point, self.direction)
    
    
    def fixzeros(self, e:int|float=1e-8) -> Ray:
        assert isinstance(e, (int, float)), "e not a real number"
        self.point.fixzeros(e)
        self.direction.fixzeros(e)
        return self
    
    
    # Indicates whether a point P lies along the ray. If so, the second 
    # returned value is the parameter of the point.
    def ispoint(self, P:PVector, e:int|float=1e-8, r:int=10) -> tuple[bool, float]:
        assert P is not None, "P is missing"
        assert isinstance(P, PVector), "P not a PVector"
        assert isinstance(e, (int, float)), "e not a real number"
        assert isinstance(r, int), "r not an integer number"
        
        # Get the normalized vector between the start point of the ray and P. 
        # Fix its zeros
        NP = (P - self.point).normalize().fixzeros(e)
        
        # Get the normalized direction vector of the ray. Fix its zeros.
        ND = self.direction.copy().normalize().fixzeros(e)
        
        dot = round(NP.dot(ND), r)
        t = None
        
        # If the absolute value of the dot product is not equal to 1 then 
        # return false since P does not lie along the ray (if it is -1 the it 
        # lies in the opposite direction with respecto to the direction vector)
        if abs(dot) != 1.0:
            return False, t
        
        # Calculate the parameter for the point along the ray. Try different
        # coordinates if the denominator is zero for any of them
        if abs(self.direction.x) >= e:
            t = (P.x - self.point.x) / self.direction.x
        elif abs(self.direction.y) >= e:
            t = (P.y - self.point.y) / self.direction.y
        elif abs(self.direction.z) >= e:
            t = (P.z - self.point.z) / self.direction.z
        else:
            assert False, "Unexpected result"
        
        return True, t
    
    
    def normalizedirection(self) -> Ray:
        self.direction.normalize()
        return self
    
    
    def set(self, px:int|float, py:int|float, pz:int|float, dx:int|float, dy:int|float, dz:int|float) -> Ray:
        assert isinstance(px, (int, float)), "px not a real number"
        assert isinstance(py, (int, float)), "py not a real number"
        assert isinstance(pz, (int, float)), "pz not a real number"
        assert isinstance(dx, (int, float)), "dx not a real number"
        assert isinstance(dy, (int, float)), "dy not a real number"
        assert isinstance(dz, (int, float)), "dz not a real number"
        self.setpoint(px, py, pz)
        self.setdirection(dx, dy, dz)
        return self
    
    
    def setdirection(self, x:int|float, y:int|float, z:int|float) -> Ray:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.direction.set(x, y, z)
        return self
    
    
    def setpoint(self, x:int|float, y:int|float, z:int|float) -> Ray:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.point.set(x, y, z)
        return self
