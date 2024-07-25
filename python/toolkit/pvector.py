# -*- coding: utf-8 -*-

from __future__ import annotations

# A class representing points and vectors in 3D.
class PVector:
    
    # Constructor of the class.
    def __init__(self, x:int|float=0, y:int|float=0, z:int|float=0):
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.x = x
        self.y = y
        self.z = z
    
    
    # Supports operations of the form A + B, where A and B are PVectors.
    def __add__(self, V:PVector) -> PVector:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return self.copy().add(V)
    
    
    def __eq__(self, V:PVector) -> PVector:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return self.isequal(V)
    
    
    def __getattr__(self, name:str):
        return self.__dict__[f"{name}"]
    
    
    def __hash__(self):
        return hash((self.x, self.y, self.z))
    
    
    # Supports operations of the form s * A, where A is a PVector and s is a 
    # real number.
    def __mul__(self, s:int|float) -> PVector:
        assert isinstance(s, (int, float)), "s not a real number"
        return self.copy().scale(s)
    
    
    def __ne__(self, V:PVector) -> bool:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return not self.isequal(V)
    
    # Supports operations of the form -A, where A is a PVector.
    def __neg__(self) -> PVector:
        return self.copy().scale(-1)
    
    
    def __setattr__(self, name:str, value):
        assert isinstance(value, (int, float)), "value not a real number"
        self.__dict__[f"{name}"] = value
    
    
    def __repr__(self):
        return f"PVector:(x:{self.x}, y:{self.y}, z:{self.z})"
    
    
    # Supports operations of the form A * s, where A is a PVector and s is a 
    # real number.
    def __rmul__(self, s:int|float) -> PVector:
        assert isinstance(s, (int, float)), "s not a real number"
        return self.__mul__(s)
    
    
    # Supports operations of the form A - B, where A and B are PVectors.
    def __sub__(self, V:PVector) -> PVector:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return self.copy().sub(V)
    
    
    def __str__(self):
        #return f"(x:{self.x}, y:{self.y}, z:{self.z})"
        return f"({self.x}, {self.y}, {self.z})"
    
    
    # Supports operations of the form A / s, where A is a PVector and s is a 
    # real number.
    def __truediv__(self, s:int|float) -> PVector:
        assert isinstance(s, (int, float)), "s not a real number"
        return self.__mul__(1/s)
    
    
    def abs(self) -> PVector:
        self.x = abs(self.x)
        self.y = abs(self.y)
        self.z = abs(self.z)
        return self
    
    
    def add(self, V:PVector) -> PVector:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        self.x += V.x
        self.y += V.y
        self.z += V.z
        return self
    
    
    # Returns the angle, in radians, between the current object and the given 
    # PVector.
    def angle(self, V:PVector) -> float:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        import math
        return math.acos((self.dot(V)) / (self.magnitude() * V.magnitude()))
    
    
    # Returns a clone of the object. This method calls copy.
    def clone(self) -> PVector:
        return self.copy()
    
    
    # Returns a copy of the object.
    def copy(self) -> PVector:
        return PVector(self.x, self.y, self.z)
    
    
    # Returns the cross product of the object with the given one.
    def cross(self, V:PVector) -> PVector:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return PVector(self.y * V.z - self.z * V.y,
                       self.z * V.x - self.x * V.z,
                       self.x * V.y - self.y * V.x)
    
    
    def distance(self, V:PVector) -> float:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return self.copy().sub(V).magnitude()
    
    
    # Returns the dot product of the object with the given one.
    def dot(self, V:PVector) -> float:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return (self.x * V.x) + (self.y * V.y) + (self.z * V.z)
    
    
    def fixzeros(self, e:int|float=1e-8) -> PVector:
        assert isinstance(e, (int, float)), "e not a real number"
        if abs(self.x) <= e:
            self.x = 0
        if abs(self.y) <= e:
            self.y = 0
        if abs(self.z) <= e:
            self.z = 0
        return self
    
    
    def isequal(self, V:PVector) -> bool:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        return self.x == V.x and self.y == V.y and self.z == V.z
    
    
    def iszerovector(self, e:int|float=1e-8) -> bool:
        return self.x <= e and self.y <= e and self.z <= e
    
    
    # Returns the length of the object. This function calls magnitude and 
    # returns its result.
    def length(self) -> float:
        return self.magnitude()
    
    
    # Returns the magnitude of the object.
    def magnitude(self) -> float:
        return (self.x ** 2 + self.y ** 2 + self.z ** 2) ** 0.5
    
    
    # Returns the norm of the object. This function calls magnitude and returns
    # its result.
    def norm(self) -> float:
        return self.magnitude()
    
    
    # Normalizes the object and returns the pointer.
    def normalize(self) -> PVector:
        mag = self.magnitude()
        return self.scale(1 / mag) if mag != 0 else self
    
    
    # Returns a normalized copy of the object.
    def normalized(self) -> PVector:
        return self.copy().normalize()
    
    
    # Uses the Axis-Angle rotation method (AKA. Rodrigues' rotation) to rotate 
    # the object around an axis vector K some angle of rotation. The rotation 
    # formula is:
    # (V * cos(a)) + ((K cross N) * sin(a)) + (K * (K dot V) * (1 - cos(a)))
    # Description: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    def rotate(self, K:PVector, angle:int|float) -> PVector:
        assert K is not None, "Missing K"
        assert isinstance(K, PVector), "K not a PVector"
        assert isinstance(angle, (int, float)), "angle not a real number"
        import math
        sin_a = math.sin(angle)
        cos_a = math.cos(angle)
        nK = K.normalized()
        V = self.copy()
        R = (V * cos_a).add(nK.cross(V) * sin_a).add(nK * nK.dot(V) * (1.0 - cos_a))
        self.set(R.x, R.y, R.z)
        return self
    
    
    #
    def roundvalues(self, r:int=10) -> PVector:
        assert isinstance(r, int), "r not an integer"
        assert r > 0, "r not a positive integer"
        self.x = round(self.x, r)
        self.y = round(self.y, r)
        self.z = round(self.z, r)
        return self
    
    
    # Scales the object by the given factor and returns it.
    def scale(self, s:int|float) -> PVector:
        assert isinstance(s, (int, float)), "s not a real number"
        self.x *= s
        self.y *= s
        self.z *= s
        return self
    
    
    # Returns a scaled copy of the object.
    def scaled(self, s:int|float) -> PVector:
        return self.copy().scale(s)
    
    
    def set(self, x:int|float, y:int|float, z:int|float) -> PVector:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.x = x
        self.y = y
        self.z = z
        return self
    
    
    def setfromlist(self, L:list) -> PVector:
        assert L is not None, "Missing L"
        assert isinstance(L, list), "L not a list"
        n = len(L)
        assert n == 2 or n == 3, "Unsupported L length"
        assert isinstance(L[0], (int, float)), "L[0] not a real number"
        assert isinstance(L[1], (int, float)), "L[1] not a real number"
        if n == 3:
            assert isinstance(L[2], (int, float)), "L[2] not a real number"
        self.x = L[0]
        self.y = L[1]
        self.z = L[2] if n == 3 else 0
        return self
    
    
    def setfrompvector(self, P:PVector) -> PVector:
        assert P is not None, "Missing P"
        assert isinstance(P, PVector), "P not a PVector"
        self.x = P.x
        self.y = P.y
        self.z = P.z
        return self
    
    
    # Subtracts the object by the given one.
    def sub(self, V:PVector) -> PVector:
        assert V is not None, "Missing V"
        assert isinstance(V, PVector), "V not a PVector"
        self.x -= V.x
        self.y -= V.y
        self.z -= V.z
        return self
    
    
    def todict(self) -> dict:
        return {'x':self.x, 'y':self.y, 'z':self.z}
    
    
    def tolist(self) -> list:
        return [self.x, self.y, self.z]
    
    
    def totuple(self) -> tuple:
        return (self.x, self.y, self.z)
    