# -*- coding: utf-8 -*-

from ds.vf import VF
from toolkit.pvector import PVector
import toolkit.utils as utils
import math

def torus(R, r, Rs, rs, C=PVector(), U=PVector(1), V=PVector(y=1)) -> VF:
    assert isinstance(R, (int, float)), "R not a real number"
    assert R > 0, "R is not positive"
    assert isinstance(r, (int, float)), "r not a real number"
    assert r > 0, "r is not positive"
    assert isinstance(Rs, int), "Rs not an integer number"
    assert Rs > 0, "Rs is not positive"
    assert isinstance(rs, int), "rs not an integer number"
    assert rs > 0, "rs is not positive"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    vf = VF()
    W = U.cross(V)
    ustep = utils.TWO_PI / Rs
    vstep = utils.TWO_PI / rs
    
    # Calculate the torus vertices
    for l in range(Rs):
        u = l * ustep
        sin_u = math.sin(u)
        cos_u = math.cos(u)
        T1 = C + (U * (R * cos_u)) + (V * (R * sin_u))
        T2 = (U * -sin_u) + (V * cos_u)
        T3 = T2.cross(W).normalize()
        
        for m in range(rs):
            v = m * vstep
            sin_v = math.sin(v)
            cos_v = math.cos(v)
            Q = T1 + (T3 * (r * cos_v)) + (W * (r * sin_v))
            vf.addvertexfrompvector(Q)
    
    # Define the torus faces
    for l in range(Rs):
        for m in range(rs):
            
            # Determine the vertex indices of the current face according to the
            # current vertex indices values
            if l < Rs - 1:
                if m < rs - 1:
                    
                    # This is the common case: All squares not at the end of 
                    # the ring or at the last ring. This case happens mostly of
                    # the times.
                    v0 = (l * rs) + m
                    v1 = ((l + 1) * rs) + m
                    v2 = v1 + 1
                    v3 = v0 + 1
                else:
                    
                    # This is the last square in a minor ring (not the last 
                    # one). This case happens once per minor ring.
                    v0 = (l * rs) + m
                    v1 = ((l + 1) * rs) + m
                    v2 = v0 + 1
                    v3 = v0 - rs + 1
            else:
                if m < rs - 1:
                    
                    # This is a square (not the last) in the last minor ring. 
                    # It uses the first vertices generated for the torus. This 
                    # case happens mostly in the last minor ring.
                    v0 = (l * rs) + m
                    v1 = m
                    v2 = v1 + 1
                    v3 = v0 + 1
                else:
                    
                    # This is the very last square in the last minor ring. This
                    # case happens one time.
                    v0 = (l * rs) + m
                    v1 = m
                    v2 = 0
                    v3 = v0 - rs + 1
            
            vf.addface([v0, v1, v2, v3])
    
    return vf