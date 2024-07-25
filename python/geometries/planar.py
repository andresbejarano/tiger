# -*- coding: utf-8 -*-

from ds.vf import VF
from toolkit.pvector import PVector
import toolkit.utils as utils
import math
import types


def barrelvault(length:int|float, radius:int|float, ls:int, rs:int, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(length, (int, float)), "length not a real number"
    assert length > 0, "length is not positive"
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    assert isinstance(ls, int), "ls not an integer number"
    assert ls > 0, "ls is not positive"
    assert isinstance(rs, int), "rs not an integer number"
    assert rs > 0, "rs is not positive"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    
    vf = VF()
    nU = U.normalized()
    nV = V.normalized()
    nW = nU.cross(nV)
    halflength = length / 2
    lengthstep = length / ls
    anglestep = utils.PI / rs
    x = halflength
    
    for i in range(ls + 1):
        angle = utils.PI
        for j in range(rs + 1):
            V = (x * nU) - (radius * nW)
            V.rotate(nU, angle)
            vf.addvertexfrompvector(C + V)
            angle -= anglestep
        x -= lengthstep
    
    for i in range(ls):
        for j in range(rs):
            v0 = (i * (rs + 1)) + j
            v1 = v0 + 1
            v2 = ((i + 1) * (rs + 1)) + j + 1
            v3 = v2 - 1
            vf.addface([v0, v1, v2, v3])
    
    return vf



def ellipticparaboloid(a:int|float, b:int|float, width:int|float, height:int|float, ws:int, hs:int, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(a, (int, float)), "a not a real number"
    assert isinstance(b, (int, float)), "b not a real number"
    assert isinstance(width, (int, float)), "width not a real number"
    assert width > 0, "width not a positive number"
    assert isinstance(height, (int, float)), "height not a real number"
    assert height > 0, "height not a positive number"
    assert isinstance(ws, int), "ws not an integer number"
    assert ws > 0, "ws not a positive number"
    assert isinstance(hs, int), "hs not an integer number"
    assert hs > 0, "hs not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    return grid(width, height, ws, hs, lambda x, y : (x ** 2 / a ** 2) + (y ** 2 / b ** 2), C, U, V)




def equilateraltriangle(length:int|float=1, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(length, (int, float)), "length not a real number"
    assert length > 0, "length not a positive real number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    vf = VF()
    halflength = length / 2
    height = (length * (3 ** 0.5)) / 2
    vf.addvertexfrompvector(C + (halflength * U) - ((height / 3) * V))
    vf.addvertexfrompvector(C + ((2 * height) / 3) * V)
    vf.addvertexfrompvector(C + (-halflength * U) - ((height / 3) * V))
    vf.addface([0, 1, 2])
    return vf


def grid(width:int|float, height:int|float, ws:int, hs:int, func:types.FunctionType, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(width, (int, float)), "width not a real number"
    assert width > 0, "width not a positive number"
    assert isinstance(height, (int, float)), "height not a real number"
    assert height > 0, "height not a positive number"
    assert isinstance(ws, int), "ws not an integer number"
    assert ws > 0, "ws not a positive number"
    assert isinstance(hs, int), "hs not an integer number"
    assert hs > 0, "hs not a positive number"
    assert func is not None, "Missing func"
    assert isinstance(func, types.FunctionType), "func not a function"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    vf = VF()
    halfwidth = width / 2
    halfheight = height / 2
    widthstep = width / ws
    heightstep = height / hs
    N = U.cross(V)
    
    for i in range(hs + 1):
        for j in range(ws + 1):
            u = -halfwidth + (j * widthstep)
            v = -halfheight + (i * heightstep)
            w = func(u, v)
            vf.addvertexfrompvector(C + (u * U) + (v * V) + (w * N))
    
    for i in range(hs):
        for j in range(ws):
            v0 = j + (i * (ws + 1))
            v1 = v0 + 1
            v2 = (j + 1) + ((i + 1) * (ws + 1))
            v3 = v2 - 1
            vf.addface([v0, v1, v2, v3])
            
    return vf


def monkeysaddle(width:int|float, height:int|float, ws:int, hs:int, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(width, (int, float)), "width not a real number"
    assert width > 0, "width not a positive number"
    assert isinstance(height, (int, float)), "height not a real number"
    assert height > 0, "height not a positive number"
    assert isinstance(ws, int), "ws not an integer number"
    assert ws > 0, "ws not a positive number"
    assert isinstance(hs, int), "hs not an integer number"
    assert hs > 0, "hs not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    return grid(width, height, ws, hs, lambda x, y : (x ** 3) - (3 * x * y ** 2), C, U, V)


def regularpolygon(sides:int, sidelength:int|float=1, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(sides, int), "sides not an integer number"
    assert sides > 2, "sides must be greater than 3"
    assert isinstance(sidelength, (int, float)), "sidelength not a real number"
    assert sidelength > 0, "sidelength not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    vf = VF()
    angle = 0
    anglestep = utils.TWO_PI / sides
    radius = sidelength / (2 * math.cos(utils.PI / sides))
    
    for i in range(sides):
        P = radius * (math.cos(angle) * U + math.sin(angle) * V)
        vf.addvertexfrompvector(P)
        angle += anglestep
    
    vf.addface([i for i in range(sides)])
    return vf


def saddle(width:int|float, height:int|float, ws:int, hs:int, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(width, (int, float)), "width not a real number"
    assert width > 0, "width not a positive number"
    assert isinstance(height, (int, float)), "height not a real number"
    assert height > 0, "height not a positive number"
    assert isinstance(ws, int), "ws not an integer number"
    assert ws > 0, "ws not a positive number"
    assert isinstance(hs, int), "hs not an integer number"
    assert hs > 0, "hs not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    return grid(width, height, ws, hs, lambda x, y : (x ** 2) - (y ** 2), C, U, V)


def square(length:int|float=1, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(length, (int, float)), "length not a real number"
    assert length > 0, "length not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    vf = VF()
    halflength = length / 2
    hU = U * halflength
    hV = V * halflength
    vf.addvertexfrompvector(C +  hU + hV)
    vf.addvertexfrompvector(C + -hU + hV)
    vf.addvertexfrompvector(C + -hU - hV)
    vf.addvertexfrompvector(C +  hU - hV)
    vf.addface([0, 1, 2, 3])
    return vf


def squaredgrid(width:int|float, height:int|float, ws:int, hs:int, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(width, (int, float)), "width not a real number"
    assert width > 0, "width not a positive number"
    assert isinstance(height, (int, float)), "height not a real number"
    assert height > 0, "height not a positive number"
    assert isinstance(ws, int), "ws not an integer number"
    assert ws > 0, "ws not a positive number"
    assert isinstance(hs, int), "hs not an integer number"
    assert hs > 0, "hs not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    return grid(width, height, ws, hs, lambda x, y : 0, C, U, V)


def wave(width:int|float, height:int|float, ws:int, hs:int, C:PVector=PVector(), U:PVector=PVector(1), V:PVector=PVector(y=1)) -> VF:
    assert isinstance(width, (int, float)), "width not a real number"
    assert width > 0, "width not a positive number"
    assert isinstance(height, (int, float)), "height not a real number"
    assert height > 0, "height not a positive number"
    assert isinstance(ws, int), "ws not an integer number"
    assert ws > 0, "ws not a positive number"
    assert isinstance(hs, int), "hs not an integer number"
    assert hs > 0, "hs not a positive number"
    assert C is not None, "Missing C"
    assert isinstance(C, PVector), "C not a PVector"
    assert U is not None, "Missing U"
    assert isinstance(U, PVector), "U not a PVector"
    assert V is not None, "Missing V"
    assert isinstance(V, PVector), "V not a PVector"
    return grid(width, height, ws, hs, lambda x, y : math.sin(x) * math.cos(y), C, U, V)
    