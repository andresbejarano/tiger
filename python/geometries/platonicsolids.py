# -*- coding: utf-8 -*-

from ds.vf import VF


def tetrahedron(radius:int|float=1) -> VF:
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    c = 1.0 / (2.0 ** 0.5)
    vf = VF()
    vf.addvertex(1.0,  0.0, -c)
    vf.addvertex(-1.0,  0.0, -c)
    vf.addvertex(0.0,  1.0,  c)
    vf.addvertex(0.0, -1.0,  c)
    vf.addface(( 0, 1, 2 ))
    vf.addface(( 0, 2, 3 ))
    vf.addface(( 1, 0, 3 ))
    vf.addface(( 1, 3, 2 ))
    vf.normalizevertices()
    vf.scale(radius)
    return vf


def hexahedron(radius:int|float=1) -> VF:
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    vf = VF()
    vf.addvertex(1.0, 1.0, 1.0)
    vf.addvertex(1.0, 1.0, -1.0)
    vf.addvertex(1.0, -1.0, 1.0)
    vf.addvertex(1.0, -1.0, -1.0)
    vf.addvertex(-1.0, 1.0, 1.0)
    vf.addvertex(-1.0, 1.0, -1.0)
    vf.addvertex(-1.0, -1.0, 1.0)
    vf.addvertex(-1.0, -1.0, -1.0)
    vf.addface(( 0, 4, 6, 2 ))
    vf.addface(( 1, 3, 7, 5 ))
    vf.addface(( 0, 2, 3, 1 ))
    vf.addface(( 6, 4, 5, 7 ))
    vf.addface(( 2, 6, 7, 3 ))
    vf.addface(( 4, 0, 1, 5 ))
    vf.normalizevertices()
    vf.scale(radius)
    return vf


def octahedron(radius:int|float=1) -> VF:
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    vf = VF()
    vf.addvertex(1.0, 0.0, 0.0)
    vf.addvertex(-1.0, 0.0, 0.0)
    vf.addvertex(0.0, 1.0, 0.0)
    vf.addvertex(0.0, -1.0, 0.0)
    vf.addvertex(0.0, 0.0, 1.0)
    vf.addvertex(0.0, 0.0, -1.0)
    vf.addface(( 0, 2, 4 ))
    vf.addface(( 2, 1, 4 ))
    vf.addface(( 1, 3, 4 ))
    vf.addface(( 3, 0, 4 ))
    vf.addface(( 5, 0, 3 ))
    vf.addface(( 5, 2, 0 ))
    vf.addface(( 5, 1, 2 ))
    vf.addface(( 5, 3, 1 ))
    vf.normalizevertices()
    vf.scale(radius)
    return vf


def dodecahedron(radius:int|float=1) -> VF:
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    c = ((5.0 ** 0.5) - 1.0) / 2.0
    c1 = 1.0 + c
    c2 = 1.0 - c * c
    vf = VF()
    vf.addvertex(1.0, 1.0, 1.0)
    vf.addvertex(1.0, 1.0, -1.0)
    vf.addvertex(1.0, -1.0, 1.0)
    vf.addvertex(1.0, -1.0, -1.0)
    vf.addvertex(-1.0, 1.0, 1.0)
    vf.addvertex(-1.0, 1.0, -1.0)
    vf.addvertex(-1.0, -1.0, 1.0)
    vf.addvertex(-1.0, -1.0, -1.0)
    vf.addvertex(0.0, c1, c2)
    vf.addvertex(0.0, c1, -c2)
    vf.addvertex(0.0, -c1, c2)
    vf.addvertex(0.0, -c1, -c2)
    vf.addvertex(c1, c2, 0.0)
    vf.addvertex(c1, -c2, 0.0)
    vf.addvertex(-c1, c2, 0.0)
    vf.addvertex(-c1, -c2, 0.0)
    vf.addvertex(c2, 0.0, c1)
    vf.addvertex(c2, 0.0, -c1)
    vf.addvertex(-c2, 0.0, c1)
    vf.addvertex(-c2, 0.0, -c1)
    vf.addface(( 2, 16, 18, 6, 10 ))
    vf.addface(( 6, 18, 4, 14, 15 ))
    vf.addface(( 4, 18, 16, 0, 8 ))
    vf.addface(( 0, 16, 2, 13, 12 ))
    vf.addface(( 1, 9, 8, 0, 12 ))
    vf.addface(( 5, 14, 4, 8, 9 ))
    vf.addface(( 7, 15, 14, 5, 19 ))
    vf.addface(( 11, 10, 6, 15, 7 ))
    vf.addface(( 3, 13, 2, 10, 11 ))
    vf.addface(( 17, 1, 12, 13, 3 ))
    vf.addface(( 17, 3, 11, 7, 19 ))
    vf.addface(( 5, 9, 1, 17, 19 ))
    vf.normalizevertices()
    vf.scale(radius)
    return vf


def icosahedron(radius:int|float=1) -> VF:
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    c = (1.0 + (5.0 ** 0.5)) / 2.0
    vf = VF()
    vf.addvertex(0.0, 1.0, c)
    vf.addvertex(0.0, 1.0, -c)
    vf.addvertex(0.0, -1.0, c)
    vf.addvertex(0.0, -1.0, -c)
    vf.addvertex(1.0, c, 0.0)
    vf.addvertex(1.0, -c, 0.0)
    vf.addvertex(-1.0, c, 0.0)
    vf.addvertex(-1.0, -c, 0.0)
    vf.addvertex(c, 0.0, 1.0)
    vf.addvertex(c, 0.0, -1.0)
    vf.addvertex(-c, 0.0, 1.0)
    vf.addvertex(-c, 0.0, -1.0)
    vf.addface(( 2, 0, 10 ))
    vf.addface(( 2, 8, 0 ))
    vf.addface(( 8, 4, 0 ))
    vf.addface(( 4, 6, 0 ))
    vf.addface(( 6, 10, 0 ))
    vf.addface(( 10, 7, 2 ))
    vf.addface(( 7, 5, 2 ))
    vf.addface(( 5, 8, 2 ))
    vf.addface(( 9, 8, 5 ))
    vf.addface(( 9, 4, 8 ))
    vf.addface(( 1, 4, 9 ))
    vf.addface(( 1, 6, 4 ))
    vf.addface(( 11, 6, 1 ))
    vf.addface(( 11, 10, 6 ))
    vf.addface(( 11, 7, 10 ))
    vf.addface(( 3, 5, 7 ))
    vf.addface(( 9, 5, 3 ))
    vf.addface(( 3, 7, 11 ))
    vf.addface(( 1, 3, 11 ))
    vf.addface(( 3, 1, 9 ))
    vf.normalizevertices()
    vf.scale(radius)
    return vf


functions = {
    'tetrahedron': tetrahedron, 
    'hexahedron': hexahedron, 
    'octahedron': octahedron, 
    'dodecahedron': dodecahedron, 
    'icosahedron': icosahedron, 
    'cube': hexahedron}


def byname(name:str, radius:int|float=1) -> VF:
    assert name is not None, "Missing name"
    assert isinstance(name, str), "name not a string"
    assert isinstance(radius, (int, float)), "radius not a real number"
    assert radius > 0, "radius is not positive"
    func = functions[name.lower()]
    if func:
        return func(radius)
    assert False, "Unexpected behavior"
    
