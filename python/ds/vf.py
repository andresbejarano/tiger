# -*- coding: utf-8 -*-

from __future__ import annotations
from toolkit.pvector import PVector

class VF:
    
    def __init__(self):
        self.V = list()
        self.F = list()
    
    
    def __assertIndices(self, indices:list|tuple):
        assert indices is not None, "Missing indices"
        assert isinstance(indices, (list, tuple)), "indices not a list or tuple"
        assert len(indices) >= 3, "indices length is not greater than or equal to 3"
        nV = len(self.V)
        for i in indices:
            assert isinstance(i, int), "indices value not an integer number"
            assert i >= 0 and i < nV, "indices value not in range"
    
    
    def __repr__(self):
        return f"VF:(V:{self.V}, F:{self.F})"
    
    
    def __str__(self):
        return f"(V:{self.V},\nF:{self.F})"
    
    
    def addface(self, indices:list|tuple) -> int:
        assert indices is not None, "Missing indices"
        assert isinstance(indices, (list, tuple)), "indices not a list or tuple"
        assert len(indices) >= 3, "indices length is not greater than or equal to 3"
        self.__assertIndices(indices)
        self.F.append([i for i in indices])
        return len(self.F) - 1
    
    
    def addvertex(self, x:int|float, y:int|float, z:int|float=0) -> int:
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        self.V.append([x, y, z])
        return len(self.V) - 1
    
    
    def addvertexfromlist(self, L:list) -> int:
        assert L is not None, "Missing L"
        n = len(L)
        assert n == 2 or n == 3, "Length of list neither 2 nor 3"
        return self.addvertex(L[0], L[1], L[2] if n == 3 else 0)
    
    
    def addvertexfrompvector(self, P:PVector) -> int:
        assert P is not None, "Missing P"
        assert isinstance(P, PVector), "P not a PVector"
        return self.addvertex(P.x, P.y, P.z)
    
    
    # Returns the centroid of the face at the given index.
    def centroid(self) -> PVector:
        C = PVector()
        P = PVector()
        nvertices = len(self.V)
        
        for v in self.V:
            P.setfromlist(v)
            C.add(P)
        
        return C.scale(1 / nvertices)
    
    
    def clone(self) -> VF:
        return self.copy()
    
    
    def copy(self) -> VF:
        vf = VF()
        for v in self.V:
            vf.addvertex(v[0], v[1], v[2])
        for f in self.F:
            vf.addface([i for i in f]) 
        return vf
    
    
    # Returns the number of edges of the geometric object.
    def countedges(self) -> int:
        return len(self.edges())
    
    
    # Returns the number of faces of the geometric object.
    def countfaces(self) -> int:
        return len(self.F)
    
    
    # Returns the number of vertices of the geometric object.
    def countvertices(self) -> int:
        return len(self.V)
    
    
    # Returns a list with the edges of the geometric object. These edges are
    # the tuples with the vertex indices required to draw the line segments
    # representing the edges of the geometry.
    def edges(self) -> set:
        edges = set()
        
        # Traverse through the faces.
        nF = self.countfaces()
        for fIdx in range(nF):
            
            nV = len(self.F[fIdx]) - 1
            for vIdx in range(nV):
                i = self.F[fIdx][vIdx]
                j = self.F[fIdx][vIdx + 1]
                edge = (min(i, j), max(i, j))
                edges.add(edge)
            
            # Check the last edge of the current face. That is the line segment
            # defined by the first and last vertices of the face.
            i = self.F[fIdx][0]
            j = self.F[fIdx][nV]
            edge = (min(i, j), max(i, j))
            edges.add(edge)
        
        return edges
    
    
    # Returns the centroid of the face at the given index.
    def facecentroid(self, idx:int) -> PVector:
        assert isinstance(idx, int), "idx not an integer"
        assert 0 <= idx and idx < len(self.F), "idx out of bounds"
        C = PVector()
        P = PVector()
        nsides = len(self.F[idx])
        
        for i in range(nsides):
            P.setfromlist(self.V[self.F[idx][i]])
            C.add(P)
        
        return C.scale(1 / nsides)
    
    
    # Flips (i.e., reverts) the face indices of the geometry.
    def flip(self) -> VF:
        for i in range(len(self.F)):
            self.F[i] = [self.F[i][j] for j in range(len(self.F[i])-1, -1, -1)]
        return self
    
    
    # Returns the normal vector of the face at the given index.
    def normal(self, idx:int) -> PVector:
        assert isinstance(idx, int), "idx not an integer"
        assert 0 <= idx and idx < len(self.F), "idx out of bounds"
        
        N = PVector()
        nsides = len(self.F[idx])
        Pprev = PVector()
        Pcurr = PVector()
        Pnext = PVector()
        
        # Don't forget the normal around the first vertex of the face
        Pprev.setfromlist(self.V[self.F[idx][nsides - 1]])
        Pcurr.setfromlist(self.V[self.F[idx][0]])
        Pnext.setfromlist(self.V[self.F[idx][1]])
        T1 = Pcurr - Pprev
        T2 = Pnext - Pcurr
        N.add(T1.cross(T2))
        
        # Traverse through the indices of the given face, starting with index 1
        for i in range(1, nsides - 1):
            Pprev.setfromlist(self.V[self.F[idx][i - 1]])
            Pcurr.setfromlist(self.V[self.F[idx][i]])
            Pnext.setfromlist(self.V[self.F[idx][i + 1]])
            T1 = Pcurr - Pprev
            T2 = Pnext - Pcurr
            N.add(T1.cross(T2))
        
        # Don't forget the normal around the last vertex of the face
        Pprev.setfromlist(self.V[self.F[idx][nsides - 2]])
        Pcurr.setfromlist(self.V[self.F[idx][nsides - 1]])
        Pnext.setfromlist(self.V[self.F[idx][0]])
        T1 = Pcurr - Pprev
        T2 = Pnext - Pcurr
        N.add(T1.cross(T2))
        
        if nsides > 1:
            N.scale(1 / nsides)
        
        return N
    
    
    # 
    def normalizevertices(self) -> VF:
        for i in range(len(self.V)):
            norm = ((self.V[i][0] ** 2) + (self.V[i][1] ** 2) + (self.V[i][2] ** 2)) ** 0.5
            if norm > 0:
                self.V[i][0] /= norm
                self.V[i][1] /= norm
                self.V[i][2] /= norm
        return self
    
    
    def scale(self, s:int|float) -> VF:
        assert isinstance(s, (int, float)), "s not a real number"
        for i in range(len(self.V)):
            self.V[i][0] *= s
            self.V[i][1] *= s
            self.V[i][2] *= s
        return self
