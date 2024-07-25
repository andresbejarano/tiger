# -*- coding: utf-8 -*-

from __future__ import annotations
from ds.vf import VF
from toolkit.pvector import PVector
from toolkit.plane import Plane

# The class representing a vertex in a DCEL.
class Vertex:
    
    # Constructor of the class.
    def __init__(self, x:int|float=0, y:int|float=0, z:int|float=0, halfedge:Halfedge=None):
        assert isinstance(x, (int, float)), "x not a real number"
        assert isinstance(y, (int, float)), "y not a real number"
        assert isinstance(z, (int, float)), "z not a real number"
        if halfedge is not None:
            assert isinstance(halfedge, Halfedge), "halfedge not a Halfedge"
        self.coords = PVector(x, y, z)
        self.halfedge = halfedge
    
    
    def checkconsistency(self):
        assert hasattr(self, 'coords'), "No coords attribute"
        assert self.coords is not None, "Missing coords"
        assert isinstance(self.coords, PVector), "coords not a PVector"
        assert hasattr(self, 'halfedge'), "No halfedge attribute"
        assert self.halfedge is not None, "Missing halfedge"
        assert isinstance(self.halfedge, Halfedge), "halfedge not a Halfedge"
        assert self == self.halfedge.start, "Inconsistent start vertex"
    
    
    # Checks whether the vertex is internal or not. A vertex is internal if its
    # respective star cannot be found due to missing connections between half
    # edges.
    def isinternal(self) -> bool:
        assert hasattr(self, 'halfedge'), "vertex has no halfedge attribute"
        h = self.halfedge
        
        while True:
            
            assert h.twin is not None, "Missing halfedge,twin"
            if h.twin.next is None:
                return False
            
            h = h.twin.next
            if h == self.halfedge:
                return True



# The class representing a face in a DCEL.
class Face:
    
    # Constructor of the class.
    def __init__(self, halfedge:Halfedge=None):
        if halfedge is not None:
            assert isinstance(halfedge, Halfedge), "halfedge not a Halfedge"
        self.halfedge = halfedge
    
    
    # Calculates the area of the face.
    def area(self) -> float:
        sumareas = 0
        A = self.halfedge.start.coords.copy()
        currenthalfedge = self.halfedge.next
        
        while True:
            
            # Get the vectors from A to the end points of the half edge
            AB = currenthalfedge.start.coords - A
            AC = currenthalfedge.twin.start.coords - A
            
            # Calculate the area of triangle between A and the end points of 
            # the current half edge
            sumareas += AB.cross(AC).magnitude()
            
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == self.halfedge.previous:
                break
        
        return sumareas / 2
    
    
    # Calculates and returns the centroid (i.e., face vertices' arithemtic 
    # mean) of the face.
    def centroid(self) -> PVector:
        assert self.halfedge is not None, "Missing halfedge"
        assert isinstance(self.halfedge, Halfedge), "Face halfedge not a Halfedge"
        C = PVector()
        nsides = 0
        currenthalfedge = self.halfedge
        
        while True:
            nsides += 1
            C.add(currenthalfedge.start.coords)
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == self.halfedge:
                break
        
        if nsides > 1:
            C.scale(1 / nsides)
        
        return C
    
    
    def checkconsistency(self):
        assert hasattr(self, 'halfedge'), "No halfedge attribute"
        assert self.halfedge is not None, "Missing halfedge"
        assert isinstance(self.halfedge, Halfedge), "halfedge not a Halfedge"
        currenthalfedge = self.halfedge
        while True:
            assert self == currenthalfedge.face, "Inconsistent face"
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == self.halfedge:
                break
    
    
    # Checks whether the face has all of its neighbors. That is, the faces 
    # incident to the twin half edges of the face exist.
    def hasallneighbors(self) -> bool:
        currenthalfedge = self.halfedge
        while True:
            
            if currenthalfedge.twin.face is None:
                return False
            
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == self.halfedge:
                break
        return True
    
    
    def hasevennumberofsides(self) -> bool:
        return self.numedges() % 2 == 0
    
    
    def isatboundary(self) -> bool:
        return not self.hasallneighbors()
    
    
    # Checks if the vertices incident to the face are in the given plane.
    def iscoplanar(self, plane:Plane, e:int|float=1e-8, r:int=10) -> bool:
        assert plane is not None, "Missing plane"
        assert isinstance(plane, Plane), "plane not a Plane"
        assert isinstance(e, (int, float)), "e not a real number"
        assert isinstance(r, int), "r not an integer number"
        currenthalfedge = self.halfedge
        
        while True:
            
            if not plane.ispoint(currenthalfedge.start.coords, e, r):
                return False
            
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == self.halfedge:
                break
        
        return True
    
    
    # Calculates the normal vector of the face.
    def normal(self) -> PVector:
        assert self.halfedge is not None, "Missing halfedge"
        assert isinstance(self.halfedge, Halfedge), "Face halfedge not a Halfedge"
        N = PVector()
        nsides = 0
        P = self.halfedge.start.coords.copy()
        currenthalfedge = self.halfedge.next
        
        while True:
            T1 = currenthalfedge.start.coords - P
            T2 = currenthalfedge.twin.start.coords - P
            N += T1.cross(T2)
            nsides += 1
            
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == self.halfedge.previous:
                break
        
        if nsides > 1:
            N.scale(1 / nsides)
        
        return N
    
    
    # Returns the number of sides (i.e., incident half edges) of the face.
    def numedges(self) -> int:
        assert hasattr(self, 'halfedge') and self.halfedge is not None, 'Missing halfedge'
        count = 0
        currentHalfedge = self.halfedge
        
        while True:
            count += 1
            currentHalfedge = currentHalfedge.next
            if currentHalfedge == self.halfedge:
                break
            
        return count
    
    
    # Returns the number of triangles required to represent the face.
    def numtriangles(self) -> int:
        return self.numedges() - 2
    
    
    # Returns the number of incident vertices of the face. It is the same as 
    # the number of sides of the face.
    def numvertices(self) -> int:
        return self.numedges()
    
    
    # Returns the geometry of the face as a VF object.
    def tovf(self) -> VF:
        vf = VF()
        
        nvertices = 0
        h = self.halfedge
        while True:
            
            vf.addvertexfrompvector(h.start.coords)
            nvertices += 1
            
            h = h.next
            if h == self.halfedge:
                break
        
        vf.addface([i for i in range(nvertices)])
        
        return vf
        


# The class representing a half edge in a DCEL.
class Halfedge:
    
    # Constructor of the class.
    def __init__(self):
        self.start = None
        self.previous = None
        self.next = None
        self.twin = None
        self.face = None
    
    
    def checkconsistency(self):
        assert hasattr(self, 'start'), "No start attribute"
        assert self.start is not None, "Missing start vertex"
        assert isinstance(self.start, Vertex), "start not a Vertex"
        assert hasattr(self, 'twin'), "No twin attribute"
        assert self.twin is not None, "Missing twin half edge"
        assert isinstance(self.twin, Halfedge), "twin not a Halfedge"
        assert self == self.twin.twin, "Itself is not twin's twin"
        assert self.start == self.twin.twin.start, "Start is not twin's twin start"
        assert hasattr(self, 'face'), "No face attribute"
        assert hasattr(self, 'previous'), "No previous attribute"
        assert hasattr(self, 'next'), "No next attribute"
        if self.face is not None:
            assert isinstance(self.face, Face), "face not a Face"
            assert self == self.previous.next, "Itself is not previous' next"
            assert self == self.next.previous, "Itself is not next's previous"
            assert self.face == self.previous.face, "Different face than previous'"
            assert self.face == self.next.face, "Different face than next's"
    
    
    # Returns the length of the half edge.
    def length(self):
        return self.start.coords.distance(self.twin.start.coords)
    
    
    # Returns the midpoint of the half edge.
    def midpoint(self) -> PVector:
        assert hasattr(self, 'start'), "No start attribute"
        assert hasattr(self, 'twin'), "No twin attribute"
        assert hasattr(self.twin, 'start'), "No twin.start attribute"
        return self.start.coords.copy().add(self.twin.start.coords).scale(0.5)
    
    
    # Returns the normal vector of the half edge. The normal vector is the 
    # average of the normal vectors from the incident face and the twin's face.
    def normal(self) -> PVector:
        assert hasattr(self, 'face'), "No face attribute"
        assert hasattr(self.twin, 'face'), "No face attribute"
        N1 = self.face.normal() if self.face is not None else PVector()
        n = 1 if self.face is not None else 0
        N2 = self.twin.face.normal() if self.twin.face is not None else PVector()
        n += 1 if self.twin.face is not None else 0
        assert n > 0, "Unexpected behavior with the number of incident faces"
        return N1.add(N2).scale(1 / n)



# The class representing a DCEL.
class DCEL:
    
    # Constructor of the class.
    def __init__(self, vf):
        assert vf is not None, "Missing vf"
        assert isinstance(vf, VF), "vf not a VF"
        self.set(vf)
    
    
    # Returns the area of the geometry.
    def area(self) -> float:
        assert hasattr(self, 'F'), "No F attribute"
        assert self.F is not None, "Missing F"
        sumareas = 0
        for f in self.F:
            sumareas += f.area()
        return sumareas
    
    
    # Checks that all faces have an even number of sides/edges.
    def arefacesevensided(self) -> bool:
        for f in self.F:
            if f.numedges() % 2 != 0:
                return False
        return True
    
    
    def clear(self):
        if hasattr(self, 'V') and self.V is not None:
            del self.V
        if hasattr(self, 'F') and self.F is not None:
            del self.F
        if hasattr(self, 'H') and self.H is not None:
            del self.H
        
        #self.V = list()
        #self.F = list()
        #self.H = list()
    
    
    # WARNING: This function won't work for certain geometries.
    def closeloops(self):
        assert hasattr(self, 'H') and self.H is not None, 'No half edges'
        
        # Traverse through the half edges
        nH = len(self.H)
        for itH in range(nH):
            
            if self.H[itH].face is not None:
                continue
            
            for jtH in range(itH + 1, nH):
                
                if self.H[jtH] is not None:
                    continue
                
                if (self.H[itH].twin.start == self.H[jtH].start 
                    and self.H[jtH].next is None 
                    and self.H[jtH].previous is None):
                    self.H[itH].next = self.H[jtH]
                    self.H[jtH].previous = self.H[itH]
                
                elif (self.H[itH].start == self.H[jtH].twin.start and self.H[itH].previous is None and self.H[jtH].next is None):
                    self.H[itH].previous = self.H[jtH]
                    self.H[jtH].next = self.H[itH]
    
    
    # Returns a VF object with the dual of the geometry.
    # Note: This funtion only considers the dual face of the internal vertices 
    # of the geometry. Vertices with no star (i.e., not internal) won't have 
    # the corresponding dual face in the returned VF object.
    def dual(self, e=1e-8) -> VF:
        assert isinstance(e, (int, float)), "e not a real number"
        
        vf = VF()

        # A dictionary for mapping face centroids to their vf indices
        fToI = {}
        
        # Traverse through the vertices of the geometry. Calculate the dual 
        # face indices only for internal vertices
        for v in self.V:
            
            # Skip vertices at the boundary of the geometry
            if not v.isinternal():
                continue
            
            indices = []
            h = v.halfedge
            
            # Traverse through the star of half edges incident to the current
            # vertex, and use their face information to determine the indices
            # for the respective dual face
            while True:
                
                fIdx = -1
                if h.face in fToI:
                    fIdx = fToI[h.face]
                else:
                    fIdx = vf.addvertexfrompvector(h.face.centroid())
                    fToI[h.face] = fIdx
                
                indices.insert(0, fIdx)
                h = h.twin.next
                if h == v.halfedge:
                    break
            
            vf.addface(indices)
        
        return vf
    
    
    # Returns a VF object with the geometry of the subdivided faces by 
    # quadrilaterals. Each original face is quadrangulated using the centroid, 
    # edge midpoints, and an incident vertex.
    def facemidpointsubdivision(self, e:int|float=1e-8) -> VF:
        assert isinstance(e, (int, float)), "e not a real number"
        
        vf = VF()
        vToI = {}
        eToI = {}
        
        # Traverse through the faces of the geometry
        for f in self.F:
            
            # Calculate the centroid of the current face, add it to the vf, and
            # store its index in the vToI dictionary
            C = f.centroid().fixzeros(e)
            v0 = vf.addvertexfrompvector(C)
            vToI[str(C)] = v0
            
            currenthalfedge = f.halfedge
            
            # Get the index for the start end point of the previous half edge
            v1key = str(currenthalfedge.previous.start.coords)
            if v1key in vToI:
                v1 = vToI[v1key]
            else:
                v1 = vf.addvertexfrompvector(currenthalfedge.previous.start.coords)
                vToI[v1key] = v1
            
            # Get the index for the other end point of the previous half edge
            v2key = str(currenthalfedge.start.coords)
            if v2key in vToI:
                v2 = vToI[v2key]
            else:
                v2 = vf.addvertexfrompvector(currenthalfedge.start.coords)
                vToI[v2key] = v2
            
            # Get the index of the midpoint of the previous half edge
            e0key = (min(v1, v2), max(v1, v2))
            if e0key in eToI:
                v4 = eToI[e0key]
            else:
                P = (currenthalfedge.previous.start.coords + currenthalfedge.start.coords)
                v4 = vf.addvertexfrompvector(P.scale(0.5).fixzeros(e))
                eToI[e0key] = v4
            
            # Traverse through the half edges incident to the current face
            while True:
                
                # Get the other end point of the current half edge
                v3key = str(currenthalfedge.twin.start.coords)
                if v3key in vToI:
                    v3 = vToI[v3key]
                else:
                    v3 = vf.addvertexfrompvector(currenthalfedge.twin.start.coords)
                    vToI[v3key] = v3
                
                # Get the index of the midpoint of the current half edge
                e1key = (min(v2, v3), max(v2, v3))
                if e1key in eToI:
                    v5 = eToI[e1key]
                else:
                    P = (currenthalfedge.start.coords + currenthalfedge.twin.start.coords)
                    v5 = vf.addvertexfrompvector(P.scale(0.5).fixzeros(e))
                    eToI[e1key] = v5
                
                # Add the triangular face to the vf
                vf.addface([v0, v4, v2, v5])
                
                v1 = v2
                v2 = v3
                v4 = v5
                e0key = e1key
                currenthalfedge = currenthalfedge.next
                if currenthalfedge == f.halfedge:
                    break
        
        return vf
    
    
    # Returns a VF object with the geometry of the subdivided faces by 
    # triangles. Each original face is triangulated using the centroid and the
    # endpoints of each edge.
    def facetriangulatesubdivision(self, e:int|float=1e-8) -> VF:
        assert isinstance(e, (int, float)), "e not a real number"
        vf = VF()
        
        # A dictionary to store the index per vertex
        vToI = {}
        
        # Traverse through the faces of the geometry
        for f in self.F:
            
            # Calculate the centroid of the current face, add it to the vf, and
            # store its index in the vToI dictionary
            C = f.centroid().fixzeros(e)
            v0 = vf.addvertexfrompvector(C)
            vToI[str(C)] = v0
            
            # Check the start vertex of the incident half edge is in the vToI
            # dictionary (which means it is already in the vf object). If so, 
            # use its index for the new face. Otherwise, insert it into the vf 
            # and the vToI dictionary.
            currenthalfedge = f.halfedge
            v1key = str(currenthalfedge.start.coords)
            if v1key in vToI:
                v1 = vToI[v1key]
            else:
                v1 = vf.addvertexfrompvector(currenthalfedge.start.coords)
                vToI[v1key] = v1
            
            # Traverse through the half edges incident to the current face
            while True:
                
                # Check the start vertex of the tein half edge is in the vToI
                # dictionary (which means it is already in the vf object). If 
                # so, use its index for the new face. Otherwise, insert it into
                # the vf and the vToI dictionary.
                v2key = str(currenthalfedge.twin.start.coords)
                if v2key in vToI:
                    v2 = vToI[v2key]
                else:
                    v2 = vf.addvertexfrompvector(currenthalfedge.twin.start.coords)
                    vToI[v2key] = v2
                
                # Add the triangular face to the vf
                vf.addface([v0, v1, v2])
                
                v1 = v2
                currenthalfedge = currenthalfedge.next
                if currenthalfedge == f.halfedge:
                    break
        
        return vf
    
    
    def faceuniformsubdivision(self, e:int|float=1e-8) -> VF:
        vf = VF()
        
        # A dictionary to map every vertex to their index
        vToI = {}
        
        # A dictionary to map every halfedge to its midpoint
        eToM = {}
        
        # Traverse through the half edges of the domain. For each half edge, 
        # add its start vertex coordinates to the vf object. Also, add its 
        # midpoint to the vf object.
        for h in self.H:
            
            if h.start.coords not in vToI:
                vToI[h.start.coords] = vf.addvertexfrompvector(h.start.coords)
            
            if h not in eToM:
                M = h.midpoint()
                eToM[h] = M
                eToM[h.twin] = M
                vToI[M] = vf.addvertexfrompvector(M)
        
        # Traverse through the faces of the domain and generate the new face
        # indices
        for f in self.F:
            
            # Initialize a list for storing the midpoint indices. These correspond
            # to the indices of the central face of the subdivision
            centralface = []
            
            h = f.halfedge
            while True:
                
                assert h.previous in eToM, "whoops #1"
                assert h.next in eToM, "whoops #2"
                assert h.start.coords in vToI, "whoops #3"
                
                v0 = vToI[eToM[h.previous]]
                v1 = vToI[h.start.coords]
                v2 = vToI[eToM[h]]
                vf.addface([v0, v1, v2])
                centralface.append(v2)
                
                h = h.next
                if h == f.halfedge:
                    break
        
            vf.addface(centralface)
        
        return vf
    
    
    def flip(self) -> VF:
        vf = VF()
        vToI = {}
        for v in self.V:
            vToI[v] = vf.addvertexfrompvector(v.coords)
        
        for f in self.F:
            h = f.halfedge
            indices = []
            while True:
                
                indices.insert(0, vToI[h.start])
                
                h = h.next
                if h == f.halfedge:
                    break
            
            vf.addface(indices)
        
        return vf
         
    
    
    def normalizevertices(self, radius:int|float=1, e:int|float=1e-8) -> DCEL:
        assert isinstance(radius, (int, float)), "radius not a real number"
        assert isinstance(e, (int, float)), "e not a real number"
        for v in self.V:
            v.coords.normalize().scale(radius)
        return self
    
    
    # Returns the number of faces of the geometry.
    def numfaces(self) -> int:
        assert hasattr(self, 'F'), "No F attribute"
        assert self.F is not None, "Missing F"
        return len(self.F)
    
    
    # Returns the number of half edges of the geometry.
    def numhalfedges(self) -> int:
        assert hasattr(self, 'H'), "No H attribute"
        assert self.H is not None, 'Missing H'
        return len(self.H)
    
    
    # Returns the number of triangles required to represent the faces of the
    # geometry.
    def numtriangles(self) -> int:
        assert hasattr(self, 'F'), "No F attribute"
        assert self.F is not None, 'Missing F'
        count = 0
        for f in self.F:
            count += f.numtriangles()
        return count
    
    
    # Returns the number of vertices of the geometry.
    def numvertices(self) -> int:
        assert hasattr(self, 'V'), 'No V attribute'
        assert self.V is not None, 'Mssing V'
        return len(self.V)
    
    
    # Scales the vertices of the geometry with respect of a point C.
    def scale(self, s:int|float, C:PVector=PVector()):
        assert isinstance(s, (int, float)), "s not a real number"
        assert C is not None, "Missing C"
        assert isinstance(C, PVector), "C not a PVector"
        
        for v in self.V:
            P = C + ((v.coords - C) * s)
            v.coords.set(P.x, P.y, P.z)
    
    
    # Sets the geometry of the DCEL based on the content of the given VF.
    def set(self, vf:VF):
        assert vf is not None, "Given VF object is none"
        assert isinstance(vf, VF), "vf not a VF"
        self.clear()
        
        # Traverse through the VF vertices and create the DCEL vertices
        nV = vf.countvertices()
        self.V = [0] * nV
        for i in range(nV):
            self.V[i] = Vertex(vf.V[i][0], vf.V[i][1], vf.V[i][2])
        
        # Initialize a dictionary to be used as a map for the half edges to be 
        # found by their start and end vertices indexes
        vToH = {}
        
        # Traverse through the faces
        nF = vf.countfaces()
        self.F = [0] * nF
        self.H = []
        for i in range(nF):
            self.F[i] = Face()
            
            # Traverse through the vertex indices for the current face
            indices = vf.F[i]
            nIndices = len(indices)
            for j in range(nIndices):
                
                # Determine the indices of the previous, current, and next
                # vertices
                prevIndex = indices[nIndices - 1] if j == 0 else indices[j - 1]
                currIndex = indices[j]
                nextIndex = indices[0] if j == nIndices - 1 else indices[j + 1]
                
                previous = None
                next = None
                
                prevKey = (prevIndex, currIndex)
                twinKey = (currIndex, prevIndex)
                
                if prevKey not in vToH:
                    previous = Halfedge()
                    previous.start = self.V[prevIndex]
                    previous.face = self.F[i]
                    twin = Halfedge()
                    twin.start = self.V[currIndex]
                    twin.twin = previous
                    previous.twin = twin
                    
                    if self.V[prevIndex].halfedge is None:
                        self.V[prevIndex].halfedge = previous
                    
                    if self.V[currIndex].halfedge is None:
                        self.V[currIndex].halfedge = twin
                    
                    vToH[prevKey] = previous
                    vToH[twinKey] = twin
                    
                    self.H.append(previous)
                    self.H.append(twin)
                
                else:
                    previous = vToH[prevKey]
                    
                    assert previous.start.halfedge is not None, "Whoopsie #1"
                    assert previous.twin.start.halfedge is not None, "Whoopsie #2"
                    
                    if previous.face is None:
                        previous.face = self.F[i]
                    
                #if self.F[i].halfedge is None:
                #    self.F[i].halfedge = previous
                
                nextKey = (currIndex, nextIndex)
                twinKey = (nextIndex, currIndex)
                
                if nextKey not in vToH:
                    next = Halfedge()
                    next.start = self.V[currIndex]
                    next.face = self.F[i]
                    twin = Halfedge()
                    twin.start = self.V[nextIndex]
                    twin.twin = next
                    next.twin = twin
                    
                    if self.V[currIndex].halfedge is None:
                        self.V[currIndex].halfedge = next
                    
                    if self.V[nextIndex].halfedge is None:
                        self.V[nextIndex].halfedge = twin
                    
                    vToH[nextKey] = next
                    vToH[twinKey] = twin
                    
                    self.H.append(next)
                    self.H.append(twin)
                
                else:
                    next = vToH[nextKey]
                    
                    assert next.start.halfedge is not None, "Whoopsie #3"
                    assert next.twin.start.halfedge is not None, "Whoopsie #4"
                    
                    if next.face is None:
                        next.face = self.F[i]
                    
                if self.F[i].halfedge is None:
                    self.F[i].halfedge = next
                
                previous.next = next
                next.previous = previous
        
        self.closeloops()
        
        # Check the geometric and incidency information is consistent
        for v in self.V:
            v.checkconsistency()
        for f in self.F:
            f.checkconsistency()
        for h in self.H:
            h.checkconsistency()
    
    
    def tovf(self) -> VF:
        vf = VF()
        vToI = {}
        for v in self.V:
            vToI[v.coords] = vf.addvertexfrompvector(v.coords)
        
        for f in self.F:
            indices = []
            
            h = f.halfedge
            while True:
                indices.append(vToI[h.start.coords])
                h = h.next
                if h == f.halfedge:
                    break
            
            vf.addface(indices)
        
        return vf
    