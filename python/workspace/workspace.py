# -*- coding: utf-8 -*-

from ds.vf import VF
from ds.dcel import DCEL, Halfedge, Face
from toolkit.plane import Plane
from toolkit.pvector import PVector
from toolkit.ray import Ray
from toolkit import utils

# The class wrapping the elements to create a TIC:
# - The geometric domain
# - Generation arguments
# - The Assembly
# - The interfaces between blocks (Future)
# - Analysis attributes (Future)
class Workspace:
    
    
    # Constructor of the class.
    def __init__(self):
        
        # The DCEL representing the geometric domain (i.e., tessellation or 3D
        # surface) after which we will design the TIC.
        self.domain = None
        
        # A dictionary that maps every face in the domain to its interlocking 
        # block
        self.assembly = None
        
        # A dictionary that maps every face in the domain to its height value.
        self.faceheights = None
        
        # A dictionary that maps every half edge of the domain to its 
        # direction. Let h be a half edge. if halfedgedirections[h] = d, then
        # halfedgedirections[h.twin] = -d.
        self.halfedgedirections = None
        
        # A dictionary that maps every half edge of the domain to its midpoint.
        self.halfedgemidpoints = None
        
        # A dictionary that maps every half edge of the domain to its incident
        # plane.
        self.halfedgeplanes = None
        
        # A dictionary that maps every half edge of the domain to its rotation
        # angle. 
        self.halfedgerotationangles = None
        
        # A dictionary that maps every half edge of the domain to its vector.
        self.halfedgevectors = None
    
    
    # Clears the content of the workspace
    def clearall(self):
        if hasattr(self, 'domain') and self.domain is not None:
            del self.domain
            self.domain = None
        
        if hasattr(self, 'assembly') and self.assembly is not None:
            del self.assembly
            self.assembly = None
        
        #if hasattr(self, 'interfaces') and self.interfaces is not None:
        #    del self.interfaces
        #    self.interfaces = None
        
        if hasattr(self, 'faceheights') and self.faceheights is not None:
            del self.faceheights
            self.faceheights = None
        
        if hasattr(self, 'halfedgedirections') and self.halfedgedirections is not None:
            del self.halfedgedirections
            self.halfedgedirections = None
        
        if hasattr(self, 'halfedgemidpoints') and self.halfedgemidpoints is not None:
            del self.halfedgemidpoints
            self.halfedgemidpoints = None
        
        if hasattr(self, 'halfedgeplanes') and self.halfedgeplanes is not None:
            del self.halfedgeplanes
            self.halfedgeplanes = None
        
        if hasattr(self, 'halfedgerotationangles') and self.halfedgerotationangles is not None:
            del self.halfedgerotationangles
            self.halfedgerotationangles = None
        
        if hasattr(self, 'halfedgevectors') and self.halfedgevectors is not None:
            del self.halfedgevectors
            self.halfedgevectors = None
    
    
    # Calculates the interlocking blocks for each face of the geometric domain.
    # Note: halfedgeplanes must exist.
    def calculateblocks(self, e:int|float=1e-8, r:int=10):
        assert isinstance(e, (int, float)), "e not a real number"
        assert isinstance(r, int), "r not an integer number"
        
        self.assembly = {}
        
        # Traverse through the faces of the geometric domain and calculate the
        # respective block
        for f in self.domain.F:
            
            # A list to store the block vertices
            blockvertices = []
            
            allvertices = True
            currenthalfedge = f.halfedge
            while True:
                
                # Get the intersection of the current three planes
                P0 = self.halfedgeplanes[currenthalfedge.previous]
                P1 = self.halfedgeplanes[currenthalfedge]
                P2 = self.halfedgeplanes[currenthalfedge.next]
                exists, point = utils.threeplanesintersection(P0, P1, P2, e, r)
                blockvertices.append(point if exists else None)
                allvertices = allvertices and exists
                
                currenthalfedge = currenthalfedge.next
                if currenthalfedge == f.halfedge:
                    break
            
            # If not all vertices were calculated then return the block vertices 
            # without the face indices
            if not allvertices:
                return blockvertices, None
            
            # A list to store the block faces
            blockfaces = []
            
            nvertices = len(blockvertices)
            direction = self.halfedgedirections[f.halfedge]
            
            for i in range(nvertices):
                
                # Calculate the indices of the previous and next vertices for the 
                # current face
                previndex = nvertices - 1 if i == 0 else i - 1
                nextindex = 0 if i == nvertices - 1 else i + 1
                
                if direction > 0:
                    blockfaces.append([i, previndex, nextindex])
                else:
                    blockfaces.append([i, nextindex, previndex])
                
                direction *= -1
            
            # If the block has more than four vertices then the top and bottom 
            # faces have to be defined. Such faces are polygons with half number of
            # sides with respect to the original face from the geometric domain
            if nvertices > 4:
                
                topface = [j for j in range(0 if direction > 0 else 1, nvertices, 2)]
                #bottomface = [j for j in range(nvertices - 1 if direction > 0 else nvertices - 2, -1, -2) if j > 1]
                
                bottomface = []
                for j in range(nvertices - 1 if direction > 0 else nvertices - 2, -1, -2):
                    bottomface.append(j)
                    if j <= 1:
                        break
                
                blockfaces.append(topface)
                blockfaces.append(bottomface)
            
            self.assembly[f] = [blockvertices, blockfaces]
    
    
    # Calculates the interlocking blocks for the given face of the geometric 
    # domain.
    # Note: halfedgeplanes and assembly must exist.
    def calculatefaceblock(self, f:Face, e:int|float=1e-8, r:int=10) -> tuple[list, list]:
        assert f is not None, "Missing f"
        assert isinstance(f, Face), "f not a Face"
        assert isinstance(e, (int, float)), "e not a real number"
        assert isinstance(r, int), "r not an integer number"
        assert self.halfedgeplanes is not None, "Missing halfedgeplanes"
        assert self.assembly is not None, "Missing assembly"
        
        # A list to store the block vertices
        blockvertices = []
        
        allvertices = True
        currenthalfedge = f.halfedge
        while True:
            
            # Get the intersection of the current three planes
            P0 = self.halfedgeplanes[currenthalfedge.previous]
            P1 = self.halfedgeplanes[currenthalfedge]
            P2 = self.halfedgeplanes[currenthalfedge.next]
            exists, point = utils.threeplanesintersection(P0, P1, P2, e, r)
            blockvertices.append(point if exists else None)
            allvertices = allvertices and exists
            
            currenthalfedge = currenthalfedge.next
            if currenthalfedge == f.halfedge:
                break
        
        # If not all vertices were calculated then return the block vertices 
        # without the face indices
        if not allvertices:
            return blockvertices, None
        
        # A list to store the block faces
        blockfaces = []
        
        nvertices = len(blockvertices)
        direction = self.halfedgedirections[f.halfedge]
        
        for i in range(nvertices):
            
            # Calculate the indices of the previous and next vertices for the 
            # current face
            previndex = nvertices - 1 if i == 0 else i - 1
            nextindex = 0 if i == nvertices - 1 else i + 1
            
            if direction > 0:
                blockfaces.append([i, previndex, nextindex])
            else:
                blockfaces.append([i, nextindex, previndex])
            
            direction *= -1
        
        # If the block has more than four vertices then the top and bottom 
        # faces have to be defined. Such faces are polygons with half number of
        # sides with respect to the original face from the geometric domain
        if nvertices > 4:
            
            topface = [j for j in range(0 if dir > 0 else 1, nvertices, 2)]
            bottomface = [j for j in range(nvertices - 1 if dir > 0 else nvertices - 2, -1, -2) if j > 1]
            blockfaces.append(topface)
            blockfaces.append(bottomface)
        
        self.assembly[f] = [blockvertices, blockfaces]
        
    
    # 
    def rotateedgevectors(self):
        assert self.halfedgedirections is not None, "Missing halfedgedirections"
        assert self.halfedgerotationangles is not None, "Missing halfedgerotationangles"
        assert self.halfedgevectors is not None, "Missing halfedgevectors"
        
        # A set for storing the visited half edges
        visited = set()
        
        for h in self.domain.H:
            if h in visited:
                continue
            
            # Calculate the direction vector of the half edge
            K = h.twin.start.coords - h.start.coords
            
            # Get the rotation angle of the half edge. Since the Axis-Angle
            # rotation method follows the right hand rule, we must toggle the
            # angle sign to adjust it to the intended rotation direction.
            assert self.halfedgedirections[h] != 0, "Not valid half edge direction value"
            direction = -1 if self.halfedgedirections[h] > 0 else 1
            angle = direction * self.halfedgerotationangles[h]
            
            # Rotate the vector associated to the half edge
            self.halfedgevectors[h].rotate(K, angle)
            
            # Add h and its twin to the set of visited halfedges since they
            # share the memory address to the rotated vector
            visited.add(h)
            visited.add(h.twin)
    
    
    # Sets the content of the geometric domain. Calling this function clears 
    # the content of the workspace.
    def setdomain(self, vf:VF):
        assert vf is not None, "Missing vf"
        assert isinstance(vf, VF), "vf not a VF"
        self.clearall()
        self.domain = DCEL(vf)
    
    
    # Sets the direction values for the given half edge and its twin. The given
    # half edge gets direction value d and its twin gets -d. A positive 
    # direction value indicates it goes inwards the incident face. A negative 
    # direction value indicates it goes outwards the incident face.
    # Note: halfedgedirections must exist
    def setedgedirection(self, h:Halfedge, d:int=1):
        assert h is not None, "Missing h"
        assert isinstance(h, Halfedge), "h not a Halfedge"
        assert isinstance(d, int), "d not an integer"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.halfedgedirections is not None, "Missing halfedgedirections"
        self.halfedgedirections[h] = d
        self.halfedgedirections[h.twin] = -d
    
    
    # Sets the direction values on the half edges of the geometric domain. A 
    # positive direction value indicates it goes inwards the incident face. A 
    # negative direction value indicates it goes outwards the incident face.
    def setedgedirections(self, d:int=1):
        assert isinstance(d, int), "d not an integer"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        self.halfedgedirections = {}
        for h in self.domain.H:
            if h in self.halfedgedirections:
                continue
            self.setedgedirectionsfromhalfedge(h, d)
    
    
    # Sets the direction values on the half edges of the geometric domain 
    # starting on the given half edge. A positive direction value indicates it 
    # goes inwards the incident face. A negative direction value indicates it 
    # goes outwards the incident face.
    def setedgedirectionsfromhalfedge(self, halfedge:Halfedge, d:int):
        assert halfedge is not None, "Missing halfedge"
        assert isinstance(halfedge, Halfedge), "halfedge not a Halfedge"
        assert isinstance(d, int), "d not an integer"
        from ds.queue import Queue
        Q = Queue()
        Q.enqueue((halfedge, d))
        while not Q.isempty():
            h, hd = Q.dequeue()
            if h in self.halfedgedirections:
                continue
            self.halfedgedirections[h] = hd
            if h.twin is not None and h.twin not in self.halfedgedirections:
                Q.enqueue((h.twin, -hd))
            if h.previous is not None and h.previous not in self.halfedgedirections:
                Q.enqueue((h.previous, -hd))
            if h.next is not None and h.next not in self.halfedgedirections:
                Q.enqueue((h.next, -hd))
    
    
    # Sets the dictionary with the half edge midpoints. Note that a half edge
    # and its twin share the reference of the same midpoint.
    def setedgemidpoints(self, e:int|float=1e-8):
        assert isinstance(e, (int, float)), "e not a real number"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        self.halfedgemidpoints = {}
        for h in self.domain.H:
            if h in self.halfedgemidpoints:
                continue
            midpoint = h.midpoint().fixzeros(e)
            self.halfedgemidpoints[h] = midpoint
            self.halfedgemidpoints[h.twin] = midpoint
    
    
    # Sets the planes incident to the edges of the geometric domain.
    # Note: halfedgevectors and halfedgemidpoints must already exist.
    def setedgeplanes(self):
        assert self.halfedgemidpoints is not None, "Missing halfedgemidpoints"
        assert self.halfedgevectors is not None, "Missing halfedgevectors"
        self.halfedgeplanes = {}
        for h in self.domain.H:
            if h in self.halfedgeplanes:
                continue
            plane = Plane(self.halfedgevectors[h], self.halfedgemidpoints[h])
            self.halfedgeplanes[h] = plane
            self.halfedgeplanes[h.twin] = plane
    
    
    # Sets the rotation angle for the given half edge and its twin in the 
    # geometric domain.
    # Note: halfedgerotationangles must already exist.
    def setedgerotationangle(self, h:Halfedge, angle:int|float):
        assert h is not None, "Missing halfedge"
        assert isinstance(h, Halfedge), "halfedge not a Halfedge"
        assert isinstance(angle, (int, float)), "angle not a real number"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.halfedgedirections is not None, "Missing edge directions"
        assert self.halfedgerotationangles is not None, "Missing edge angles"
        self.halfedgerotationangles[h] = angle
        self.halfedgerotationangles[h.twin] = angle
    
    
    # Sets the rotation angle for all half edges in the geometric domain.
    def setedgerotationangles(self, angle:int|float):
        assert isinstance(angle, (int, float)), "angle not a real number"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        self.halfedgerotationangles = {}
        for h in self.domain.H:
            if h in self.halfedgerotationangles:
                continue
            self.halfedgerotationangles[h] = angle
            self.halfedgerotationangles[h.twin] = angle
    
    
    # Note: faceheights must already exist.
    def setedgerotationanglesfromfaceheights(self, e:int|float=1e-8):
        assert isinstance(e, (int, float)), "e not a real number"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert hasattr(self, 'faceheights'), "No face heights attribute"
        self.halfedgerotationangles = {}
        
        # Calculate the top and bottom values for each face in the geometric
        # domain
        facepoints = {}
        for f in self.domain.F:
            C = f.centroid().fixzeros(e)
            N = f.normal().fixzeros(e)
            facepoints[f] = {'T':C + N * self.faceheights[f]['top'], 
                             'B':C - N * self.faceheights[f]['bottom']}
        
        # TODO: Check Tessellation::SetEdgeRotatedVectorsUsingTopBottomSectionPoints
        
        
    
    
    # Sets the vector for the given half edge and its twin in the geometric 
    # domain. Note that a half edge and its twin share the reference of the 
    # same vector.
    # Note: halfedgevectors must already exist.
    def setedgevector(self, h:Halfedge, vector:PVector):
        assert h is not None, "Missing halfedge"
        assert isinstance(h, Halfedge), "halfedge not a Halfedge"
        assert vector is not None, "Missing vector"
        assert isinstance(vector, PVector), "vector not a PVector"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.halfedgevectors is not None, "Missing normals"
        self.halfedgevectors[h] = vector
        self.halfedgevectors[h.twin] = vector
    
    
    # Sets the vector for all half edges in the geometric domain using their 
    # respective normals. Note that a half edge and its twin share the 
    # reference of the same vector.
    def setedgevectorsfromnormals(self, normalize:bool=True):
        assert isinstance(normalize, bool), "normalize not a bool"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        self.halfedgevectors = {}
        for h in self.domain.H:
            if h in self.halfedgevectors:
                continue
            normal = h.normal().normalize() if normalize else h.normal()
            self.halfedgevectors[h] = normal
            self.halfedgevectors[h.twin] = normal
    
    
    # Sets the bottom height for the given face in the geometric domain.
    # Note: faceheights and faceheights[f] must already exist.
    def setfacebottomheight(self, f:Face, bottom:int|float):
        assert f is not None, "Missing face"
        assert isinstance(f, Face), "f not a Face"
        assert isinstance(bottom, (int, float)), "bottom not a real number"
        assert bottom > 0, "bottom must be positive"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.faceheights is not None, "Missing face heights"
        self.faceheights[f]['bottom'] = bottom
    
    
    # Sets the top and bottom heights for the given face in the geometric 
    # domain.
    # Note: faceheights must already exist.
    def setfaceheight(self, f:Face, top:int|float, bottom:int|float):
        assert f is not None, "Missing face"
        assert isinstance(f, Face), "f not a Face"
        assert isinstance(top, (int, float)), "top not a real number"
        assert top > 0, "top must be positive"
        assert isinstance(bottom, (int, float)), "bottom not a real number"
        assert bottom > 0, "bottom must be positive"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.faceheights is not None, "Missing face heights"
        self.faceheights[f] = {'top':top, 'bottom':bottom}
    
    
    # Sets the top and bottom heights for all faces in the geometric domain.
    def setfaceheights(self, top:int|float, bottom:int|float):
        assert isinstance(top, (int, float)), "top not a real number"
        assert top > 0, "top must be positive"
        assert isinstance(bottom, (int, float)), "bottom not a real number"
        assert bottom > 0, "bottom must be positive"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        self.faceheights = {}
        for f in self.domain.F:
            self.faceheights[f] = {'top':top, 'bottom':bottom}
    
    
    # Sets the top height for the given face in the geometric domain.
    # Note: faceheights and faceheights[f] must already exist.
    def setfacetopheight(self, f:Face, top:int|float):
        assert f is not None, "Missing face"
        assert isinstance(f, Face), "f not a Face"
        assert isinstance(top, (int, float)), "top not a real number"
        assert top > 0, "top must be positive"
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.faceheights is not None, "Missing face heights"
        self.faceheights[f]['top'] = top
    
    
    # Rotates the edge vectors and calculates their rotation angle using the
    # face heights according to the Height-Bisection method. Rotation angles 
    # are calculated between the normal vector of each halfedge and the 
    # respective bisector vector.
    # NOTE: halfedgedirections and faceheights must already exist.
    def setheightbisectionelements(self):
        assert hasattr(self, 'domain'), "No domain attribute"
        assert self.domain is not None, "Missing domain"
        assert self.faceheights is not None, "Missing face heights"
        self.halfedgerotationangles = {}
        self.halfedgevectors = {}
        facecentroid = {}
        facenormal = {}
        facetop = {}
        facebottom = {}
        for h in self.domain.H:
            if h in self.halfedgevectors or h.face is None:
                continue
            A = h.start.coords
            B = h.twin.start.coords
            if h.face not in facecentroid:
                facecentroid[h.face] = h.face.centroid()
                facenormal[h.face] = h.face.normal().normalize()
                facetop[h.face] = facecentroid[h.face] + facenormal[h.face] * self.faceheights[h.face]['top']
                facebottom[h.face] = facecentroid[h.face] - facenormal[h.face] * self.faceheights[h.face]['bottom']
            assert self.halfedgedirections[h] != 0, "Not valid half edge direction value"
            D = facebottom[h.face] if self.halfedgedirections[h] > 0 else facetop[h.face]
            
            V = None
            
            # If the twin half edge is incident to a face then consider the 
            # respective point to calculate the rotated vector and find the 
            # bisector. Otherwise, the current face is at the borderline of the
            # geometric domain. So, use the points A, B, and D to determine the
            # rotated vector for the current half edge.
            if h.twin.face is not None:
                if h.twin.face not in facecentroid:
                    facecentroid[h.twin.face] = h.twin.face.centroid()
                    facenormal[h.twin.face] = h.twin.face.normal().normalize()
                    facetop[h.twin.face] = facecentroid[h.twin.face] + facenormal[h.twin.face] * self.faceheights[h.twin.face]['top']
                    facebottom[h.twin.face] = facecentroid[h.twin.face] - facenormal[h.twin.face] * self.faceheights[h.twin.face]['bottom']
                E = facetop[h.twin.face] if self.halfedgedirections[h] > 0 else facebottom[h.twin.face]
                N1 = (A - D).cross(B - D).normalize()
                N2 = (B - E).cross(A - E).normalize()
                V = N1.add(N2).scale(0.5)
            else:
                V = (A - D).cross(B - D).normalize()
            
            self.halfedgevectors[h] = V
            self.halfedgevectors[h.twin] = V
            angle = V.angle(h.normal())
            self.halfedgerotationangles[h] = angle
            self.halfedgerotationangles[h.twin] = angle
    
    
    def setvertexrays(self, fixdirection=True):
        assert self.halfedgeplanes is not None, "Missing halfedgeplanes"
        self.vertexrays = {}
        for face in self.domain.F:
            
            self.vertexrays[face] = {}
            facenormal = face.normal()
            
            h = face.halfedge
            while True:
                
                N = self.halfedgeplanes[h.previous].normal.cross(self.halfedgeplanes[h].normal)
                R = Ray(h.start.coords, N.normalize())
                
                if fixdirection and R.direction.dot(facenormal) < 0:
                    R.direction *= -1
                
                self.vertexrays[face][h.start] = R
                
                
                
                h = h.next
                if h == face.halfedge:
                    break
    
    