# -*- coding: utf-8 -*-

from ds.vf import VF
from workspace.workspace import Workspace
from toolkit.plane import Plane
from toolkit.ray import Ray
import toolkit.utils as utils

# 
# Notes:
# - Positive angles will rotate towards inside the face, while negative angles 
#   will rotate outwards the face.
# - This version of the evolve function adds the evolved points and the edges 
#   between them and the vertices of the seed polygon to the allpoints and
#   alledges objects.
def evolve(seed:VF, angles:list, height:int|float, allpoints:dict, alledges:set, allfaces:set, e:int|float=1e-8, r:int=10):
    assert seed is not None, "Missing seed"
    assert isinstance(seed, VF), "seed is not a VF"
    assert angles is not None, "Missing angles"
    assert isinstance(angles, list), "angles not a list"
    assert isinstance(height, (int, float)), "height not a real number"
    assert allpoints is not None, "Missing allpoints"
    assert isinstance(allpoints, dict), "allpoints not a dictionary"
    assert alledges is not None, "Missing alledges"
    assert isinstance(alledges, set), "alledges not a set"
    assert isinstance(e, (int, float)), "e not a real number"
    
    W = Workspace()
    W.setdomain(seed)
    
    W.halfedgedirections = {}
    W.halfedgerotationangles = {}
    i = 0
    h = W.domain.F[0].halfedge
    
    while True:
        
        # Check if the current vertex of the seed
        v = h.start.coords.clone().fixzeros(e)
        if v not in allpoints:
            allpoints[v] = len(allpoints)
        
        direction = 1 if angles[i] >= 0 else -1
        W.setedgedirection(h, direction)
        W.setedgerotationangle(h, abs(angles[i]))
        i += 1
        
        # Move to the next half edge in the face. Break the while loop if we 
        # return to the first halfe dge
        h = h.next
        if h == W.domain.F[0].halfedge:
            break
    
    # Set the default elements required to calculate the rotated edge planes
    W.setedgemidpoints()
    W.setedgevectorsfromnormals()
    W.rotateedgevectors()
    W.setedgeplanes()
    
    # Set the face and evolution planes
    normal = W.domain.F[0].normal().normalize()
    centroid = W.domain.F[0].centroid()
    faceplane = Plane(normal, centroid)
    evolutionplane = Plane(normal, centroid + normal * height)
    
    # Calculate the rays representing the intersection between two adjacent
    # planes
    vertexrays = {}
    h = W.domain.F[0].halfedge
    while True:
        N = W.halfedgeplanes[h.previous].normal.cross(W.halfedgeplanes[h].normal)
        vertexrays[h.start] = Ray(h.start.coords, N.normalize())
        
        # It might be possible that the direction vector of the ray points at
        # a different direction than the normal of the evolution plane. In 
        # that case, flip the direction of the ray.
        if vertexrays[h.start].direction.dot(evolutionplane.normal) < 0:
           vertexrays[h.start].direction *= -1 
        
        # Move to the next half edge in the face. Break the while loop if we 
        # return to the first halfe dge
        h = h.next
        if h == W.domain.F[0].halfedge:
            break
    
    
    
    # Initialize an empty dictionary to store the evolution points. We use a 
    # dictionary since it preserves the insertion order and avoid storing 
    # replicated points. Points are stored as keys
    points = {}
    h = W.domain.F[0].halfedge
    while True:
        
        # Get the intersection (if any) between the rays incident to the end
        # points of the current half edge
        R1 = vertexrays[h.start]
        R2 = vertexrays[h.twin.start]
        
        v0idx = allpoints[h.start.coords.clone().fixzeros(e)]
        v1idx = allpoints[h.twin.start.coords.clone().fixzeros(e)]
        
        result, Q, s, t = utils.raysintersection(R1, R2, e)
        Q.fixzeros(e).roundvalues(r)
        
        # If there is an intersection of the ray then it fall into one of the 
        # following cases:
        # 1) Above the evolution plane: In this case, we need the intersection
        # of the rays with the evolution plane.
        # 2) At the evolution plane: In this case, we only need the 
        # intersection point.
        # 3) Below the evolution plane and above the face plane: In this case,
        # we only need the intersection plane.
        # 4) Below the face plane: In this case, we calculate the intersection
        # of the ray with the evolution plane.
        if result:
            
            # Check case 1:
            if evolutionplane.pointlocation(Q, e) > 0:
                result1, t1 = utils.planerayintersection(evolutionplane, R1, e)
                result2, t2 = utils.planerayintersection(evolutionplane, R2, e)
                assert result1, "Whoops"
                assert result2, "Whoops"
                
                assert t1 > 0, "Whoops"
                assert t2 > 0, "Whoops"
                if t1 > 0 and t2 > 0:
                    Q1 = R1.at(t1).fixzeros(e).roundvalues(r)
                    Q2 = R2.at(t2).fixzeros(e).roundvalues(r)
                    
                    idx1 = -1
                    if Q1 not in allpoints:
                        idx1 = len(allpoints)
                        allpoints[Q1] = idx1
                        points[Q1] = None
                    else:
                        idx1 = allpoints[Q1]
                    
                    idx2 = -1
                    if Q2 not in allpoints:
                        idx2 = len(allpoints)
                        allpoints[Q2] = idx2
                        points[Q2] = None
                    else:
                        idx2 = allpoints[Q2]
                    
                    # Define the face for this case
                    allfaces.add((v0idx, v1idx, idx2, idx1))
                    
                    # Define the edges for this case
                    alledges.add((min(v0idx, v1idx), max(v0idx, v1idx)))
                    alledges.add((min(v1idx, idx2), max(v1idx, idx2)))
                    alledges.add((min(idx2, idx1), max(idx2, idx1)))
                    alledges.add((min(idx1, v0idx), max(idx1, v0idx)))
            
            # Check case 2
            elif evolutionplane.pointlocation(Q, e) == 0:
                
                idx = -1
                if Q not in allpoints:
                    idx = len(allpoints)
                    allpoints[Q] = idx
                    points[Q] = None
                else:
                    idx = allpoints[Q]
                
                # Define the face for this case
                allfaces.add((v0idx, v1idx, idx))
                
                # Define the edges for this case
                alledges.add((min(v0idx, v1idx), max(v0idx, v1idx)))
                alledges.add((min(v1idx, idx), max(v1idx, idx)))
                alledges.add((min(idx, v0idx), max(idx, v0idx)))
                    
            
            # Check case 3
            elif evolutionplane.pointlocation(Q, e) < 0 and faceplane.pointlocation(Q, e) > 0:
                
                idx = -1
                if Q not in allpoints:
                    idx = len(allpoints)
                    allpoints[Q] = idx
                    points[Q] = None
                else:
                    idx = allpoints[Q]
                
                # Define the face for this case
                allfaces.add((v0idx, v1idx, idx))
                
                # Define the edges for this case
                alledges.add((min(v0idx, v1idx), max(v0idx, v1idx)))
                alledges.add((min(v1idx, idx), max(v1idx, idx)))
                alledges.add((min(idx, v0idx), max(idx, v0idx)))
            
            # Weird case that should not happen
            elif faceplane.pointlocation(Q, e) == 0:
                assert False, "Why though?"
            
            # Check case 4:
            elif faceplane.pointlocation(Q, e) < 0:
                result1, t1 = utils.planerayintersection(evolutionplane, R1, e)
                result2, t2 = utils.planerayintersection(evolutionplane, R2, e)
                assert result1, "Whoops"
                assert result2, "Whoops"
                
                assert t1 > 0, "Whoops"
                assert t2 > 0, "Whoops"
                if t1 > 0 and t2 > 0:
                    Q1 = R1.at(t1).fixzeros(e).roundvalues(r)
                    Q2 = R2.at(t2).fixzeros(e).roundvalues(r)
                    
                    idx1 = -1
                    if Q1 not in allpoints:
                        idx1 = len(allpoints)
                        allpoints[Q1] = idx1
                        points[Q1] = None
                    else:
                        idx1 = allpoints[Q1]
                    
                    idx2 = -1
                    if Q2 not in allpoints:
                        idx2 = len(allpoints)
                        allpoints[Q2] = idx2
                        points[Q2] = None
                    else:
                        idx2 = allpoints[Q2]
                    
                    # Define the face for this case
                    allfaces.add((v0idx, v1idx, idx2, idx1))
                    
                    # Define the edges for this case
                    alledges.add((min(v0idx, v1idx), max(v0idx, v1idx)))
                    alledges.add((min(v1idx, idx2), max(v1idx, idx2)))
                    alledges.add((min(idx2, idx1), max(idx2, idx1)))
                    alledges.add((min(idx1, v0idx), max(idx1, v0idx)))
        
        # Move to the next half edge in the face. Break the while loop if we 
        # return to the first halfe dge
        h = h.next
        if h == W.domain.F[0].halfedge:
            break
    
    return points