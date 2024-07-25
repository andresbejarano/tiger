# -*- coding: utf-8 -*-

from ds.dcel import DCEL
import geometries.planar as planar
import geometries.platonicsolids as platonic
import geometries.closed as closed
#G = platonic.tetrahedron()
#G = platonic.hexahedron()
#G = platonic.octahedron()
#G = platonic.dodecahedron()
#G = platonic.icosahedron()

#G = planar.barrelvault(2, 2, 25, 15)
#G = planar.equilateraltriangle()
#G = planar.square()
#G = planar.squaredgrid(3, 3, 6, 6)
#G = planar.saddle(2, 2, 10, 10)
#G = planar.monkeysaddle(2, 2, 20, 20)
#G = planar.wave(5, 5, 20, 20)
#G = planar.regularpolygon(6, 3)
#G = planar.ellipticparaboloid(3, 3, 3, 3, 10, 10)

#G = closed.torus(1, 0.5, 30, 20)

#D = DCEL(G)
#G = D.dual()
#D = DCEL(G)



S = planar.regularpolygon(5)
print(S)

S.flip()
print(S)

S.flip()
print(S)


#import apps.pyqt_vtk.helpers.load_vtkactor_with_dcel as tiger
#tiger.testactor(D)