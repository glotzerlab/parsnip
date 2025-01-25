# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

import freud
import gsd.hoomd
import numpy as np
import rowan
from coxeter.families import ArchimedeanFamily, JohnsonFamily, PlatonicFamily

from parsnip import CifFile

gbf = JohnsonFamily.get_shape("Gyrobifastigium")
gbf.volume = 4
oct = PlatonicFamily.get_shape("Octahedron")
oct.volume = 4

ttet = ArchimedeanFamily.get_shape("Truncated Tetrahedron")
ttet.volume = 4

# cif = CifFile("../aflow_cif_db/AFLOW/A_tI4_141_a.cif")
# cif = CifFile("../aflow_cif_db/AFLOW/A_cI2_229_a.cif")
cif = CifFile("../aflow_cif_db/AFLOW/A_hR2_166_c.alpha-As.cif")
N = 1
box, pos = freud.data.UnitCell(cif.box, cif.build_unit_cell()).generate_system(N)
print(pos)


# q1 = rowan.from_axis_angle([0,1,0],np.pi)
# q1 = rowan.from_matrix(rowan.to_matrix(rowan.from_axis_angle([0,1,0],np.pi/2)) @ rowan.to_matrix(rowan.from_axis_angle([0,1,0], np.pi/2)))
# q2 = rowan.from_matrix(rowan.to_matrix(rowan.from_axis_angle([0,1,0],np.pi/2)) @ rowan.to_matrix(rowan.from_axis_angle([0,1,0], np.pi)))
# q2 = (rowan.from_axis_angle([0,0,1],np.pi/2))
#

# q1 = rowan.from_matrix(rowan.to_matrix(rowan.from_axis_angle([0,1,0.5],np.pi/2)) @ rowan.to_matrix(rowan.from_axis_angle([0,0,1], np.pi)))
# q2 = rowan.from_matrix(rowan.to_matrix(rowan.from_axis_angle([0,1,0.5],np.pi/2)) @ rowan.to_matrix(rowan.from_axis_angle([0,1,0], np.pi)))


q1 = rowan.from_matrix(
    # rowan.to_matrix(rowan.from_axis_angle([0,1,0],np.pi/2)) @
    # rowan.to_matrix(rowan.from_axis_angle([0,0,1],np.pi/2)) @
    rowan.to_matrix(rowan.from_axis_angle([0, 0, 1], np.pi))
)
q2 = rowan.from_matrix(
    rowan.to_matrix(rowan.from_axis_angle([0, 1, 0], np.pi / 2))
    @
    # rowan.to_matrix(rowan.from_axis_angle([0,0,1],np.pi/2)) @
    rowan.to_matrix(rowan.from_axis_angle([0, 1, 0], np.pi))
)
# q2 = rowan.from_matrix(rowan.to_matrix(rowan.from_axis_angle([0,1,0],np.pi/2)) @ rowan.to_matrix(rowan.from_axis_angle([0,1,0], np.pi)))

q1 = [1, 0, 0, 0]
q2 = [np.cos(np.pi / 4), np.sin(np.pi / 4), 0, 0]
# q1 = rowan.from_matrix(
#     rowan.to_matrix(rowan.from_axis_angle([0,0,1], np.pi/3))
# )
# q2 = rowan.from_matrix(
#     rowan.to_matrix(rowan.from_axis_angle([0,1,0],np.pi)) @
#     rowan.to_matrix(rowan.from_axis_angle([0,0,1], np.pi/3))
# )

print("")
frame = gsd.hoomd.Frame()
frame.particles.N = len(pos)
frame.particles.position = pos
# frame.particles.orientation = [q1, q1, q2, q2]
# frame.particles.orientation = [q1, q2] * (len(pos) // 2)
frame.particles.orientation = [q2, q1] * (len(pos) // 2)
# frame.particles.orientation = [q2, q2, q2, q1, q1, q1]
# frame.particles.orientation = [[1,0,0,0]] * len(pos)
frame.configuration.box = [*box.L, box.xy, box.xz, box.yz]
frame.particles.type_shapes = [ttet.gsd_shape_spec]

with gsd.hoomd.open("lattice.gsd", "w") as f:
    f.append(frame)
