import subprocess
from .make_mesh import mesh_grid, mesh_stitch, mesh_gear, mesh_convex_polygon, mesh_ring
from .save_stl import save_stl
from .gears import bevel_herringbone_gear_data
from .math import rad

from math import cos, sin, pi


def regular_polygon(n, r):
    return [(r*cos(2*pi*i/n), r*sin(2*pi*i/n), 0) for i in range(n)]


#tooth_west, tooth_east = bevel_herringbone_gear_data(
#    modul=2,
#    tooth_number=21,
#    partial_cone_angle=rad(50),
#    tooth_width=10,
#    pressure_angle=rad(20),
#    helix_angle=rad(40),
#    tooth_step=16,
#    helix_step=10,
#)
mesh = []

inner = regular_polygon(500, 1)
outer = regular_polygon(3, 2.4)
mesh_ring(mesh, outer, inner)


save_stl(
    mesh,
    "./a.stl",
)
import os
os.system("fstl a.stl")