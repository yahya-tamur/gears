from .make_mesh import mesh_gear
from .save_stl import save_stl
from .gears import bevel_herringbone_gear_data
from .math import rad

from math import cos, sin, pi, atan

from os import system

def regular_polygon(n, r):
    return [(r * cos(2 * pi * i / n), r * sin(2 * pi * i / n), 0) for i in range(n)]

axis_angle = rad(80)
pinion_teeth=21
gear_teeth=32

delta_gear = atan(sin(axis_angle) / (pinion_teeth / gear_teeth + cos(axis_angle)))
delta_pinion = atan(sin(axis_angle) / (gear_teeth / pinion_teeth + cos(axis_angle)))

tooth_west, tooth_east = bevel_herringbone_gear_data(
    modul=2,
    tooth_number=pinion_teeth,
    partial_cone_angle=delta_pinion,
    tooth_width=15,
    pressure_angle=rad(20),
    helix_angle=rad(50),
    tooth_step=16,
    helix_step=10,
)
mesh = []

# inner = regular_polygon(500, 1)
# outer = regular_polygon(3, 2.4)
# mesh_ring(mesh, outer, inner)
mesh_gear(mesh, tooth_west, tooth_east, 21, 10,  bore=0.2, bore_steps=20)

save_stl(
    mesh,
    "./a.stl",
)

system("fstl a.stl")


tooth_west, tooth_east = bevel_herringbone_gear_data(
    modul=2,
    tooth_number=gear_teeth,
    partial_cone_angle=delta_gear,
    tooth_width=15,
    pressure_angle=rad(20),
    helix_angle=-rad(50),
    tooth_step=16,
    helix_step=10,
)
mesh = []

# inner = regular_polygon(500, 1)
# outer = regular_polygon(3, 2.4)
# mesh_ring(mesh, outer, inner)
mesh_gear(mesh, tooth_west, tooth_east, 21, 10,  bore=0.2, bore_steps=20)

save_stl(
    mesh,
    "./b.stl",
)

system("fstl b.stl")
