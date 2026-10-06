from .mesh import Mesh
from .gears import bevel_herringbone_gear_data
from .math import rad


from os import system


tooth_west, tooth_east = bevel_herringbone_gear_data(
    modul=1,
    tooth_number=13,
    partial_cone_angle=rad(20),
    tooth_width=4,
    pressure_angle=rad(20),
    helix_angle=-rad(50),
    tooth_step=16,
    helix_step=10,
)
mesh = Mesh()
mesh.add_gear(tooth_west, tooth_east, 13, 10, 0, 20)
mesh.save_stl("./a.stl")

system("fstl a.stl")
