from gears import bevel_gear
from make_stl import make_stl

gear = bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, bore=1, tooth_width=10, pressure_angle = 20, helix_angle=20);

print(len(gear))

make_stl(gear, "tooth.stl")
