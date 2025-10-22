from gears import bevel_gear, rad
from make_stl import make_stl

gear = bevel_gear(modul=1, tooth_number=20, partial_cone_angle=rad(20), tooth_width=10, pressure_angle = rad(20), helix_angle=rad(20));

print(len(gear))

make_stl(gear, "tooth.stl")
