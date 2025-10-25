from gears import bevel_gear
from make_stl import make_stl

from gears_internal import bevel_herringbone_gear_pair_assembly, rad


if __name__ == "__main__":
    #gear = bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, bore=1, tooth_width=10, pressure_angle = 20, helix_angle=20);
    gear = bevel_herringbone_gear_pair_assembly(modul=1, gear_teeth=20, pinion_teeth=35, axis_angle=rad(30), gear_bore=1, pinion_bore=1, tooth_width=10, pressure_angle = rad(20), helix_angle=rad(20), tooth_step=16, flat_step=10, together_built=True);

    print(len(gear))

    make_stl(gear, "tooth.stl")

import os

def disp(mesh, viewer="fstl", filename="a.stl"):
    make_stl(mesh, filename)
    os.system(f"{viewer} {filename}")
