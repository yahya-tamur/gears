from gears import bevel_gear
from make_stl import make_stl

from gears_internal import bevel_herringbone_gear_data, rad


if __name__ == "__main__":
    #gear = bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, bore=1, tooth_width=10, pressure_angle = 20, helix_angle=20);
    gear = bevel_herringbone_gear_data(modul=1, tooth_number=20, partial_cone_angle=rad(20), bore=1, tooth_width=10, pressure_angle = rad(20), helix_angle=rad(20), tooth_step=16);

    print(len(gear))

    make_stl(gear, "tooth.stl")

import os

def disp(mesh, viewer="fstl", filename="a.stl"):
    make_stl(mesh, filename)
    os.system(f"{viewer} {filename}")
