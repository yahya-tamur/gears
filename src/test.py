from .lib import bevel_gear

from os import system

bevel_gear(
    modul=2,
    tooth_number=20,
    partial_cone_angle=40,
    tooth_width=5,
    bore=3,
    # pressure_angle=20,
    # helix_angle=10,
    # tooth_steps=16,
    # flat_steps=5,
    # helix_steps=10,
    # bore_steps=20,
    # filename="a.stl"
)

system("fstl a.stl")
