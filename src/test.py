from sys import argv

from .make_comparison import run_openscad
from .lib import make_gear_pair

argv = set(argv[1:])
print(argv)

make_gear_pair(
    modul=2,
    gear_teeth=64,
    pinion_teeth=16,
    axis_angle=30,
    tooth_width=8,
    bore=0.5,
    helix_angle=-20,
    pressure_angle=20,
    bore_steps=30,
    together_built=True,
    filename="./a1.stl",
)

run_openscad(
    "bevel_herringbone_gear_pair("
    + "modul=2,"
    + "gear_teeth=64,"
    + "pinion_teeth=16,"
    + "tooth_width=8,"
    + "axis_angle=30,"
    + "gear_bore=0.5,"
    + "pinion_bore=0.5,"
    + "helix_angle=-20,"
    + "pressure_angle=20,"
    + "together_built=true"
    + ")",
    "./b.stl",
)
