from .lib import bevel_gear


bevel_gear(
    2,
    21,
    0,
    20,
    bore=1,
    pressure_angle=20,
    helix_angle=10,
    tooth_steps=16,
    flat_steps=5,
    helix_steps=10,
    bore_steps=20,
    da_factor=1,
    filename="a.stl",
)
bevel_gear(
    2,
    21,
    0.01,
    tooth_width=20,
    bore=1,
    pressure_angle=20,
    helix_angle=10,
    tooth_steps=16,
    flat_steps=5,
    helix_steps=10,
    bore_steps=20,
    da_factor=1,
    filename="b.stl",
)
