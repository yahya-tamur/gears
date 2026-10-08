from os import system

from .lib import make_gear_pair, make_gear, make_ring_gear, make_planetary_gear
make_planetary_gear(
    modul=2,
    sun_teeth=64,
    planet_teeth=16,
    width=30,
    number_planets=8,
    rim_width=5,
    bore=10,
    helix_angle=10,
    tooth_steps=8,
    flat_steps=2,
    helix_steps=20,
    bore_steps=50,
    ring_shortening_factor=0.6,
    together_built=True,
)
system("fstl a.stl")
