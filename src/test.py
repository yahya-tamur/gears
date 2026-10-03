from .lib import disp, planetary_gear


mesh = planetary_gear(
    modul=2,
    sun_teeth=48,
    planet_teeth=16,
    number_planets=0,
    width=20,
    rim_width=4,
    helix_angle=20,
)
disp(mesh)
