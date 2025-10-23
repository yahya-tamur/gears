# There's probably a more canonical way to do this in python.

# this is meant to be something you can import from another file.

# Only functions in this file are in degrees


from gears_internal import rad, bevel_gear_assembly

def bevel_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle = 20, helix_angle=0, tooth_step=50, flat_step=20):
    return bevel_gear_assembly(modul, tooth_number, rad(partial_cone_angle), tooth_width, bore, rad(pressure_angle), rad(helix_angle), tooth_step, flat_step)

