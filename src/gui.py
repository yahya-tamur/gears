from .gui_lib import run_gui
from .lib import bevel_herringbone_gear_pair, planetary_gear
from sys import argv


def run_bevel_gear_pair_gui():
    bevel_gear_pair_parameters = [
        ("modul", 2.0, "Tooth Size", "float"),
        ("gear_teeth", 40, "Number of Gear Teeth", "int"),
        ("pinion_teeth", 22, "Number of Pinion Teeth", "int"),
        ("tooth_width", 25.0, "Width of Gear at Teeth", "float"),
        ("axis_angle", 70.0, "Axis Angle", "float"),
        ("gear_bore", 0.0, "Gear Bore Radius", "float"),
        ("pinion_bore", 0.0, "Pinion Bore Radius", "float"),
        ("pressure_angle", 20.0, "Pressure Angle", "float"),
        ("helix_angle", 40.0, "Herringbone Angle", "float"),
        ("together_built", "True", "Assemble Model", "bool"),
        ("tooth_step", 16, "Teeth Resolution", "int"),
    ]

    run_gui(
        "Bevel Gear Pair Generator",
        bevel_gear_pair_parameters,
        bevel_herringbone_gear_pair,
    )


def run_planetary_gear_gui():
    planetary_gear_parameters = [
        ("modul", 2.0, "Tooth Size", "float"),
        ("sun_teeth", 16, "Number of Sun Teeth", "int"),
        ("planet_teeth", 16, "Number of Planet Teeth", "int"),
        ("width", 30.0, "Width", "float"),
        ("number_planets", 4, "Number of Planets", "int"),
        ("rim_width", 5.0, "Ring Gear Width", "float"),
        ("sun_bore", 2.0, "Sun Gear Bore", "float"),
        ("planet_bore", 1.0, "Planet Gear Bore", "float"),
        ("pressure_angle", 20, "Pressure Angle", "float"),
        ("helix_angle", 30, "Herringbone Angle", "float"),
        ("together_built", "True", "Assemble Model", "bool"),
        ("tooth_step", 16, "Teeth Resolution", "int"),
        ("ring_shortening_factor", 1, "Ring Teeth Shortening", "float"),
    ]

    run_gui("Planetary Gear Generator", planetary_gear_parameters, planetary_gear)


if __name__ == "__main__":
    match argv[1]:
        case "planetary":
            run_planetary_gear_gui()
        case "bevel_pair":
            run_bevel_gear_pair_gui()
        case _:
            print()
            print("Run with `python -m src.gui <mode>`")
            print("Supported modes: 'planetary', 'bevel_pair'")
            print()
