from .gui_lib import run_gui
from .lib import make_gear_pair, make_planetary_gear
from sys import argv


def run_gear_pair_gui():
    parameters = [
        ("modul", 2.0, "Tooth Size", "float"),
        ("gear_teeth", 40, "Number of Gear Teeth", "int"),
        ("pinion_teeth", 22, "Number of Pinion Teeth", "int"),
        ("axis_angle", 30.0, "Axis Angle", "float"),
        ("tooth_width", 25.0, "Width at Teeth", "float"),
        ("bore", 1.0, "Bore Diameter", "float"),
        ("pressure_angle", 20.0, "Pressure Angle", "float"),
        ("helix_angle", 40.0, "Helix Angle", "float"),
        ("together_built", "True", "Assemble Model", "bool"),
        ("tooth_steps", 16, "Teeth Steps", "int"),
        ("helix_steps", 5, "Helix Steps", "int"),
        ("flat_steps", 3, "Flat Steps", "int"),
        ("bore_steps", 10, "Bore Steps", "int"),
        ("da_factor", 1, "Addendum Factor", "float"),
    ]

    run_gui(
        "Gear Pair Generator",
        parameters,
        make_gear_pair,
    )


def run_planetary_gear_gui():
    planetary_gear_parameters = [
        ("modul", 2.0, "Tooth Size", "float"),
        ("sun_teeth", 16, "Number of Sun Teeth", "int"),
        ("planet_teeth", 16, "Number of Planet Teeth", "int"),
        ("number_planets", 4, "Number of Planets", "int"),
        ("width", 30.0, "Width", "float"),
        ("rim_width", 5.0, "Ring Width", "float"),
        ("bore", 2.0, "Bore Diameter", "float"),
        ("pressure_angle", 20, "Pressure Angle", "float"),
        ("helix_angle", 30, "Helix Angle", "float"),
        ("ring_shortening_factor", 0.3, "Ring Shortening Factor", "float"),
        ("together_built", "True", "Assemble Model", "bool"),
        ("tooth_steps", 16, "Teeth Steps", "int"),
        ("helix_steps", 5, "Helix Steps", "int"),
        ("flat_steps", 3, "Flat Steps", "int"),
        ("bore_steps", 20, "Ring Steps", "int"),
        ("da_factor", 1, "Addendum Factor", "float"),
    ]

    run_gui("Planetary Gear Generator", planetary_gear_parameters, make_planetary_gear)


if __name__ == "__main__":
    match argv[1]:
        case "planetary":
            run_planetary_gear_gui()
        case "gear_pair":
            run_gear_pair_gui()
        case _:
            print()
            print("Run with `python -m src.gui <mode>`")
            print("Supported modes: 'planetary', 'gear_pair'")
            print()
