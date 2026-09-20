from tkinter import *
from tkinter import ttk, font

from tkinter import filedialog

root = Tk()

# 1. Fetch all default Tkinter fonts and update their sizes
default_font = font.nametofont("TkDefaultFont")  # Used for Labels, Buttons, etc.
text_font = font.nametofont("TkTextFont")        # Used for Entry, Text widgets
menu_font = font.nametofont("TkMenuFont")        # Used for Menus

# 2. Configure a new larger size (e.g., 16)
for f in (default_font, text_font, menu_font):
    f.configure(size=16)


frm = ttk.Frame(root, padding=10)
frm.grid()

ttk.Label(frm, text="Bevel Herringbone Gear Pair Generator").grid(column=0, row=0, columnspan=2, sticky="ew")
frm.grid_columnconfigure(0, minsize=300)

filename = "./a.stl"

filename_label = ttk.Label(frm, text=filename)
filename_label.grid(column=1, row=1)

def pick_file():
    # Hide the main Tkinter root window if you only want the dialog
    # root.withdraw() 
    global filename
    
    new_filename = filedialog.asksaveasfilename(
        title="Save As",
        initialfile="a.stl",
        initialdir=".",
        filetypes=(
            ("Stl File", "*.stl"),
        )
    )
    
    if new_filename:
        filename_label.config(text=new_filename)
        filename = new_filename
    else:
        print("No file was selected.")

ttk.Button(frm, text="Save As", command=pick_file).grid(column=0,row=1)

line = 2

def make_entry(text, default):
    global line
    ttk.Label(frm, text=text).grid(column=0, row=line)
    entry = ttk.Entry(frm)
    entry.grid(column=1, row=line)
    entry.insert(0, default)
    line += 1
    return entry.get


get_viewer = make_entry("stl viewer", "fstl")
get_modul = make_entry("modul", 2)
get_gear_teeth = make_entry("Gear Teeth", 40)
get_pinion_teeth = make_entry("Pinion Teeth", 22)
get_tooth_width = make_entry("Tooth Width", 25)
get_axis_angle = make_entry("Axis Angle", 70)
get_gear_bore = make_entry("Gear Bore", 0)
get_pinion_bore = make_entry("Pinion Bore", 0)
get_pressure_angle = make_entry("Pressure Angle", 20)
get_helix_angle = make_entry("Helix Angle", 40)
get_together_built = make_entry("Build Together", "True")
get_tooth_step = make_entry("Tooth Step", 16)
get_flat_step = make_entry("Flat Step", 10)

import subprocess
from gears import *
from make_stl import make_stl
def create(open_viewer=True):
    global filename

    mesh = bevel_herringbone_gear_pair( \
        modul = float(get_modul()), \
        gear_teeth = int(get_gear_teeth()), \
        pinion_teeth = int(get_pinion_teeth()), \
        tooth_width = float(get_tooth_width()), \
        axis_angle = float(get_axis_angle()), \
        gear_bore = float(get_gear_bore()), \
        pinion_bore = float(get_pinion_bore()), \
        pressure_angle = float(get_pressure_angle()), \
        helix_angle = float(get_helix_angle()), \
        together_built = get_together_built().lower() in ("true", "t"), \
        tooth_step = int(get_tooth_step()), \
        flat_step = int(get_flat_step()), \
    )
    viewer = get_viewer()
    make_stl(mesh, filename)
    if open_viewer:
        subprocess.Popen([viewer, filename])

def update():
    create(open_viewer=False)

ttk.Button(frm, text="Create", command=create).grid(column=0, row=line, sticky="ew")
ttk.Button(frm, text="Update", command=update).grid(column=1, row=line, sticky="ew")

#ttk.Button(frm, text="Quit", command=root.destroy).grid(column=1, row=0)
root.mainloop()
input()