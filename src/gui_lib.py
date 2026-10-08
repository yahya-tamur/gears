from tkinter import Tk, ttk, font, filedialog
import subprocess


# parameters = (name of variable, default, text, type)[]
def run_gui(title, parameters, callback):

    root = Tk()
    root.title(title)

    default_font = font.nametofont("TkDefaultFont")
    text_font = font.nametofont("TkTextFont")
    menu_font = font.nametofont("TkMenuFont")
    for f in (default_font, text_font, menu_font):
        f.configure(size=16)

    frm = ttk.Frame(root, padding=10)
    frm.grid()

    frm.grid_columnconfigure(0, minsize=300)
    frm.grid_columnconfigure(0, minsize=300, weight=0)

    filename = "./a.stl"

    filename_label = ttk.Label(frm, text=filename)
    filename_label.grid(column=1, row=0)

    def pick_file():
        nonlocal filename

        new_filename = filedialog.asksaveasfilename(
            title="Save As",
            initialfile="a.stl",
            initialdir=".",
            filetypes=(("Stl File", "*.stl"),),
        )

        if new_filename:
            filename = new_filename
            if len(new_filename) > 23:
                filename_label.config(text="..." + new_filename[-20:])
            else:
                filename_label.config(text=new_filename)

    ttk.Button(frm, text="Save As", command=pick_file).grid(column=0, row=0)

    line = 1

    def make_entry(text, default):
        nonlocal line
        ttk.Label(frm, text=text).grid(column=0, row=line)
        entry = ttk.Entry(frm)
        entry.grid(column=1, row=line)
        entry.insert(0, default)
        line += 1
        return entry.get

    get_viewer = make_entry("STL Viewer", "fstl")
    get_param = {
        name: make_entry(text, default) for (name, default, text, _) in parameters
    }

    def string_to_bool(s):
        if s.lower() in ("true", "t"):
            return True
        if s.lower() in ("false", "f"):
            return False
        raise ValueError(f"Invalid boolean value: {s}")

    def eval_with_type(inp, typ):
        match typ:
            case "int":
                return int(inp)
            case "float":
                return float(inp)
            case "bool":
                return string_to_bool(inp)
            case "string":
                return inp
            case _:
                raise ValueError(f"Unsupported type {typ}")

    def create(open_viewer=True):
        nonlocal filename, get_param, eval_with_type
        params = {
            name: eval_with_type(get_param[name](), typ)
            for (name, _, _, typ) in parameters
        }
        params["filename"] = filename

        callback(**params)

        viewer = get_viewer()
        if open_viewer:
            subprocess.Popen([viewer, filename])

    def update():
        create(open_viewer=False)

    ttk.Button(frm, text="Create", command=create).grid(column=0, row=line, sticky="ew")
    ttk.Button(frm, text="Update", command=update).grid(column=1, row=line, sticky="ew")

    root.mainloop()
