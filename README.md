This is an effort to port the openscad code to a python script creating the stl
directly. This creates a nicer stl file.

You can look at `main.py` for example usage and the original README for
explanations of all the variables. The `disp` function in `main.py` can be used
in the python repl to iterate through different gear designs quickly: (You need
some lightweight stl viewer, fstl by default)

```
from main import *
from gears import *
disp(bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, tooth_width=2, helix_angle=20))
```

to do:

* n/a

nice to have:

* go over numerical errors

* port other gear designs?
