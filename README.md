This is an effort to port the openscad code to a python script.

The main reason I started to do this is that the openscad stl creation
process is complicated and error-prone. Since it supports differences,
unions, etc. of 3d volumes, it has a hard time converting these into
stl files, which are just sets of triangles in 3d space.

The original 'gears' script works by creating rectangular prisms,
then stretching them and rotating them to make the gear shapes.

This program by comparison, calculates the correct positions of vertices,
and saves them directly into stl files.

Compare the number of errors below, or check out the differences in the
'comparison' folder.

![screenshot](./screenshot.png)

This should make it easier to work with these in modeling software.
Also, this program is significantly faster. Generating the comparison
images took `0.228346 seconds` with this project and `235.214264 seconds`
with the original.

The easiest way to use this program is by clicking on the scripts under
`bin`. These just call `src/gui.py` with the corresponding argument. You
can also import the code to create the meshes from `lib.py` and write them
into stl files using `make_stl` from `make_stl.py`.

to do:


* print example planetary gear

* windows bins may not work idk

* update readme

* make a post

nice to have:
* more polygons for linear_extrusion!
    - bevel gears may be more difficult.
* I adapted bevel gears which use spherical coordinates to make flat gears (called spur
    gears in original), but the library contains math to do this using polar instead of spherical coordinates. So, I should probably just port that code instead.
* tkinter gui for some other functions?
* port other gear designs?