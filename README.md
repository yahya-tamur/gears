This is an effort to port the openscad code to a python script.

The main reason I started to do this is that the openscad stl creation
process is complicated and error-prone. Since it supports differences,
unions, etc. of 3d volumes ("CSG"), it has a hard time converting these
into stl files, which are just sets of triangles in 3d space.

This program by comparison, calculates the correct positions of vertices,
and saves them directly into stl files.

Compare the number of errors below, or check out the generated stl's in the
'comparison' folder.

before:

![before](./openscad_bevel_herringbone_gear_analysis.png)

after:

![before](./bevel_herringbone_gear_analysis.png)

This should make it easier to work with these in modeling software.
Also, this program is significantly faster. Generating the comparison
images took `0.228346 seconds` with this project and `235.214264 seconds`
with the original.

The easiest way to use this program is by clicking on the scripts under
`bin`. These just call `src/gui.py` with the corresponding argument. You
can also import the code to create the meshes from `lib.py` and write them
into stl files using `make_stl` from `make_stl.py`.

![gui example](./gui_example.png)

to do:
* more polygons for linear_extrusion!
    - better final step of mesh_ring

    - adapt gears:
        - make gear_pair

        - make ring_gear
        - make spur_gear
        - make planetary_gear

        - delete everything else
        - make lib call the correct stuff

    - make sure openscad comparison uses updated version
    - check gear pairs and examples
  
    - finish readme, screenshots, make post
    
    - is filename thing fixed?

nice to have:
* tkinter gui for some other functions?
    * auto generate these??
* port other gear designs?
