import struct

# input = array of triangles
# each triangle is three points or a normal vector, then 3 points.

# !! If no normal vector is provided, this program will set it to zero.
# And I don't provide normal vectors for the gears.

# This doesn't cause an issue in any program I've tried.
# You can change this behavior to calculate a normal if one isn't provided,
# However, it's hard to make sure the calculated normal vector points out, not in.

# Since I couldn't test to make sure every normal points out,
# I decided to not calculate.

# Something to look at if you have an issue.

def make_stl(triangles, filename):
    with open(filename, 'wb') as f:
        f.write(bytearray(80))
        f.write(struct.pack('<i', len(triangles)))

        for tri in triangles:
            if len(tri) == 3:
                for _ in range(3):
                    f.write(struct.pack('<f', 0))



            for p in tri:
                for pi in p:
                    f.write(struct.pack('<f', pi))
            f.write(bytearray(2))

cube = [ \
[[-1,0,0], [0,0,0], [0,1,0], [0,0,1]], \
[[-1,0,0], [0,1,1], [0,1,0], [0,0,1]], \
[[1,0,0], [1,0,0], [1,1,0], [1,0,1]], \
[[1,0,0], [1,1,1], [1,1,0], [1,0,1]], \

[[0,-1,0], [0,0,0], [1,0,0], [0,0,1]], \
[[0,-1,0], [1,0,1], [1,0,0], [0,0,1]], \
[[0,1,0], [0,1,0], [1,1,0], [0,1,1]], \
[[0,1,0], [1,1,1], [1,1,0], [0,1,1]], \

[[0,0,-1], [0,0,0], [1,0,0], [0,1,0]], \
[[0,0,-1], [1,1,0], [1,0,0], [0,1,0]], \
[[0,0,1], [0,0,1], [1,0,1], [0,1,1]], \
[[0,0,1], [1,1,1], [1,0,1], [0,1,1]]]

make_stl(cube, 'cube.stl')
