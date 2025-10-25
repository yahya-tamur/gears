import struct

from math import sqrt # for normal vector

# The two functions below are included for convenience, make_stl can put any
# array of triangles (optionally including a normal vector) into an stl file.


# In case the direction of the normal matters, both functions below create a
# consistent direction, which can reversed:
# triangulate_polyhedron, using the 'reverse' argument, and
# triangulate_prism, by exchanging the top and bottom arguments.

# This can be used to stitch together any two sequences of points, not just the
# sides of a prism. 'closed' means the first and last points will be connected.

def triangulate_prism(top, bottom, closed=True):
    ans = []

    if closed:
        ans.append([top[-1],bottom[-1], top[0]]) 
        ans.append([top[0],bottom[-1], bottom[0]]) 

    for i in range(len(top)-1):
        ans.append([top[i], bottom[i], top[i+1]])
        ans.append([top[i+1], bottom[i], bottom[i+1]])
    return ans

# All triangles have a common center, which, if not provided, is the first
# element of the list.

def triangulate_polyhedron(p, center=None, reverse=False):
    start = 0
    ans = []
    if center is None:
        center = p[0]
        start = 1
    else:
        if reverse:
            ans.append([center, p[0], p[-1]])
        else:
            ans.append([center, p[-1], p[0]])


    for i in range(start, len(p)-1):
        if reverse:
            ans.append([center, p[i+1], p[i]])
        else:
            ans.append([center, p[i], p[i+1]])

    return ans

def normal(t):
    (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = t

    d1x, d1y, d1z = x1 - x0, y1 - y0, z1 - z0
    d2x, d2y, d2z = x2 - x0, y2 - y0, z2 - z0

    ans_x, ans_y, ans_z = d1y*d2z - d1z*d2y, - d1x*d2z + d1z*d2x, d1x*d2y - d1y*d2x
    
    n = sqrt(ans_x*ans_x + ans_y*ans_y + ans_z*ans_z)

    return (ans_x/n, ans_y/n, ans_z/n)

def make_stl(triangles, filename):
    with open(filename, 'wb') as f:
        f.write(bytearray(80))
        f.write(struct.pack('<i', len(triangles)))

        for tri in triangles:
            if len(tri) == 3:
                for ni in normal(tri):
                    f.write(struct.pack('<f', ni))



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
