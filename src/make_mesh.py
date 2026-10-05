from .math import center, interpolate_line, rotate, distance, diameter, segment_point_distance
from math import pi

# 'closed' means the first and last points will be connected.


def mesh_prism(top, bottom, closed=True, make_center=True):
    ans = []

    if closed:
        ans += triangulate_polyhedron([top[-1], top[0], bottom[0], bottom[-1]])

    for i in range(len(top) - 1):
        ans += triangulate_polyhedron([top[i], top[i + 1], bottom[i + 1], bottom[i]])

    return ans

def mesh_polygon(mesh, p, reverse=False):

    center_point = center(p)

    for i in range(-1, len(p) - 1):
        if reverse:
            mesh.append([center_point, p[i + 1], p[i]])
        else:
            mesh.append([center_point, p[i], p[i + 1]])


# factor = 0 -> point,
# factor = 1 -> center 
# factor = -1 -> center + 2*(point - center)
def pull(point, center, factor):

    # point + factor*(center - point)
    # (1 - factor)point + factor*center
    a, b, c = point
    aa, bb, cc = center
    
    return ((1-factor)*a + factor*aa, (1-factor)*b + factor*bb, (1-factor)*c + factor*cc)

#  direction = 1 for inside, -1 for outside
# this will try to move by max_delta, but only move by the
# max diameter of the smallest segment
# if full_move is True, if the max diameter of the smallest segment is
# more than max_delta, instead of moving by max_delta, it will fail.
# 1 + (direction)max delta x = most scaling possible
# min_delta: if scaling 
def scale_reduce(mesh, line, max_delta, direction, reverse=False):
    center_point = center(line)
    segments = []
    line.append(line[0])
    i = 0

    while i < len(line):
        segments.append( line[i:i+4])
        i += 3

    delta = max_delta
    for segment in segments:
        c = center(segment)
        if len(segment) > 1: # ?? and distance(c, center_point) > 0:
            d = diameter(segment)
            delta = min(delta, d/distance(c, center_point))

    line_ = []
    for i in range(len(segments)):
        c = pull(center(segments[i]), center_point, delta*direction)
        for j in range(len(segments[i])-1):
            if reverse ^ (direction == -1):
                mesh.append([segments[i][j], c, segments[i][j+1]])
            else:
                mesh.append([segments[i][j], segments[i][j+1], c])
        line_.append(c)
    for i in range(len(line_)):
        if reverse ^ (direction == -1):
            mesh.append([line_[i], line[i*3], line_[i-1]])
        else:
            mesh.append([line_[i], line_[i-1], line[i*3]])
    return line_, delta

# nicely triangulates convex polygon made by line (set of coplanar points)
# reduces large number of points
def mesh_convex_polygon(mesh, line, reverse=False):
    if reverse:
        line = line[::-1]
    center_point = center(line)
    print(center_point)
    while len(line) > 6:
        pass
        segments = []
        line.append(line[0])
        i = 0
        while i < len(line):
            segments.append( line[i:i+4])
            i += 3
        factor = 0.3
        for segment in segments:
            c = center(segment)
            print("c", c)
            if len(segment) > 1: # ?? and distance(c, center_point) > 0:
                d = diameter(segment)
                factor = min(factor, d/distance(c, center_point))

        line_ = []
        for i in range(len(segments)):
            c = pull(center(segments[i]), center_point, factor)
            print("pulled",center(segments[i]), center_point, c)
            for j in range(len(segments[i])-1):
                mesh.append([segments[i][j], segments[i][j+1], c])
            line_.append(c)
        for i in range(len(line_)):
            mesh.append([line_[i], line_[i-1], line[i*3]])
        line = line_
    mesh_polygon(mesh, line)

def mesh_ring(mesh, outer, inner):
    # centers should be equal...
    center_point = center(outer+inner)
    def r_inner():
        return max(distance(p, center_point) for p in inner)
    def r_outer():
        return min(segment_point_distance(outer[i], outer[i+1], center_point) for i in range(-1, len(outer)-1))
    def ring_distance():
        return (r_outer() - r_inner())
    print("RD", ring_distance())
    def point_distance(ring):
        return max(distance(ring[i-1], ring[i]) for i in range(-1, len(ring)))
    inner_extent = r_inner() + 0.4*ring_distance()
    outer_extent = r_outer() - 0.4*ring_distance()

    while 2*point_distance(inner) < inner_extent - r_inner():
        inner, _delta = scale_reduce(mesh, inner, inner_extent/r_inner() - 1, -1)

    while  2*point_distance(outer) < r_outer() - outer_extent:
        outer, _delta = scale_reduce(mesh, outer, 1 - outer_extent/r_outer(), 1)

    # find closest outer point to inner[0]
    _, k = min(((distance(inner[0], outer[k]), k) for k in range(len(outer))))

    i = -1
    j = k - len(outer)
    print(len(inner), k)
    while i < len(inner) and j < k:
        print(i, j)
        if distance(inner[i+1-len(inner)], outer[j]) < distance(inner[i], outer[j+1]):
            mesh.append((inner[i], outer[j], inner[i+1-len(inner)]))
            i += 1
        else:
            mesh.append((inner[i], outer[j], outer[j+1]))
            j += 1

    # only one of the following will run based on how previous loop ended
    while i < len(inner):
        print("stitch A", j == k, i, j)
        mesh.append((inner[i], outer[k], inner[i+1 - len(inner)]))
        i += 1
    while j < k:
        print("stitch B", i == len(inner), i, j)
        mesh.append((inner[0], outer[j], outer[j+1]))
        j += 1



    #while outer_extent > point_distance(outer):
    #    outer, factor = scale_reduce(mesh, outer, point_distance(outer), 1)
    #    outer_extent -= factor


## grid = i by j array of points
def mesh_grid(mesh, grid, reverse=False):
    for i in range(len(grid)-1):
        for j in range(len(grid[0])-1):
            if reverse:
                mesh.append([grid[i][j], grid[i][j+1], grid[i+1][j]])
                mesh.append([grid[i+1][j], grid[i][j+1], grid[i+1][j+1]])
            else:                
                mesh.append([grid[i][j], grid[i+1][j], grid[i][j+1]])
                mesh.append([grid[i+1][j], grid[i+1][j+1], grid[i][j+1]])

def mesh_stitch(mesh, first_line, second_line, n):
    print(first_line[:3], second_line[:3])
    mesh_grid(mesh, [interpolate_line(v, w, n) for (v, w) in zip(first_line, second_line)])

def mesh_gear(mesh, tooth_west, tooth_east, tooth_number, flat_step, bore):

    teeth_west = [[line.copy() for line in tooth_west] for _ in range(tooth_number)]
    teeth_east = [[line.copy() for line in tooth_east] for _ in range(tooth_number)]

    for i in range(1, tooth_number):
        for line in teeth_west[i]:
            rotate([0, 0, 2*pi*i/tooth_number], line)
        for line in teeth_east[i]:
            rotate([0, 0, 2*pi*i/tooth_number], line)

    for tooth_west in teeth_west:
        mesh_grid(mesh,tooth_west)
    
    for tooth_east in teeth_east:
        mesh_grid(mesh,tooth_east, reverse=True)

    for i in range(-1, tooth_number-1):
        mesh_stitch(mesh, [line[0] for line in teeth_east[i+1]], [line[0] for line in teeth_west[i]], flat_step )
        mesh_stitch(mesh, [line[-1] for line in teeth_west[i]], [line[-1] for line in teeth_east[i]], flat_step )
        mesh_convex_polygon(mesh,
            teeth_west[i][0] +
            interpolate_line(teeth_west[i][0][-1], teeth_east[i][0][-1], flat_step, endpoints=False) +
            teeth_east[i][0][::-1] +
            interpolate_line(teeth_east[i][0][0], teeth_west[i][0][0], flat_step, endpoints=False) 
            )
        mesh_convex_polygon(mesh,
            teeth_west[i][-1] +
            interpolate_line(teeth_west[i][-1][-1], teeth_east[i][-1][-1], flat_step, endpoints=False) +
            teeth_east[i][-1][::-1] +
            interpolate_line(teeth_east[i][-1][0], teeth_west[i][-1][0], flat_step, endpoints=False),
            reverse=True
        )
        #mesh_polygon(mesh,
        #    teeth_west[i][-1] +
        #    interpolate_line(teeth_west[i][-1][-1], teeth_east[i][-1][-1], flat_step, endpoints=False) +
        #    teeth_east[i][-1][::-1],
        #    reverse=True)