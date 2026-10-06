from .math import (
    center,
    interpolate_line,
    rotate,
    distance,
    diameter,
    segment_point_distance,
    gram_schmidt3,
)
from math import pi, sin, cos

# 'closed' means the first and last points will be connected.


#def mesh_prism(top, bottom, closed=True, make_center=True):
#    ans = []

#    if closed:
#        ans += triangulate_polyhedron([top[-1], top[0], bottom[0], bottom[-1]])

#    for i in range(len(top) - 1):
#        ans += triangulate_polyhedron([top[i], top[i + 1], bottom[i + 1], bottom[i]])

#    return ans


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

    return (
        (1 - factor) * a + factor * aa,
        (1 - factor) * b + factor * bb,
        (1 - factor) * c + factor * cc,
    )


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
        segments.append(line[i : i + 4])
        i += 3

    delta = max_delta
    for segment in segments:
        c = center(segment)
        if len(segment) > 1:  # ?? and distance(c, center_point) > 0:
            d = diameter(segment)
            delta = min(delta, d / distance(c, center_point))

    line_ = []
    for i in range(len(segments)):
        c = pull(center(segments[i]), center_point, delta * direction)
        for j in range(len(segments[i]) - 1):
            if reverse ^ (direction == -1):
                mesh.append([segments[i][j], c, segments[i][j + 1]])
            else:
                mesh.append([segments[i][j], segments[i][j + 1], c])
        line_.append(c)
    for i in range(len(line_)):
        if reverse ^ (direction == -1):
            mesh.append([line_[i], line[i * 3], line_[i - 1]])
        else:
            mesh.append([line_[i], line_[i - 1], line[i * 3]])
    line.pop()
    return line_


# nicely triangulates convex polygon made by line (set of coplanar points)
# reduces large number of points
def mesh_convex_polygon(mesh, line, reverse=False):
    if reverse:
        line = line[::-1]
    while len(line) > 6:
        line = scale_reduce(mesh, line, 0.3, 1, reverse=False)
    mesh_polygon(mesh, line)


def mesh_ring(mesh, outer, inner, reverse=False):
    # centers should be equal...
    center_point = center(outer + inner)

    def r_inner():
        return max(distance(p, center_point) for p in inner)

    def r_outer():
        return min(
            segment_point_distance(outer[i], outer[i + 1], center_point)
            for i in range(-1, len(outer) - 1)
        )

    def ring_distance():
        return r_outer() - r_inner()

    def point_distance(ring):
        return max(distance(ring[i - 1], ring[i]) for i in range(-1, len(ring)))

    inner_extent = r_inner() + 0.4 * ring_distance()
    outer_extent = r_outer() - 0.4 * ring_distance()

    while 2 * point_distance(inner) < inner_extent - r_inner() and point_distance(inner)/r_inner() < 1/6:
        inner = scale_reduce(
            mesh, inner, inner_extent / r_inner() - 1, -1, reverse=reverse
        )

    while 2 * point_distance(outer) < r_outer() - outer_extent:
        outer = scale_reduce(
            mesh, outer, 1 - outer_extent / r_outer(), 1, reverse=reverse
        )

    # find closest outer point to inner[0]
    def closest(point, ls):
        _, ix = min((distance(point, ls[i]), i) for i in range(len(ls)))
        return ix

    i, j = 0, 0

    while (j != closest(inner[i], outer)) or (i != closest(outer[j], inner)):
        j = closest(inner[i], outer)
        i = closest(outer[j], inner)

    ii, jj = i, j
    i -= len(inner)
    j -= len(outer)
    a = 0
    while i < ii and j < jj:
        a += 1
        if distance(inner[i + 1], outer[j]) < distance(inner[i], outer[j + 1]):
            if reverse:
                mesh.append((inner[i], inner[i + 1], outer[j]))
            else:
                mesh.append((inner[i], outer[j], inner[i + 1]))
            i += 1
        else:
            if reverse:
                mesh.append((inner[i], outer[j + 1], outer[j]))
            else:
                mesh.append((inner[i], outer[j], outer[j + 1]))
            j += 1

    while i < ii:
        if reverse:
            mesh.append((inner[i], inner[i + 1], outer[jj]))
        else:
            mesh.append((inner[i], outer[jj], inner[i + 1]))
        i += 1
    while j < jj:
        if reverse:
            mesh.append((inner[ii], outer[j + 1], outer[j]))
        else:
            mesh.append((inner[ii], outer[j], outer[j + 1]))
        j += 1


## grid = i by j array of points
def mesh_grid(mesh, grid, reverse=False):
    for i in range(len(grid) - 1):
        for j in range(len(grid[0]) - 1):
            if reverse:
                mesh.append([grid[i][j], grid[i][j + 1], grid[i + 1][j]])
                mesh.append([grid[i + 1][j], grid[i][j + 1], grid[i + 1][j + 1]])
            else:
                mesh.append([grid[i][j], grid[i + 1][j], grid[i][j + 1]])
                mesh.append([grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]])


def mesh_stitch(mesh, first_line, second_line, n):
    grid = [interpolate_line(v, w, n) for (v, w) in zip(first_line, second_line)]
    mesh_grid(
        mesh, [interpolate_line(v, w, n) for (v, w) in zip(first_line, second_line)]
    )
    return grid[0], grid[-1]


def mesh_gear(mesh, tooth_west, tooth_east, tooth_number, flat_step, bore, bore_steps):

    teeth_west = [[line.copy() for line in tooth_west] for _ in range(tooth_number)]
    teeth_east = [[line.copy() for line in tooth_east] for _ in range(tooth_number)]

    for i in range(1, tooth_number):
        for line in teeth_west[i]:
            rotate([0, 0, 2 * pi * i / tooth_number], line)
        for line in teeth_east[i]:
            rotate([0, 0, 2 * pi * i / tooth_number], line)

    for tooth_west in teeth_west:
        mesh_grid(mesh, tooth_west)

    for tooth_east in teeth_east:
        mesh_grid(mesh, tooth_east, reverse=True)

    top_edge = []
    bottom_edge = []

    for i in range(-1, tooth_number - 1):
        bottom_between, top_between = mesh_stitch(
            mesh,
            [line[0] for line in teeth_east[i + 1]],
            [line[0] for line in teeth_west[i]],
            flat_step,
        )
        bottom_tip, top_tip = mesh_stitch(
            mesh,
            [line[-1] for line in teeth_west[i]],
            [line[-1] for line in teeth_east[i]],
            flat_step,
        )
        bottom_base = interpolate_line(
            teeth_east[i][0][0], teeth_west[i][0][0], flat_step, endpoints=False
        )
        top_base = interpolate_line(
            teeth_east[i][-1][0], teeth_west[i][-1][0], flat_step, endpoints=False
        )
        mesh_convex_polygon(
            mesh,
            teeth_west[i][-1] + top_tip[1:-1] + teeth_east[i][-1][::-1] + top_base,
            reverse=True,
        )
        mesh_convex_polygon(
            mesh,
            teeth_west[i][0] + bottom_tip[1:-1] + teeth_east[i][0][::-1] + bottom_base,
        )
        top_edge += top_base + top_between[::-1]
        bottom_edge += bottom_base + bottom_between[::-1]

    if bore == 0:
        mesh_convex_polygon(mesh, top_edge)
        mesh_convex_polygon(mesh, bottom_edge, reverse=True)
    else:
        top_1, top_2, top_3 = center(top_edge)
        bottom_1, bottom_2, bottom_3 = center(bottom_edge)
        (v1, v2, v3), (w1, w2, w3) = gram_schmidt3(
            (top_1 - bottom_1, top_2 - bottom_2, top_3 - bottom_3),
            (top_edge[0][0] - top_1, top_edge[0][1] - top_2, top_edge[0][2] - top_3),
            (
                top_edge[1][0] - top_edge[0][0],
                top_edge[1][1] - top_edge[0][1],
                top_edge[2][2] - top_edge[0][2],
            ),
        )
        top_bore = []
        bottom_bore = []
        for i in range(bore_steps):
            t = i * 2 * pi / bore_steps
            top_bore.append(
                (
                    top_1 + bore * cos(t) * v1 + bore * sin(t) * w1,
                    top_2 + bore * cos(t) * v2 + bore * sin(t) * w2,
                    top_3 + bore * cos(t) * v3 + bore * sin(t) * w3,
                )
            )
            bottom_bore.append(
                (
                    bottom_1 + bore * cos(t) * v1 + bore * sin(t) * w1,
                    bottom_2 + bore * cos(t) * v2 + bore * sin(t) * w2,
                    bottom_3 + bore * cos(t) * v3 + bore * sin(t) * w3,
                )
            )

        mesh_ring(mesh, top_edge, top_bore)
        mesh_ring(mesh, bottom_edge, bottom_bore, reverse=True)

        top_bore.append(top_bore[0])
        bottom_bore.append(bottom_bore[0])

        mesh_stitch(mesh, top_bore, bottom_bore, len(teeth_west[0]))
