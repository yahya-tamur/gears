from math import sqrt


def plus3(v, w):
    (a, b, c) = v
    (d, e, f) = w
    return (a + d, b + e, c + f)


def minus3(v, w):
    (a, b, c) = v
    (d, e, f) = w
    return (a - d, b - e, c - f)


def scale(k, v):
    (a, b, c) = v
    return (k * a, k * b, k * c)


def dot3(v, w):
    (a, b, c) = v
    (d, e, f) = w
    return a * d + b * e + c * f


def norm3(v):
    return sqrt(dot3(v, v))


def unit(v):
    (a, b, c) = v
    n = norm3(v)
    return scale(1 / n, v)


def cross(v, w):
    (a, b, c) = v
    (d, e, f) = w
    return (b * f - c * e, c * d - a * f, a * e - b * d)


def distance(v, w):
    return norm3(minus3(v, w))


# line must have at least two points
def diameter(line):
    return max(
        max(distance(line[i], line[j]) for i in range(j)) for j in range(1, len(line))
    )


def segment_point_distance(p, q, center):

    c_minus_p = minus3(center, p)
    q_minus_p = minus3(q, p)

    t = dot3(c_minus_p, q_minus_p) / dot3(q_minus_p, q_minus_p)
    if t <= 0:
        return distance(p, center)
    if t >= 1:
        return distance(q, center)

    w = plus3(p, scale(t, q_minus_p))

    return distance(w, center)


def gram_schmidt3(v, w, z):
    v = unit(v)

    v1, v2, v3 = v
    w1, w2, w3 = w

    c = dot3(v, w)

    w = minus3(w, scale(c, v))
    w = unit(w)

    z_ = cross(v, w)

    if (dot3(z, z_)) < 0:
        z_ = scale(-1, z_)

    return (w, z_)


## next time you read this apply dot, scale, cross methods to methods above please
