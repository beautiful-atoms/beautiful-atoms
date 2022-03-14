import numpy as np
from math import pi


def ellipse(n, w, h):
    profile = np.zeros((n, 3))
    t = np.linspace(0, 1, n, endpoint=False)
    a = t*2*pi + pi/4
    profile[:, 0] = np.cos(a)*w
    profile[:, 2] = np.sin(a)*h
    return profile


def rectangle(w, h):
    profile = np.array([
        [-1*w,  0, -1*h],
        [-1*w,  0, 1*h],
        [1*w,  0, 1*h],
        [1*w,  0, -1*h],
    ])

    return profile


def build_mesh(vertices, normals, sides, profile, scales=None):
    """
    todo: improve performance
    """
    # calc side interp
    # tstart = time()
    n = len(vertices)
    m = len(profile)
    forwards = np.cross(normals, sides)
    rotations = np.zeros((n, 3, 3))
    sides = sides.reshape(n, 1, 3)
    forwards = forwards.reshape(n, 1, 3)
    normals = normals.reshape(n, 1, 3)
    # build new axis for a given vertices
    newAxis = np.concatenate((sides, forwards, normals), axis=1)
    # rotate and scale the profile
    rotations = np.linalg.inv(newAxis)
    if scales is not None:
        rotations = rotations*scales[:, None]
    vertices = np.tile(vertices, (m, 1))
    #
    for i in range(m):
        vertices[i*n:(i + 1)*n] += np.dot(rotations, profile[i])
    #
    faces = np.zeros((m*(n - 1), 4), dtype=int)
    for i in range(m - 1):
        faces[i*(n-1):(i+1)*(n-1), 0] = np.arange(i*n, (i+1)*n-1)
        faces[i*(n-1):(i+1)*(n-1), 1] = np.arange(i*n + 1, (i+1)*n)
        faces[i*(n-1):(i+1)*(n-1), 2] = np.arange((i+1)*n + 1, (i+2)*n)
        faces[i*(n-1):(i+1)*(n-1), 3] = np.arange((i+1)*n, (i+2)*n - 1)
    i = m - 1
    faces[i*(n-1):(i+1)*(n-1), 0] = np.arange(i*n, (i+1)*n-1)
    faces[i*(n-1):(i+1)*(n-1), 1] = np.arange(i*n + 1, (i+1)*n)
    faces[i*(n-1):(i+1)*(n-1), 2] = np.arange(1, n)
    faces[i*(n-1):(i+1)*(n-1), 3] = np.arange(0, n - 1)
    # cap
    faces = faces.tolist()
    faces.append([i*n for i in range(m)])
    faces.append([(i + 1)*n - 1 for i in range(m)])
    # print('build mesh: %s'%(time() - tstart))
    return vertices, faces


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # profile = ellipse(16, 1, 1)
    profile = rectangle(1, 1)
    plt.plot(profile[:, 0], profile[:, 1])
    plt.show()
