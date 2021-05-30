import getfem as gf
import numpy as np
import pyvista as pv

pv.set_plot_theme("document")

# Yong modulus (MPa)
E = 10000.0
# Poisson ratio
nu = 0.0
# First Lame coefficient (MPa)
clambda = E * nu / ((1 + nu) * (1 - 2 * nu))
# Second Lame coefficient (MPa)
cmu = E / (2 * (1 + nu))
# Initially imperfect (mm)
alphas = [0.0, 0.02, 0.2]
# X length (mm)
b = 1.0
# Y length (mm)
h = 1.0
# Z length (mm)
L = 20.0
# Degree of the finite element methods
elements_degree = 2
# Displacements of the TOP boundary
displacements = np.linspace(0.0, 2.0, 100)

meshs = []
p = pv.Plotter(shape=(1, len(alphas)))
for i, alpha in enumerate(alphas):

    x = np.linspace(0.0, b, 4)
    y = np.linspace(0.0, h, 4)
    z = np.linspace(0.0, L, 80)

    mesh = gf.Mesh("cartesian", x, y, z)
    pts = mesh.pts()
    pts += np.array(
        [0.0 * pts[2], alpha * (1.0 - np.cos(np.pi * pts[2] / (L / 2.0))), 0.0 * pts[2]]
    )
    mesh.set_pts(pts)
    meshs.append(mesh)
    mesh.export_to_vtk("mesh" + str(i) + ".vtk", "ascii")

for i, alpha in enumerate(alphas):

    p.subplot(0, i)
    m = pv.read("mesh" + str(i) + ".vtk")
    p.add_mesh(m, show_edges=True)
    p.camera.zoom(2)
    p.show_grid()

p.show(screenshot="mesh.png", window_size=[1200, 1400])

for mesh in meshs:

    fb1 = mesh.outer_faces_with_direction([0.0, 0.0, 1.0], 0.01)
    fb2 = mesh.outer_faces_with_direction([0.0, 0.0, -1.0], 0.01)

    TOP_BOUND = 1
    BOTTOM_BOUND = 2

    mesh.set_region(TOP_BOUND, fb1)
    mesh.set_region(BOTTOM_BOUND, fb2)

mfus = []
mims = []
for mesh in meshs:

    mfu = gf.MeshFem(mesh, 3)
    mfu.set_classical_fem(elements_degree)
    mfus.append(mfu)

    mim = gf.MeshIm(mesh, pow(elements_degree, 2))
    mims.append(mim)

mds = []
for mfu, mim in zip(mfus, mims):

    md = gf.Model("real")
    md.add_fem_variable("u", mfu)
    mds.append(md)

for md, mfu, mim in zip(mds, mfus, mims):

    md.add_initialized_data("params", [clambda, cmu])
    lawname = "SaintVenant Kirchhoff"
    md.add_finite_strain_elasticity_brick(mim, lawname, "u", "params")

for md, mfu, mim in zip(mds, mfus, mims):

    md.add_initialized_data("r1", [0.0, 0.0, -2.0])
    md.add_initialized_data("r2", [0.0, 0.0, 0.0])
    md.add_initialized_data("H1", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    md.add_initialized_data("H2", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    md.add_generalized_Dirichlet_condition_with_multipliers(
        mim, "u", mfu, TOP_BOUND, "r1", "H1"
    )
    md.add_generalized_Dirichlet_condition_with_multipliers(
        mim, "u", mfu, BOTTOM_BOUND, "r2", "H2"
    )

for i, (md, mfu) in enumerate(zip(mds, mfus)):

    md.solve(
        "max_res",
        1e-9,
        "max_iter",
        100,
        "noisy",
        "lsearch",
        "simplest",
        "alpha min",
        0.8,
    )

    U = md.variable("u")
    mfu.export_to_vtk(
        "displacement" + str(i) + ".vtk", "ascii", mfu, U, "Displacements"
    )


p = pv.Plotter(shape=(1, len(alphas)))
for i, alpha in enumerate(alphas):

    p.subplot(0, i)
    d = pv.read("displacement" + str(i) + ".vtk")
    d.set_active_vectors("Displacements")
    p.add_mesh(d.warp_by_vector(factor=1.00), show_edges=True, color="white")
    p.camera_position = "yz"
    p.camera.enable_parallel_projection()
    p.camera.zoom(1.75)
    p.show_grid()

p.show(screenshot="displacement.png", window_size=[1200, 1400])
