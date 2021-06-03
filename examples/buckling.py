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
# Linear
linear = False
# Force
# forces = np.arange(0.0, 100.0 + 1.0, 1.0)
# forces = np.array([0.0, 10.0, 20.0, 30.0, 4.0, 50.0, 60.0, 70.0, 80.0, 90.0 , 100.0])
forces = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0 , 100.0])

mesh1s = []
p = pv.Plotter(shape=(1, len(alphas)))
for i, alpha in enumerate(alphas):

    x = np.linspace(0.0, b, 4)
    y = np.linspace(0.0, h, 4)
    z = np.linspace(0.0, L, 80)

    mesh1 = gf.Mesh("cartesian", x, y, z)
    pts = mesh1.pts()
    pts += np.array(
        [0.0 * pts[2], alpha * (1.0 - np.cos(np.pi * pts[2] / (L / 2.0))), 0.0 * pts[2]]
    )
    mesh1.set_pts(pts)
    # mesh1.merge(gf.Mesh("cartesian", x, y, np.array([L, L + 1.0])))

    mesh1s.append(mesh1)
    mesh1.export_to_vtk("mesh1" + str(i) + ".vtk", "ascii")

for i, alpha in enumerate(alphas):

    p.subplot(0, i)
    m = pv.read("mesh1" + str(i) + ".vtk")
    p.add_mesh(m, show_edges=True)
    p.camera.zoom(2)
    p.show_grid()

p.show(screenshot="mesh1.png", window_size=[1200, 1400])

for mesh1 in mesh1s:

    # P = mesh1.pts()
    # c1 = (P[2, :] > L - 1e-6)
    # pid1 = np.compress(c1, list(range(0, mesh1.nbpts())))
    # fb1 = mesh1.faces_from_pid(pid1)
    fb1 = mesh1.outer_faces_with_direction([0.0, 0.0, 1.0], 0.01)
    fb2 = mesh1.outer_faces_with_direction([0.0, 0.0, -1.0], 0.01)

    TOP_BOUND = 1
    BOTTOM_BOUND = 2

    mesh1.set_region(TOP_BOUND, fb1)
    mesh1.set_region(BOTTOM_BOUND, fb2)

mfus = []
mims = []
for mesh1 in mesh1s:

    mfu = gf.MeshFem(mesh1, 3)
    mfu.set_classical_fem(elements_degree)
    mfus.append(mfu)

    mim = gf.MeshIm(mesh1, pow(elements_degree, 2))
    mims.append(mim)

mds = []
for mfu, mim in zip(mfus, mims):

    md = gf.Model("real")
    md.add_fem_variable("u", mfu)
    mds.append(md)

for md, mfu, mim in zip(mds, mfus, mims):

    if linear:
        md.add_initialized_data('cmu1', cmu)
        md.add_initialized_data('clambda1', clambda)
        md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'clambda1', 'cmu1')
        md.add_initialized_data('cmu2', cmu*1000.0)
        md.add_initialized_data('clambda2', clambda*1000.0)
        md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'clambda2', 'cmu2', TOP_BOUND)
    else:
        md.add_initialized_data("params", [clambda, cmu])
        md.add_initialized_data("lambda", clambda*1000.0)
        md.add_initialized_data("mu", cmu*1000.0)
        lawname = "SaintVenant Kirchhoff"
        md.add_finite_strain_elasticity_brick(mim, lawname, "u", "params")
        md.add_isotropic_linearized_elasticity_brick(mim, "u", "lambda", "mu", TOP_BOUND)

for md, mfu, mim in zip(mds, mfus, mims):

    md.add_initialized_data("r1", [0.0, 0.0, 0.0])
    md.add_initialized_data("r2", [0.0, 0.0, 0.0])
    md.add_initialized_data("H1", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    md.add_initialized_data("H2", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    md.add_generalized_Dirichlet_condition_with_multipliers(
        mim, "u", mfu, TOP_BOUND, "r1", "H1"
    )
    md.add_generalized_Dirichlet_condition_with_multipliers(
        mim, "u", mfu, BOTTOM_BOUND, "r2", "H2"
    )

curves = []
for i, (md, mfu, mim) in enumerate(zip(mds, mfus, mims)):
    # md.add_initialized_data("force", [0.0, 0.0, 0.0])

    Fs = []
    displacements = []
    # for j, force in enumerate(forces):
    for j in range(100):
        # md.set_variable("force", [0.0, 0.0, -force])
        md.add_source_term_brick(mim, "u", "[0.0, 0.0, -1.0]", TOP_BOUND)
        iter_number = md.solve(
            "max_res",
            1e-6,
            "max_iter",
            200,
            "noisy",
            "lsearch",
            "simplest",
            "alpha min",
            0.8,
        )

        if iter_number == 200:
           break
        displacement = (((md.variable("u"))[mfu.dof_on_region(TOP_BOUND)]).reshape(-1, 3))[0, 2]
        if np.abs(displacement) > 20.0:
          break

        displacements.append(displacement)
        Fs.append(j*1.0)
        U = md.variable("u")
        mfu.export_to_vtk(
            "displacement" + str(i) + "-" + "{:0=3}".format(j) + ".vtk", "ascii", mfu, U, "Displacements"
        )

    curves.append([Fs, displacements])


# p = pv.Plotter(shape=(1, len(alphas)))
# for i, alpha in enumerate(alphas):
# 
#     p.subplot(0, i)
#     d = pv.read("displacement" + str(i) + ".vtk")
#     d.set_active_vectors("Displacements")
#     p.add_mesh(d.warp_by_vector(factor=1.00), show_edges=True, color="white")
#     p.camera_position = "yz"
#     p.camera.enable_parallel_projection()
#     p.camera.zoom(1.75)
#     p.show_grid()
# 
# p.show(screenshot="displacement.png", window_size=[1200, 1400])
