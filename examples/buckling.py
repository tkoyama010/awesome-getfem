import getfem as gf
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt


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
alphas = [0.0, 0.01, 0.02, 0.1, 0.2]
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
forces = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
# Augmentation parameter for the augmented Lagrangian
gamma0 = 1.0 / E

mesh1s = []
mesh2s = []
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

    mesh2 = gf.Mesh("cartesian", x, y, np.array([L, L + 1.0]))
    mesh2s.append(mesh2)
    mesh2.export_to_vtk("mesh2" + str(i) + ".vtk", "ascii")

del i
del alpha

for i, alpha in enumerate(alphas):

    p.subplot(0, i)
    m = pv.read("mesh1" + str(i) + ".vtk")
    p.add_mesh(m, show_edges=True)
    p.camera.zoom(2)
    p.show_grid()

p.show(screenshot="mesh1.png", window_size=[1200, 1400])

del i
del alpha

for mesh1, mesh2 in zip(mesh1s, mesh2s):

    TOP_BOUND = 1
    BOTTOM_BOUND = 2
    CONTACT_BOUND = 3

    # P = mesh1.pts()
    # c1 = (P[2, :] > L - 1e-6)
    # pid1 = np.compress(c1, list(range(0, mesh1.nbpts())))
    # fb11 = mesh1.faces_from_pid(pid1)
    fb11 = mesh1.outer_faces_with_direction([0.0, 0.0, 1.0], 0.01)
    fb12 = mesh1.outer_faces_with_direction([0.0, 0.0, -1.0], 0.01)

    mesh1.set_region(CONTACT_BOUND, fb11)
    mesh1.set_region(BOTTOM_BOUND, fb12)

    fb21 = mesh2.outer_faces_with_direction([0.0, 0.0, 1.0], 0.01)
    fb22 = mesh2.outer_faces_with_direction([0.0, 0.0, -1.0], 0.01)

    mesh2.set_region(TOP_BOUND, fb21)
    mesh2.set_region(CONTACT_BOUND, fb22)

del mesh1
del mesh2

mfu1s = []
mim1s = []
mflambda_Cs = []
mflambdas = []
for mesh1, mesh2 in zip(mesh1s, mesh2s):

    mfu1 = gf.MeshFem(mesh1, 3)
    mfu1.set_classical_fem(elements_degree)
    mfu1s.append(mfu1)

    mim1 = gf.MeshIm(mesh1, pow(elements_degree, 2))
    mim1s.append(mim1)

    mflambda = gf.MeshFem(mesh1, 3)
    mflambda.set_classical_fem(elements_degree - 1)
    mflambdas.append(mflambda)

    mflambda_C = gf.MeshFem(mesh1, 1)
    mflambda_C.set_classical_fem(elements_degree - 1)
    mflambda_Cs.append(mflambda_C)

del mesh1
del mesh2
del mfu1

mds = []
for mfu1 in mfu1s:

    md = gf.Model("real")
    md.add_fem_variable("u1", mfu1)
    mds.append(md)

del mfu1

for md, mim1 in zip(mds, mim1s):

    if linear:
        md.add_initialized_data("cmu1", cmu)
        md.add_initialized_data("clambda1", clambda)
        md.add_isotropic_linearized_elasticity_brick(mim1, "u1", "clambda1", "cmu1")
    else:
        md.add_initialized_data("params", [clambda, cmu])
        lawname = "SaintVenant Kirchhoff"
        md.add_finite_strain_elasticity_brick(mim1, lawname, "u1", "params")

del md
del mim1

for md, mfu1, mim1, mflambda in zip(mds, mfu1s, mim1s, mflambdas):

    md.add_initialized_data("r11", [0.0, 0.0, 0.0])
    md.add_initialized_data("r12", [0.0, 0.0, 0.0])
    md.add_initialized_data("H11", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    md.add_initialized_data("H12", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    md.add_generalized_Dirichlet_condition_with_multipliers(
        mim1, "u1", mfu1, CONTACT_BOUND, "r11", "H11"
    )
    md.add_generalized_Dirichlet_condition_with_multipliers(
        mim1, "u1", mfu1, BOTTOM_BOUND, "r12", "H12"
    )

    md.add_variable("alpha_D", 1)
    md.add_filtered_fem_variable('lambda_D', mflambda, CONTACT_BOUND)
    md.add_initialized_data('F', [0.0])
    md.add_linear_term(mim1, "-lambda_D.Test_u1 + (alpha_D*[0;0;1]-u1).Test_lambda_D"
                             " + (lambda_D.[0;0;1]+F)*Test_alpha_D", CONTACT_BOUND)

del md
del mfu1
del mim1

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([0.0, -2.0])
ax.set_ylim([0.0, 100.0])

for i, (md, mfu1) in enumerate(zip(mds, mfu1s)):

    Fs = []
    alpha_Ds = []
    # for j, force in enumerate(forces):
    for j in range(100):
        md.set_variable("F", [j*1.0])
        iter_number = md.solve(
            "max_res",
            1e-9,
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
        alpha_D = md.variable("alpha_D")[0]
        if np.abs(alpha_D) > 20.0:
            break

        alpha_Ds.append(alpha_D)
        Fs.append(j * 1.0)
        U1 = md.variable("u1")
        mfu1.export_to_vtk(
            "mfu1-" + str(i) + "-" + "{:0=3}".format(j) + ".vtk",
            "ascii",
            mfu1,
            U1,
            "Displacements",
        )

    ax.plot(alpha_Ds, Fs)

plt.show()
fig.savefig("buckling.png")

del i
del md
del mfu1


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
