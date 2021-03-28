"""
Compare the FEM results and analytical results in Hertz problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run a Hertz problem and compare the results.

"""
###############################################################################
# # A finite element analysis of elastic contact problem
#
# In this example of a deformable ``cylinder 1`` enters in contact with a
# deformable  ``cylinder 2``. We use here python interface, translate this
# program for another interface or in C++ is easy (see the previous example).
#
# ![Kontakt_paralleler_Zylinder](https://upload.wikimedia.org/wikipedia/commons/e/ef/Kontakt_paralleler_Zylinder.jpg)

# ## The problem setting
#
# Let $\Omega^{1} \subset \mathbb{R}^{2}$ be the reference of a 2D cylinder 1
# and $\Omega^{2} \subset \mathbb{R}^{2}$ the reference configuration of a
# deformable cylinder. We consider small deformation of these two bodies
# (linearized elasticity) and the contact between them.

###############################################################################
# ## Building the program
#
# Let us begin by loading Getfem and fixing the parameters of the problem


import getfem as gf
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

pv.set_plot_theme("document")

# Yong Modulus (MPa)
E = 200000.0
# Poisson ratio
nu = 0.3
# First Lame coefficient (MPa)
clambda = E * nu / ((1 + nu) * (1 - 2 * nu))
# Second Lame coefficient (MPa)
cmu = E / (2 * (1 + nu))
# Lame coefficient for Plane stress (MPa)
clambdastar = 2 * clambda * cmu / (clambda + 2 * cmu)
# Degree of the finite element methods
elements_degree = 2
# Force at the top boundary (N/mm)
F = -200.0 / 10
# Augmentation parameter for the augmented Lagrangian
gamma0 = 1.0 / E


###############################################################################
#
# We consider that the radius of the two cylinder is 5mm. We load the mesh of
# the cylinder using the load of a mesh from a GetFEM ascii mesh file (see the
# documentation of the Mesh object in the python interface).
# !gmsh hertz.mesh -f msh2 -save -o hertz.msh
mesh = gf.Mesh("import", "gmsh", "hertz.msh")
mesh.translate([0.0, 5.0])
P = mesh.pts()

# cvid1
c1 = P[1, :] >= 0.0
pid1 = np.compress(c1, list(range(0, mesh.nbpts())))
cvid1 = mesh.cvid_from_pid(pid1)

# cvid2
c2 = P[1, :] <= 0.0
pid2 = np.compress(c2, list(range(0, mesh.nbpts())))
cvid2 = mesh.cvid_from_pid(pid2)

# Approximate mesh size
h = 1.0

mesh1 = gf.Mesh("clone", mesh)
mesh1.del_convex(cvid2)

mesh2 = gf.Mesh("clone", mesh)
mesh2.del_convex(cvid1)

mesh3 = gf.Mesh("clone", mesh2)
mesh3.translate([0.0, 10.0])
mesh1.merge(mesh3)
mesh1.export_to_vtk("mesh1.vtk", "ascii")

mesh4 = gf.Mesh("clone", mesh1)
mesh4.translate([0.0, -10.0])
mesh2.merge(mesh4)
mesh2.export_to_vtk("mesh2.vtk", "ascii")

X = np.arange(0.0, 5.0 + 0.1, 0.1)
Y = np.arange(0.0, 2.0 + 0.1, 0.1)

mesh5 = gf.Mesh("cartesian", X, Y)
mesh5.translate([0.0, 10.0])
mesh5.export_to_vtk("mesh5.vtk", "ascii")

mesh6 = gf.Mesh("cartesian", X, Y)
mesh6.translate([0.0, -12.0])
mesh6.export_to_vtk("mesh6.vtk", "ascii")

###############################################################################
#
# The result is the following


m1 = pv.read("mesh1.vtk")
m2 = pv.read("mesh2.vtk")
m5 = pv.read("mesh5.vtk")
m6 = pv.read("mesh6.vtk")
p = pv.Plotter(shape=(1, 1))
p.subplot(0, 0)
p.add_mesh(m1, show_edges=True)
p.add_mesh(m2, show_edges=True)
p.add_mesh(m5, show_edges=True)
p.add_mesh(m6, show_edges=True)
p.show_grid()
p.show(screenshot="mesh.png", window_size=[1200, 1400], cpos="xy")


###############################################################################
#
# ## Boundary selection

# We have to select the different parts of the boundary where we will set some
# boundary conditions, namely the boundary of the rim (in order to apply a
# force and the fact that the rim is rigid), the contact boundary of the wheel
# and the bottom boundary of the foundation that we will assume clamped.

fb11 = mesh1.outer_faces_with_direction([+0.0, +1.0], np.pi / 4.0)
fb12 = mesh1.outer_faces_with_direction([+0.0, -1.0], np.pi / 4.0)
fb13 = mesh1.outer_faces_with_direction([-1.0, +0.0], np.pi / 4.0)
fb21 = mesh2.outer_faces_with_direction([+0.0, +1.0], np.pi / 4.0)
fb22 = mesh2.outer_faces_with_direction([+0.0, -1.0], np.pi / 4.0)
fb23 = mesh2.outer_faces_with_direction([-1.0, +0.0], np.pi / 4.0)
fb51 = mesh5.outer_faces_with_direction([+0.0, +1.0], np.pi / 4.0)
fb52 = mesh5.outer_faces_with_direction([+0.0, -1.0], np.pi / 4.0)
fb53 = mesh5.outer_faces_with_direction([-1.0, +0.0], np.pi / 4.0)
fb61 = mesh6.outer_faces_with_direction([+0.0, +1.0], np.pi / 4.0)
fb62 = mesh6.outer_faces_with_direction([+0.0, -1.0], np.pi / 4.0)
fb63 = mesh6.outer_faces_with_direction([-1.0, +0.0], np.pi / 4.0)

TOP_BOUND = 1
BOTTOM_BOUND = 2
LEFT_BOUND = 3

mesh1.set_region(TOP_BOUND, fb11)
mesh2.set_region(TOP_BOUND, fb21)
mesh5.set_region(TOP_BOUND, fb51)
mesh6.set_region(TOP_BOUND, fb61)
mesh1.set_region(BOTTOM_BOUND, fb12)
mesh2.set_region(BOTTOM_BOUND, fb22)
mesh5.set_region(BOTTOM_BOUND, fb52)
mesh6.set_region(BOTTOM_BOUND, fb62)
mesh1.set_region(LEFT_BOUND, fb13)
mesh2.set_region(LEFT_BOUND, fb23)
mesh5.set_region(LEFT_BOUND, fb53)
mesh6.set_region(LEFT_BOUND, fb63)


###############################################################################
#
# Note that the command `mesh1.outer_faces_with_direction([0., -1.0], np.pi/6)`
# allows to select all the faces having a unit outward normal having an angle
# less or equal to `np.pi/6` with the vector `[0., -1.0]`.
#
# ## Definition of finite elements methods and integration method
#
# We define mfu1, mfu2 two finite element methods which will approximate the
# displacements in the `cylinder1` and `cylinder2`. `mflambda_C1` is to
# approximate the contact multiplier (contact pressure) and `mfvm1`, `mfvm2`
# will be used to interpolate the Von Mises stresses of the wheel and the
# foundation for post-processing. `mim1`, `mim2` are two integration methods
# on the `cylinder1` and the `cylinder2`.


mfu1 = gf.MeshFem(mesh1, 2)
mfu1.set_classical_fem(elements_degree)
mfd1 = gf.MeshFem(mesh1, 1)
mfd1.set_classical_fem(elements_degree)
mflambda_C1 = gf.MeshFem(mesh1, 1)
mflambda_C1.set_classical_fem(elements_degree - 1)
mfu2 = gf.MeshFem(mesh2, 2)
mfu2.set_classical_fem(elements_degree)
mfd2 = gf.MeshFem(mesh2, 1)
mfd2.set_classical_fem(elements_degree)
mflambda_C2 = gf.MeshFem(mesh2, 1)
mflambda_C2.set_classical_fem(elements_degree)
mfvm1 = gf.MeshFem(mesh1, 1)
mfvm1.set_classical_discontinuous_fem(elements_degree)
mfvm2 = gf.MeshFem(mesh2, 1)
mfvm2.set_classical_discontinuous_fem(elements_degree)
mim1 = gf.MeshIm(mesh1, pow(elements_degree, 2))
mim2 = gf.MeshIm(mesh2, pow(elements_degree, 2))

mfu5 = gf.MeshFem(mesh5, 2)
mfu5.set_classical_fem(elements_degree)
mfd5 = gf.MeshFem(mesh5, 1)
mfd5.set_classical_fem(elements_degree)
mfu6 = gf.MeshFem(mesh6, 2)
mfu6.set_classical_fem(elements_degree)
mfd6 = gf.MeshFem(mesh6, 1)
mfd6.set_classical_fem(elements_degree)
mfvm5 = gf.MeshFem(mesh5, 1)
mfvm5.set_classical_discontinuous_fem(elements_degree)
mfvm6 = gf.MeshFem(mesh6, 1)
mfvm6.set_classical_discontinuous_fem(elements_degree)
mim5 = gf.MeshIm(mesh5, pow(elements_degree, 2))
mim6 = gf.MeshIm(mesh6, pow(elements_degree, 2))

###############################################################################
#
# ## Model definition
#
# We use a real model and declare the two variables which will represent the
# displacements:


md = gf.Model("real")
md.add_fem_variable("u1", mfu1)
md.add_fem_variable("u2", mfu2)
md.add_fem_variable("u5", mfu5)
md.add_fem_variable("u6", mfu6)

###############################################################################
#
# ## Linearized elasticity bricks
#
# We add the Lame coefficients as data of the model and add a linearized
# elasticity brick for the wheel and the foundation:

md.add_initialized_data("cmu", [cmu])
md.add_initialized_data("clambdastar", [clambdastar])
md.add_initialized_data("E1", [E])
md.add_initialized_data("E2", [E * 1000.0])
md.add_initialized_data("nu", [nu])
md.add_isotropic_linearized_elasticity_brick_pstrain(mim1, "u1", "E1", "nu")
md.add_isotropic_linearized_elasticity_brick_pstrain(mim2, "u2", "E1", "nu")
md.add_isotropic_linearized_elasticity_brick_pstrain(mim5, "u5", "E2", "nu")
md.add_isotropic_linearized_elasticity_brick_pstrain(mim6, "u6", "E2", "nu")

###############################################################################
#
# ## Clamped condition at the bottom boundary
#
# We prescribed the displacement at bottom face of the foundation to vanish,
# for instance with a multiplier with the add of the following brick:

md.add_initialized_data("r", [0, 0])
md.add_initialized_data("H1", [[1, 0], [0, 0]])
md.add_initialized_data("H2", [[0, 0], [0, 1]])
md.add_initialized_data("F", F)

md.add_generalized_Dirichlet_condition_with_multipliers(
    mim1, "u1", mfu1, LEFT_BOUND, "r", "H1"
)
md.add_generalized_Dirichlet_condition_with_multipliers(
    mim2, "u2", mfu2, LEFT_BOUND, "r", "H1"
)
md.add_source_term_brick(mim5, "u5", "[0.0, F]", TOP_BOUND)
md.add_generalized_Dirichlet_condition_with_multipliers(
    mim5, "u5", mfu5, LEFT_BOUND, "r", "H1"
)
md.add_generalized_Dirichlet_condition_with_multipliers(
    mim6, "u6", mfu6, BOTTOM_BOUND, "r", "H2"
)
md.add_generalized_Dirichlet_condition_with_multipliers(
    mim6, "u6", mfu6, LEFT_BOUND, "r", "H1"
)

###############################################################################
#

# ## Contact condition (use of interpolate transformations)

# Now, let us see how to prescribed the contact condition between the two
# structures. It is possible to use predefined bricks (see [Small sliding
# contact with friction bricks](https://getfem.readthedocs.io/en/latest/userdoc/model_contact_friction.html#ud-model-contact-friction)
# for small deformation/small sliding contact and [Large sliding/large deformation contact with friction bricks](https://getfem.readthedocs.io/en/latest/userdoc/model_contact_friction_large_sliding.html#ud-model-contact-friction-large)
# for large deformation/large sliding contact).
# However, we will see here how to directly prescribe a contact condition using
# an augmented Lagrangian formulation and the interpolate transformations.
#
# For small deformation contact, the correspondence between points of one
# contact surface to the other have to be described on the reference
# configuration and is not evolving, which is of course simpler but is an
# approximation.
#
# We consider that the contact boundary of the wheel is the slave one and we
# have to describe the transformation from the contact boundary of the wheel to
# the contact boundary of the foundation. This is quite simple here, since the
# contact boundary of the foundation corresponds to a vanishing vertical
# coordinate. So we define the transformation
#
# $$
# X \longmapsto (X(1), -X(2))
# $$
#
# where $X$ is the vector of coordinates of the point. We add this
# transformation to the model with the command

md.add_interpolate_transformation_from_expression(
    "Proj12", mesh1, mesh2, "[X(1);-X(2)-0.001]"
)
md.add_interpolate_transformation_from_expression(
    "Proj15", mesh1, mesh5, "[X(1);+10.000]"
)
md.add_interpolate_transformation_from_expression(
    "Proj26", mesh2, mesh6, "[X(1);-10.000]"
)

###############################################################################
#

# As a consequence, it will be possible to use this transformation, from the
# mesh of the wheel to the mesh of the foundation, into GWFL expressions.
# Notes that this is here a very simple constant expression.
# More complex expressions depending on the data or even the variables of the
# model can be used.
# If the expression of a transformation depends on the variable of the model,
# the tangent linear system will automatically takes into account this
# dependence (see [Interpolate transformations](https://getfem.readthedocs.io/en/latest/userdoc/gasm_high.html#ud-gasm-high-transf) for more details).
# Note also that transformation corresponding to a large sliding contact and
# automatically searching for the correspondence between contact boundaries
# exist in GetFEM (see [Integral contact brick with raytrace](https://getfem.readthedocs.io/en/latest/userdoc/model_contact_friction_large_sliding.html#ud-model-contact-friction-large-hlgav)).
#
# Using the defined transformation, we can write an integral contact condition
# using an augmented Lagrangian formulation (see [Small sliding contact with friction bricks](https://getfem.readthedocs.io/en/latest/userdoc/model_contact_friction.html#ud-model-contact-friction)
# for more details).
# The corresponding term (to be added to the rest of the weak formulation)
# reads:
#
# $$
#   \cdots +  \int_{\Gamma_c} \lambda_N(X) (\delta_{u^1}(X)-\delta_{u^2}(\Pi(X)))\cdot n d\Gamma \\
#   -   \int_{\Gamma_c} \left(\lambda_N(X) + \left(\lambda_N(X) + \dfrac{1}{h_T\gamma_0}((X + u^1(X))\cdot n - (\Pi(X) - u^2(\Pi(X)))\cdot n\right)_-\right)\delta_{\lambda_N}(X) d\Gamma = 0 ~~~~ \forall \delta_{\lambda_N}, \forall \delta_{u^1}, \forall \delta_{u^2},
# $$
#
# where $\Gamma_c$ is the slave contact boundary, $\lambda_N$ is the contact
# multiplier (contact pressure), $h_T$ is the radius of the element, $\Pi$ is
# the transformation, $n$ is the outward normal vector to the master contact
# boundary (here $n = (0,1)$), $\gamma_0$ is an augmentation parameter, $(\cdot)_-:I\hspace{-0.2em}R\rightarrow I\hspace{-0.2em}R_+$
# is the negative part and $\delta_{\lambda_N}, \delta_{u^1}, \delta_{u^2}$ are
# the test  functions corresponding to $\lambda_N, u^1, u^2$, respectively.
#
# Using GWFL, the contact condition can be added by:

md.add_initialized_data("gamma0", [gamma0])
md.add_filtered_fem_variable("lambda15", mflambda_C1, TOP_BOUND)
md.add_macro("n15", "[0;-1]")
md.add_nonlinear_term(
    mim1,
    "lambda15*(Test_u1.n15)-lambda15*(Interpolate(Test_u5,Proj15).n15)",
    TOP_BOUND,
)
md.add_nonlinear_term(
    mim1,
    "-(gamma0*element_size)"
    "*(lambda15 + neg_part(lambda15+(1/(gamma0*element_size))"
    "*((u1-Interpolate(u5,Proj15)+X-Interpolate(X,Proj15)).n15)))*Test_lambda15",
    TOP_BOUND,
)
md.add_filtered_fem_variable("lambda12", mflambda_C1, BOTTOM_BOUND)
md.add_macro("n12", "[0;+1]")
md.add_nonlinear_term(
    mim1,
    "lambda12*(Test_u1.n12)-lambda12*(Interpolate(Test_u2,Proj12).n12)",
    BOTTOM_BOUND,
)
md.add_nonlinear_term(
    mim1,
    "-(gamma0*element_size)"
    "*(lambda12 + neg_part(lambda12+(1/(gamma0*element_size))"
    "*((u1-Interpolate(u2,Proj12)+X-Interpolate(X,Proj12)).n12)))*Test_lambda12",
    BOTTOM_BOUND,
)
md.add_filtered_fem_variable("lambda26", mflambda_C2, BOTTOM_BOUND)
md.add_nonlinear_term(
    mim2,
    "lambda26*(Test_u2.n12)-lambda26*(Interpolate(Test_u6,Proj26).n12)",
    BOTTOM_BOUND,
)
md.add_nonlinear_term(
    mim2,
    "-(gamma0*element_size)"
    "*(lambda26 + neg_part(lambda26+(1/(gamma0*element_size))"
    "*((u2-Interpolate(u6,Proj26)+X-Interpolate(X,Proj26)).n12)))*Test_lambda26",
    BOTTOM_BOUND,
)


###############################################################################
#

# ## Prescribing the rigidity of the rim and the vertical force
# We have now to prescribe the rigidity of the rim. This is a non-standard
# condition, since we do not know a priori what will be the vertical
# displacement of the rim. We can use an additional unknown for that vertical
# displacement. We need a multiplier to prescribe the displacement on the rim
# boundary:

# ## Model solve
# We can now solve our problem with:

###############################################################################
#


md.solve("max_res", 1e-5, "max_iter", 100, "noisy")


###############################################################################
# Compute Cauchy stress tensor with direction Y-axis in global cooridnate.
#

U1 = md.variable("u1")

Grad_u1 = gf.compute_gradient(mfu1, U1, mfd1)
Div_u1 = np.trace(Grad_u1).reshape(-1, 1)
I = np.eye(2)
sigma = clambda * np.outer(Div_u1, I).reshape(2, 2, -1) + cmu * (
    Grad_u1 + Grad_u1.transpose(1, 0, 2)
)
n = np.array([0.0, 1.0])
t1 = sigma[1, 0] * n[0] + sigma[1, 1] * n[1]

U2 = md.variable("u2")

Grad_u2 = gf.compute_gradient(mfu2, U2, mfd2)
Div_u2 = np.trace(Grad_u2).reshape(-1, 1)
I = np.eye(2)
sigma = clambda * np.outer(Div_u2, I).reshape(2, 2, -1) + cmu * (
    Grad_u2 + Grad_u2.transpose(1, 0, 2)
)
n = np.array([0.0, 1.0])
t2 = sigma[1, 0] * n[0] + sigma[1, 1] * n[1]

U5 = md.variable("u5")

Grad_u5 = gf.compute_gradient(mfu5, U5, mfd5)
Div_u5 = np.trace(Grad_u5).reshape(-1, 1)
I = np.eye(2)
sigma = clambda * np.outer(Div_u5, I).reshape(2, 2, -1) + cmu * (
    Grad_u5 + Grad_u5.transpose(1, 0, 2)
)
n = np.array([0.0, 1.0])
t5 = sigma[1, 0] * n[0] + sigma[1, 1] * n[1]

U6 = md.variable("u6")

Grad_u6 = gf.compute_gradient(mfu6, U6, mfd6)
Div_u6 = np.trace(Grad_u6).reshape(-1, 1)
I = np.eye(2)
sigma = clambda * np.outer(Div_u6, I).reshape(2, 2, -1) + cmu * (
    Grad_u6 + Grad_u6.transpose(1, 0, 2)
)
n = np.array([0.0, 1.0])
t6 = sigma[1, 0] * n[0] + sigma[1, 1] * n[1]


###############################################################################
# Note that in some configuration, it is preferable to use a more basic line
# search than the default one:
# md.solve(
#     "max_res", 1e-9, "max_iter", 100, "noisy", "lsearch", "simplest", "alpha min", 0.8
# )
#
# ## Export the solution
#
# Now the code to export the solution with the VonMises stress:


VM1 = md.compute_isotropic_linearized_Von_Mises_or_Tresca(
    "u1", "clambdastar", "cmu", mfvm1
)
VM2 = md.compute_isotropic_linearized_Von_Mises_or_Tresca(
    "u2", "clambdastar", "cmu", mfvm2
)
VM5 = md.compute_isotropic_linearized_Von_Mises_or_Tresca(
    "u5", "clambdastar", "cmu", mfvm5
)
VM6 = md.compute_isotropic_linearized_Von_Mises_or_Tresca(
    "u6", "clambdastar", "cmu", mfvm6
)

mfvm1.export_to_vtk(
    "displacement_with_von_mises1.vtk",
    "ascii",
    mfvm1,
    VM1,
    "Von Mises Stresses",
    mfu1,
    U1,
    "Displacements",
)

mfd1.export_to_vtk(
    "stress1.vtk", "ascii", mfd1, t1, "Stress",
)

mfvm2.export_to_vtk(
    "displacement_with_von_mises2.vtk",
    "ascii",
    mfvm2,
    VM2,
    "Von Mises Stresses",
    mfu2,
    U2,
    "Displacements",
)

mfd2.export_to_vtk(
    "stress2.vtk", "ascii", mfd2, t2, "Stress",
)

mfvm5.export_to_vtk(
    "displacement_with_von_mises5.vtk",
    "ascii",
    mfvm5,
    VM5,
    "Von Mises Stresses",
    mfu5,
    U5,
    "Displacements",
)

mfd5.export_to_vtk(
    "stress5.vtk", "ascii", mfd5, t5, "Stress",
)

mfvm6.export_to_vtk(
    "displacement_with_von_mises6.vtk",
    "ascii",
    mfvm6,
    VM6,
    "Von Mises Stresses",
    mfu6,
    U6,
    "Displacements",
)

mfd6.export_to_vtk(
    "stress6.vtk", "ascii", mfd6, t6, "Stress",
)


###############################################################################
# You can view solutions with pyvista:
#

m1 = pv.read("mesh1.vtk")
m2 = pv.read("mesh2.vtk")
m5 = pv.read("mesh5.vtk")
m6 = pv.read("mesh6.vtk")
d1 = pv.read("displacement_with_von_mises1.vtk")
d2 = pv.read("displacement_with_von_mises2.vtk")
d5 = pv.read("displacement_with_von_mises5.vtk")
d6 = pv.read("displacement_with_von_mises6.vtk")
s1 = pv.read("stress1.vtk")
s2 = pv.read("stress2.vtk")
s5 = pv.read("stress5.vtk")
s6 = pv.read("stress6.vtk")
p = pv.Plotter(shape=(1, 3))

p.subplot(0, 0)
p.add_text("Mesh")
p.add_mesh(m1, color="white", show_edges=True)
p.add_mesh(m2, color="white", show_edges=True)
p.add_mesh(m5, color="white", show_edges=True)
p.add_mesh(m6, color="white", show_edges=True)
p.show_grid()

p.subplot(0, 1)
p.add_text("Von Mises Stresses")
cmap = plt.cm.get_cmap("rainbow", 20)
p.add_mesh(d1, clim=[0.0, 500.0], cmap=cmap)
p.add_mesh(d2, clim=[0.0, 500.0], cmap=cmap)
p.add_mesh(d5, clim=[0.0, 500.0], cmap=cmap)
p.add_mesh(d6, clim=[0.0, 500.0], cmap=cmap)
p.show_grid()

p.subplot(0, 2)
p.add_text("Surface traction (Y direction)")
cmap = plt.cm.get_cmap("rainbow", 20)
p.add_mesh(s1, clim=[0.0, 500.0], cmap=cmap)
p.add_mesh(s2, clim=[0.0, 500.0], cmap=cmap)
p.add_mesh(s5, clim=[0.0, 500.0], cmap=cmap)
p.add_mesh(s6, clim=[0.0, 500.0], cmap=cmap)

p.show(screenshot="contour.png", window_size=[1200, 1200], cpos="xy")


###############################################################################
# Plot the values of a dataset over a line through that dataset
#

# Make two points to construct the line between
a = [0.000, 5.000, 0.000]
b = [0.000, 0.000, 0.000]

# Run the filter and produce a line plot
sampled = d1.sample_over_line(a, b)
values = sampled.get_array("Displacements")
distance = sampled["Distance"]
plt.plot(distance, values[:, 0])
plt.show()
plt.savefig("displacements.png")

# Make two points to construct the line between

center = [0.0, 5.0, 0.0]
normal = [0.0, 0.0, 1.0]
polar = [0.0, -(5.0 - 0.001), 0.0]
angle = 180.0
sampled = s1.sample_over_circular_arc2(center, normal=normal, polar=polar)
values = sampled.get_array("Stress")
distance = sampled["Distance"]
plt.plot(distance, values)

center = [0.0, -5.0, 0.0]
normal = [0.0, 0.0, -1.0]
polar = [0.0, 5.0 - 0.001, 0.0]
angle = 180.0
sampled = s2.sample_over_circular_arc2(center, normal=normal, polar=polar)
values = sampled.get_array("Stress")
distance = sampled["Distance"]
plt.plot(distance, values)

plt.show()

a = [0.000, 10.000, 0.000]
b = [5.000, 10.000, 0.000]
sampled = s5.sample_over_line(a, b)
values = sampled.get_array("Stress")
distance = sampled["Distance"]
plt.plot(distance, values)

a = [0.000, -10.000, 0.000]
b = [5.000, -10.000, 0.000]
sampled = s6.sample_over_line(a, b)
values = sampled.get_array("Stress")
distance = sampled["Distance"]
plt.plot(distance, values)
plt.show()
plt.savefig("surface_traction.png")
