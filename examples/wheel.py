"""
.. _wheel_example:

"""
import getfem as gf
import numpy as np

import pyvista as pv

###############################################################################
#

E = 21e6  # Young Modulus (N/cm^2)
nu = 0.3  # Poisson ratio
clambda = E * nu / ((1 + nu) * (1 - 2 * nu))  # First Lame coefficient (N/cm^2)
cmu = E / (2 * (1 + nu))  # Second Lame coefficient (N/cm^2)
clambdastar = 2 * clambda * cmu / (clambda + 2 * cmu)  # Lame coefficient for Plane stress (N/cm^2)
applied_force = 1e7  # Force at the hole boundary (N)

h = 1  # Approximate mesh size
elements_degree = 2  # Degree of the finite element methods
gamma0 = 1.0 / E
# Augmentation parameter for the augmented Lagrangian


###############################################################################
#


mo1 = gf.MesherObject("ball", [0.0, 15.0], 15.0)
mo2 = gf.MesherObject("ball", [0.0, 15.0], 8.0)
mo3 = gf.MesherObject("set minus", mo1, mo2)
gf.util("trace level", 2)  # No trace for mesh generation
mesh1 = gf.Mesh("generate", mo3, h, 2)
mesh1.export_to_vtk("mesh1.vtk")

mo4 = gf.MesherObject("ball", [0.0, -15.0], 15.0)
mo5 = gf.MesherObject("ball", [0.0, -15.0], 8.0)
mo6 = gf.MesherObject("set minus", mo4, mo5)
gf.util("trace level", 2)  # No trace for mesh generation
mesh2 = gf.Mesh("generate", mo6, h, 2)
mesh2.export_to_vtk("mesh2.vtk")


###############################################################################
#



m1 = pv.read("mesh1.vtk")
m2 = pv.read("mesh2.vtk")
p = pv.Plotter(shape=(1, 1))
p.subplot(0, 0)
p.add_mesh(m1, show_edges=True)
p.add_text("Neumann\ncondition", position=[525, 650])
p.add_mesh(m2, show_edges=True)
p.add_text("Dirichlet\ncondition", position=[525, 300])
p.show_grid()
p.show(screenshot="mesh.png", window_size=[1200, 1400], cpos="xy")

###############################################################################
#

fb1 = mesh1.outer_faces_in_box([-8.1, 6.9], [8.1, 23.1])  # Boundary of the hole
fb2 = mesh1.outer_faces_with_direction([0.0, -1.0], np.pi / 4.5)  # Contact boundary of the wheel
fb3 = mesh2.outer_faces_in_box([-8.1, -6.9], [8.1, -23.1])  # Boundary of the hole
fb4 = mesh2.outer_faces_with_direction([0.0, 1.0], np.pi / 4.5)  # Contact boundary of the wheel

HOLE_BOUND = 1
CONTACT_BOUND = 2
BOTTOM_BOUND = 3

mesh1.set_region(HOLE_BOUND, fb1)
mesh1.set_region(CONTACT_BOUND, fb2)
mesh1.region_subtract(CONTACT_BOUND, HOLE_BOUND)
mesh2.set_region(HOLE_BOUND, fb3)
mesh2.set_region(CONTACT_BOUND, fb4)
mesh2.region_subtract(CONTACT_BOUND, HOLE_BOUND)


###############################################################################
#

mfu1 = gf.MeshFem(mesh1, 2)
mfu1.set_classical_fem(elements_degree)
mflambda = gf.MeshFem(mesh1, 2)
mflambda.set_classical_fem(elements_degree - 1)
mflambda_C = gf.MeshFem(mesh1, 1)
mflambda_C.set_classical_fem(elements_degree - 1)
mfu2 = gf.MeshFem(mesh2, 2)
mfu2.set_classical_fem(elements_degree)
mfvm1 = gf.MeshFem(mesh1, 1)
mfvm1.set_classical_discontinuous_fem(elements_degree)
mfvm2 = gf.MeshFem(mesh2, 1)
mfvm2.set_classical_discontinuous_fem(elements_degree)
mim1 = gf.MeshIm(mesh1, pow(elements_degree, 2))
mim2 = gf.MeshIm(mesh2, pow(elements_degree, 2))


###############################################################################
#

md = gf.Model("real")
md.add_fem_variable("u1", mfu1)  # Displacement of the structure 1
md.add_fem_variable("u2", mfu2)  # Displacement of the structure 2

###############################################################################
#

md.add_initialized_data("cmu", [cmu])
md.add_initialized_data("clambdastar", [clambdastar])
md.add_isotropic_linearized_elasticity_brick(mim1, "u1", "clambdastar", "cmu")
md.add_isotropic_linearized_elasticity_brick(mim2, "u2", "clambdastar", "cmu")

###############################################################################
#

md.add_Dirichlet_condition_with_multipliers(mim2, "u2", elements_degree - 1, HOLE_BOUND)

###############################################################################
#

md.add_interpolate_transformation_from_expression("Proj1", mesh1, mesh2, "[X(1);-X(2)]")

###############################################################################
#

md.add_initialized_data("gamma0", [gamma0])
md.add_filtered_fem_variable("lambda1", mflambda_C, CONTACT_BOUND)
md.add_nonlinear_term(mim1, "lambda1*(Test_u1.[0;1])-lambda1*(Interpolate(Test_u2,Proj1).[0;1])", CONTACT_BOUND)
md.add_nonlinear_term(mim1, "-(gamma0*element_size)*(lambda1 + neg_part(lambda1+(1/(gamma0*element_size))" "*((u1-Interpolate(u2,Proj1)+X-Interpolate(X,Proj1)).[0;1])))*Test_lambda1", CONTACT_BOUND)

###############################################################################
#

md.add_variable("alpha_D", 1)

###############################################################################
#

md.add_filtered_fem_variable("lambda_D", mflambda, HOLE_BOUND)

###############################################################################
#

md.add_initialized_data("F", [applied_force / (8 * 2 * np.pi)])
md.add_linear_term(mim1, "-lambda_D.Test_u1 + (alpha_D*[0;1]-u1).Test_lambda_D" " + (lambda_D.[0;1]+F)*Test_alpha_D", HOLE_BOUND)

###############################################################################
#

md.add_linear_term(mim1, "1E-6*alpha_D*Test_alpha_D")

###############################################################################
#

md.solve("max_res", 1e-9, "max_iter", 100, "noisy")

###############################################################################
#

md.solve("max_res", 1e-9, "max_iter", 100, "noisy", "lsearch", "simplest", "alpha min", 0.8)

###############################################################################
#

U1 = md.variable("u1")
U2 = md.variable("u2")
VM1 = md.compute_isotropic_linearized_Von_Mises_or_Tresca("u1", "clambdastar", "cmu", mfvm1)
VM2 = md.compute_isotropic_linearized_Von_Mises_or_Tresca("u2", "clambdastar", "cmu", mfvm2)

mfvm1.export_to_vtk("displacement_with_von_mises1.vtk", mfvm1, VM1, "Von Mises Stresses", mfu1, U1, "Displacements")

mfvm2.export_to_vtk("displacement_with_von_mises2.vtk", mfvm2, VM2, "Von Mises Stresses", mfu2, U2, "Displacements")

###############################################################################
# You can view solutions with pyvista:
#

m1 = pv.read("displacement_with_von_mises1.vtk")
m2 = pv.read("displacement_with_von_mises2.vtk")
c1 = m1.contour()
c2 = m2.contour()
p = pv.Plotter(shape=(1, 1))

p.subplot(0, 0)
p.add_text("Von Mises Stresses")
p.add_mesh(c1, color="black", line_width=2)
p.add_mesh(m1, show_edges=False, opacity=0.85)
p.add_mesh(c2, color="black", line_width=2)
p.add_mesh(m2, show_edges=False, opacity=0.85)
p.show_grid()

p.show(screenshot="von_mises.png", window_size=[1200, 900], cpos="xy")
