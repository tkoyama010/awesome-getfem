"""
.. _ref_plane_stress_concentration:

2D Plane Stress Concentration Analysis
--------------------------------------

This tutorial shows how you can use GetFEM to determine and
verify the "stress concentration factor" when modeling using 2D plane
elements and then verify this using 3D elements.

First, start GetFEM as a service.
"""
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np

import getfem as gf
import pyvista as pv

pv.set_plot_theme("document")


###############################################################################
# Element Type and Material Properties
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This example will use plain stress elements as a thin plate can be
# modeled with 3d plane elements and a thickness is provided.
#
# This example will use SI units.

epsilon = 0.001  # thickness of plate (m)

E = 210e9  # Elastic moduli in Pa (kg/(m*s**2))
nu = 0.3  # Poisson's Ratio
F = 100000000.0  # Force density at the right boundary(Pa/m3)

elements_degree = 2  # Degree of the finite element methods

###############################################################################
# Geometry
# ~~~~~~~~
# Create a rectangular area with the hole in the middle.  To correctly
# approximate an infinite plate, the maximum stress must occur far
# away from the edges of the plate.  A length to width factor can
# approximate this.

length = 0.4
width = 0.1

ratio = 0.3  # diameter/width
diameter = width * ratio
radius = diameter * 0.5


# create the rectangle
rect_anum = gf.MesherObject("rectangle", [0.0, 0.0], [length, width / 2])

# create a circle in the middle of the rectangle
circ_anum = gf.MesherObject("ball", [length / 2, width / 2], radius)

# Note how GetFEM parses the output and returns the area numbers
# created by each command.  This can be used to execute a boolean
# operation on these areas to cut the circle out of the rectangle.
plate_with_hole_anum = gf.MesherObject("set minus", rect_anum, circ_anum)

###############################################################################
# Meshing
# ~~~~~~~
# Mesh the plate using an approximate mesh size.
# ensure there are at 50 elements around the hole
hole_esize = np.pi * diameter / 50  # 0.0002
plate_esize = 0.01

mesh = gf.Mesh("generate", plate_with_hole_anum, hole_esize)
mesh.export_to_vtk("mesh.vtk")

m = pv.read("mesh.vtk")
p = pv.Plotter(shape=(1, 1))
p.subplot(0, 0)
p.add_mesh(m)
p.show_grid()
p.show(screenshot="mesh.png", cpos="xy")

###############################################################################
# Boundary Conditions
# ~~~~~~~~~~~~~~~~~~~
# Fix the left-hand side of the plate in the X direction and set a
# force of 1 kN in the positive X direction.
#

mfu = gf.MeshFem(mesh, 2)
mfu.set_classical_fem(elements_degree)

mfd = gf.MeshFem(mesh, 1)
mfd.set_classical_fem(elements_degree)

md = gf.Model("real")
md.add_fem_variable("u", mfu)

mim = gf.MeshIm(mesh, elements_degree * 2)

HOLE_BOUND = 1
fb1 = mesh.outer_faces_in_box(
    [length / 2 - radius * 1.01, width / 2 - radius * 1.01],
    [length / 2 + radius * 1.01, width / 2 + radius * 1.01],
)
mesh.set_region(HOLE_BOUND, fb1)

# Fix the left-hand side.
LEFT_BOUND = 2
fb2 = mesh.outer_faces_with_direction([-1.0, 0.0], 0.01)

mesh.set_region(LEFT_BOUND, fb2)
mesh.region_subtract(LEFT_BOUND, HOLE_BOUND)

md.add_initialized_data("r2", [0, 0])
md.add_initialized_data("H2", [[1, 0], [0, 0]])
md.add_generalized_Dirichlet_condition_with_multipliers(
    mim, "u", mfu, LEFT_BOUND, "r2", "H2"
)

# Fix nodes on the top side of the plate in the Y
# direction.  Otherwise, the mesh would be allowed to move in the y
# direction and would be an improperly constrained mesh.
TOP_BOUND = 3
fb3 = mesh.outer_faces_with_direction([0.0, 1.0], 0.01)

mesh.set_region(TOP_BOUND, fb3)
mesh.region_subtract(TOP_BOUND, HOLE_BOUND)

md.add_initialized_data("r3", [0, 0])
md.add_initialized_data("H3", [[0, 0], [0, 1]])
md.add_generalized_Dirichlet_condition_with_multipliers(
    mim, "u", mfu, TOP_BOUND, "r3", "H3"
)

# Apply a force on the right-hand side of the plate.  For this
# example, we select the nodes at the right-most side of the plate.
RIGHT_BOUND = 4
fb4 = mesh.outer_faces_with_direction([1.0, 0.0], 0.01)

mesh.set_region(RIGHT_BOUND, fb4)
mesh.region_subtract(RIGHT_BOUND, HOLE_BOUND)

# TODO: Verify that only the nodes at length have been selected:

# Next, couple the DOF for these nodes.  This lets us provide a force
# to one node that will be spread throughout all nodes in this coupled
# set.
# Select a single node in this set and apply a force to it
# We use "R" to re-select from the current node group
md.add_initialized_data("width", [width])
md.add_initialized_data("plate_esize", [plate_esize])
md.add_initialized_data("F", [F])
md.add_source_term_brick(mim, "u", "[F/width*plate_esize, 0]", RIGHT_BOUND)

# TODO: finally, be sure to select all nodes again to solve the entire solution


###############################################################################
# Solve the Static Problem
# ~~~~~~~~~~~~~~~~~~~~~~~~
# define a plain stress element type with thickness
md.add_initialized_data("E", [E])
md.add_initialized_data("nu", [nu])
md.add_isotropic_linearized_elasticity_brick_pstress(mim, "u", "E", "nu")
# Solve the static analysis
md.solve("max_res", 1e-9, "max_iter", 100, "noisy")

###############################################################################
# Post-Processing
# ~~~~~~~~~~~~~~~
# The static result can be post-processed outside of GetFEM using ``pyvista``.
# This example shows how to extract the von Mises stress and plot it using the
# ``pyvista`` result reader.

# grab the result from the ``getfem`` instance

# mfvm = gf.MeshFem(mesh, 1)
# mfvm.set_classical_discontinuous_fem(elements_degree)
# von_mises = md.compute_isotropic_linearized_Von_Mises_pstress("u", "E", "nu", mfvm)
# mfvm.export_to_vtk("von_mises.vtk", mfvm, von_mises, "Von Mises Stresses")

U = md.variable("u")
mfu.export_to_vtk("displacement.vtk", mfu, U, "Displacements")

Grad_u = gf.compute_gradient(mfu, U, mfd)
clambda = E * nu / ((1 + nu) * (1 - 2 * nu))
cmu = E / (2 * (1 + nu))
clambdastar = 2 * clambda * cmu / (clambda + 2 * cmu)
sigmaxx = clambdastar * (Grad_u[0, 0] + Grad_u[1, 1]) + 2.0 * cmu * Grad_u[0, 0]
sigmayy = clambdastar * (Grad_u[0, 0] + Grad_u[1, 1]) + 2.0 * cmu * Grad_u[1, 1]
sigmazz = 0.0
sigmaxy = cmu * (Grad_u[0, 1] + Grad_u[1, 0])
sigmayx = cmu * (Grad_u[1, 0] + Grad_u[0, 1])
sigmayz = 0.0
sigmazx = 0.0

von_mises = np.sqrt(
    0.5
    * (
        (sigmaxx - sigmayy) ** 2
        + (sigmayy - sigmazz) ** 2
        + (sigmazz - sigmaxx) ** 2
        + 6.0 * (sigmaxy ** 2 + sigmayz ** 2 + sigmazx ** 2)
    )
)

# Must use nanmax as stress is not computed at mid-side nodes
max_stress = np.nanmax(von_mises)

v = pv.read("von_mises.vtk")
p = pv.Plotter(shape=(1, 1))
p.subplot(0, 0)
cmap = plt.cm.get_cmap("rainbow", 10)
p.add_mesh(v, cmap=cmap)
p.show_grid()
p.show(screenshot="von_mises.png", cpos="xy")

mfd.export_to_vtk("sigmaxx.vtk", mfd, sigmaxx, "Sigmaxx")
mfd.export_to_vtk("sigmayy.vtk", mfd, sigmayy, "Sigmayy")
mfd.export_to_vtk("sigmaxy.vtk", mfd, sigmaxy, "Sigmaxy")
mfd.export_to_vtk("sigmayx.vtk", mfd, sigmayx, "Sigmayx")
mfd.export_to_vtk("von_mises.vtk", mfd, von_mises, "Von_Mises")

s1 = pv.read("sigmaxx.vtk")
s2 = pv.read("sigmayy.vtk")
s3 = pv.read("sigmaxy.vtk")
s4 = pv.read("sigmayx.vtk")
s5 = pv.read("von_mises.vtk")

a = [0.0, width / 2, 0.0]
b = [length, width / 2, 0.0]

fig = plt.figure()

ax = fig.add_subplot(511)
ax.set_ylabel("Sigmaxx")
sampled = s1.sample_over_line(a, b)
values = sampled.get_array("Sigmaxx")
position = sampled.points[:, 0]
ax.set_ylim([-20000000.0, 20000000.0])
ax.plot(position, values)

ax = fig.add_subplot(512)
ax.set_ylabel("Sigmayy")
sampled = s2.sample_over_line(a, b)
values = sampled.get_array("Sigmayy")
position = sampled.points[:, 0]
ax.set_ylim([-20000000.0, 20000000.0])
ax.plot(position, values)

ax = fig.add_subplot(513)
ax.set_ylabel("Sigmaxy")
sampled = s3.sample_over_line(a, b)
values = sampled.get_array("Sigmaxy")
position = sampled.points[:, 0]
ax.set_ylim([-20000000.0, 20000000.0])
ax.plot(position, values)

ax = fig.add_subplot(514)
ax.set_ylabel("Sigmayx")
sampled = s4.sample_over_line(a, b)
values = sampled.get_array("Sigmayx")
position = sampled.points[:, 0]
ax.set_ylim([-20000000.0, 20000000.0])
ax.plot(position, values)

ax = fig.add_subplot(515)
ax.set_ylabel("Von Mises")
sampled = s5.sample_over_line(a, b)
values = sampled.get_array("Von_Mises")
position = sampled.points[:, 0]
ax.set_ylim([-20000000.0, 20000000.0])
ax.plot(position, values)

plt.show()

###############################################################################
# Compute the Stress Concentration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The stress concentration :math:`K_t` is the ratio of the maximum
# stress at the hole to the far-field stress, or the mean cross
# sectional stress at a point far from the hole.  Analytically, this
# can be computed with:
#
# :math:`\sigma_{nom} = \frac{F}{wt}`
#
# Where
#
# - :math:`F` is the force
# - :math:`w` is the width of the plate
# - :math:`t` is the thickness of the plate.
#
# Experimentally, this is computed by taking the mean of the nodes at
# the right-most side of the plate.

# We use nanmean here because mid-side nodes have no stress
mask = mfd.basic_dof_nodes()[0, :] == length
far_field_stress = np.nanmean(von_mises[mask])
print("Far field von mises stress: %e" % far_field_stress)
# Which almost exactly equals the analytical value of 10000000.0 Pa

###############################################################################
# Since the expected nominal stress across the cross section of the
# hole will increase as the size of the hole increases, regardless of
# the stress concentration, the stress must be adjusted to arrive at
# the correct stress.  This stress is adjusted by the ratio of the
# width over the modified cross section width.
adj = width / (width - diameter)
stress_adj = far_field_stress * adj

# The stress concentration is then simply the maximum stress divided
# by the adjusted far-field stress.
stress_con = max_stress / stress_adj
print("Stress Concentration: %.2f" % stress_con)


###############################################################################
# Batch Analysis
# ~~~~~~~~~~~~~~


###############################################################################
# Run the batch and record the stress concentration


###############################################################################
# Analytical Comparison
# ~~~~~~~~~~~~~~~~~~~~~
# Stress concentrations are often obtained by referencing tablular
# results or polynominal fits for a variety of geometries.  According
# to Peterson's Stress Concentration Factors (ISBN 0470048247), the analytical
# equation for a hole in a thin plate in uniaxial tension:
#
# :math:`k_t = 3 - 3.14\frac{d}{h} + 3.667\left(\frac{d}{h}\right)^2 - 1.527\left(\frac{d}{h}\right)^3`
#
# Where:
#
# - :math:`k_t` is the stress concentration
# - :math:`d` is the diameter of the circle
# - :math:`h` is the height of the plate
#
# As shown in the following plot, GetFEM matches the known tabular
# result for this geometry remarkably well using plain stress elements.
# The fit to the results may vary depending on the ratio between the
# height and width of the plate.

# where ratio is (d/h)
