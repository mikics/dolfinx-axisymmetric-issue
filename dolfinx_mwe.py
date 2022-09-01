import sys
from functools import partial

import numpy as np
from scipy.special import jv, jvp

from dolfinx import fem
from dolfinx.io import VTXWriter
from dolfinx.io.gmshio import model_to_mesh
from ufl import (FacetNormal, FiniteElement, Measure, MixedElement,
                 SpatialCoordinate, TestFunction, TrialFunction,
                 as_vector, cross, inner, lhs, rhs, sqrt
                 )
import gmsh
from mpi4py import MPI
from petsc4py import PETSc

degree = 2


def background_field_rz(theta, k0, m, x):

    rho, z = x[0], x[1]

    a_r = (np.cos(theta) * np.exp(1j * k0 * z * np.cos(theta))
           * (1j)**(-m + 1) * jvp(m, k0 * rho * np.sin(theta), 1))

    a_z = (np.sin(theta) * np.exp(1j * k0 * z * np.cos(theta))
           * (1j)**-m * jv(m, k0 * rho * np.sin(theta)))

    return (a_r, a_z)


def background_field_p(theta, k0, m, x):

    rho, z = x[0], x[1]

    a_p = (np.cos(theta) / (k0 * rho * np.sin(theta))
           * np.exp(1j * k0 * z * np.cos(theta)) * m
           * (1j)**(-m) * jv(m, k0 * rho * np.sin(theta)))

    return a_p


def curl_axis(a, m, rho):

    curl_r = -a[2].dx(1) - 1j * m / rho * a[1]
    curl_z = a[2] / rho + a[2].dx(0) + 1j * m / rho * a[0]
    curl_p = a[0].dx(1) - a[1].dx(0)

    return as_vector((curl_r, curl_z, curl_p))


um = 1
nm = um * 10**-3

radius_sph = 0.025 * um
radius_dom = 1 * um

in_sph_size = 2 * nm
on_sph_size = 0.8 * nm
scatt_size = 40 * nm

au_tag = 1
bkg_tag = 2
scatt_tag = 3

gmsh.initialize(sys.argv)

gmsh.model.add("geometry")

gmsh.model.occ.addCircle(0, 0, 0, radius_sph * 0.5,
                         angle1=-np.pi / 2, angle2=np.pi / 2, tag=1)
gmsh.model.occ.addCircle(0, 0, 0, radius_sph,
                         angle1=-np.pi / 2, angle2=np.pi / 2, tag=2)
gmsh.model.occ.addCircle(0, 0, 0, radius_dom,
                         angle1=-np.pi / 2, angle2=np.pi / 2, tag=3)

gmsh.model.occ.addLine(6, 4, tag=4)
gmsh.model.occ.addLine(4, 2, tag=5)
gmsh.model.occ.addLine(2, 1, tag=6)
gmsh.model.occ.addLine(1, 3, tag=7)
gmsh.model.occ.addLine(3, 5, tag=8)

gmsh.model.occ.addCurveLoop([6, 1], tag=1)
gmsh.model.occ.addPlaneSurface([1], tag=1)
gmsh.model.occ.addCurveLoop([7, 2, 5, -1], tag=2)
gmsh.model.occ.addPlaneSurface([2], tag=2)
gmsh.model.occ.addCurveLoop([4, -3, 8, 2], tag=3)
gmsh.model.occ.addPlaneSurface([3], tag=3)

gmsh.model.occ.synchronize()

gmsh.model.addPhysicalGroup(2, [1, 2], tag=au_tag)
gmsh.model.addPhysicalGroup(2, [3], tag=bkg_tag)
gmsh.model.addPhysicalGroup(1, [3], tag=scatt_tag)

gmsh.model.mesh.setSize([(0, 1)], size=in_sph_size)
gmsh.model.mesh.setSize([(0, 2)], size=in_sph_size)
gmsh.model.mesh.setSize([(0, 3)], size=on_sph_size)
gmsh.model.mesh.setSize([(0, 4)], size=on_sph_size)
gmsh.model.mesh.setSize([(0, 5)], size=scatt_size)
gmsh.model.mesh.setSize([(0, 6)], size=scatt_size)

gmsh.model.mesh.generate(2)

model = gmsh.model

domain, cell_tags, facet_tags = model_to_mesh(
    model, MPI.COMM_WORLD, 0, gdim=2)

gmsh.finalize()

dx = Measure("dx", domain, subdomain_data=cell_tags,
             metadata={'quadrature_degree': 40})

ds = Measure("ds", domain, subdomain_data=facet_tags)

dAbs = dx(au_tag)
dBkg = dx(bkg_tag)
dDom = dx((au_tag, bkg_tag))
dSbc = ds(scatt_tag)

n = FacetNormal(domain)
n_3d = as_vector((n[0], n[1], 0))

m = 1

theta = 45 * np.pi / 180
wl0 = 0.4 * um
k0 = 2 * np.pi / wl0

eps_au = -1.0782 + 1j * 5.8089

curl_el = FiniteElement("N1curl", domain.ufl_cell(), degree)
lagr_el = FiniteElement("Lagrange", domain.ufl_cell(), degree)
V = fem.FunctionSpace(domain, MixedElement([curl_el, lagr_el]))

rho, z = SpatialCoordinate(domain)
r = sqrt(rho**2 + z**2)

Eh_m = fem.Function(V)

Es_m = TrialFunction(V)
v_m = TestFunction(V)

Eb_m = fem.Function(V)
f_rz = partial(background_field_rz, theta, k0, m)
f_p = partial(background_field_p, theta, k0, m)
Eb_m.sub(0).interpolate(f_rz)
Eb_m.sub(1).interpolate(f_p)

curl_Es_m = curl_axis(Es_m, m, rho)
curl_v_m = curl_axis(v_m, m, rho)

curl_term = inner(curl_Es_m, curl_v_m) * rho * dDom
eps_term_1 = - eps_au * k0 ** 2 * inner(Es_m, v_m) * rho * dAbs
eps_term_2 = - k0 ** 2 * inner(Es_m, v_m) * rho * dBkg
field_term = - k0 ** 2 * (eps_au - 1) * inner(Eb_m, v_m) * rho * dAbs
sbc_term = - (1j * k0 + 1 / r) \
    * inner(cross(Es_m, n_3d), cross(v_m, n_3d)) * rho * dSbc

F = curl_term + eps_term_1 + eps_term_2 + field_term + sbc_term

a, L = lhs(F), rhs(F)

problem = fem.petsc.LinearProblem(a, L, bcs=[], petsc_options={
                                  "ksp_type": "preonly", "pc_type": "lu"})
Esh_m = problem.solve()

Esh_rz_m, Esh_p_m = Esh_m.split()

Eh_m.x.array[:] = Eb_m.x.array[:] + Esh_m.x.array[:]

Q = 2 * np.pi * eps_au.imag * k0 * (inner(Eh_m, Eh_m))

q_abs_fenics_proc = (fem.assemble_scalar(
    fem.form(Q * rho * dAbs))).real
q_abs_fenics = domain.comm.allreduce(q_abs_fenics_proc, op=MPI.SUM)

with VTXWriter(domain.comm, f"sols/dolfinx/Es_p_m{m}_deg{degree}.bp", Esh_p_m) as f:
    f.write(0.0)

q_abs_analyt = 0.0004908126796809729

err_abs = np.abs(q_abs_analyt - q_abs_fenics) / q_abs_analyt

print()
print(f"The analytical absorption efficiency is {q_abs_analyt}")
print(f"The numerical absorption efficiency is {q_abs_fenics}")
print(f"The error is {err_abs*100}%")
print()

assert err_abs < 1e-2
