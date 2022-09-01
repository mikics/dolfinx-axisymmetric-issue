import os
import sys

import gmsh
import meshio
import numpy as np
from dolfin import (MPI, Constant, Expression,
                    FacetNormal, FiniteElement, Function, FunctionSpace,
                    Measure, Mesh, MeshValueCollection, MixedElement,
                    PETScLUSolver, SpatialCoordinate,
                    TestFunctions, TrialFunction,
                    XDMFFile, as_vector,
                    assemble, cpp, cross, dot, inner, lhs,
                    parameters, project, rhs, split, sqrt)

degree = 3


def background_field_r(m):

    B = "cos(theta)"

    if m == 0:

        Jm_prime = "(-jn(1, k0*x[0]*sin(theta)))"

    else:

        Jm_prime = "1/2*(jn(m-1, k0*x[0]*sin(theta)) - jn(m+1, k0*x[0]*sin(theta)))"

    if (m+4) % 4 == 0:

        exp_real = "(-sin(k0*x[1]*cos(theta)))"
        exp_imag = "cos(k0*x[1]*cos(theta))"

    elif (m+3) % 4 == 0:

        exp_real = "cos(k0*x[1]*cos(theta))"
        exp_imag = "sin(k0*x[1]*cos(theta))"

    elif (m+2) % 4 == 0:

        exp_real = "sin(k0*x[1]*cos(theta))"
        exp_imag = "(-cos(k0*x[1]*cos(theta)))"

    elif (m+1) % 4 == 0:

        exp_real = "(-cos(k0*x[1]*cos(theta)))"
        exp_imag = "(-sin(k0*x[1]*cos(theta)))"

    real = B + "*" + exp_real + "*" + Jm_prime
    imag = B + "*" + exp_imag + "*" + Jm_prime

    return real, imag


def background_field_z(m):

    A = "sin(theta)"

    Jm = "jn(m, k0*x[0]*sin(theta))"

    if (m+4) % 4 == 0:

        exp_real = "cos(k0*x[1]*cos(theta))"
        exp_imag = "sin(k0*x[1]*cos(theta))"

    elif (m+3) % 4 == 0:

        exp_real = "sin(k0*x[1]*cos(theta))"
        exp_imag = "(-cos(k0*x[1]*cos(theta)))"

    elif (m+2) % 4 == 0:

        exp_real = "(-cos(k0*x[1]*cos(theta)))"
        exp_imag = "(-sin(k0*x[1]*cos(theta)))"

    elif (m+1) % 4 == 0:

        exp_real = "(-sin(k0*x[1]*cos(theta)))"
        exp_imag = "cos(k0*x[1]*cos(theta))"

    real = A + "*" + exp_real + "*" + Jm
    imag = A + "*" + exp_imag + "*" + Jm

    return real, imag


def background_field_p(m):

    B = "cos(theta)"
    u = "(k0*x[0]*sin(theta))"
    Jm = "jn(m, k0*x[0]*sin(theta))"
    coeff = "m" + "*" + B
    if (m+4) % 4 == 0:

        exp_real = "cos(k0*x[1]*cos(theta))"
        exp_imag = "sin(k0*x[1]*cos(theta))"

    elif (m+3) % 4 == 0:

        exp_real = "sin(k0*x[1]*cos(theta))"
        exp_imag = "(-cos(k0*x[1]*cos(theta)))"

    elif (m+2) % 4 == 0:

        exp_real = "(-cos(k0*x[1]*cos(theta)))"
        exp_imag = "(-sin(k0*x[1]*cos(theta)))"

    elif (m+1) % 4 == 0:

        exp_real = "(-sin(k0*x[1]*cos(theta)))"
        exp_imag = "cos(k0*x[1]*cos(theta))"

    real = coeff + "*" + exp_real + "*" + Jm + "/" + u
    imag = coeff + "*" + exp_imag + "*" + Jm + "/" + u

    return real, imag


def curl_r(vector_real, vector_imag, m, rho):

    real = -vector_real[2].dx(1) + m/rho*vector_imag[1]
    imag = -vector_imag[2].dx(1) - m/rho*vector_real[1]

    return real, imag


def curl_z(vector_real, vector_imag, m, rho):

    real = +vector_real[2]/rho + vector_real[2].dx(0) - m/rho*vector_imag[0]
    imag = +vector_imag[2]/rho + vector_imag[2].dx(0) + m/rho*vector_real[0]

    return real, imag


def curl_p(vector_real, vector_imag):

    real = +vector_real[0].dx(1) - vector_real[1].dx(0)
    imag = +vector_imag[0].dx(1) - vector_imag[1].dx(0)
    return real, imag


def curl_axis(vector_real, vector_imag, m, rho):

    curlr_real, curlr_imag = curl_r(vector_real, vector_imag, m, rho)
    curlz_real, curlz_imag = curl_z(vector_real, vector_imag, m, rho)
    curlp_real, curlp_imag = curl_p(vector_real, vector_imag)

    real = as_vector((curlr_real, curlz_real, curlp_real))
    imag = as_vector((curlr_imag, curlz_imag, curlp_imag))

    return real, imag


um = 1
nm = um * 10**-3

radius_sph = 0.025*um
radius_dom = 1*um

in_sph_size = 2*nm
on_sph_size = 0.8*nm
scatt_size = 40*nm

au_tag = 1
bkg_tag = 2
scatt_tag = 3

dim = 2

gmsh.initialize(sys.argv)

gmsh.model.add("geometry")

gmsh.model.occ.addCircle(0, 0, 0, radius_sph*0.5,
                         angle1=-np.pi/2, angle2=np.pi/2, tag=1)
gmsh.model.occ.addCircle(0, 0, 0, radius_sph,
                         angle1=-np.pi/2, angle2=np.pi/2, tag=2)
gmsh.model.occ.addCircle(0, 0, 0, radius_dom,
                         angle1=-np.pi/2, angle2=np.pi/2, tag=3)

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
gmsh.write("mesh.msh")

gmsh.finalize()

if dim == 3:

    cell_str = "tetra"
    facet_str = "triangle"
    edge_str = "line"

if dim == 2:

    cell_str = "triangle"
    facet_str = "line"

if MPI.comm_world.rank == 0:

    msh = meshio.read("mesh.msh")

    cells = np.vstack(np.array([cells.data for cells in msh.cells
                                if cells.type == cell_str]))

    total_mesh = meshio.Mesh(points=msh.points,
                             cells=[(cell_str, cells)])
    total_mesh.prune_z_0()

    cells_data = msh.cell_data_dict["gmsh:physical"][cell_str]

    cells_mesh = meshio.Mesh(points=msh.points,
                             cells=[(cell_str, cells)],
                             cell_data={"name_to_read": [cells_data]})

    cells_mesh.prune_z_0()
    facet_cells = np.vstack(np.array([cells.data for cells in msh.cells
                                      if cells.type == facet_str]))

    facet_data = msh.cell_data_dict["gmsh:physical"][facet_str]

    facet_mesh = meshio.Mesh(points=msh.points,
                             cells=[(facet_str, facet_cells)],
                             cell_data={"name_to_read": [facet_data]})

    facet_mesh.prune_z_0()

    meshio.write("total.xdmf", total_mesh)
    meshio.write("facets.xdmf", facet_mesh)
    meshio.write("cells.xdmf", cells_mesh)

parameters["ghost_mode"] = "shared_vertex"

comm = MPI.comm_world
mesh = Mesh(comm)

with XDMFFile(comm, "cells.xdmf") as infile:
    infile.read(mesh)

mvc_facets = MeshValueCollection("size_t", mesh, dim-1)

with XDMFFile(comm, "facets.xdmf") as infile:
    infile.read(mvc_facets, "name_to_read")

mf_facets = cpp.mesh.MeshFunctionSizet(mesh, mvc_facets)

mvc_cells = MeshValueCollection("size_t", mesh, dim)

with XDMFFile(comm, "cells.xdmf") as infile:
    infile.read(mvc_cells, "name_to_read")

mf_cells = cpp.mesh.MeshFunctionSizet(mesh, mvc_cells)

dx = Measure("dx", domain=mesh, subdomain_data=mf_cells)
ds = Measure("ds", domain=mesh, subdomain_data=mf_facets)

dAbs = dx(au_tag)
dBkg = dx(bkg_tag)
dDom = dx(au_tag) + dx(bkg_tag)
dSbc = ds(scatt_tag)

n = FacetNormal(mesh)
n_3d = as_vector((n[0], n[1], 0))

m = 1

theta = 45 * np.pi / 180
wl0 = 0.4 * um
k0 = 2 * np.pi / wl0

reps = -1.0782
ieps = +5.8089

curl_el = FiniteElement('N1curl', mesh.ufl_cell(), degree)
lagr_el = FiniteElement('Lagrange', mesh.ufl_cell(), degree)
element = MixedElement([curl_el, curl_el, lagr_el, lagr_el])

V = FunctionSpace(mesh, element)

curl_space = FunctionSpace(mesh, "N1curl", degree)
lagr_space = FunctionSpace(mesh, "Lagrange", degree)

rEs_rz_m, iEs_rz_m, rEs_p_m, iEs_p_m = split(TrialFunction(V))
rv_rz_m, iv_rz_m, rv_p_m, iv_p_m = TestFunctions(V)

rho, z = SpatialCoordinate(mesh)
r = sqrt(rho**2 + z**2)

rEb_r_m_str, iEb_r_m_str = background_field_r(m)

rEb_z_m_str, iEb_z_m_str = background_field_z(m)

rEb_p_m_str, iEb_p_m_str = background_field_p(m)

rEb_rz_m = Expression(
    (rEb_r_m_str, rEb_z_m_str),
    m=m, k0=k0, theta=theta, degree=degree, domain=mesh)
iEb_rz_m = Expression(
    (iEb_r_m_str, iEb_z_m_str),
    m=m, k0=k0, theta=theta, degree=degree, domain=mesh)
rEb_p_m = Expression(
    (rEb_p_m_str),
    m=m, k0=k0, theta=theta, degree=degree, domain=mesh)
iEb_p_m = Expression(
    (iEb_p_m_str),
    m=m, k0=k0, theta=theta, degree=degree, domain=mesh)

rEb_m = as_vector((rEb_rz_m[0], rEb_rz_m[1], rEb_p_m))
iEb_m = as_vector((iEb_rz_m[0], iEb_rz_m[1], iEb_p_m))

rEs_m = as_vector((rEs_rz_m[0], rEs_rz_m[1], rEs_p_m))
iEs_m = as_vector((iEs_rz_m[0], iEs_rz_m[1], iEs_p_m))

rv_m = as_vector((rv_rz_m[0], rv_rz_m[1], rv_p_m))
iv_m = as_vector((iv_rz_m[0], iv_rz_m[1], iv_p_m))

rcurlEs_m, icurlEs_m = curl_axis(rEs_m, iEs_m, m, rho)
rcurlv_m, icurlv_m = curl_axis(rv_m, iv_m, m, rho)

curl_term = dot(rcurlEs_m, rcurlv_m)*rho*dDom \
    + dot(icurlEs_m, icurlv_m)*rho*dDom

eps_term_1 = -reps * k0 ** 2 * dot(
    rEs_m, rv_m) * rho * dAbs + ieps * k0 ** 2 * dot(
    iEs_m, rv_m) * rho * dAbs - reps * k0 ** 2 * dot(
    iEs_m, iv_m) * rho * dAbs - ieps * k0 ** 2 * dot(
    rEs_m, iv_m) * rho * dAbs

eps_term_2 = - k0 ** 2 * dot(
    rEs_m, rv_m) * rho * dBkg - k0 ** 2 * dot(
    iEs_m, iv_m) * rho * dBkg

field_term = (
    -(reps - 1) * k0 ** 2 * dot(rEb_m, rv_m) + ieps * k0 ** 2
    * dot(iEb_m, rv_m)) * rho * dAbs + (
    -(reps - 1) * k0 ** 2 * dot(iEb_m, iv_m) - ieps * k0 ** 2
    * dot(rEb_m, iv_m)) * rho * dAbs

sbc_term = inner(
    Constant(k0) * cross(iEs_m, n_3d),
    cross(rv_m, n_3d)) * rho * dSbc + inner(
    -Constant(k0) * cross(rEs_m, n_3d),
    cross(iv_m, n_3d)) * rho * dSbc - inner(
    1/r * cross(rEs_m, n_3d),
    cross(rv_m, n_3d)) * rho * dSbc - inner(
    1/r * cross(iEs_m, n_3d),
    cross(iv_m, n_3d)) * rho * dSbc

F = curl_term + eps_term_1 + eps_term_2 + field_term + sbc_term

a, L = lhs(F), rhs(F)
A = assemble(a)
b = assemble(L)

sol = PETScLUSolver(method="mumps")

x = Function(V)

sol.solve(A, x.vector(), b)
rEsh_rz_m, iEsh_rz_m, rEsh_p_m, iEsh_p_m = x.split()

rEsh_m = as_vector(
    (rEsh_rz_m[0],
     rEsh_rz_m[1],
     rEsh_p_m))
iEsh_m = as_vector(
    (iEsh_rz_m[0],
     iEsh_rz_m[1],
     iEsh_p_m))

rEh_rz_m = project(rEb_rz_m + rEsh_rz_m, curl_space)
iEh_rz_m = project(iEb_rz_m + iEsh_rz_m, curl_space)
rEh_p_m = project(rEb_p_m + rEsh_p_m, lagr_space)
iEh_p_m = project(iEb_p_m + iEsh_p_m, lagr_space)

rEh_m = as_vector((rEh_rz_m[0], rEh_rz_m[1], rEh_p_m))
iEh_m = as_vector((iEh_rz_m[0], iEh_rz_m[1], iEh_p_m))

Q = 2*np.pi*ieps*(dot(rEh_m,
                      rEh_m)+dot(iEh_m, iEh_m))*k0

rEb_p_m_proj = project(rEb_p_m, lagr_space)
iEb_p_m_proj = project(iEb_p_m, lagr_space)

fout = XDMFFile(MPI.comm_world, f"sols/legacy/rEs_p_m{m}_deg{degree}.xdmf")
fout.write_checkpoint(rEsh_p_m, "rEs_p_m", 0,
                      XDMFFile.Encoding.HDF5, False)
fout.close()

fout = XDMFFile(MPI.comm_world, f"sols/legacy/iEs_p_m{m}_deg{degree}.xdmf")
fout.write_checkpoint(iEsh_p_m, "iEs_p_m", 0,
                      XDMFFile.Encoding.HDF5, False)
fout.close()

q_abs_fenics = assemble(Q*rho*dAbs)
q_abs_analyt = 0.0004908126796809729

err_abs = np.abs(q_abs_analyt - q_abs_fenics) / q_abs_analyt

print()
print(f"The analytical absorption efficiency is {q_abs_analyt}")
print(f"The numerical absorption efficiency is {q_abs_fenics}")
print(f"The error is {err_abs*100}%")
print()

assert err_abs < 1e-2
