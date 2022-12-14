# DOLFINx issue with axisymmetric electromagnetic problems

This repository has been made for showing a DOLFINx issue when solving
axisymmetric Maxwell's equations. In particular, the problem seems to be related
to the `degree` of the finite elements. When using `degree = 2`, the solution
calculated by DOLFINx looks correct, as confirmed by a calculation of the
absorption efficiency, which is pretty close to the analytical value. However,
for `degree = 3` the solution presents some artifacts, which results in
an absorption efficiency not close to its analytical value.  

The image here below compares the solution in legacy DOLFIN with
the solutions in DOLFINx. The DOLFIN version of the code have been
extensively tested against analytical results, and its output it is therefore
reliable, and does not present any inconsistency using different `degree` for
the finite elements.

![image](comparison.png)

## Minimal working examples

In order to understand the issue, I have added the file `dolfinx_mwe.py`, which
is a minimal working example showing the issue when passing from `degree = 2` to
`degree = 3`. The file `legacy_mwe.py` contains the same problem solved in
legacy DOLFIN for a comparison.

## Problem

The problem solved  by `dolfinx_mwe.py` and `legacy_mwe.py`
is the electromagnetic scattering of a plane wave
by a metallic sphere. Due to the symmetry of the problem,
we can exploit the expansion of the field in cylindrical harmonics to simplify
the 3D problem into few 2D problems, corresponding to the multiple cylindrical
harmonics that propagate independently. For this reason,
the reference coordinate system is the cylindrical one, $(\rho, z, \phi)$, and the
scattered electric field is expressed in terms of its components $(E_{\rho}, E_z, E_{\phi})$.
In the examples, $E_\rho$ and $E_z$ are discretized with the first-kind Nedelec elements,
while $E_\phi$ is discretized with Lagrange elements.
The minimal working examples only solve the problem for the harmonic number `m = 1`,
since for `m = 0` DOLFINx and legacy DOLFIN seem to agree, whatever the `degree`.
Scattering boundary conditions are used for making the boundary transparent to
outgoing waves.

Here below a quick rundown of the DOLFIN/DOLFINx functions and their corresponding
mathematical representation. The implementation of such functions in legacy DOLFIN is
necessarily more complicated due to the lack of support of complex numbers.


- `curl_r` $\rightarrow\left(-\frac{\partial a_{\phi}}{\partial z}-i \frac{m}{\rho}a_{z}\right)$
- `curl_z` $\rightarrow\left(\frac{a_{\phi}}{\rho}+\frac{\partial a_{\phi}}{\partial \rho}+i \frac{m}{\rho} a_{\rho}\right)$
- `curl_p` $\rightarrow\left(\frac{\partial a_{\rho}}{\partial z}-\frac{\partial a_{z}}{\partial \rho}\right)$
- `background_field_r` $\rightarrow\left(\cos \theta e^{i k_0 z \cos \theta} i^{-m+1} J_{m}^{\prime}\left(k_{0} \rho \sin\theta\right)\right)$
- `background_field_z` $\rightarrow\left( \sin \theta e^{i k_0 z \cos \theta}i^{-m} J_{m}\left(k_0 \rho \sin \theta\right)\right)$
- `background_field_p` $\rightarrow\left( \frac{\cos \theta}{k_0 \rho \sin \theta}e^{i k_0 z \cos \theta} i^{-m} J_{m}\left(k_0 \rho \sin \theta\right)\right)$
- $J_m(x)\rightarrow m$-th order Bessel function
- `curl_term` $\rightarrow \int_{\Omega_{dom}}-(\nabla \times \mathbf{E}^{(m)}_s)\cdot (\nabla \times \bar{\mathbf{v}}^{(m)})~dx$
- `eps_term_1` $\rightarrow \int_{\Omega_{abs}}-\varepsilon_r k_0^{2}\mathbf{E}^{(m)}_s \cdot \bar{\mathbf{v}}^{(m)} ~dx$ 
- `eps_term_2` $\rightarrow \int_{\Omega_{bkg}}-k_0^{2}\mathbf{E}^{(m)}_s \cdot \bar{\mathbf{v}}^{(m)} ~dx$
- `field_term` $\rightarrow \int_{\Omega_{abs}}-k_0^{2}\left(\varepsilon_r - 1\right)\mathbf{E}^{(m)}_b \cdot \bar{\mathbf{v}}^{(m)}  ~dx$
- `sbc_term` $\rightarrow \int_{\partial\Omega_{sbc}}-(jk_0 + \frac{1}{r})(\mathbf{E}^{(m)}_s\times\mathbf{n}) \cdot (\bar{\mathbf{v}}^{(m)}\times\mathbf{n}) ~ds$
