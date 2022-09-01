# DOLFINx issue with axisymmetric electromagnetic problems

This repository was made for showing a DOLFINx issue when solving
axisymmetric Maxwell's equations. In particular, the problem arises when
the `degree` of the finite elements in DOLFINx is increased from `degree = 2`
to `degree=3`. The image here below compares the solution in legacy DOLFIN with
the solutions in DOLFINx, and it clearly shows that for `degree=3` there are
some artifacts in the DOLFINx solution, which do not appear in legacy DOLFIN
and in DOLFINx for `degree=2`. The DOLFIN version of the code have been
extensively tested against analytical results, and therefore its output is
reliable.

## Minimal working examples

The file `dolfinx_mwe.py` contains a minimal working example showing the issue.
The file `legacy_mwe.py` contains the same problem solved in legacy DOLFIN for
a comparison.

## Problem

The problem being solved is the calculation of the scattered electric field
from a metallic sphere hit by a plane wave. Due to the symmetry of the problem,
we can exploit the expansion of the field in cylindrical harmonics to simplify
the 3D problems into few 2D problems, corresponding to the multiple cylindrical
harmonics propagating independently.

$$
\begin{align}
\sum_{m}\int_{\Omega_{cs}}&-(\nabla \times \mathbf{E}^{(m)}_s)
\cdot (\nabla \times \bar{\mathbf{v}}^{(m)})+\varepsilon_{r} k_{0}^{2}
\mathbf{E}^{(m)}_s \cdot \bar{\mathbf{v}}^{(m)}
+k_{0}^{2}\left(\varepsilon_{r}
-\varepsilon_b\right)\mathbf{E}^{(m)}_b \cdot \bar{\mathbf{v}}^{(m)}\\
&+\left(\boldsymbol{\mu}^{-1}_{pml} \nabla \times \mathbf{E}^{(m)}_s
\right)\cdot \nabla \times \bar{\mathbf{v}}^{(m)}-k_{0}^{2}
\left(\boldsymbol{\varepsilon}_{pml} \mathbf{E}^{(m)}_s \right)\cdot
\bar{\mathbf{v}}^{(m)}~ \rho d\rho dz =0
\end{align}
$$

$$
\begin{align}
\mathbf{E}^{(m)}_b = &\hat{\rho} \left(E_{0} \cos \theta
e^{i k z \cos \theta} i^{-m+1} J_{m}^{\prime}\left(k_{0} \rho \sin
\theta\right)\right)\\
+&\hat{z} \left(E_{0} \sin \theta e^{i k z \cos \theta}i^{-m} J_{m}
\left(k \rho \sin \theta\right)\right)\\
+&\hat{\phi} \left(\frac{E_{0} \cos \theta}{k \rho \sin \theta}
e^{i k z \cos \theta} i^{-m} J_{m}\left(k \rho \sin \theta\right)\right)
\end{align}
$$

$$
\begin{align}
\left(\nabla \times \mathbf{a}^{(m)}\right) = &\left[\hat{\rho}
\left(-\frac{\partial a_{\phi}^{(m)}}{\partial z}
-i \frac{m}{\rho} a_{z}^{(m)}\right)+\\ \hat{\phi}
\left(\frac{\partial a_{\rho}^{(m)}}{\partial z}
-\frac{\partial a_{z}^{(m)}}{\partial \rho}\right)+\right.\\
&\left.+\hat{z}\left(\frac{a_{\phi}^{(m)}}{\rho}
+\frac{\partial a_{\phi}^{(m)}}{\partial \rho}
+i \frac{m}{\rho} a_{\rho}^{(m)}\right)\right]
\end{align}
$$
