# RKPM (Reproducing kernel particle methods)

This code implements an RKPM-based interpolation extension method proposed for immersed boundary-type frameworks. The method proposed for immersed boundary-type frame works. This method based on the conservation properties inherrent to RKPM, enforcing that function evaluations are performed on Lagrangian markers distributed over the embedded geometry tather than on the surrounding Eulerian grid points. It is compatible with the curently mainstream spatial discretization solver.

Main Reference:

[1] Pinelli, A., Naqavi, I., Favier, J., & Pasquetti, R. (2010). Immersed-boundary methods for general finite-difference and finite-volume Navier–Stokes solvers. Journal of Computational Physics, 229(24), 9073-9091. https://doi.org/10.1016/j.jcp.2010.08.021

[2] W.K. Liu, S. Jun, Y.F. Zhang, Reproducing kernel particle methods, Int. J. Numer. Methods Fluids 20 (8) (1995) 1081–1106.

[3] A. M. Roma, C. S. Peskin, and M. J. Berger, "An adaptive version of the immersed boundary method," Journal of Computational Physics, vol. 153, no. 2, pp. 509–534, 1999. doi: 10.1006/jcph.1999.6293.
