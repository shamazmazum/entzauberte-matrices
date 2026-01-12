# Changelog

## Version 0.1

* Incompatible change: Specialized LAPACK wrappers are now total. This means
  that you have to check returned values for non-nullness. The condition
  `lapack-error` is removed. The outer level API functions still can signal an
  error if the input is incompatible (like two matrices having incompatible
  shapes or different element types).
* Incompatible change: The functions `add` and `scale` now use BLAS. Hence,
  `add` now computes αA + B.
* Improvement: A new function `sub` was added
* Improvement: All functions now work with arrays of type `(array TYPE)`, not
  just with simple arrays. Returned values are still simple arrays.
* Improvement: Add `dot` and `norm`. These function work only with real-valued
  arrays as for now.
* Improvement: Add helper functions `vstack`, `hstack`, `reshape` and
  `reshape-unsafe`.

## Version 0.1-rc1

* Matrix multiplication with `mult`.
* Matrix/vector addition and scaling (`add`, `scale`).
* Eigenvalue and eigenvectors (`eig`, `eig-self-adjoint`).
* SVD decomposition (`svd`).
* Determinant (`det`).
* Inversion (`invert`).
* Solver (`solve`).
