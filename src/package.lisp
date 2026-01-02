(defpackage entzauberte-matrices
  (:use #:cl #:cffi)
  (:export #:lapack-error
           #:mult #:add #:scale
           #:det #:inversions #:invert
           #:eig #:eig-self-adjoint
           #:set-num-threads))
