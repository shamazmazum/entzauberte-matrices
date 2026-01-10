(defpackage entzauberte-matrices
  (:use #:cl #:cffi)
  (:export #:lapack-error
           #:dot #:norm
           #:mult #:add #:sub #:scale
           #:det #:inversions #:invert
           #:eig #:eig-self-adjoint #:svd
           #:solve #:transpose
           #:set-num-threads))
