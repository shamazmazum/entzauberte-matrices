(defpackage entzauberte-matrices
  (:use #:cl #:cffi)
  (:export #:lapack-error
           #:mult #:add #:scale
           #:det #:inversions #:invert
           #:set-num-threads))
