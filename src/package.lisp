(defpackage entzauberte-matrices
  (:use #:cl #:cffi)
  (:shadow #:map)
  (:export #:lapack-error
           #:reshape #:reshape-unsafe
           #:vector->column #:vector->row
           #:vector->column-unsafe #:vector->row-unsafe
           #:map #:zip-with
           #:row #:column #:vstack #:hstack
           #:dot #:norm
           #:@ #:mult #:add #:sub #:scale
           #:det #:inversions #:invert
           #:eig #:eig-self-adjoint #:svd
           #:solve #:transpose
           #:set-num-threads))
