(in-package :entzauberte-matrices)

(define-foreign-library openblas
  (:unix  (:or "libopenblas.so"
               "libopenblas.so.0"))
  (t (:default "libopenblas")))

(use-foreign-library openblas)

;; Good helper macro
(defmacro with-array-pointers (bindings &body body)
  (reduce
   (lambda (binding acc)
     (destructuring-bind (ptr array)
         binding
       `(with-pointer-to-vector-data (,ptr ,array)
          ,acc)))
   bindings
   :initial-value `(progn ,@body)
   :from-end t))

;; Utility function to tell OpenMP not to use all available CPU resources
(cffi:defcfun ("omp_set_num_threads" set-num-threads) :void
  (num-threads :int))

;; Useful type
(deftype mat-or-vec (type) `(or (simple-array ,type 2) (simple-array ,type 1)))
