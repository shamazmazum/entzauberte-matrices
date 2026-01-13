(in-package :entzauberte-matrices)

(macrolet ((def-scale (name)
             (multiple-value-bind (lisp-name fortran-name)
                 (wrapper-names name)
               `(defcfun (,lisp-name ,fortran-name) :void
                  (n    (:pointer :int))
                  (sa   :pointer)
                  (sx   :pointer)
                  (incx (:pointer :int))))))
  (def-scale sscal)
  (def-scale dscal)
  (def-scale cscal)
  (def-scale zscal))

(macrolet ((def-scale (name low-level-fn lisp-type foreign-type)
             (let ((complexp (listp lisp-type)))
               `(defun ,name (v s)
                  (let ((v (copy-array v))
                        (n (array-total-size v)))
                    (with-foreign-objects ((nptr    :int)
                                           (saptr   ,foreign-type ,(if complexp 2 1))
                                           (incxptr :int))
                      (setf (mem-ref nptr    :int) n
                            (mem-ref incxptr :int) 1
                            ,@(if complexp
                                  `((mem-aref saptr ,foreign-type 0)
                                    (realpart s)
                                    (mem-aref saptr ,foreign-type 1)
                                    (imagpart s))
                                  `((mem-ref saptr ,foreign-type) s)))
                      (with-array-pointers ((vptr v))
                        (,low-level-fn nptr saptr vptr incxptr)))
                      v)))))
  (def-scale %scale-rs %sscal single-float :float)
  (def-scale %scale-rd %dscal double-float :double)
  (def-scale %scale-cs %cscal (complex single-float) :float)
  (def-scale %scale-cd %zscal (complex double-float) :double))

;; See add.lisp
(sb-c:defknown %scale ((mat-or-vec *) number)
    (smat-or-svec *)
    (sb-c:foldable sb-c:flushable)
  :overwrite-fndb-silently t
  :derive-type #'first-arg-array-type-deriver)

(def-op-specializers %scale (v s)
  ((          single-float . %scale-rs)
   (          double-float . %scale-rd)
   ((complex single-float) . %scale-cs)
   ((complex double-float) . %scale-cd))
  "Cannot scale a vector: Unsupported array type")

(serapeum:-> scale ((mat-or-vec *) number)
             (values (smat-or-svec *) &optional))
(declaim (inline scale))
(defun scale (m s)
  "Compute \\(s m\\) where \\(s\\) is a scalar and \\(m\\) is a matrix
or a vector."
  (unless (typep s (array-element-type m))
    (error "Cannot scale a matrix: Incompatible types"))
  (%scale m s))
