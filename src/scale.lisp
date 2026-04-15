(in-package :entzauberte-matrices)

(macrolet ((def-scale (name scalar-type)
             (multiple-value-bind (lisp-name fortran-name)
                 (capi-wrapper-names name :blas)
               `(defcfun (,lisp-name ,fortran-name) :void
                  (n    blas-int)
                  (sa   ,scalar-type)
                  (sx   :pointer)
                  (incx blas-int)))))
  (def-scale sscal :float)
  (def-scale dscal :double)
  (def-scale cscal :pointer)
  (def-scale zscal :pointer))

(macrolet ((def-scale (name low-level-fn lisp-type &optional foreign-type)
             (let ((complexp foreign-type))
               `(progn
                  (serapeum:-> ,name ((mat-or-vec ,lisp-type) ,lisp-type)
                               (values (smat-or-svec ,lisp-type) &optional))
                  (defun ,name (v s)
                    (let ((v (copy-array v))
                          (n (array-total-size v)))
                      ,(let ((call
                               `(with-array-pointers ((vptr v))
                                  (,low-level-fn
                                   n ,(if complexp `sptr `s) vptr 1))))
                         (if complexp
                             `(with-foreign-object (sptr ,foreign-type 2)
                                (setf (mem-aref sptr ,foreign-type 0)
                                      (realpart s)
                                      (mem-aref sptr ,foreign-type 1)
                                      (imagpart s))
                                ,call)
                             call))
                      v))))))
  (def-scale %scale-rs %sscal single-float)
  (def-scale %scale-rd %dscal double-float)
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
