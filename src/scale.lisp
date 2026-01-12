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
               `(progn
                  (serapeum:-> ,name ((mat-or-vec ,lisp-type) ,lisp-type)
                               (values (smat-or-svec ,lisp-type) &optional))
                  (defun ,name (v s)
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
                      v))))))
  (def-scale scale-rs %sscal single-float :float)
  (def-scale scale-rd %dscal double-float :double)
  (def-scale scale-cs %cscal (complex single-float) :float)
  (def-scale scale-cd %zscal (complex double-float) :double))

(serapeum:-> scale ((mat-or-vec *) number)
             (values (smat-or-svec *) &optional))
(declaim (inline scale))
(defun scale (m s)
  "Compute \\(s m\\) where \\(s\\) is a scalar and \\(m\\) is a matrix
or a vector."
  (unless (typep s (array-element-type m))
    (error "Cannot scale a matrix: Incompatible types"))
  (funcall
   (cond
     ((eq (array-element-type m) 'single-float)
      #'scale-rs)
     ((eq (array-element-type m) 'double-float)
      #'scale-rd)
     ((equalp (array-element-type m) '(complex single-float))
      #'scale-cs)
     ((equalp (array-element-type m) '(complex double-float))
      #'scale-cd)
     (t
      (error "Cannot scale a matrix: Unknown element-type")))
   m s))
