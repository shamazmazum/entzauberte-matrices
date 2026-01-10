(in-package :entzauberte-matrices)

(macrolet ((def-foreign-dot (name foreign-type)
             (multiple-value-bind (lisp-name fortran-name)
                 (wrapper-names name)
               `(defcfun (,lisp-name ,fortran-name) ,foreign-type
                  (n    (:pointer :int))
                  (sx   :pointer)
                  (incx (:pointer :int))
                  (sy   :pointer)
                  (incy (:pointer :int))))))
  ;; Real
  (def-foreign-dot sdot :float)
  (def-foreign-dot ddot :double)
  ;; Complex is hard
  #+nil
  (progn
    (def-foreign-dot cdotc :float)
    (def-foreign-dot zdotc :double)))

(macrolet ((def-dot (name low-level-fn lisp-type)
             `(progn
                (serapeum:-> ,name ((vec ,lisp-type) (vec ,lisp-type))
                             (values ,lisp-type &optional))
                (defun ,name (v1 v2)
                  (let ((n (length v1)))
                    (with-foreign-objects ((nptr    :int)
                                           (incxptr :int)
                                           (incyptr :int))
                      (setf (mem-ref nptr    :int) n
                            (mem-ref incxptr :int) 1
                            (mem-ref incyptr :int) 1)
                      (with-array-pointers ((v1ptr v1)
                                            (v2ptr v2))
                        ;; The order does not matter now, but for complex
                        ;; vectors it's v2 v1. BLAS have the arguments swapper
                        ;; for some reason (it computes v1^H v2, and the dot
                        ;; product is v2^H v1).
                        (,low-level-fn nptr v2ptr incxptr v1ptr incyptr))))))))
  (def-dot dot-rs-unsafe %sdot single-float)
  (def-dot dot-rd-unsafe %ddot double-float))

(serapeum:-> dot ((vec *) (vec *))
             (values number &optional))
(declaim (inline dot))
(defun dot (v1 v2)
  "Compute the dot product \\(\\langle v_1, v_2 \\rangle = \\sum_i {v_1}_i
\\overline{{v_2}_i}\\)."
  (assert (= (length v1) (length v2)))
  (funcall
   (cond
     ((and (eq (array-element-type v1) 'single-float)
           (eq (array-element-type v2) 'single-float))
      #'dot-rs-unsafe)
     ((and (eq (array-element-type v1) 'double-float)
           (eq (array-element-type v2) 'double-float))
      #'dot-rd-unsafe)
     (t
      (error "Cannot compute the dot product: Unknown array element types")))
   v1 v2))

(serapeum:-> norm ((vec *))
             (values number &optional))
(declaim (inline norm))
(defun norm (v)
  "Compute the norm of \\(v\\), i.e. \\(\\sqrt{\\langle v, v\\rangle}\\)."
  (sqrt
   ;; Use REALPART to convert to a real. IMAGPART is always 0
   (realpart
    (dot v v))))
