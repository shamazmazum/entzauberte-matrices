(in-package :entzauberte-matrices)

(macrolet ((def-add (name scalar-type)
             (multiple-value-bind (lisp-name c-name)
                 (capi-wrapper-names name :blas)
               `(defcfun (,lisp-name ,c-name) :void
                  (n    blas-int)
                  (sa   ,scalar-type)
                  (sx   :pointer)
                  (incx blas-int)
                  (sy   :pointer)
                  (incy blas-int)))))
  (def-add saxpy :float)
  (def-add daxpy :double)
  (def-add caxpy :pointer)
  (def-add zaxpy :pointer))

(macrolet ((def-add (name low-level-fn lisp-type &optional foreign-type)
             (let ((complexp foreign-type))
               `(progn
                  (serapeum:-> ,name ((mat-or-vec ,lisp-type)
                                      (mat-or-vec ,lisp-type)
                                      ,lisp-type)
                               (values (smat-or-svec ,lisp-type) &optional))
                  (defun ,name (v1 v2 s)
                    (let ((v2 (copy-array v2))
                          (n  (array-total-size v1)))
                      ,(let ((call
                               `(with-array-pointers ((v1ptr v1)
                                                      (v2ptr v2))
                                  (,low-level-fn
                                   n ,(if complexp `sptr `s) v1ptr 1 v2ptr 1))))
                         (if complexp
                             `(with-foreign-object (sptr ,foreign-type 2)
                                (setf (mem-aref sptr ,foreign-type 0)
                                      (realpart s)
                                      (mem-aref sptr ,foreign-type 1)
                                      (imagpart s))
                                ,call)
                             call))
                      v2))))))
  (def-add add-rs-unsafe %saxpy single-float)
  (def-add add-rd-unsafe %daxpy double-float)
  (def-add add-cs-unsafe %caxpy (complex single-float) :float)
  (def-add add-cd-unsafe %zaxpy (complex double-float) :double))

;; There is ADD-UNSAFE because otherwise SBCL cannot infer dimensions
;; of the result. ADD-UNSAFE uses a stupid type deriver.
(sb-c:defknown add-unsafe ((mat-or-vec *) (mat-or-vec *) number)
    (smat-or-svec *)
    (sb-c:foldable sb-c:flushable)
  :overwrite-fndb-silently t
  :derive-type #'first-arg-array-type-deriver)

(def-op-specializers add-unsafe (v1 v2 s)
  ((          single-float . add-rs-unsafe)
   (          double-float . add-rd-unsafe)
   ((complex single-float) . add-cs-unsafe)
   ((complex double-float) . add-cd-unsafe))
  "Cannot add vectors: Unsupported array type")

(serapeum:-> add ((mat-or-vec *) (mat-or-vec *) &optional number)
             (values (smat-or-svec *) &optional))
(declaim (inline add))
(defun add (m1 m2 &optional (s (coerce 1 (array-element-type m1))))
  "Compute \\(s m_1 + m_2\\) where \\(s\\) is a scalar and \\(m_1\\)
and \\(m_2\\) are two matrices of vectors.. By default \\(s\\) is
\\(1\\)."
  (unless (equalp (array-dimensions m1)
                  (array-dimensions m2))
    (error "Cannot add matrices: Incompatible dimensions"))
  (let ((t1 (array-element-type m1))
        (t2 (array-element-type m2)))
    (unless (and (typep s t1) (equalp t1 t2))
      (error "Cannot add: incompatible element types: ~a, ~a and ~a"
             t1 t2 (type-of s))))
  (add-unsafe m1 m2 s))

(serapeum:-> sub ((mat-or-vec *) (mat-or-vec *))
             (values (smat-or-svec *) &optional))
(declaim (inline sub))
(defun sub (m1 m2)
  "Compute \\(m_1 - m_2\\). A special case of @c(add)."
  (add m2 m1 (coerce -1 (array-element-type m1))))
