(in-package :entzauberte-matrices)

;; Low-level bindings
(macrolet ((def-foreign-lu (name)
             (multiple-value-bind (lisp-name foreign-name)
                 (capi-wrapper-names name :lapack)
               `(defcfun (,lisp-name ,foreign-name) lapack-int
                  (layout lapack-order)
                  (m      lapack-int)
                  (n      lapack-int)
                  (a      :pointer)
                  (lda    lapack-int)
                  (ipiv   :pointer)))))
  (def-foreign-lu sgetrf)
  (def-foreign-lu dgetrf)
  (def-foreign-lu cgetrf)
  (def-foreign-lu zgetrf))

;; Semi-high-level functions
(macrolet ((def-lisp-lu (name low-level-fn type)
             `(progn
                (serapeum:-> ,name ((mat ,type))
                             (values (or (smat ,type) null)
                                     (or (svec (unsigned-byte 32)) null)
                                     integer &optional))
                (defun ,name (a)
                  (let* ((m (array-dimension a 0)) ; Number of rows
                         (n (array-dimension a 1)) ; Number of columns
                         (lda m)
                         (ipiv (make-array (min m n) :element-type '(unsigned-byte 32)))
                         (acopy (copy-array a)))
                    (with-array-pointers ((aptr    acopy)
                                          (ipivptr ipiv))
                      (let ((info (,low-level-fn :row-major m n aptr lda ipivptr)))
                        (if (zerop info)
                            (values acopy (fix-pivot ipiv) 0)
                            (values nil nil info)))))))))
  (def-lisp-lu %lu-rs %sgetrf single-float)
  (def-lisp-lu %lu-rd %dgetrf double-float)
  (def-lisp-lu %lu-cs %cgetrf (complex single-float))
  (def-lisp-lu %lu-cd %zgetrf (complex double-float)))

(serapeum:-> %lu ((mat *))
             (values (or (smat *) null)
                     (or (svec (unsigned-byte 32)) null)
                     integer &optional))
(declaim (inline %lu))
(defun %lu (a)
  (when (/= (array-dimension a 0)
            (array-dimension a 1))
    (error "This is not tested"))
  (cond
    ((eq (array-element-type a) 'single-float)
     (%lu-rs a))
    ((eq (array-element-type a) 'double-float)
     (%lu-rd a))
    ((equalp (array-element-type a) '(complex single-float))
     (%lu-cs a))
    ((equalp (array-element-type a) '(complex double-float))
     (%lu-cd a))
    (t
     (error "Cannot perform LU factorization: Unknown array element type"))))

;; TODO: Define high-level functions
