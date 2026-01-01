(in-package :entzauberte-matrices)

;; Low-level bindings
(macrolet ((def-foreign-lu (name)
             (let* ((name (symbol-name name))
                    (lisp-name (intern (format nil "%~a" name)))
                    (foreign-name (format nil "~a_" (string-downcase name))))
               `(defcfun (,lisp-name ,foreign-name) :void
                  (m    (:pointer :int))
                  (n    (:pointer :int))
                  (a    :pointer)
                  (lda  (:pointer :int))
                  (ipiv :pointer)
                  (info (:pointer :int))))))
  (def-foreign-lu sgetrf)
  (def-foreign-lu dgetrf)
  (def-foreign-lu cgetrf)
  (def-foreign-lu zgetrf))

;; Semi-high-level functions
(macrolet ((def-lisp-lu (name low-level-fn type)
             `(progn
                (serapeum:-> ,name ((simple-array ,type 2))
                             (values (simple-array ,type 2)
                                     (simple-array (unsigned-byte 32) 1)
                                     integer &optional))
                (defun ,name (a)
                  (let* (; Number of rows, but our matrix is row-major
                         (m (array-dimension a 1))
                         (n (array-dimension a 0)) ; Number of columns
                         (lda m)
                         (ipiv (make-array (min m n) :element-type '(unsigned-byte 32)))
                         (acopy (copy-for-ffi a)))
                    (with-foreign-objects ((mptr    :int)
                                           (nptr    :int)
                                           (ldaptr  :int)
                                           (infoptr :int))
                      (setf (mem-ref mptr   :int) m
                            (mem-ref nptr   :int) n
                            (mem-ref ldaptr :int) lda)
                      (with-array-pointers ((aptr    acopy)
                                            (ipivptr ipiv))
                        (,low-level-fn mptr nptr aptr ldaptr ipivptr infoptr))
                      (let ((info (mem-ref infoptr :int)))
                        (when (< info 0)
                          (error 'lapack-error :message "Cannot perform LU decomposition"))
                        (values acopy (fix-pivot ipiv) info))))))))
  (def-lisp-lu %lu-rs %sgetrf single-float)
  (def-lisp-lu %lu-rd %dgetrf double-float)
  (def-lisp-lu %lu-cs %cgetrf (complex single-float))
  (def-lisp-lu %lu-cd %zgetrf (complex double-float)))

(declaim (inline %lu))
(defun %lu (a)
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
