(in-package :entzauberte-matrices)

(macrolet ((def-svd-low-level (name)
             (let* ((name (symbol-name name))
                    (lisp-name (intern (format nil "%~a" name)))
                    (foreign-name (format nil "~a_" (string-downcase name))))
               `(defcfun (,lisp-name ,foreign-name) :void
                  (jobu  (:pointer :char))
                  (jobvt (:pointer :char))
                  (m     (:pointer :int))
                  (n     (:pointer :int))
                  (a     :pointer)
                  (lda   (:pointer :int))
                  (s     :pointer)
                  (u     :pointer)
                  (ldu   (:pointer :int))
                  (vt    :pointer)
                  (ldvt  (:pointer :int))
                  (work  :pointer)
                  (lwork (:pointer :int))
                  (info  (:pointer :int))))))
  (def-svd-low-level sgesvd)
  (def-svd-low-level dgesvd))

(macrolet ((def-svd (name low-level-fn lisp-type foreign-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type))
                             (values (mat ,lisp-type)
                                     (vec ,lisp-type)
                                     (mat ,lisp-type)
                                     &optional))
                (defun ,name (a)
                  (let* ((m (array-dimension a 1)) ; Number of rows in A^T
                         (n (array-dimension a 0)) ; Number of columns in A^T
                         (acopy (copy-array a))
                         (s  (make-array (min n m)  :element-type ',lisp-type))
                         (u  (make-array (list m m) :element-type ',lisp-type))
                         (vt (make-array (list n n) :element-type ',lisp-type)))
                    (with-foreign-objects ((jobuptr  :char)
                                           (jobvtptr :char)
                                           (mptr     :int)
                                           (nptr     :int)
                                           (ldaptr   :int)
                                           (lduptr   :int)
                                           (ldvtptr  :int)
                                           (infoptr  :int))
                      (flet ((check-info ()
                               (unless (zerop (mem-ref infoptr :int))
                                 (error 'lapack-error
                                        :message "Cannot compute SVD decomposition"))))
                        (with-array-pointers ((aptr  acopy)
                                              (sptr  s)
                                              (uptr  u)
                                              (vtptr vt))
                          (setf (mem-ref jobuptr :char)
                                (char-code #\A)
                                (mem-ref jobvtptr :char)
                                (char-code #\A)
                                (mem-ref mptr :int) m
                                (mem-ref nptr :int) n
                                (mem-ref ldaptr :int) m
                                (mem-ref lduptr :int) m
                                (mem-ref ldvtptr :int) n)
                          (let ((work
                                  (with-foreign-objects ((workptr  ,foreign-type)
                                                         (lworkptr :int))
                                    (setf (mem-ref lworkptr :int) -1)
                                    (,low-level-fn jobuptr jobvtptr mptr nptr aptr ldaptr
                                                   sptr uptr lduptr vtptr ldvtptr workptr
                                                   lworkptr infoptr)
                                    (check-info)
                                    (round (mem-ref workptr ,foreign-type)))))
                            (with-foreign-objects ((workptr  ,foreign-type work)
                                                   (lworkptr :int))
                              (setf (mem-ref lworkptr :int) work)
                              (,low-level-fn jobuptr jobvtptr mptr nptr aptr ldaptr
                                             sptr uptr lduptr vtptr ldvtptr workptr
                                             lworkptr infoptr)
                              (check-info))))))
                    (values vt s u))))))
  (def-svd svd-rs %sgesvd single-float :float)
  (def-svd svd-rd %dgesvd double-float :double))

(serapeum:-> svd ((mat *))
             (values (mat *) (vec *) (mat *) &optional))
(declaim (inline svd))
(defun svd (m)
  (cond
    ((eq (array-element-type m) 'single-float)
     (svd-rs m))
    ((eq (array-element-type m) 'double-float)
     (svd-rd m))
    #|
    ((equalp (array-element-type m) '(complex single-float))
     (svd-cs m))
    ((equalp (array-element-type m) '(complex double-float))
     (svd-cd m))
    |#
    (t
     (error "Cannot compute SVD decomposition: unknown array element type"))))
