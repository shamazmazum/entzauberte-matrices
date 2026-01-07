(in-package :entzauberte-matrices)

(macrolet ((def-svd-low-level (name complexp)
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
                  ,@(if complexp `((rwork :pointer)))
                  (info  (:pointer :int))))))
  (def-svd-low-level sgesvd nil)
  (def-svd-low-level dgesvd nil)
  (def-svd-low-level cgesvd t)
  (def-svd-low-level zgesvd t))

(macrolet ((def-svd (name low-level-fn lisp-type foreign-type)
             (let* ((complexp (listp lisp-type))
                    (real-type (if complexp (second lisp-type) lisp-type)))
               `(progn
                  (serapeum:-> ,name ((mat ,lisp-type) boolean)
                               (values (mat ,lisp-type)
                                       (vec ,real-type)
                                       (mat ,lisp-type)
                                       &optional))
                  (defun ,name (a compact)
                    (let* ((m (array-dimension a 1)) ; Number of rows in A^T
                           (n (array-dimension a 0)) ; Number of columns in A^T
                           (min (min n m))
                           (acopy (copy-array a))
                           (s  (make-array min  :element-type ',real-type))
                           (u  (make-array (list (if compact min m) m)
                                           :element-type ',lisp-type))
                           (vt (make-array (list n (if compact min n))
                                           :element-type ',lisp-type)))
                      (with-foreign-objects ((jobuptr  :char)
                                             (jobvtptr :char)
                                             (mptr     :int)
                                             (nptr     :int)
                                             (ldaptr   :int)
                                             (lduptr   :int)
                                             (ldvtptr  :int)
                                             ,@(if complexp
                                                   `((rworkptr ,foreign-type
                                                               (* min 5))))
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
                                  (char-code (if compact #\S #\A))
                                  (mem-ref jobvtptr :char)
                                  (char-code (if compact #\S #\A))
                                  (mem-ref mptr :int) m
                                  (mem-ref nptr :int) n
                                  (mem-ref ldaptr :int) m
                                  (mem-ref lduptr :int) m
                                  (mem-ref ldvtptr :int) (if compact min n))
                            (let ((work
                                    (with-foreign-objects ((workptr
                                                            ,foreign-type ,(if complexp 2 1))
                                                           (lworkptr :int))
                                      (setf (mem-ref lworkptr :int) -1)
                                      (,low-level-fn jobuptr jobvtptr mptr nptr aptr ldaptr
                                                     sptr uptr lduptr vtptr ldvtptr workptr
                                                     lworkptr ,@(if complexp '(rworkptr))
                                                     infoptr)
                                      (check-info)
                                      (round (mem-ref workptr ,foreign-type)))))
                              (with-foreign-objects ((workptr
                                                      ,foreign-type
                                                      ,(if complexp '(* work 2) 'work))
                                                     (lworkptr :int))
                                (setf (mem-ref lworkptr :int) work)
                                (,low-level-fn jobuptr jobvtptr mptr nptr aptr ldaptr
                                               sptr uptr lduptr vtptr ldvtptr workptr
                                               lworkptr ,@(if complexp '(rworkptr))
                                               infoptr)
                                (check-info))))))
                      (values vt s u)))))))
  (def-svd svd-rs %sgesvd single-float :float)
  (def-svd svd-rd %dgesvd double-float :double)
  (def-svd svd-cs %cgesvd (complex single-float) :float)
  (def-svd svd-cd %zgesvd (complex double-float) :double))

(serapeum:-> svd ((mat *) &key (:compact boolean))
             (values (mat *) (vec *) (mat *) &optional))
(declaim (inline svd))
(defun svd (m &key compact)
  (funcall
   (cond
     ((eq (array-element-type m) 'single-float)
      #'svd-rs)
     ((eq (array-element-type m) 'double-float)
      #'svd-rd)
     ((equalp (array-element-type m) '(complex single-float))
      #'svd-cs)
     ((equalp (array-element-type m) '(complex double-float))
      #'svd-cd)
     (t
      (error "Cannot compute SVD decomposition: unknown array element type")))
   m compact))
