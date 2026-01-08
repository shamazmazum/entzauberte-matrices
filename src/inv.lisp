(in-package :entzauberte-matrices)

(macrolet ((def-inv-low-level (name)
             (let* ((name (symbol-name name))
                    (lisp-name (intern (format nil "%~a" name)))
                    (foreign-name (format nil "~a_" (string-downcase name))))
               `(defcfun (,lisp-name ,foreign-name) :void
                  (n     (:pointer :int))
                  (a     :pointer)
                  (lda   (:pointer :int))
                  (ipiv  :pointer)
                  (work  :pointer)
                  (lwork (:pointer :int))
                  (info  (:pointer :int))))))
  (def-inv-low-level sgetri)
  (def-inv-low-level dgetri)
  (def-inv-low-level cgetri)
  (def-inv-low-level zgetri))

(macrolet ((def-inv-lisp (name lisp-type foreign-type lufn invfn)
             (let ((complexp (listp lisp-type)))
               `(progn
                  (serapeum:-> ,name ((simple-array ,lisp-type 2))
                               (values (simple-array ,lisp-type 2) &optional))
                  (defun ,name (a)
                    (let* ((n (array-dimension a 0))
                           (lda n)
                           (acopy (copy-array a)))
                      (with-foreign-objects ((mptr    :int)
                                             (nptr    :int)
                                             (ipivptr :int n)
                                             (ldaptr  :int)
                                             (infoptr :int))
                        (flet ((check-info (infoptr)
                                 (let ((info (mem-ref infoptr :int)))
                                   (when (< info 0)
                                     (error 'lapack-error
                                            :message "Cannot invert a matrix")))))
                          (setf (mem-ref mptr   :int) n
                                (mem-ref nptr   :int) n
                                (mem-ref ldaptr :int) lda)
                          (with-array-pointers ((aptr acopy))
                            (,lufn mptr nptr aptr ldaptr ipivptr infoptr)
                            (check-info infoptr)
                            (let ((work
                                    (with-foreign-objects ((workptr
                                                            ,foreign-type
                                                            ,@(if complexp '(2)))
                                                           (lworkptr :int))
                                      (setf (mem-ref lworkptr :int) -1)
                                      (,invfn nptr aptr ldaptr ipivptr workptr
                                              lworkptr infoptr)
                                      (check-info infoptr)
                                      (round (mem-ref workptr ,foreign-type)))))
                              (with-foreign-objects ((workptr
                                                      :float
                                                      ,(if complexp `(* work 2) 'work))
                                                     (lworkptr :int))
                                (setf (mem-ref lworkptr :int) work)
                                (,invfn nptr aptr ldaptr ipivptr workptr lworkptr infoptr)
                                (check-info infoptr))))))
                      acopy))))))
  (def-inv-lisp inv-rs-unsafe single-float :float  %sgetrf %sgetri)
  (def-inv-lisp inv-rd-unsafe double-float :double %dgetrf %dgetri)
  (def-inv-lisp inv-cs-unsafe (complex single-float) :float  %cgetrf %cgetri)
  (def-inv-lisp inv-cd-unsafe (complex double-float) :double %zgetrf %zgetri))

(serapeum:-> invert ((mat *))
             (values (mat *) &optional))
(declaim (inline invert))
(defun invert (m)
  "Compute \\(M^{-1}\\) using \\(LU\\) factorization."
  (unless (= (array-dimension m 0)
             (array-dimension m 1))
    (error "Cannot invert a matrix: non-square matrix"))
  (cond
    ((eq (array-element-type m) 'single-float)
     (inv-rs-unsafe m))
    ((eq (array-element-type m) 'double-float)
     (inv-rd-unsafe m))
    ((equalp (array-element-type m) '(complex single-float))
     (inv-cs-unsafe m))
    ((equalp (array-element-type m) '(complex double-float))
     (inv-cd-unsafe m))
    (t
     (error "Cannot invert a matrix: Unknown array element type"))))

