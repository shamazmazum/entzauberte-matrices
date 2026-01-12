(in-package :entzauberte-matrices)

(macrolet ((def-svd-low-level (name complexp)
             (multiple-value-bind (lisp-name foreign-name)
                 (wrapper-names name)
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
                               (values (or (smat ,lisp-type) null)
                                       (or (svec ,real-type) null)
                                       (or (smat ,lisp-type) null)
                                       integer
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
                                 (let ((info (mem-ref infoptr :int)))
                                   (unless (zerop info)
                                     (return-from ,name (values nil nil nil info))))))
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
                      (values vt s u 0)))))))
  (def-svd svd-rs %sgesvd single-float :float)
  (def-svd svd-rd %dgesvd double-float :double)
  (def-svd svd-cs %cgesvd (complex single-float) :float)
  (def-svd svd-cd %zgesvd (complex double-float) :double))

(serapeum:-> svd ((mat *) &key (:compact boolean))
             (values (or (smat *) null)
                     (or (svec *) null)
                     (or (smat *) null)
                     integer &optional))
(declaim (inline svd))
(defun svd (m &key compact)
  "Do the singular value decomposition of a matrix \\(M = U\\Sigma
V^T\\). Matrices \\(U\\), \\(\\Sigma\\) and \\(V^T\\) are returned in
three values. For a matrix \\(m\\times n\\) and \\(r = \\min(m, n)\\),
if @c(:compact) is @c(NIL) then \\(U\\) is an orthogonal \\(m \\times
m\\) matrix, \\(V^T\\) is an orthogonal \\(n \\times n\\) matrix and
\\(\\Sigma\\) is a vector of the length \\(r\\), so that in pseudocode

@begin[lang=lisp](code)
(multiple-value-bind (U Σ VT)
    (svd m)
  (approx= (mult U (mult (from-diag Σ m n) VT)) m))
@end(code)

is @c(T). If @c(:compact) is @c(T) then \\(U\\) is a semi-orthogonal
\\(m \\times r\\) matrix, \\(V^T\\) is a semi-orthogonal \\(r \\times
n\\) matrix and \\(\\Sigma\\) is as with @c(:compact) being @c(NIL), so
that in pseudocode

@begin[lang=lisp](code)
(multiple-value-bind (U Σ VT)
    (svd m)
  (let ((r (min (array-dimension m 0)
                (array-dimension m 1))))
    (approx= (mult U (mult (from-diag Σ r r) VT)) m)))
@end(code)

is @c(T). The first three returned values may be @c(NIL) if the
decomposition fails. The fourth value is the info code returned by
LAPACK."
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
