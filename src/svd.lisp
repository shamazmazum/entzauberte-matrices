(in-package :entzauberte-matrices)

(macrolet ((def-svd-low-level (name)
             (multiple-value-bind (lisp-name foreign-name)
                 (capi-wrapper-names name :lapack)
               `(defcfun (,lisp-name ,foreign-name) lapack-int
                  (layout lapack-order)
                  (jobu   :char)
                  (jobvt  :char)
                  (m      lapack-int)
                  (n      lapack-int)
                  (a      :pointer)
                  (lda    lapack-int)
                  (s      :pointer)
                  (u      :pointer)
                  (ldu    lapack-int)
                  (vt     :pointer)
                  (ldvt   lapack-int)
                  (superb :pointer)))))
  (def-svd-low-level sgesvd)
  (def-svd-low-level dgesvd)
  (def-svd-low-level cgesvd)
  (def-svd-low-level zgesvd))

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
                    (let* ((m (array-dimension a 0))
                           (n (array-dimension a 1))
                           (min (min n m))
                           (acopy (copy-array a))
                           (s  (make-array min  :element-type ',real-type))
                           (u  (make-array (list m (if compact min m))
                                           :element-type ',lisp-type))
                           (vt (make-array (list (if compact min n) n)
                                           :element-type ',lisp-type))
                           (job (char-code (if compact #\S #\A))))
                      (with-array-pointers ((aptr    acopy)
                                            (sptr    s)
                                            (uptr    u)
                                            (vtptr   vt))
                        (with-foreign-object (sprbptr ,foreign-type (1- min))
                          (let ((info (,low-level-fn :row-major
                                                     job job m n aptr
                                                     n sptr uptr
                                                     (if compact min m)
                                                     vtptr n sprbptr)))
                            (declare (type fixnum info))
                            (if (zerop info)
                                (values u s vt 0)
                                (values nil nil nil info)))))))))))
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
    (svd m :compact t)
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
