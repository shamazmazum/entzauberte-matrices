(in-package :entzauberte-matrices)

(deftype uplo () '(member :upper :lower))

(macrolet ((def-self-adjoint-low-level (name)
             (multiple-value-bind (lisp-name foreign-name)
                 (capi-wrapper-names name :lapack)
               `(defcfun (,lisp-name ,foreign-name) lapack-int
                  (layout lapack-order)
                  (jobz   :char)
                  (uplo   :char)
                  (n      lapack-int)
                  (a      :pointer)
                  (lda    lapack-int)
                  (w      :pointer)))))
  (def-self-adjoint-low-level ssyev)
  (def-self-adjoint-low-level dsyev)
  (def-self-adjoint-low-level cheev)
  (def-self-adjoint-low-level zheev))

(macrolet ((def-self-adjoint-lisp (name low-level-fn lisp-type)
             (let ((real-type (if (listp lisp-type)
                                  (second lisp-type)
                                  lisp-type)))
               `(progn
                  (serapeum:-> ,name ((mat ,lisp-type) uplo)
                               (values (or (svec ,real-type) null)
                                       (or (smat ,lisp-type) null)
                                       integer
                                       &optional))
                  (defun ,name (a where)
                    (let* ((n (array-dimension a 0))
                           (lda n)
                           (copy (copy-array a))
                           (vals (make-array n :element-type ',real-type)))
                      (with-array-pointers ((aptr copy)
                                            (wptr vals))
                        (let ((info (,low-level-fn
                                     :row-major
                                     (char-code #\V)
                                     (char-code
                                      (if (eq where :upper) #\U #\L))
                                     n aptr lda wptr)))
                          (declare (type fixnum info))
                          (if (zerop info)
                              (values vals copy 0)
                              (values nil nil info))))))))))
  (def-self-adjoint-lisp eig-self-adjoint-rs-unsafe %ssyev single-float)
  (def-self-adjoint-lisp eig-self-adjoint-rd-unsafe %dsyev double-float)
  (def-self-adjoint-lisp eig-self-adjoint-cs-unsafe %cheev (complex single-float))
  (def-self-adjoint-lisp eig-self-adjoint-cd-unsafe %zheev (complex double-float)))

(serapeum:-> eig-self-adjoint ((mat *) uplo)
             (values (or (svec *) null)
                     (or (smat *) null)
                     integer &optional))
(declaim (inline eig-self-adjoint))
(defun eig-self-adjoint (m where)
    "Compute eigenvalues and eigenvectors of a self-adjoint matrix
\\(M\\). Since version 0.3 eigenvectors are stored @b(in columns) as
usual, so that (in pseudo-code)

@begin[lang=lisp](code)
(multiple-value-bind (vals vecs)
    (eig-self-adjoint m where)
  (approx= (mult m vecs)
           (mult vecs (from-diag vals))))
@end(code)

is @c(T). Before version 0.3 eigenvectors were stored @b(in rows). The
returned values @c(vals) and @c(vecs) may be @c(NIL) if the
decomposition fails. The third returned value is the info code
returned by LAPACK."
  (assert (= (array-dimension m 0)
             (array-dimension m 1)))
  (funcall
   (cond
     ((eq (array-element-type m) 'single-float)
      #'eig-self-adjoint-rs-unsafe)
     ((eq (array-element-type m) 'double-float)
      #'eig-self-adjoint-rd-unsafe)
     ((equalp (array-element-type m) '(complex single-float))
      #'eig-self-adjoint-cs-unsafe)
     ((equalp (array-element-type m) '(complex double-float))
      #'eig-self-adjoint-cd-unsafe)
     (t
      (error "Cannot compute eigenvalues: unknown array element type")))
   m where))
