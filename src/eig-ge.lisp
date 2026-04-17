(in-package :entzauberte-matrices)

(macrolet ((def-eig-low-level (name complexp)
             (multiple-value-bind (lisp-name foreign-name)
                 (capi-wrapper-names name :lapack)
               `(defcfun (,lisp-name ,foreign-name) lapack-int
                  (layout lapack-order)
                  (jobvl :char)
                  (jobvr :char)
                  (n     lapack-int)
                  (a     :pointer)
                  (lda   lapack-int)
                  ,@(if complexp
                        `((w  :pointer))
                        `((wr :pointer)
                          (wi :pointer)))
                  (vl    :pointer)
                  (ldvl  lapack-int)
                  (vr    :pointer)
                  (ldvr  lapack-int)))))
  (def-eig-low-level sgeev nil)
  (def-eig-low-level dgeev nil)
  (def-eig-low-level cgeev t)
  (def-eig-low-level zgeev t))

;; Complex
(macrolet ((def-eig (name low-level-fn lisp-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type))
                             (values (or (svec ,lisp-type) null)
                                     (or (smat ,lisp-type) null)
                                     integer
                                     &optional))
                (defun ,name (a)
                  (let* ((n (array-dimension a 0))
                         (lda n)
                         (copy (copy-array a))
                         (vals    (make-array n :element-type ',lisp-type))
                         (vecs    (make-array (list n n) :element-type ',lisp-type)))
                    (with-array-pointers ((aptr  copy)
                                          (wptr  vals)
                                          (vrptr vecs))
                      (let ((info (,low-level-fn :row-major
                                                 (char-code #\N)
                                                 (char-code #\V)
                                                 n aptr lda wptr
                                                 (null-pointer) 1
                                                 vrptr n)))
                        (if (zerop info)
                            (values vals vecs 0)
                            (values nil nil info)))))))))
  (def-eig eig-cs-unsafe %cgeev (complex single-float))
  (def-eig eig-cd-unsafe %zgeev (complex double-float)))

(macrolet ((def-eig (name low-level-fn lisp-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type))
                             (values (or (svec (complex ,lisp-type)) null)
                                     (or (smat (complex ,lisp-type)) null)
                                     integer
                                     &optional))
                (defun ,name (a)
                  (let* ((n (array-dimension a 0))
                         (lda n)
                         (copy (copy-array a))
                         (vals    (make-array n :element-type '(complex ,lisp-type)))
                         (vals-re (make-array n :element-type ',lisp-type))
                         (vals-im (make-array n :element-type ',lisp-type))
                         (vecs-re (make-array (list n n) :element-type ',lisp-type))
                         (vecs    (make-array (list n n) :element-type '(complex ,lisp-type))))
                    (with-array-pointers ((aptr copy)
                                          (wrptr vals-re)
                                          (wiptr vals-im)
                                          (vrptr vecs-re))
                      (let ((info (,low-level-fn :row-major
                                                 (char-code #\N)
                                                 (char-code #\V)
                                                 n aptr lda wrptr wiptr
                                                 (null-pointer) 1
                                                 vrptr n)))
                        (declare (type fixnum info))
                        (cond
                          ((not (zerop info))
                           (values nil nil info))
                          (t
                           (loop for i below n do
                             (setf (aref vals i)
                                   (complex
                                    (aref vals-re i)
                                    (aref vals-im i))))
                           (loop for i below n do
                             (loop for j below n do
                               (let ((im (imagpart (aref vals j))))
                                 (setf (aref vecs i j)
                                       (cond
                                         ((zerop im)
                                          (complex (aref vecs-re i j)))
                                         ((> im 0)
                                          (complex (aref vecs-re i ( + j))
                                                   (aref vecs-re i (1+ j))))
                                         (t
                                          (complex (+ (aref vecs-re i (1- j)))
                                                   (- (aref vecs-re i ( + j))))))))))
                           (values vals vecs 0))))))))))
  (def-eig eig-rs-unsafe %sgeev single-float)
  (def-eig eig-rd-unsafe %dgeev double-float))

(serapeum:-> eig ((mat *))
             (values (or (svec *) null)
                     (or (smat *) null)
                     integer &optional))
(declaim (inline eig))
(defun eig (m)
  "Compute eigenvalues and eigenvectors of \\(M\\). Since version 0.3
eigenvectors are stored @b(in columns), so that (in pseudo-code)

@begin[lang=lisp](code)
(multiple-value-bind (vals vecs)
    (eig m)
  (approx= (mult m vecs)
           (mult vecs (from-diag vals))))
@end(code)

is @c(T). Before version 0.2 eigenvectors were stored @b(in rows). The
returned values @c(vals) and @c(vecs) may be @c(NIL) if the
decomposition fails. The third returned value is the info code
returned by LAPACK."
  (assert (= (array-dimension m 0)
             (array-dimension m 1)))
  (cond
    ((eq (array-element-type m) 'single-float)
     (eig-rs-unsafe m))
    ((eq (array-element-type m) 'double-float)
     (eig-rd-unsafe m))
    ((equalp (array-element-type m) '(complex single-float))
     (eig-cs-unsafe m))
    ((equalp (array-element-type m) '(complex double-float))
     (eig-cd-unsafe m))
    (t
     (error "Cannot compute eigenvalues: unknown array element type"))))
