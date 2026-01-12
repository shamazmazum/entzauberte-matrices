(in-package :entzauberte-matrices)

(declaim (inline transpose))
(defun transpose (m)
  "Transpose a matrix \\(M\\) (slow for big matrices, avoid this
function if possible)."
  (let ((res (make-array (reverse (array-dimensions m))
                         :element-type (array-element-type m))))
    (loop for i below (array-dimension m 0) do
      (loop for j below (array-dimension m 1) do
        (setf (aref res j i) (aref m i j))))
    res))

(macrolet ((def-solver-low-level (name)
             (multiple-value-bind (lisp-name foreign-name)
                 (wrapper-names name)
               `(defcfun (,lisp-name ,foreign-name) :void
                  (n    (:pointer :int))
                  (nrhs (:pointer :int))
                  (a    :pointer)
                  (lda  (:pointer :int))
                  (ipiv :pointer)
                  (b    :pointer)
                  (ldb  (:pointer :int))
                  (info (:pointer :int))))))
  (def-solver-low-level sgesv)
  (def-solver-low-level dgesv)
  (def-solver-low-level cgesv)
  (def-solver-low-level zgesv))

;; Here we must to transpose, not copy. Otherwise we'll solve xA=b
(macrolet ((def-solve (name low-level-fn lisp-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type)
                                    (mat ,lisp-type))
                             (values (or (smat ,lisp-type) null) integer &optional))
                (defun ,name (a b)
                  (let ((acopy (transpose a))
                        (bcopy (transpose b))
                        (n    (array-dimension a 0))
                        (nrhs (array-dimension b 1)))
                    (with-foreign-objects ((nptr    :int)
                                           (nrhsptr :int)
                                           (ldaptr  :int)
                                           (ipivptr :int n)
                                           (ldbptr  :int)
                                           (infoptr :int))
                      (flet ((check-info ()
                               (let ((info (mem-ref infoptr :int)))
                                 (unless (zerop info)
                                   (return-from ,name (values nil info))))))
                        (setf (mem-ref nptr    :int) n
                              (mem-ref nrhsptr :int) nrhs
                              (mem-ref ldaptr  :int) n
                              (mem-ref ldbptr  :int) n)
                        (with-array-pointers ((aptr acopy)
                                              (bptr bcopy))
                          (,low-level-fn
                           nptr nrhsptr aptr ldaptr ipivptr bptr ldbptr infoptr)
                          (check-info))))
                    (values (transpose bcopy) 0))))))
  (def-solve solve-rs-unsafe %sgesv single-float)
  (def-solve solve-rd-unsafe %dgesv double-float)
  (def-solve solve-cs-unsafe %cgesv (complex single-float))
  (def-solve solve-cd-unsafe %zgesv (complex double-float)))

(serapeum:-> solve ((mat *) (mat *))
             (values (or (smat *) null) integer &optional))
(declaim (inline solve))
(defun solve (a b)
  "Solve an equation \\(AX = B\\). This function uses @c(transpose)
internally. Return @c((values X infocode)) where @c(X) may be @c(NIL)
if the solver fails."
  (assert (= (array-dimension a 0)
             (array-dimension a 1)
             (array-dimension b 0)))
  (funcall
   (cond
     ((eq (array-element-type a) 'single-float)
      #'solve-rs-unsafe)
     ((eq (array-element-type a) 'double-float)
      #'solve-rd-unsafe)
     ((equalp (array-element-type a) '(complex single-float))
      #'solve-cs-unsafe)
     ((equalp (array-element-type a) '(complex double-float))
      #'solve-cd-unsafe)
     (t
      (error "Cannot solve an equation: unknown array element type")))
   a b))
