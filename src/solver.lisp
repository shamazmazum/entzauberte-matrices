(in-package :entzauberte-matrices)

(declaim (inline transpose))
(defun transpose (m)
  (let ((res (make-array (reverse (array-dimensions m))
                         :element-type (array-element-type m))))
    (loop for i below (array-dimension m 0) do
      (loop for j below (array-dimension m 1) do
        (setf (aref res j i) (aref m i j))))
    res))

(macrolet ((def-solver-low-level (name)
             (let* ((name (symbol-name name))
                    (lisp-name (intern (format nil "%~a" name)))
                    (foreign-name (format nil "~a_" (string-downcase name))))
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
                             (values (mat ,lisp-type) &optional))
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
                               (unless (zerop (mem-ref infoptr :int))
                                 (error 'lapack-error :message "Cannot solve an equation"))))
                        (setf (mem-ref nptr    :int) n
                              (mem-ref nrhsptr :int) nrhs
                              (mem-ref ldaptr  :int) n
                              (mem-ref ldbptr  :int) n)
                        (with-array-pointers ((aptr acopy)
                                              (bptr bcopy))
                          (,low-level-fn
                           nptr nrhsptr aptr ldaptr ipivptr bptr ldbptr infoptr)
                          (check-info))))
                    (transpose bcopy))))))
  (def-solve solve-rs-unsafe %sgesv single-float)
  (def-solve solve-rd-unsafe %dgesv double-float)
  (def-solve solve-cs-unsafe %cgesv (complex single-float))
  (def-solve solve-cd-unsafe %zgesv (complex double-float)))

(serapeum:-> solve ((mat *) (mat *))
             (values (mat *) &optional))
(declaim (inline solve))
(defun solve (m b)
  (assert (= (array-dimension m 0)
             (array-dimension m 1)
             (array-dimension b 0)))
  (funcall
   (cond
     ((eq (array-element-type m) 'single-float)
      #'solve-rs-unsafe)
     ((eq (array-element-type m) 'double-float)
      #'solve-rd-unsafe)
     ((equalp (array-element-type m) '(complex single-float))
      #'solve-cs-unsafe)
     ((equalp (array-element-type m) '(complex double-float))
      #'solve-cd-unsafe)
     (t
      (error "Cannot solve an equation: unknown array element type")))
   m b))
