(in-package :entzauberte-matrices)

;; This is really trivial
(serapeum:-> scale ((mat-or-vec *) number)
             (values (smat-or-svec *) &optional))
(declaim (inline scale))
(defun scale (v s)
  "Compute \\(sv\\) where \\(s\\) is a scalar and \\(v\\) is a vector
or a matrix. This function is implemented in lisp."
  (unless (typep s (array-element-type v))
    (error "The matrix and the scalar are of incompatible types."))
  (let ((result (make-array (array-dimensions v)
                            :element-type (array-element-type v))))
    (locally
        (declare (optimize (sb-c:insert-array-bounds-checks 0)))
      (loop for i below (array-total-size v) do
        (setf (row-major-aref result i)
              (* s (row-major-aref v i)))))
    result))
