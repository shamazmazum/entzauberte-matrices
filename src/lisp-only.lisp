(in-package :entzauberte-matrices)

;; BLAS function requires a copy of one of the input matrices (the
;; operation is destructive), so I think, there is no point in using
;; BLAS. Here is one (inlined) function for all types.

(serapeum:-> add ((mat-or-vec *) (mat-or-vec *)
                  &key (:c1 number) (:c2 number))
             (values (mat-or-vec *) &optional))
(declaim (inline add))
(defun add (m1 m2 &key c1 c2)
  (let ((c1 (if c1 c1 (coerce 1 (array-element-type m1))))
        (c2 (if c2 c2 (coerce 1 (array-element-type m2)))))
    (unless (and (equalp (array-element-type m1)
                         (array-element-type m2))
                 (typep c1 (array-element-type m1))
                 (typep c2 (array-element-type m2)))
      (error "Cannot add matrices: Incompatible types"))
    ;; NB: SBCL cannot elliminate a call to ARRAY-DIMENSIONS
    (unless (or (and (= (array-rank m1)
                        (array-rank m2) 2)
                     (= (array-dimension m1 0)
                        (array-dimension m2 0))
                     (= (array-dimension m1 1)
                        (array-dimension m2 1)))
                (and (= (array-rank m1)
                        (array-rank m2) 1)
                     (= (length m1) (length m2))))
      (error "Cannot add matrices: Incompatible dimensions"))
    (let ((result (make-array (array-dimensions m1)
                              :element-type (array-element-type m1))))
      ;; I hope someday this declaration wouldn't be needed.
      (locally
          (declare (optimize (sb-c:insert-array-bounds-checks 0)))
        (loop for i below (array-total-size m1) do
          (setf (row-major-aref result i)
                (+ (* c1 (row-major-aref m1 i))
                   (* c2 (row-major-aref m2 i))))))
      result)))
