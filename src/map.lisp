(in-package :entzauberte-matrices)

(serapeum:-> map (function t (mat *))
             (values (mat *) &optional))
(declaim (inline map))
(defun map (f element-type m)
  "Apply an unary operation @c(f) to a matrix @c(m)
element-wise. @c(element-type) must match the element type of the
result."
  (let ((result (make-array (array-dimensions m)
                            :element-type element-type)))
    (map-into
     (array-storage-vector result) f
     (array-storage-vector m))
    result))

(serapeum:-> zip-with (function t (mat *) (mat *))
             (values (mat *) &optional))
(declaim (inline zip-with))
(defun zip-with (f element-type x y)
  "Apply a binary operation @c(f) to matrices @c(x) and @c(y)
element-wise. @c(element-type) must match the element type of the
result."
  (unless (equal (array-dimensions x)
                 (array-dimensions y))
    (error "Matrices have incompatible dimensions"))
  (let ((result (make-array (array-dimensions x)
                            :element-type element-type)))
    (map-into
     (array-storage-vector result) f
     (array-storage-vector x)
     (array-storage-vector y))
    result))
