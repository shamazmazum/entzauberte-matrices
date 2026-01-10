(in-package :entzauberte-matrices)

(declaim (inline reshape))
(defun reshape (m shape)
  "Reshape an array."
  (assert (= (reduce #'* shape)
             (array-total-size m)))
  (let ((result (make-array shape :element-type (array-element-type m))))
    (replace (array-storage-vector result)
             (array-storage-vector m))
    result))

(declaim (inline reshape-unsafe))
(defun reshape-unsafe (m shape)
  "Reshape an array. The new array will share the storage with the argument."
  (assert (= (reduce #'* shape)
             (array-total-size m)))
  (make-array shape
              :element-type (array-element-type m)
              :displaced-to m))
