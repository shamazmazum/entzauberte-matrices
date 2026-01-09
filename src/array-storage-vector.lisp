(in-package :entzauberte-matrices)

(sb-c:defknown array-storage-vector (array) (simple-array * (*))
    (sb-c:flushable sb-c::recursive))

(defun array-storage-vector (a)
  "This is like SB-EXT:ARRAY-STORAGE-VECTOR, but works also with
displaced arrays."
  (let ((displacement (array-displacement a)))
    (if displacement
        (array-storage-vector displacement)
        (sb-ext:array-storage-vector a))))

(sb-c:defoptimizer (array-storage-vector sb-c:derive-type) ((array))
  (let ((type (sb-c::lvar-type array)))
    (when (sb-kernel:array-type-p type)
      (let ((dimensions (sb-kernel:array-type-dimensions type)))
        (sb-kernel:make-array-type
         (cond
           ((atom dimensions) dimensions)
           ((every #'integerp dimensions)
            (list (reduce #'* dimensions)))
           (t '*))
         :element-type (sb-kernel:array-type-element-type type)
         :specialized-element-type (sb-kernel:array-type-specialized-element-type type)
         :complexp nil)))))

(sb-c:deftransform array-storage-vector ((array) ((simple-array * (*))))
  'array)
