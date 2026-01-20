(in-package :entzauberte-matrices)

(sb-c:defknown array-storage-vector (array) (simple-array * (*))
    (sb-c:flushable sb-c::recursive)
  :overwrite-fndb-silently t)

(defun array-storage-vector (a)
  "This is like SB-EXT:ARRAY-STORAGE-VECTOR, but works also with
displaced arrays."
  (let ((displacement (array-displacement a)))
    (if displacement
        (array-storage-vector displacement)
        (sb-ext:array-storage-vector a))))

(sb-c:defoptimizer (array-storage-vector sb-c:derive-type) ((array))
  (let* ((type (sb-c::lvar-conservative-type array))
         (dims (array-type-dimensions type))
         (et   (sb-c::array-type-upgraded-element-type type)))
    (sb-kernel:make-array-type
     (cond
       ((atom dims) dims)
       ((every #'integerp dims)
        (list (reduce #'* dims)))
       (t '(*)))
     :element-type et
     :specialized-element-type et
     :complexp nil)))

(sb-c:deftransform array-storage-vector ((array) (simple-array))
  '(sb-ext:array-storage-vector array))

;; Good helper macro
(defmacro with-array-pointers (bindings &body body)
  (reduce
   (lambda (binding acc)
     (destructuring-bind (ptr array)
         binding
       `(with-pointer-to-vector-data (,ptr (array-storage-vector ,array))
          ,acc)))
   bindings
   :initial-value `(progn ,@body)
   :from-end t))
