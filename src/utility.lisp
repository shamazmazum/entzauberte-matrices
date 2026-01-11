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

(serapeum:-> row ((mat *) (integer 0))
             (values (svec *) &optional))
(declaim (inline row))
(defun row (m ri)
  "Get @c(ri)-th row from a matrix @c(m)."
  (let ((v (make-array (array-dimension m 1)
                       :element-type (array-element-type m)))
        (idx (array-row-major-index m ri 0)))
    (loop for i below (length v) do
      (setf (aref v i) (row-major-aref m (+ idx i))))
    v))

(serapeum:-> column ((mat *) (integer 0))
             (values (svec *) &optional))
(declaim (inline column))
(defun column (m ci)
  "Get @c(ci)-th column from a matrix @c(m)."
  (let ((v (make-array (array-dimension m 0)
                       :element-type (array-element-type m))))
    (loop for i below (length v) do
      (setf (aref v i) (aref m i ci)))
    v))

;; It's too big and ineffective to inline it. Write a simple type deriver instead.
(defun stack-type-deriver (call)
  (let ((type (second (sb-c::combination-args call))))
    (when (sb-c::constant-lvar-p type)
      (let ((type (sb-c::lvar-value type)))
        (sb-kernel:specifier-type `(simple-array ,type (* *)))))))

(sb-c:defknown (hstack vstack) (list t) (smat *)
    (sb-c:foldable sb-c:unsafely-flushable)
  :overwrite-fndb-silently t
  :derive-type #'stack-type-deriver)

(defun vstack (ms type)
  "Concatenate matrices vertically. The matrices must have the same
number of columns and the same element type @c(type)."
  (unless ms
    (error "The list is empty"))
  (let* ((m (car ms))
         (ncols (array-dimension m 1))
         (nrows (reduce #'+ ms
                        :key (lambda (m) (array-dimension m 0))
                        :initial-value 0))
         (result (make-array (list nrows ncols) :element-type type))
         (offset 0))
    (declare (type fixnum offset))
    (loop for m in ms do
      (unless (and (= (array-dimension m 1) ncols)
                   (equalp (array-element-type m) type))
        (error "Cannot concatenate matrices"))
      (loop for i below (array-total-size m) do
        (setf (row-major-aref result (+ offset i))
              (row-major-aref m i)))
      (incf offset (array-total-size m)))
    result))

(defun hstack (ms type)
  "Concatenate matrices horizontally. The matrices must have the same
number of rows and the same element type @c(type)."
  (unless ms
    (error "The list is empty"))
  (let* ((m (car ms))
         (nrows (array-dimension m 0))
         (ncols (reduce #'+ ms
                        :key (lambda (m) (array-dimension m 1))
                        :initial-value 0))
         (result (make-array (list nrows ncols) :element-type type))
         (offset 0))
    (declare (type fixnum offset))
    (loop for m in ms do
      (unless (and (= (array-dimension m 0) nrows)
                   (equalp (array-element-type m) type))
        (error "Cannot concatenate matrices"))
      (loop for i below nrows do
        (loop for j below (array-dimension m 1) do
          (setf (aref result i (+ j offset))
                (aref m i j))))
      (incf offset (array-dimension m 1)))
    result))
