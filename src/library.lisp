(in-package :entzauberte-matrices)

(define-foreign-library openblas
  (:unix  (:or "libopenblas.so"
               "libopenblas.so.0"))
  (t (:default "libopenblas")))

(use-foreign-library openblas)

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

;; Utility function to tell OpenBLAS not to use all available CPU
;; resources.
(cffi:defcfun ("goto_set_num_threads" %goto-set-num-threads) :void
  (num-threads :int))
(cffi:defcfun ("openblas_set_num_threads" %openblas-set-num-threads) :void
  (num-threads :int))
(cffi:defcfun ("omp_set_num_threads" %omp-set-num-threads) :void
  (num-threads :int))

(declaim (inline set-num-threads))
(defun set-num-threads (n)
  "Set a number of threads for OpenBLAS."
  (%openblas-set-num-threads n)
  (%goto-set-num-threads     n)
  (when (foreign-symbol-pointer "omp_set_num_threads")
    (%omp-set-num-threads    n)))

;; Useful types
(deftype smat (type) `(simple-array ,type 2))
(deftype svec (type) `(simple-array ,type 1))
(deftype smat-or-svec (type) `(or (smat ,type) (svec ,type)))

(deftype mat (type) `(array ,type 2))
(deftype vec (type) `(array ,type 1))
(deftype mat-or-vec (type) `(or (mat ,type) (vec ,type)))

;; And another useful function
(declaim (inline copy-array))
(defun copy-array (a)
  (let ((result (make-array (array-dimensions a)
                            :element-type (array-element-type a))))
    (replace (array-storage-vector result)
             (array-storage-vector a))
    result))

;; And another one
(declaim (inline map-array))
(defun map-array (f a)
  (let ((result (make-array (array-dimensions a)
                            :element-type (array-element-type a))))
    (map-into (array-storage-vector result) f (array-storage-vector a))
    result))

;; Transform LAPACK "pivot" indices to "normal" permutation matrix indices.
(serapeum:-> fix-pivot ((simple-array (unsigned-byte 32) (*)))
             (values (simple-array (unsigned-byte 32) (*)) &optional))
(defun fix-pivot (indices)
  (declare (optimize (speed 3)))
  (let* ((length (length indices))
         (new-indices (make-array length :element-type '(unsigned-byte 32))))
    (loop for i below length do
      (setf (aref new-indices i) i))
    (loop for i below length
          ;; Fortran starts indexing from 1
          for idx = (1- (aref indices i))
          when (/= idx i) do
            (rotatef (aref new-indices i)
                     (aref new-indices idx)))
    new-indices))

;; Function for generating wrapper names
(serapeum:-> wrapper-names (symbol)
             (values symbol string &optional))
(defun wrapper-names (name)
  (let ((name (symbol-name name)))
    (values
     (intern (format nil "%~a" name))
     (format nil "~a_" (string-downcase name)))))
