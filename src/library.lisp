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
       `(with-pointer-to-vector-data (,ptr (sb-ext:array-storage-vector ,array))
          ,acc)))
   bindings
   :initial-value `(progn ,@body)
   :from-end t))

;; Utility function to tell OpenMP not to use all available CPU resources
(cffi:defcfun ("omp_set_num_threads" %set-num-threads) :void
  (num-threads :int))

(declaim (inline set-num-threads))
(defun set-num-threads (n)
  "Set a number of threads for OpenBLAS."
  (%set-num-threads n))

;; Useful types
(deftype mat (type) `(simple-array ,type 2))
(deftype vec (type) `(simple-array ,type 1))
(deftype mat-or-vec (type) `(or (mat ,type) (vec ,type)))

;; And another useful function
(declaim (inline copy-array))
(defun copy-array (a)
  (let ((result (make-array (array-dimensions a)
                            :element-type (array-element-type a))))
    (replace (sb-ext:array-storage-vector result)
             (sb-ext:array-storage-vector a))
    result))

;; And another one
(declaim (inline map-array))
(defun map-array (f a)
  (let ((result (make-array (array-dimensions a)
                            :element-type (array-element-type a))))
    (map-into (sb-ext:array-storage-vector result) f (sb-ext:array-storage-vector a))
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

(define-condition lapack-error (error)
  ((message :type          string
            :initarg       :message
            :reader        error-message
            :documentation "A message")
   (info    :type          (or integer null)
            :initform      nil
            :initarg       :info
            :reader        error-info
            :documentation "A code returned by LAPACK"))
  (:report
   (lambda (c s)
     (format s "LAPACK error: ~a " (error-message c))
     (let ((info (error-info c)))
       (when info
         (format s "INFO: ~d" info)))))
  (:documentation "An error condition which is signaled with a
descriptive message when LAPACK fails to do its job."))
