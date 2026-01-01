(in-package :entzauberte-matrices)

(serapeum:-> inversions ((simple-array (unsigned-byte 32) (*)))
             (values unsigned-byte &optional))
(defun inversions (ipiv)
  "Return the number of inversions in permutation array"
  (declare (optimize (speed 3)))
  (let ((result 0)
        (length (length ipiv)))
    (loop for i below length
          for n = (aref ipiv i) do
            (loop for j from (1+ i) below length
                  for m = (aref ipiv j) do
                    (when (< m n)
                      (incf result))))
    result))


(macrolet ((def-det (name type)
             `(progn
                (serapeum:-> ,name ((simple-array ,type 2))
                             (values ,type &optional))
                (defun ,name (m)
                  (multiple-value-bind (lu ipiv)
                      (%lu m)
                    (* (expt -1 (inversions ipiv))
                       (let ((det (coerce 1 ',type)))
                         ;; SBCL cannot infer the type of DET, hence all this macrolet shit
                         (declare (type ,type det))
                         (loop for i below (array-dimension m 0) do
                           (setq det (* det (aref lu i i))))
                         det)))))))
  (def-det det-rs-unsafe single-float)
  (def-det det-rd-unsafe double-float)
  (def-det det-cs-unsafe (complex single-float))
  (def-det det-cd-unsafe (complex double-float)))

(serapeum:-> det ((mat *))
             (values number &optional))
(declaim (inline det))
(defun det (m)
  (unless (= (array-dimension m 0)
             (array-dimension m 1))
    (error "The matrix is not square"))
  (cond
    ((eq (array-element-type m) 'single-float)
     (det-rs-unsafe m))
    ((eq (array-element-type m) 'double-float)
     (det-rd-unsafe m))
    ((equalp (array-element-type m) '(complex single-float))
     (det-cs-unsafe m))
    ((equalp (array-element-type m) '(complex double-float))
     (det-cd-unsafe m))
    (t (error "Cannot compute DET: Unknown matrix element type"))))
