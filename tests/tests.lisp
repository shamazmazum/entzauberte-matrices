(in-package :entzauberte-matrices/tests)

(def-suite algebra :description "Algebraic operations")

(defun run-tests ()
  (every #'identity
         (mapcar (lambda (suite)
                   (let ((status (run suite)))
                     (explain! status)
                     (results-status status)))
                 '(algebra))))

(declaim (inline transpose))
(defun transpose (m)
  (let ((res (make-array (reverse (array-dimensions m))
                         :element-type (array-element-type m))))
    (loop for i below (array-dimension m 0) do
      (loop for j below (array-dimension m 1) do
        (setf (aref res j i) (aref m i j))))
    res))

(declaim (inline random-matrix))
(defun random-matrix (m n type)
  (make-array (list m n)
              :element-type type
              :initial-contents
              (loop repeat m collect
                    (loop repeat n collect
                          (if (listp type)
                              (complex
                               (random (coerce 1 (second type)))
                               (random (coerce 1 (second type))))
                              (random (coerce 1 type)))))))

(in-suite algebra)

(defun matrix-mul (a b)
  (assert (= (array-dimension a 1) (array-dimension b 0)))
  (let ((result (make-array (list (array-dimension a 0)
                                  (array-dimension b 1)))))
    (loop for i below (array-dimension result 0) do
      (loop for j below (array-dimension result 1) do
        (setf (aref result i j)
              (loop for k below (array-dimension a 1) sum
                    (* (aref a i k) (aref b k j))))))
    result))

(macrolet ((def-multiplication-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "MULTIPLICATION/COMPLEX-~a" (second type))
                              (format nil "MULTIPLICATION/~a" type)))))
               `(test ,name
                  (loop repeat 100
                        for n = (+ (random 100) 20)
                        for m = (+ (random 100) 20)
                        for k = (+ (random 100) 20)
                        for a  = (random-matrix n k ',type)
                        for at = (random-matrix k n ',type)
                        for b  = (random-matrix k m ',type)
                        for bt = (random-matrix m k ',type) do
                          (is-true (array-approx-p (matrix-mul a b) (em:mult a b)))
                          (is-true (array-approx-p (matrix-mul (transpose at) b)
                                                   (em:mult at b :ta t)))
                          (is-true (array-approx-p (matrix-mul a (transpose bt))
                                                   (em:mult a bt :tb t)))
                          (is-true (array-approx-p (matrix-mul (transpose at) (transpose bt))
                                   (em:mult at bt :ta t :tb t))))))))
  (def-multiplication-test single-float)
  (def-multiplication-test double-float)
  (def-multiplication-test (complex single-float))
  (def-multiplication-test (complex double-float)))
