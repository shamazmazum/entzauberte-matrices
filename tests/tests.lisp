(in-package :entzauberte-matrices/tests)

(def-suite algebra :description "Algebraic operations")
(def-suite det     :description "Determinant")
(def-suite inv     :description "Inversion")

(defun run-tests ()
  (em:set-num-threads 4)
  (every #'identity
         (mapcar (lambda (suite)
                   (let ((status (run suite)))
                     (explain! status)
                     (results-status status)))
                 '(algebra det))))

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

(declaim (inline random-vector))
(defun random-vector (n type)
  (make-array n
              :element-type type
              :initial-contents
              (loop repeat n
                    collect
                    (if (listp type)
                        (complex
                         (random (coerce 1 (second type)))
                         (random (coerce 1 (second type))))
                        (random (coerce 1 type))))))

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

(defun add (a b)
  (assert (equalp (array-dimensions a) (array-dimensions b)))
  (let ((result (make-array (array-dimensions a))))
    (loop for i below (array-total-size a) do
          (setf (row-major-aref result i)
                (+ (row-major-aref a i)
                   (row-major-aref b i))))
    result))

(defun scale (v s)
  (let ((result (make-array (array-dimensions v))))
    (loop for i below (array-total-size v) do
      (setf (row-major-aref result i)
            (* (row-major-aref v i) s)))
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

(test add-matrices
  (loop repeat 100
        for n = (+ (random 100) 20)
        for m = (+ (random 100) 20)
        for m1 = (random-matrix m n 'single-float)
        for m2 = (random-matrix m n 'single-float) do
          (is-true (array-approx-p (add m1 m2) (em:add m1 m2)))))

(test add-vectors
  (loop repeat 100
        for n = (+ (random 100) 20)
        for v1 = (random-vector n 'single-float)
        for v2 = (random-vector n 'single-float) do
          (is-true (array-approx-p (add v1 v2) (em:add v1 v2)))))

(test scale-matrix
  (loop repeat 100
        for n = (+ (random 100) 20)
        for m = (+ (random 100) 20)
        for mat = (random-matrix m n 'single-float)
        for s = (random 1.0) do
          (is-true (array-approx-p (scale mat s) (em:scale mat s)))))

(test scale-vector
  (loop repeat 100
        for n = (+ (random 100) 20)
        for v = (random-vector n 'single-float)
        for s = (random 1.0) do
          (is-true (array-approx-p (scale v s) (em:scale v s)))))

(in-suite det)

(defun permutations (n)
  (let (perm)
    (alexandria:map-permutations
     (lambda (p) (push p perm))
     (loop for i below n collect i))
    perm))

(defun det (m)
  (assert (= (array-dimension m 0)
             (array-dimension m 1)))
  (let ((length (array-dimension m 0)))
    (loop for %perm in (permutations length)
          for perm = (make-array length
                                :element-type '(unsigned-byte 32)
                                :initial-contents %perm)
          sum
          (* (expt -1 (em:inversions perm))
             (let ((d 1))
               (loop for i below length do
                 (setq d (* d (aref m i (aref perm i)))))
               d)))))

(macrolet ((def-det-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "DET/COMPLEX-~a" (second type))
                              (format nil "DET/~a" type)))))
               `(test ,name
                  (loop repeat 400
                        for n = (+ (random 6) 2)
                        for a = (random-matrix n n ',type) do
                          (is (approxp (det a) (em:det a)
                                       :rtol
                                       (/ (coerce 100
                                                  ',(if (listp type)
                                                        (second type)
                                                        type))))))))))
  (def-det-test single-float)
  (def-det-test double-float)
  (def-det-test (complex single-float))
  (def-det-test (complex double-float)))

(in-suite inv)

(declaim (inline make-identity))
(defun make-identity (n type)
  (let ((id (make-array (list n n)
                        :element-type type
                        :initial-element (coerce 0 type))))
    (loop for i below n do
      (setf (aref id i i) (coerce 1 type)))
    id))

(macrolet ((def-inv-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "INV/COMPLEX-~a" (second type))
                              (format nil "INV/~a" type)))))
               `(test ,name
                  (loop repeat 400
                        for n   = (+ (random 6) 2)
                        for a   = (random-matrix n n ',type)
                        for inv = (em:invert a)
                        for id  = (em:mult a inv) do
                          (is-true (array-approx-p id (make-identity n ',type))))))))
  (def-inv-test single-float)
  (def-inv-test double-float)
  (def-inv-test (complex single-float))
  (def-inv-test (complex double-float)))
