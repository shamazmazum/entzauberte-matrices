(in-package :entzauberte-matrices/tests)

(def-suite dot     :description "Dot product")
(def-suite algebra :description "Algebraic operations")
(def-suite det     :description "Determinant")
(def-suite inv     :description "Inversion")
(def-suite eig     :description "Eigenvalue decomposition")
(def-suite svd     :description "SVD decomposition")
(def-suite solver  :description "Solver")

(defun run-tests ()
  (em:set-num-threads 4)
  (every #'identity
         (mapcar (lambda (suite)
                   (let ((status (run suite)))
                     (explain! status)
                     (results-status status)))
                 '(dot algebra det inv eig svd solver))))

(declaim (inline conjugate-transpose))
(defun conjugate-transpose (m)
  (let ((res (make-array (reverse (array-dimensions m))
                         :element-type (array-element-type m))))
    (loop for i below (array-dimension m 0) do
      (loop for j below (array-dimension m 1) do
        (setf (aref res j i) (conjugate (aref m i j)))))
    res))

(declaim (inline random-number))
(defun random-number (type &optional imagpart-zero)
  (if (listp type)
      (complex
       (random (coerce 1 (second type)))
       (if imagpart-zero
           (coerce 0 (second type))
           (random (coerce 1 (second type)))))
      (random (coerce 1 type))))

(declaim (inline random-matrix))
(defun random-matrix (m n type)
  (make-array (list m n)
              :element-type type
              :initial-contents
              (loop repeat m collect
                    (loop repeat n collect (random-number type)))))

(declaim (inline random-self-adjoint))
(defun random-self-adjoint (n type)
  (let ((result (make-array (list n n)
                            :element-type type)))
    (loop for i below n do
      (loop for j from i below n do
        (if (<= i j)
            (let ((v (random-number type (= i j))))
              (setf (aref result i j) v
                    (aref result j i) (conjugate v))))))
    result))

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

(in-suite dot)

(declaim (inline dot))
(defun dot (v1 v2)
  (assert (= (length v1) (length v2)))
  (loop for i below (length v1) sum
        (* (aref v1 i) (conjugate (aref v2 i)))))

(macrolet ((def-dot-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "DOT/COMPLEX-~a" (second type))
                              (format nil "DOT/~a" type)))))
               `(test ,name
                  (loop repeat 1000
                        for n  = (+ (random 100) 20)
                        for v1 = (random-vector n ',type)
                        for v2 = (random-vector n ',type) do
                          (is (approxp (dot v1 v2) (em:dot v1 v2))))))))
  (def-dot-test single-float)
  (def-dot-test double-float))

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

(defun add (a b s)
  (assert (equalp (array-dimensions a) (array-dimensions b)))
  (let ((result (make-array (array-dimensions a))))
    (loop for i below (array-total-size a) do
          (setf (row-major-aref result i)
                (+ (* (row-major-aref a i) s)
                   (* (row-major-aref b i)))))
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
                          (is-true (array-approx-p (matrix-mul (em:transpose at) b)
                                                   (em:mult at b :ta t)))
                          (is-true (array-approx-p (matrix-mul a (em:transpose bt))
                                                   (em:mult a bt :tb t)))
                          (is-true (array-approx-p (matrix-mul (em:transpose at)
                                                               (em:transpose bt))
                                   (em:mult at bt :ta t :tb t))))))))
  (def-multiplication-test single-float)
  (def-multiplication-test double-float)
  (def-multiplication-test (complex single-float))
  (def-multiplication-test (complex double-float)))

(macrolet ((def-addition-test (type)
             (let* ((complexp (listp type))
                    (name (intern
                           (if complexp
                               (format nil "ADDITION/COMPLEX-~a" (second type))
                               (format nil "ADDITION/~a" type)))))
               `(test ,name
                  (loop repeat 100
                        for n = (+ (random 100) 20)
                        for m = (+ (random 100) 20)
                        for m1 = (random-matrix m n ',type)
                        for m2 = (random-matrix m n ',type)
                        for s  = (random-number ',type) do
                          (is-true (array-approx-p (add m1 m2 s) (em:add m1 m2 s))))))))
  (def-addition-test single-float)
  (def-addition-test double-float)
  (def-addition-test (complex single-float))
  (def-addition-test (complex double-float)))

(macrolet ((def-scale-test (type)
             (let* ((complexp (listp type))
                    (name (intern
                           (if complexp
                               (format nil "SCALE/COMPLEX-~a" (second type))
                               (format nil "SCALE/~a" type)))))
               `(test ,name
                  (loop repeat 100
                        for n = (+ (random 100) 20)
                        for m = (+ (random 100) 20)
                        for mat = (random-matrix m n ',type)
                        for s = (random-number ',type) do
                          (is-true (array-approx-p (scale mat s) (em:scale mat s))))))))
  (def-scale-test single-float)
  (def-scale-test double-float)
  (def-scale-test (complex single-float))
  (def-scale-test (complex double-float)))

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
                  (loop repeat 2000
                        for n  = (+ (random 6) 2)
                        for a  = (random-matrix n n ',type)
                        for d1 = (det a)
                        for d2 = (em:det a)
                        when (> (abs d2) 0.1) do
                          (is (approxp d1 d2)))))))
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
                  (loop repeat 2000
                        for n   = (+ (random 30) 2)
                        for a   = (random-matrix n n ',type)
                        for det = (em:det a)
                        when (> (abs det) 0.1) do
                        (let* ((inv (em:invert a))
                               (id  (em:mult a inv)))
                          (is-true (array-approx-p id (make-identity n ',type)))))))))
  (def-inv-test single-float)
  (def-inv-test double-float)
  (def-inv-test (complex single-float))
  (def-inv-test (complex double-float)))

(in-suite eig)

(defun multiply-eig (%t Λ)
  (let ((result (make-array (array-dimensions %t)
                            :element-type (array-element-type %t))))
    (loop for i below (array-dimension result 0) do
      (loop for j below (array-dimension result 1) do
        (setf (aref result i j)
              (* (aref %t i j)
                 (aref Λ i)))))
    result))

(defun convert-to-complex (a)
  (let ((type (array-element-type a)))
    (if (listp type) a
        (let ((result (make-array (array-dimensions a)
                                  :element-type (list 'complex type))))
          (map-into (sb-ext:array-storage-vector result)
                    #'complex
                    (sb-ext:array-storage-vector a))
          result))))

(macrolet ((def-eig-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "EIG-HERM/COMPLEX-~a" (second type))
                              (format nil "EIG-SYM/~a" type)))))
               `(test ,name
                  (loop repeat 400
                        for n = (+ (random 100) 2)
                        for a = (random-self-adjoint n ',type) do
                          (flet ((check (Λ %t)
                                   (is-true (array-approx-p
                                             (em:mult %t a)
                                             (multiply-eig %t Λ)))))
                            (multiple-value-call #'check
                              (em:eig-self-adjoint a :upper))
                            (multiple-value-call #'check
                              (em:eig-self-adjoint a :lower))))))))
  (def-eig-test single-float)
  (def-eig-test double-float)
  (def-eig-test (complex single-float))
  (def-eig-test (complex double-float)))

(macrolet ((def-eig-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "EIG/COMPLEX-~a" (second type))
                              (format nil "EIG/~a" type)))))
               `(test ,name
                  (loop repeat 400
                        for n = (+ (random 100) 2)
                        for a = (random-matrix n n ',type) do
                          (multiple-value-bind (Λ %T)
                              (em:eig a)
                            (is-true (array-approx-p
                                      (em:mult %T (convert-to-complex a))
                                      (multiply-eig %T Λ)))))))))
  (def-eig-test single-float)
  (def-eig-test double-float)
  (def-eig-test (complex single-float))
  (def-eig-test (complex double-float)))

(in-suite svd)

(declaim (inline from-diag))
(defun from-diag (d m n)
  (assert (= (length d) (min m n)))
  (let* ((length (length d))
         (type (array-element-type d))
         (result (make-array (list m n)
                             :element-type type
                             :initial-element (coerce 0 type))))
    (loop for i below length do
      (setf (aref result i i) (aref d i)))
    result))

(macrolet ((def-svd-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "SVD/COMPLEX-~a" (second type))
                              (format nil "SVD/~a" type))))
                   (complexp (listp type)))
               `(test ,name
                  (loop repeat 400
                        for m   = (+ (random 100) 2)
                        for n   = (+ (random 100) 2)
                        for min = (min m n)
                        for a   = (random-matrix m n ',type) do
                          (flet ((%test (compactp u s vt)
                                   (let ((s (if compactp
                                                (from-diag s min min)
                                                (from-diag s m n))))
                                     (is-true (array-approx-p
                                               (em:mult (conjugate-transpose u) u)
                                               (make-identity (if compactp min m) ',type)))
                                     (is-true (array-approx-p
                                               (em:mult vt (conjugate-transpose vt))
                                               (make-identity (if compactp min n) ',type)))
                                     (is-true (array-approx-p
                                               (em:mult
                                                (em:mult
                                                 u ,(if complexp '(convert-to-complex s) 's))
                                                vt)
                                               a)))))
                            (multiple-value-call #'%test nil (em:svd a :compact nil))
                            (multiple-value-call #'%test   t (em:svd a :compact t))))))))
  (def-svd-test single-float)
  (def-svd-test double-float)
  (def-svd-test (complex single-float))
  (def-svd-test (complex double-float)))

(in-suite solver)

(macrolet ((def-solver-test (type)
             (let ((name (intern
                          (if (listp type)
                              (format nil "SOLVER/COMPLEX-~a" (second type))
                              (format nil "SOLVER/~a" type)))))
               `(test ,name
                  (loop repeat 2000
                        for n    = (+ (random 30) 2)
                        for cols = (+ (random 30) 2)
                        for a    = (random-matrix n n ',type)
                        for det  = (em:det a)
                        when (> (abs det) 1f-1) do
                          (let* ((b  (random-matrix n cols ',type))
                                 (x  (em:solve a b))
                                 (%b (em:mult a x)))
                            (is-true (array-approx-p b %b))))))))
  (def-solver-test single-float)
  (def-solver-test double-float)
  (def-solver-test (complex single-float))
  (def-solver-test (complex double-float)))
