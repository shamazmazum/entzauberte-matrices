(in-package :entzauberte-matrices)

(deftype transpose-legacy ()
    "Control matrix transposition in @c(mult). Values @c(T) and
@c(:trans) mean that the corresponding matrix is
transposed. @c(:no-trans) means the matrix remains as is,
@c(:conj-trans) means Hermitean conjugate and @c(:conj-no-trans) means
that elements of the matrix become complex adjoints."
  '(or (eql t) transpose))

(deftype transpose ()
  "Same as @c(transpose-legacy) but without @c(T)."
  '(member t :no-trans :trans :conj-trans :conj-no-trans))

(serapeum:-> transpose-required-p (transpose)
             (values boolean &optional))
(declaim (inline transpose-required-p))
(defun transpose-required-p (transpose)
  (if (member transpose '(:trans :conj-trans)) t))

(serapeum:-> transpose-to-keyword (transpose-legacy)
             (values transpose &optional))
(declaim (inline transpose-to-keyword))
(defun transpose-to-keyword (transpose)
  (if (eq transpose t) :trans transpose))

(macrolet ((define-foreign-mult (name foreign-type)
             (multiple-value-bind (lisp-name fortran-name)
                 (capi-wrapper-names name :blas)
               `(defcfun (,lisp-name ,fortran-name) :void
                  (order  cblas-order)
                  (transa cblas-transpose)
                  (transb cblas-transpose)
                  (m      blas-int)
                  (n      blas-int)
                  (k      blas-int)
                  (alpha  ,foreign-type)
                  (a       :pointer)
                  (lda    blas-int)
                  (b      :pointer)
                  (ldb    blas-int)
                  (beta   ,foreign-type)
                  (c      :pointer)
                  (ldc    blas-int)))))
  ;; Real
  (define-foreign-mult sgemm :float)
  (define-foreign-mult dgemm :double)
  ;; Complex
  (define-foreign-mult cgemm (:pointer :float))
  (define-foreign-mult zgemm (:pointer :double)))

(macrolet ((define-lisp-mult (high-level-name low-level-name lisp-type &optional foreign-type)
             (let* ((complexp foreign-type)
                    (float-type (if complexp (second lisp-type) lisp-type)))
               `(progn
                  (serapeum:-> ,high-level-name ((mat ,lisp-type)
                                                 (mat ,lisp-type)
                                                 transpose transpose ,lisp-type)
                               (values (smat ,lisp-type) &optional))
                  (defun ,high-level-name (a b ta tb scale)
                    (let* ((tap (transpose-required-p ta))
                           (tbp (transpose-required-p tb))
                           (c (make-array (list (array-dimension a (if tap 1 0))
                                                (array-dimension b (if tbp 0 1)))
                                          :element-type ',lisp-type)))
                      (with-array-pointers ((aptr a)
                                            (bptr b)
                                            (cptr c))
                        (with-foreign-objects ,(if complexp
                                                   `((alphaptr ,foreign-type 2)
                                                     (betaptr  ,foreign-type 2)))
                          (let ((m (array-dimension a (if tap 1 0)))
                                (n (array-dimension b (if tbp 0 1)))
                                (k (array-dimension a (if tap 0 1)))
                                (lda (array-dimension a 1))
                                (ldb (array-dimension b 1))
                                (ldc (array-dimension c 1))
                                (alpha ,(if complexp `alphaptr `scale))
                                (beta  ,(if complexp `betaptr  `(coerce 0 ',float-type))))
                            ,@(when complexp
                                `((setf (mem-aref alphaptr ,foreign-type 0)
                                        (realpart scale)
                                        (mem-aref alphaptr ,foreign-type 1)
                                        (imagpart scale)
                                        (mem-aref betaptr ,foreign-type 0)
                                        (coerce 0 ',float-type)
                                        (mem-aref betaptr ,foreign-type 1)
                                        (coerce 0 ',float-type))))
                            (,low-level-name
                             :row-major
                             ta tb m n k alpha aptr lda bptr ldb beta cptr ldc)))
                        c)))))))
  (define-lisp-mult mult-rs-unsafe %sgemm single-float)
  (define-lisp-mult mult-rd-unsafe %dgemm double-float)
  (define-lisp-mult mult-cs-unsafe %cgemm (complex single-float) :float)
  (define-lisp-mult mult-cd-unsafe %zgemm (complex double-float) :double))

(serapeum:-> mult ((mat *) (mat *)
                   &key (:ta transpose-legacy) (:tb transpose-legacy) (:scale number))
             (values (smat *) &optional))
(declaim (inline mult))
(defun mult (a b &key (ta :no-trans) (tb :no-trans) (scale (coerce 1 (array-element-type a))))
  "Compute a product of two matrices \\(s A^* B^*\\) where \\(A^*\\)
and \\(B^*\\) are derived from \\(A\\) and \\(B\\) accordingly to
@c(:ta) and @c(:tb) arguments of type @c(transpose) and \\(s\\) is a
scalar @c(scale)."
  (let* ((ta  (transpose-to-keyword ta))
         (tb  (transpose-to-keyword tb))
         (tap (transpose-required-p ta))
         (tbp (transpose-required-p tb)))
    (unless (= (array-dimension a (if tap 0 1))
               (array-dimension b (if tbp 1 0)))
      (error "Cannot multiply: matrices have incompatible dimensions"))
    (let ((t1 (array-element-type a))
          (t2 (array-element-type b)))
      (funcall
       (cond
         ((and (eq t1 'single-float)
               (eq t2 'single-float))
          #'mult-rs-unsafe)
         ((and (eq t1 'double-float)
               (eq t2 'double-float))
          #'mult-rd-unsafe)
         ((and (equalp t1 '(complex single-float))
               (equalp t2 '(complex single-float)))
          #'mult-cs-unsafe)
         ((and (equalp t1 '(complex double-float))
               (equalp t2 '(complex double-float)))
          #'mult-cd-unsafe)
         (t (error "Cannot multiply: incompatible element types: ~a and ~a" t1 t2)))
       a b ta tb scale))))

(defun @ (&rest ms)
  "Left associative multiplication of matrices."
  (reduce #'mult ms))

(defun @-expander (ms)
  (reduce
   (lambda (acc m)
     `(mult ,acc ,m))
   ms))

(define-compiler-macro @ (&rest ms)
  (@-expander ms))
