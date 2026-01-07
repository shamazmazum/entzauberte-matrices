(in-package :entzauberte-matrices)

(macrolet ((def-eig-low-level (name complexp)
             (let* ((name (symbol-name name))
                    (lisp-name (intern (format nil "%~a" name)))
                    (foreign-name (format nil "~a_" (string-downcase name))))
               `(defcfun (,lisp-name ,foreign-name) :void
                  (jobvl (:pointer :char))
                  (jobvr (:pointer :char))
                  (n     (:pointer :int))
                  (a     :pointer)
                  (lda   (:pointer :int))
                  ,@(if complexp
                        `((w  :pointer))
                        `((wr :pointer)
                          (wi :pointer)))
                  (vl    :pointer)
                  (ldvl  (:pointer :int))
                  (vr    :pointer)
                  (ldvr  (:pointer :int))
                  (work  :pointer)
                  (lwork (:pointer :int))
                  ,@(if complexp `((rwork :pointer)))
                  (info  (:pointer :int))))))
  (def-eig-low-level sgeev nil)
  (def-eig-low-level dgeev nil)
  (def-eig-low-level cgeev t)
  (def-eig-low-level zgeev t))

;; Real
(macrolet ((def-eig (name low-level-fn lisp-type foreign-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type))
                             (values (vec ,lisp-type)
                                     (mat ,lisp-type)
                                     &optional))
                (defun ,name (a)
                  (let* ((n (array-dimension a 0))
                         (lda n)
                         (copy (copy-array a))
                         (vals    (make-array n :element-type ',lisp-type))
                         (vecs    (make-array (list n n) :element-type ',lisp-type)))
                    (with-foreign-objects ((jobvlptr :char)
                                           (jobvrptr :char)
                                           (nptr     :int)
                                           (ldaptr   :int)
                                           (ldvlptr  :int)
                                           (ldvrptr  :int)
                                           (rworkptr ,foreign-type (* n 2))
                                           (infoptr  :int))
                      (flet ((check-info ()
                               (unless (zerop (mem-ref infoptr :int))
                                 (error 'lapack-error
                                        :message "Cannot compute eigenvalues"))))
                        (setf (mem-ref jobvlptr :char)
                              ;; To be consistent with EIG-SELF-ADJOINT
                              (char-code #\N)
                              (mem-ref jobvrptr :char)
                              (char-code #\V)
                              (mem-ref nptr :int) n
                              (mem-ref ldaptr :int) lda
                              ;; FIXME ldvl >= 1
                              (mem-ref ldvlptr :int) 1
                              (mem-ref ldvrptr :int) n)
                        (with-array-pointers ((aptr  copy)
                                              (wptr  vals)
                                              (vrptr vecs))
                          (let ((work
                                  (with-foreign-objects ((workptr  ,foreign-type 2)
                                                         (lworkptr :int))
                                    (setf (mem-ref lworkptr :int) -1)
                                    (,low-level-fn jobvlptr jobvrptr nptr aptr ldaptr
                                                   wptr (null-pointer) ldvlptr vrptr
                                                   ldvrptr workptr lworkptr rworkptr
                                                   infoptr)
                                    (check-info)
                                    (round (mem-ref workptr ,foreign-type)))))
                            (with-foreign-objects ((workptr  ,foreign-type (* work 2))
                                                   (lworkptr :int))
                              (setf (mem-ref lworkptr :int) work)
                              (,low-level-fn jobvlptr jobvrptr nptr aptr ldaptr
                                                   wptr (null-pointer) ldvlptr vrptr
                                                   ldvrptr workptr lworkptr rworkptr
                                                   infoptr)
                              (check-info))))))
                    (values vals vecs))))))
  (def-eig eig-cs-unsafe %cgeev (complex single-float) :float)
  (def-eig eig-cd-unsafe %zgeev (complex double-float) :double))

;; Complex
(macrolet ((def-eig (name low-level-fn lisp-type foreign-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type))
                             (values (vec (complex ,lisp-type))
                                     (mat (complex ,lisp-type))
                                     &optional))
                (defun ,name (a)
                  (let* ((n (array-dimension a 0))
                         (lda n)
                         (copy (copy-array a))
                         (vals    (make-array n :element-type '(complex ,lisp-type)))
                         (vals-re (make-array n :element-type ',lisp-type))
                         (vals-im (make-array n :element-type ',lisp-type))
                         (vecs-re (make-array (list n n) :element-type ',lisp-type))
                         (vecs    (make-array (list n n) :element-type '(complex ,lisp-type))))
                    (with-foreign-objects ((jobvlptr :char)
                                           (jobvrptr :char)
                                           (nptr     :int)
                                           (ldaptr   :int)
                                           (ldvlptr  :int)
                                           (ldvrptr  :int)
                                           (infoptr  :int))
                      (flet ((check-info ()
                               (unless (zerop (mem-ref infoptr :int))
                                 (error 'lapack-error
                                        :message "Cannot compute eigenvalues"))))
                        (setf (mem-ref jobvlptr :char)
                              ;; To be consistent with EIG-SELF-ADJOINT
                              (char-code #\N)
                              (mem-ref jobvrptr :char)
                              (char-code #\V)
                              (mem-ref nptr :int) n
                              (mem-ref ldaptr :int) lda
                              ;; FIXME ldvl >= 1
                              (mem-ref ldvlptr :int) 1
                              (mem-ref ldvrptr :int) n)
                        (with-array-pointers ((aptr copy)
                                              (wrptr vals-re)
                                              (wiptr vals-im)
                                              (vrptr vecs-re))
                          (let ((work
                                  (with-foreign-objects ((workptr  ,foreign-type)
                                                         (lworkptr :int))
                                    (setf (mem-ref lworkptr :int) -1)
                                    (,low-level-fn jobvlptr jobvrptr nptr aptr ldaptr
                                                   wrptr wiptr (null-pointer)
                                                   ldvlptr vrptr ldvrptr workptr lworkptr
                                                   infoptr)
                                    (check-info)
                                    (round (mem-ref workptr ,foreign-type)))))
                            (with-foreign-objects ((workptr  ,foreign-type work)
                                                   (lworkptr :int))
                              (setf (mem-ref lworkptr :int) work)
                              (,low-level-fn jobvlptr jobvrptr nptr aptr ldaptr
                                             wrptr wiptr (null-pointer)
                                             ldvlptr vrptr ldvrptr workptr lworkptr infoptr)
                              (check-info))))))
                    (loop for i below n do
                      (setf (aref vals i)
                            (complex
                             (aref vals-re i)
                             (aref vals-im i))))
                    (loop for i below n
                          for im = (aref vals-im i) do
                      (cond
                        ((zerop im)
                         (loop for j below n do
                           (setf (aref vecs i j)
                                 (complex (aref vecs-re i j)))))
                        ((> im 0)
                         (loop for j below n do
                           (setf (aref vecs i j)
                                 (complex (aref vecs-re ( + i) j)
                                          (aref vecs-re (1+ i) j)))))
                        ((< im 0)
                         (loop for j below n do
                           (setf (aref vecs i j)
                                 (complex (+ (aref vecs-re (1- i) j))
                                          (- (aref vecs-re ( + i) j))))))))
                    (values vals vecs))))))
  (def-eig eig-rs-unsafe %sgeev single-float :float)
  (def-eig eig-rd-unsafe %dgeev double-float :double))

(serapeum:-> eig ((mat *))
             (values (vec *) (mat *) &optional))
(declaim (inline eig))
(defun eig (m)
  (assert (= (array-dimension m 0)
             (array-dimension m 1)))
  (cond
    ((eq (array-element-type m) 'single-float)
     (eig-rs-unsafe m))
     ((eq (array-element-type m) 'double-float)
      (eig-rd-unsafe m))
     ((equalp (array-element-type m) '(complex single-float))
      (eig-cs-unsafe m))
     ((equalp (array-element-type m) '(complex double-float))
      (eig-cd-unsafe m))
     (t
      (error "Cannot compute eigenvalues: unknown array element type"))))
