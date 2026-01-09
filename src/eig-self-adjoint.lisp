(in-package :entzauberte-matrices)

(deftype uplo () '(member :upper :lower))

(macrolet ((def-self-adjoint-low-level (name complexp)
             (multiple-value-bind (lisp-name foreign-name)
                 (wrapper-names name)
               `(defcfun (,lisp-name ,foreign-name) :void
                  (jobz  (:pointer :char))
                  (uplo  (:pointer :char))
                  (n     (:pointer :int))
                  (a     :pointer)
                  (lda   (:pointer :int))
                  (w     :pointer)
                  (work  :pointer)
                  (lwork (:pointer :int))
                  ,@(if complexp '((rwork :pointer)))
                  (info  (:pointer :int))))))
  (def-self-adjoint-low-level ssyev nil)
  (def-self-adjoint-low-level dsyev nil)
  (def-self-adjoint-low-level cheev t)
  (def-self-adjoint-low-level zheev t))

;; TODO: Do not forget to document that the eigenvalues are in ROWS
;; (so-called left eignedvectors)!
(macrolet ((def-self-adjoint-lisp (name low-level-fn lisp-type foreign-type)
             `(progn
                (serapeum:-> ,name ((mat ,lisp-type) uplo)
                             (values (svec ,lisp-type)
                                     (smat ,lisp-type)
                                     &optional))
                (defun ,name (a where)
                  (let* ((n (array-dimension a 0))
                         (lda n)
                         (copy (copy-array a))
                         (vals (make-array n :element-type ',lisp-type)))
                    (with-array-pointers ((aptr copy)
                                          (wptr vals))
                      (with-foreign-objects ((jobzptr :char)
                                             (uploptr :char)
                                             (nptr    :int)
                                             (ldaptr  :int)
                                             (infoptr :int))
                        (flet ((check-info ()
                                 (let ((info (mem-ref infoptr :int)))
                                   (unless (zerop info)
                                     (error 'lapack-error
                                            :message "Cannot compute eigenvalues"
                                            :info info)))))
                          (setf (mem-ref jobzptr :char)
                                (char-code #\V)
                                (mem-ref uploptr :char)
                                (char-code
                                 (if (eq where :upper) #\L #\U))
                                (mem-ref nptr   :int) n
                                (mem-ref ldaptr :int) lda)
                          (let ((work
                                  (with-foreign-objects ((workptr  ,foreign-type)
                                                         (lworkptr :int))
                                    (setf (mem-ref lworkptr :int) -1)
                                    (,low-level-fn jobzptr uploptr nptr aptr ldaptr wptr
                                                   workptr lworkptr infoptr)
                                    (check-info)
                                    (round (mem-ref workptr ,foreign-type)))))
                            (with-foreign-objects ((workptr  ,foreign-type work)
                                                   (lworkptr :int))
                              (setf (mem-ref lworkptr :int) work)
                              (,low-level-fn jobzptr uploptr nptr aptr ldaptr wptr
                                             workptr lworkptr infoptr)
                              (check-info))))))
                    (values vals copy))))))
  (def-self-adjoint-lisp eig-self-adjoint-rs-unsafe %ssyev single-float :float)
  (def-self-adjoint-lisp eig-self-adjoint-rd-unsafe %dsyev double-float :double))


(macrolet ((def-self-adjoint-lisp (name low-level-fn lisp-type foreign-type)
             (let ((real-type (second lisp-type)))
               `(progn
                  (serapeum:-> ,name ((mat ,lisp-type) uplo)
                               (values (svec ,real-type)
                                       (smat ,lisp-type)
                                       &optional))
                  (defun ,name (a where)
                    (let* ((n (array-dimension a 0))
                           (lda n)
                           ;; Make conjugate transpose column-major matrix
                           (copy (map-array #'conjugate a))
                           (vals (make-array n :element-type ',real-type)))
                      (with-array-pointers ((aptr copy)
                                            (wptr vals))
                        (with-foreign-objects ((jobzptr  :char)
                                               (uploptr  :char)
                                               (nptr     :int)
                                               (ldaptr   :int)
                                               (infoptr  :int)
                                               (rworkptr ,foreign-type (- (* 3 n) 2)))
                          (flet ((check-info ()
                                   (let ((info (mem-ref infoptr :int)))
                                     (unless (zerop info)
                                       (error 'lapack-error
                                              :message "Cannot compute eigenvalues"
                                              :info info)))))
                            (setf (mem-ref jobzptr :char)
                                  (char-code #\V)
                                  (mem-ref uploptr :char)
                                  (char-code
                                   (if (eq where :upper) #\L #\U))
                                  (mem-ref nptr   :int) n
                                  (mem-ref ldaptr :int) lda)
                            (let ((work
                                    (with-foreign-objects ((workptr  ,foreign-type 2)
                                                           (lworkptr :int))
                                      (setf (mem-ref lworkptr :int) -1)
                                      (,low-level-fn jobzptr uploptr nptr aptr ldaptr wptr
                                                     workptr lworkptr rworkptr infoptr)
                                      (check-info)
                                      (round (mem-ref workptr ,foreign-type)))))
                              (with-foreign-objects ((workptr  ,foreign-type (* work 2))
                                                     (lworkptr :int))
                                (setf (mem-ref lworkptr :int) work)
                                (,low-level-fn jobzptr uploptr nptr aptr ldaptr wptr
                                               workptr lworkptr rworkptr infoptr)
                                (check-info))))))
                      ;; Conjugate back
                      (values vals (map-array #'conjugate copy))))))))
  (def-self-adjoint-lisp eig-self-adjoint-cs-unsafe %cheev (complex single-float) :float)
  (def-self-adjoint-lisp eig-self-adjoint-cd-unsafe %zheev (complex double-float) :double))

(serapeum:-> eig-self-adjoint ((mat *) uplo)
             (values (svec *) (smat *) &optional))
(declaim (inline eig-self-adjoint))
(defun eig-self-adjoint (m where)
  "Compute eigenvalues (the first returned value) and eigenvectors
(the second returned value) of a self-adjoint matrix
\\(M\\). Eigenvectors are stored @b(in rows), so that (in pseudo-code)

@begin[lang=lisp](code)
(multiple-value-bind (vals vecs)
    (eig m)
  (approx= (mult vecs m)
           (mult (from-diag vals) vecs)))
@end(code)

is @c(T)."
  (assert (= (array-dimension m 0)
             (array-dimension m 1)))
  (funcall
   (cond
     ((eq (array-element-type m) 'single-float)
      #'eig-self-adjoint-rs-unsafe)
     ((eq (array-element-type m) 'double-float)
      #'eig-self-adjoint-rd-unsafe)
     ((equalp (array-element-type m) '(complex single-float))
      #'eig-self-adjoint-cs-unsafe)
     ((equalp (array-element-type m) '(complex double-float))
      #'eig-self-adjoint-cd-unsafe)
     (t
      (error "Cannot compute eigenvalues: unknown array element type")))
   m where))
