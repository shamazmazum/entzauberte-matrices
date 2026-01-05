(defsystem :entzauberte-matrices
    :version "0.1"
    :author "Vasily Postnicov <shamaz.mazum@gmail.com>"
    :description "Fast wrapper around LAPACK"
    :license "2-clause BSD"
    :serial t
    :pathname "src"
    :components ((:file "package")
                 (:file "library")
                 (:file "lisp-only")
                 (:file "mult")
                 (:file "lu")
                 (:file "det")
                 (:file "inv")
                 (:file "eig-self-adjoint")
                 (:file "eig-ge"))
    :depends-on (:serapeum :cffi)
    :in-order-to ((test-op (load-op "entzauberte-matrices/tests")))
    :perform (test-op (op system)
                      (declare (ignore op system))
                      (uiop:symbol-call :entzauberte-matrices/tests '#:run-tests)))

(defsystem :entzauberte-matrices/tests
    :version "0.1"
    :author "Vasily Postnicov <shamaz.mazum@gmail.com>"
    :license "2-clause BSD"
    :pathname "tests"
    :components ((:file "package")
                 (:file "tests" :depends-on ("package")))
    :depends-on (:entzauberte-matrices
                 :alexandria
                 :fiveam
                 :approx))
