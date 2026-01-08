(defun do-all()
  (ql:quickload :entzauberte-matrices/tests)
  (uiop:quit
   (if (uiop:call-function "entzauberte-matrices/tests:run-tests")
       0 1)))

(do-all)
