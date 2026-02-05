(defun do-all()
  (handler-case
      (asdf:load-system :entzauberte-matrices/tests)
    (error ()
      (uiop:quit 1)))
  (uiop:quit
   (if (uiop:call-function "entzauberte-matrices/tests:run-tests")
       0 1)))

(do-all)
