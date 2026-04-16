(in-package :entzauberte-matrices)

#+bsd
(cc-flags "-I/usr/local/include")
(include "cblas.h")
(include "lapacke.h")

(ctype blas-int   "blasint")
(ctype lapack-int "lapack_int")

(cenum cblas-order
  ((:row-major "CblasRowMajor"))
  ((:col-major "CblasColMajor")))

(constantenum lapack-order
  ((:row-major "LAPACK_ROW_MAJOR"))
  ((:col-major "LAPACK_COL_MAJOR")))

(cenum cblas-transpose
  ((:no-trans      "CblasNoTrans"))
  ((:trans         "CblasTrans"))
  ((:conj-trans    "CblasConjTrans"))
  ((:conj-no-trans "CblasConjNoTrans")))

(cenum cblas-uplo
  ((:upper "CblasUpper"))
  ((:lower "CblasLower")))

(cenum cblas-diag
  ((:non-unit "CblasNonUnit"))
  ((:unit     "CblasUnit")))

(cenum cblas-side
  ((:left  "CblasLeft"))
  ((:right "CblasRight")))
