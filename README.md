# Entzauberte-matrices
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://shamazmazum.github.io/entzauberte-matrices)

This is a tiny wrapper around LAPACK which uses OpenBLAS.

**For Ubuntu users:** Ubuntu provides OpenBLAS without LAPACKE functions. On
this linux distribution `liblapacke.so` is used which seems to be a slow
reference implementation from netlib.
