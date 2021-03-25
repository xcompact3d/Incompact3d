;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((f90-mode
  (fortan-comment-indent-style . 'relative)
  (flycheck-gfortran-language-standard . "f2003")
  (flycheck-gfortran-args . ("-fcray-pointer" "-cpp"))
  (flycheck-fortran-gfortran-executable . "mpif90")))
