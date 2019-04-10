((nil . ((eval . (setq flycheck-fortran-gfortran-executable "mpif90"))
				 (eval . (setf flycheck-fortran-args "-fcray-pointer -cpp"))
				 (eval . (setq flycheck-gfortran-include-path
											 ;; Find this file and use it as the project root directory.
											 (list (file-name-directory
															(let ((d (dir-locals-find-file ".")))
																(if (stringp d)
																		d
																	(car d)))))))
				 (eval . (setq flycheck-gfortran-language-standard "f2003")))))
