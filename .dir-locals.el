((nil . ((eval . (setq flycheck-gfortran-include-path
											 (list (file-name-directory
															(let ((d (dir-locals-find-file ".")))
																(if (stringp d)
																		d
																	(car d)))))))
				 (eval . (setq flycheck-gfortran-language-standard "f2003")))))
