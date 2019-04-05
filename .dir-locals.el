((nil . ((eval . (setq flycheck-gfortran-include-path
			(file-name-directory
			 (let ((d (dir-locals-find-file ".")))
			   (if (stringp d)
			       d
			     (car d)))))))))
