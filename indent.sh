#!/bin/bash

#
#        FILE: indent.sh
#      AUTHOR: Paul Bartholomew
# DESCRIPTION: Use emacs to indent files
#

# decomp2d/
for f in decomp2d/*.f90
do
	echo "Indenting ${f}"
	emacs -batch ${f} --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
done
for f in decomp2d/*.inc
do
	echo "Indenting ${f}"
	emacs -batch ${f} --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
done

# src/
for f in src/*.f90
do
	echo "Indenting ${f}"
	emacs -batch ${f} --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
done

# Make sure
# no .eq., use ==
# no .ne., use /=
# no .ge., use >=
# no .le., use <=
# no .gt., use >
# no .le., use <
for f in decomp2d/*f90 decomp2d/*.inc src/*f90
do
	sed -i 's!\.eq\.! == !g' $f
	sed -i 's!\.ne\.! \/= !g' $f
	sed -i 's!\.ge\.! >= !g' $f
	sed -i 's!\.le\.! <= !g' $f
	sed -i 's!\.gt\.! > !g' $f
	sed -i 's!\.lt\.! < !g' $f
done
