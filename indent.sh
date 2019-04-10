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
