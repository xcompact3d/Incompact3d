Guidelines for Developers
==========================

Keep coding consistent!

--> Fortran 1995 or above should be used (so, no “.gt.” but “>”, for instance)

--> Indentation using fprettify (3 spaces)

--> Avoid “tabs” by any means, only spaces

--> Use a maximum of 80 characters per line

--> No withespace at the end of the lines

--> As Fortran is not case-sensitive, do not use capital letters (“function” rather than “Function”, because an error found by the compiler will return “function” in both cases

--> Avoid goto

--> Prefer predefined values, like zero, one, two, to 0.0_mytype, 1.0_mytype, 2.0_mytype

--> Variables should all be initialised before use

--> All allocated arrays should be de-allocated

--> All variable names should have a meaning in English (apart from i, j, k, l)
