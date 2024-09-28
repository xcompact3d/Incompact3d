# ABL
set(case "Wind-Turbine")
set(example_files "")
list(APPEND example_files "NREL-5MW_ALM")
list(APPEND example_files "precursor-succesor-ADM")
AddExample("${examples_dir}" "${case}" "${example_files}")


