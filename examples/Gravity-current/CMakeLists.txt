# Gravity-current
set(case "Gravity-current")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_gravitycurrent_3DRe2236_LR.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


