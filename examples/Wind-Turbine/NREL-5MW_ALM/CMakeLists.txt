# Wind-Turbine: NREL
set(case "Wind-Turbine-NREL-5MW")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "NREL-5MW.i3d")
list(APPEND example_files "NREL-5MW.turb")
AddExample("${examples_dir}" "${case}" "${example_files}")
add_subdirectory(NRELResources)
