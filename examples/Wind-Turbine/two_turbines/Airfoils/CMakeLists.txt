set(dir "Wind-Turbine-two-turbines/Airfoils")
file(GLOB_RECURSE example_files "${CMAKE_CURRENT_SOURCE_DIR}/*.air")
AddExample("${examples_dir}" "${dir}" "${example_files}")
