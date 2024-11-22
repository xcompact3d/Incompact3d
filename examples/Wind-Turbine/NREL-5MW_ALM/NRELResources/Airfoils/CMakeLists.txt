set(dir "Wind-Turbine-NREL-5MW/NRELResources/Airfoils")
file(GLOB_RECURSE example_files "${CMAKE_CURRENT_SOURCE_DIR}/*.air")
AddExample("${examples_dir}" "${dir}" "${example_files}")
