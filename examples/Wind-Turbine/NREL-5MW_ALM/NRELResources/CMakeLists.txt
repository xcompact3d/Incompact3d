set(dir "Wind-Turbine-NREL-5MW/NRELResources")
file(GLOB_RECURSE example_files "${CMAKE_CURRENT_SOURCE_DIR}/*.al")
AddExample("${examples_dir}" "${dir}" "${example_files}")
add_subdirectory(Airfoils)
