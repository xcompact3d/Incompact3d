macro(AddExample example_dir case list_files)
  install(DIRECTORY DESTINATION ${example_dir}/${case})
  set(local_list "")
  list(APPEND local_list ${list_files})
  foreach(ff IN LISTS local_list)
    install(FILES ${ff} DESTINATION ${example_dir}/${case})
  endforeach()
endmacro()

