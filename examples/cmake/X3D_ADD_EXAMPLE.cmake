macro(AddExample root_dir case list_files)
  install(DIRECTORY DESTINATION ${root_dir}/${case})
  set(local_list "")
  list(APPEND local_list ${list_files})
  foreach(ff IN LISTS local_list)
    install(FILES ${ff} DESTINATION ${root_dir}/${case})
  endforeach()
endmacro()

