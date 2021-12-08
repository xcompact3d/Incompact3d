# Add the x3d unit testing library
add_subdirectory(lib)

# Add x3d verification (derive and filter subroutines)
add_subdirectory(verif)

# Add x3d mesh convergence
add_subdirectory(convergence)

# Add 2D TGV test case
add_subdirectory(Taylor-Green-Vortex-2D)
