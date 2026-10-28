file(REMOVE_RECURSE
  "lib/Linux/libtb2-archive.a"
  "lib/Linux/libtb2-archive.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/tb2-archive.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
