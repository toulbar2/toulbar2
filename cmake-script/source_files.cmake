
# list of files
# new files need to be added to this list
# you can also define your own list and add it to the whole list

#########################
# file list used for library building
########################


file(
  GLOB
  LIBTB2FILES
  ${My_Source}/applis/*.*pp
  ${My_Source}/core/*.*pp
  ${My_Source}/cpd/*.*pp
  ${My_Source}/globals/tb2graph.*pp
  ${My_Source}/globals/tb2globalconstr.*pp
  ${My_Source}/globals/tb2flowbasedconstr.*pp
  ${My_Source}/globals/tb2alldiffconstr.*pp
  ${My_Source}/globals/tb2globalcardinalityconstr.*pp
  ${My_Source}/globals/tb2sameconstr.*pp
  ${My_Source}/globals/tb2regularflowconstr.*pp
  ${My_Source}/globals/tb2dpglobalconstr.*pp
  ${My_Source}/globals/tb2regulardpconstr.*pp
  ${My_Source}/globals/tb2amongconstr.*pp
  ${My_Source}/globals/tb2grammarconstr.*pp
  ${My_Source}/globals/tb2grammarutils.*pp
  ${My_Source}/globals/tb2maxconstr.*pp
  ${My_Source}/globals/tb2treeconstr.*pp
  #                        ${My_Source}/globals/tb2linearconstr.*pp
  #                        ${My_Source}/globals/tb2lpsconstr.*pp
  #                        ${My_Source}/globals/tb2mipsolver.*pp
  ${My_Source}/incop/*.cpp 
  ${My_Source}/search/*.*pp
  ${My_Source}/utils/*.*pp
  ${My_Source}/vns/*.*pp 
  ${My_Source}/ToulbarVersion.*pp
  )

set(source_files ${LIBTB2FILES} ${My_Source}/tb2main.cpp)
set(pysource_files ${LIBTB2FILES} ${My_Source}/pytoulbar2.cpp)

MESSAGE(STATUS "-------------\n")
MESSAGE(STATUS "toulbar2 source files = ${source_files}\n")
MESSAGE(STATUS "-------------\n")

