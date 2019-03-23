# list of files
# new files need to be added in this list

##################################
# file list used for libs building
##################################

file(GLOB LIBTB2FILE
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
##  ${My_Source}/globals/tb2linearconstr.*pp
##  ${My_Source}/globals/tb2lpsconstr.*pp
##  ${My_Source}/globals/tb2mipsolver.*pp
    ${My_Source}/incop/*.h
    ${My_Source}/incop/*.cpp 
    ${My_Source}/search/*.*pp
    ${My_Source}/utils/*.*pp
    ${My_Source}/vns/*.*pp 
    ${My_Source}/ToulbarVersion.*pp
)

file(GLOB LIBCTB2FILE
    ${My_Source}/toulbar2_cstub.cpp
)

################
# toulbar2 files
################
file(GLOB source_files
    ${LIBTB2FILE}
    ${My_Source}/tb2main.cpp
)

MESSAGE(STATUS "------------------------------\n")
MESSAGE(STATUS "toulbar file = ${source_files}\n")
MESSAGE(STATUS "------------------------------\n")


