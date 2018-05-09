
# liste of file 
# new file news to be add in this list
# you can also define your own list and add it to the wall list
#  for exmape : set ( ${WALLFILE} ${my_file_list_to_add})

        file(
                        GLOB
                        source_files
                        ${My_Source}/applis/*.cpp
                        ${My_Source}/core/*.cpp
                        ${My_Source}/globals/tb2graph.cpp
                        ${My_Source}/globals/tb2globalconstr.cpp
                        ${My_Source}/globals/tb2flowbasedconstr.cpp
                        ${My_Source}/globals/tb2alldiffconstr.cpp
                        ${My_Source}/globals/tb2globalcardinalityconstr.cpp
                        ${My_Source}/globals/tb2sameconstr.cpp
                        ${My_Source}/globals/tb2regularflowconstr.cpp
                        ${My_Source}/globals/tb2dpglobalconstr.cpp
                        ${My_Source}/globals/tb2regulardpconstr.cpp
                        ${My_Source}/globals/tb2amongconstr.cpp
                        ${My_Source}/globals/tb2grammarconstr.cpp
                        ${My_Source}/globals/tb2grammarutils.cpp
                        ${My_Source}/globals/tb2maxconstr.cpp
                        ${My_Source}/globals/tb2treeconstr.cpp
##                        ${My_Source}/globals/tb2linearconstr.cpp
##                        ${My_Source}/globals/tb2lpsconstr.cpp
##                        ${My_Source}/globals/tb2mipsolver.cpp
                        ${My_Source}/incop/*.cpp
                        ${My_Source}/search/*.cpp
                        ${My_Source}/utils/*.cpp
                        ${My_Source}/vns/*.cpp
                        ${My_Source}/tb2main.cpp
                        ${My_Source}/ToulbarVersion.cpp
        )

MESSAGE(STATUS "-------------\n")
MESSAGE(STATUS "toulbar file = ${source_files}\n")
MESSAGE(STATUS "-------------\n")

#########################
# file list used for lib building
########################


        file(
                        GLOB
                        LIBTB2FILE
                        ${My_Source}/applis/*.*pp
                        ${My_Source}/core/*.*pp
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
##                        ${My_Source}/globals/tb2linearconstr.*pp
##                        ${My_Source}/globals/tb2lpsconstr.*pp
##                        ${My_Source}/globals/tb2mipsolver.*pp
                        ${My_Source}/incop/*.h
                        ${My_Source}/incop/*.cpp 
                        ${My_Source}/search/*.*pp
                        ${My_Source}/utils/*.*pp
                        ${My_Source}/vns/*.*pp 
                        ${My_Source}/ToulbarVersion.*pp
        )


