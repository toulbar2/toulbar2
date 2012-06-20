
# liste of file 
# new file news to be add in this list
# you can also define your own list and add it to the wall list
#  for exmape : set ( ${WALLFILE} ${my_file_list_to_add})

        file(
                        GLOB
                        source_files
                        ${My_Source}/tb2abstractconstr.cpp
                        ${My_Source}/tb2arithmetic.cpp
			${My_Source}/tb2automaton.cpp
                        ${My_Source}/tb2bep.cpp
                        ${My_Source}/tb2binconstr.cpp
                        ${My_Source}/tb2btd.cpp
                        ${My_Source}/tb2btlist.cpp
                        ${My_Source}/tb2btqueue.cpp
                        ${My_Source}/tb2clusters.cpp
                        ${My_Source}/tb2constraint.cpp
                        ${My_Source}/tb2graph.cpp
                        ${My_Source}/tb2globaldecomposable.cpp
                        ${My_Source}/tb2globalconstr.cpp
                        ${My_Source}/tb2flowbasedconstr.cpp
                        ${My_Source}/tb2alldiffconstr.cpp
                        ${My_Source}/tb2globalcardinalityconstr.cpp
                        ${My_Source}/tb2sameconstr.cpp
                        ${My_Source}/tb2regularconstr.cpp
                        ${My_Source}/tb2domain.cpp
                        ${My_Source}/tb2enumvar.cpp
                        ${My_Source}/tb2intervar.cpp
                        ${My_Source}/tb2main.cpp
                        ${My_Source}/tb2naryconstr.cpp
                        ${My_Source}/tb2pedigree.cpp
                        ${My_Source}/tb2haplotype.cpp
                        ${My_Source}/tb2queue.cpp
                        ${My_Source}/tb2randomgen.cpp
                        ${My_Source}/tb2reader.cpp
                        ${My_Source}/tb2solver.cpp
                        ${My_Source}/tb2system.cpp
                        ${My_Source}/tb2ternaryconstr.cpp
                        ${My_Source}/tb2vac.cpp
                        ${My_Source}/tb2vacutils.cpp
                        ${My_Source}/tb2variable.cpp
                        ${My_Source}/tb2wcsp.cpp
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
                        ${My_Source}/tb2abstractconstr.*pp 
                        ${My_Source}/tb2arithmetic.*pp 
                        ${My_Source}/tb2bep.*pp 
                        ${My_Source}/tb2binconstr.*pp 
                        ${My_Source}/tb2btlist.*pp
                        ${My_Source}/tb2btqueue.*pp
                        ${My_Source}/tb2clusters.*pp
                        ${My_Source}/tb2constraint.*pp
                        ${My_Source}/tb2graph.*pp
                        ${My_Source}/tb2globalconstr.*pp
                        ${My_Source}/tb2flowbasedconstr.*pp
                        ${My_Source}/tb2alldiffconstr.*pp
                        ${My_Source}/tb2globalcardinalityconstr.*pp
                        ${My_Source}/tb2sameconstr.*pp
                        ${My_Source}/tb2regularconstr.*pp
                        ${My_Source}/tb2domain.*pp
                        ${My_Source}/tb2enumvar.*pp
                        ${My_Source}/tb2intervar.*pp
                        ${My_Source}/tb2naryconstr.*pp
                        ${My_Source}/tb2pedigree.*pp
                        ${My_Source}/tb2haplotype.cpp
                        ${My_Source}/tb2queue.*pp
                        ${My_Source}/tb2randomgen.*pp
                        ${My_Source}/tb2reader.*pp
                        ${My_Source}/tb2system.*pp
                        ${My_Source}/tb2ternaryconstr.*pp
                        ${My_Source}/tb2vac.*pp
                        ${My_Source}/tb2vacutils.*pp
                        ${My_Source}/tb2variable.*pp
                        ${My_Source}/tb2wcsp.*pp
        )


