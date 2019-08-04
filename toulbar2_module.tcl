#%Module1.0#################################################################

# description : module to load environnment on SLURM cluster
# kad version 1.0
# usage: module load -f path_to/my_own_module
# To unload one module: module unload bioinfo/bowtie2-2.2.9
# To unload all module and specific variable: module purge

#limit coredumpsize 0
module purge
module load compiler/gcc-7.2.0
module load mpi/openmpi-2.1.2
module load compiler/cmake-3.12.3

setenv BOOST_ROOT /tools/libraries/Boost/boost_1_70_0_openmpi-2.1.2
setenv Boost_INCLUDE_DIR /tools/libraries/Boost/boost_1_70_0_openmpi-2.1.2/include
setenv Boost_LIBRARY_DIR /tools/libraries/Boost/boost_1_70_0_openmpi-2.1.2/lib

prepend-path	 LD_LIBRARY_PATH  /tools/libraries/Boost/boost_1_70_0_openmpi-2.1.2/lib

#export BOOST_ROOT=${chemin}
#export Boost_LIBRARY_DIRS=${chemin}/lib/
#export Boost_INCLUDE_DIRS=${chemin}/include/
#setenv		 MPIHOME /tools/cluster/mpi/openmpi/2.1.2/gcc-4.8.5 
#setenv		 MPI openmpi 

#prepend-path	 PATH /tools/cluster/mpi/openmpi/2.1.2/gcc-4.8.5/bin 
#prepend-path	 LD_LIBRARY_PATH /tools/cluster/mpi/openmpi/2.1.2/gcc-4.8.5/lib 
#prepend-path	 MANPATH /tools/cluster/mpi/openmpi/2.1.2/gcc-4.8.5/share/man 


#The Boost C++ Libraries were successfully built!
#The following directory should be added to compiler include paths:
#    /home/abeldjilali/Bureau/boost_1_70_0
#The following directory should be added to linker library paths:
#    /home/abeldjilali/Bureau/boost_1_70_0/stage/lib

