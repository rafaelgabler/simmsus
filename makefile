#-----------------------------------------------------#
# 	COMPILATION INSTRUCTIONS FOR SIMMSUS 	      #
#                                                     #
#               Written initially by       	      #
#               Rafael Gabler Gontijo                 #
#                                                     #
# Goal: to compile a given version of SIMMSUS	      #
#-----------------------------------------------------#
#
# Code source files (Linux)
#
SRC-LNX = subroutines.f90 variables.f90 input.f90  \
          main.f90 statistics.f90 simmsus.f90  \
        
 OBJ-LNX = $(SRC-LNX:.f90=.o)

# Compiler definition and flags

ENGINE-LNX = ifort
FLAGS-LNX = -m64 -O2 #-qopenmp

# Executable file generation (Linux)

simmsus.ex : $(SRC-LNX) 
	$(ENGINE-LNX) $(FLAGS-LNX) -o simmsus.ex $(SRC-LNX)

# Cleaning command
clean :
	rm simmsus.ex
	
	
