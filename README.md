First, you need to compile the .c file, with the following comman : 

R CMD SHLIB lqmm.so lqmm.c init.c

This will create a .so file which is then loaded using dyn.load in R.

The file fonctions.R contains all fonctions needed to run the method ( called with lqmmtrue ) . fonctions.R has to be called in the beginning of the code with source().

The file HLQMM_simu_conv.R shows the coding of simulation from the article presenting HLQMM.
