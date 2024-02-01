First, you need to compile the .c file, with the following command : 

R CMD SHLIB lqmm.so lqmm_test.c initè_test.c

This will create a .so file which is then loaded using dyn.load in R.

The file fonctions.R contains all functions needed to run the method ( called with lqmmtrue ) . fonctions.R has to be called at the beginning of the code with source().

The file HLQMM_simu_conv.R shows the simulation coding from the article presenting HLQMM.
