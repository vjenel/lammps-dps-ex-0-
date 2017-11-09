An example of an implementation of a pair_style in lammps for the potential:
U(r)=A*exp(-b*r)-C/r^6+D/r^12 + U_q(r) [ with U_q(r) = qi*qj/r) ]  (1)

The pair-style is defined as: 
pair_style exp612/cut/coul/long
and the pair coefficients are declared as:
pair_coeff i j A b C D 

The modified parts of the original code (starting from the files: “pair_lj_cut_coul_cut.cpp” and “pair_lj_cut_coul_cut.h”) were commented 
and the new added/edited code [files “pair_exp612_cut_coul_long.cpp” and “pair_exp612_cut_coul_long.h”] was labeled as “added” (to be easier to follow). 

To use it: copy the two files “pair_exp612_cut_coul_long.cpp” and “pair_exp612_cut_coul_long.h” in the src of lammps and recompile. 

I tested it on a smaller system(s) from refs DOI:10.1021/jp301399b (given in the tables from SI) and DOI 10.1021/acs.jpclett.5b02378 with modifications for testing the 14 reductions. 
The input file in.lammps.overlay runs the system with overlay pairstyle and it requires no modification to lammps code.
The input file in.lammps will use the modified code for the implemented pairstyle exp612/cut/coul/long .
