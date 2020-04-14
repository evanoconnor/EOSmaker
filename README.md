!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!                   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!   EOSMaker v1.0   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Copyright Evan O'Connor and C.D. Ott, October 2011 !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!                   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!                   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DISCLAIMER: Please note that the routines and tables provided here
come with absolutely no warranty and we are unable to guarantee that
we will be able to provide support or help if you run into problems
integrating them with your simulation code. If you decide to use the
provided EOS routines and data in published work, it is YOUR
responsibility to check their physical correctness and consistency.


 If you have any questions or have discovered a bug in our
routines/tables, please e-mail us at evanoc@tapir.caltech.edu and/or
cott@tapir.caltech.edu.  The routines are finicky, especially with
never before tested EOS tables, expect runtime errors, hopefully the
error messages will be sufficient to track down issues. In particular,
the fix_compositions.F90 routine.  There is no good way to fix the
compositions besides doing the original calculation at each grid
point.

These routines are tested with the eos2.tab from H. Shen et al. 2011
and the NL3eosb1.03.dat from G. Shen et al. 2011.  You must copy the
read_table_GShen_NL3.F90 file to read_table.F90 and the
eos_table_GShen_NL3.F90 file to eos_table.F90 to compile are create
the table for G. Shen et al. 2011, the make file will yell at you if
you do not do this.

The layout of the program in detailed in the comments of
eos_table.F90.  This file and possibly read_eos_table.F90 should be
the only files you need to change, unless you add something like
exotic particles or something.

You must have HDF5 compiled with the _same_ compiler you compile the
EOSmaker with.  This usually means downloading the source from
http://www.hdfgroup.org/HDF5/release/obtain5.html, configuring with
your version of:

./configure --enable-fortran FC=ifort --prefix=/Users/evanoc/opt/hdf5-current-ifort12

and then 

make
make install

the HDF5DIR variable in make.inc would then be set to /Users/evanoc/opt/hdf5-current-ifort12

After this, a simple make should create the eos_table executables in
the main directory.

