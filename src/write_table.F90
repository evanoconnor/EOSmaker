!-*-f90-*-
subroutine write_table(h5filename,timestamp)

  use hdf5
  use eos_table_module
  implicit none

  integer(HID_T) file_id,dset_id,dspace_id,aspace_id,attr_id
  integer(HID_T) atype_id
  integer(SIZE_T) attrlen

  integer(HSIZE_T) dims1(1), dims3(3)
  integer rank,error
  integer i

  character(len=256) h5filename

  real*8 logenergy(nlogrho,nlogtemp,nye)
  real*8 logpress(nlogrho,nlogtemp,nye)
  real*8 munul(nlogrho,nlogtemp,nye)
  real*8 muhat(nlogrho,nlogtemp,nye)
  real*8 xa(nlogrho,nlogtemp,nye)
  real*8 xn(nlogrho,nlogtemp,nye)
  real*8 xp(nlogrho,nlogtemp,nye)
  real*8 xh(nlogrho,nlogtemp,nye)
  real*8 xd(nlogrho,nlogtemp,nye)
  real*8 xt(nlogrho,nlogtemp,nye)
  real*8 x3he(nlogrho,nlogtemp,nye)
  real*8 x4li(nlogrho,nlogtemp,nye)
  real*8 xL(nlogrho,nlogtemp,nye)
  
  real*8 norm

  real*8 :: timestamp

  logpress(:,:,:) = eos_table(:,:,:,ipress)
  logenergy(:,:,:) = eos_table(:,:,:,ienergy)
  munul(:,:,:) = eos_table(:,:,:,imue) - eos_table(:,:,:,imun) &
       + eos_table(:,:,:,imup)

  muhat(:,:,:) = eos_table(:,:,:,imun) - eos_table(:,:,:,imup)

  xa = eos_table(:,:,:,ixa)
  xh = eos_table(:,:,:,ixh)
  xn = eos_table(:,:,:,ixn)
  xp = eos_table(:,:,:,ixp)
  if (Hempel) then
     xd = eos_table(:,:,:,ixd)
     xt = eos_table(:,:,:,ixt)
     x3he = eos_table(:,:,:,ix3he)
     x4li = eos_table(:,:,:,ix4li)
  endif
  if (Hyperon) then
     xL = eos_table(:,:,:,ixL)
  endif

  write(*,*) "Writing HDF5 EOS Table: ", trim(adjustl(h5filename)), "at time: ", timestamp

  call h5open_f(error)
  call h5fcreate_f(trim(adjustl(h5filename)), H5F_ACC_TRUNC_F, file_id, error)

  ! scalars go next
  rank = 1
  dims1(1) = 1
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "timestamp", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, timestamp, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  rank = 1
  dims1(1) = 1
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "pointsrho", H5T_NATIVE_INTEGER, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nlogrho, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)


  rank = 1
  dims1(1) = 1
  energy_shift = energy_shift/massn_cgs*mev_to_erg
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "energy_shift", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, energy_shift, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)


  rank = 1
  dims1(1) = 1
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "pointstemp", H5T_NATIVE_INTEGER, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nlogtemp, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  rank = 1
  dims1(1) = 1
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "pointsye", H5T_NATIVE_INTEGER, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nye, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

! 1D arrays
  
  rank=1
  dims1(1) = nlogrho
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "logrho", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, logrho, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  rank=1
  dims1(1) = nlogtemp
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "logtemp", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, logtemp, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  rank=1
  dims1(1) = nye
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "ye", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ye, &
       dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

! 3D arrays
  rank=3
  dims3(1) = nlogrho
  dims3(2) = nlogtemp
  dims3(3) = nye

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "logpress", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, logpress, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "logenergy", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, logenergy, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

#if 0
  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "epscgs", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, logenergy, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
#endif

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "entropy", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
       eos_table(:,:,:,ientropy),dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "dedt", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
       dedt,dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "mu_e", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_table(:,:,:,imue), &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "mu_p", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_table(:,:,:,imup), &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "mu_n", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_table(:,:,:,imun), &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "munu", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, munul, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "muhat", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, muhat, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "dpdrhoe", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dpdrho, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "dpderho", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dpde, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "gamma", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gamma, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "cs2", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cs2, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Xn", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xn, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Xp", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xp, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Xh", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xh, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Xa", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xa, &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)


  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Abar", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_table(:,:,:,iabar), &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Zbar", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_table(:,:,:,izbar), &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  call h5screate_simple_f(rank, dims3, dspace_id, error)
  call h5dcreate_f(file_id, "Meff", H5T_NATIVE_DOUBLE, &
       dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_table(:,:,:,imasseff), &
       dims3, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)

  if (Hempel) then
     call h5screate_simple_f(rank, dims3, dspace_id, error)
     call h5dcreate_f(file_id, "Xd", H5T_NATIVE_DOUBLE, &
          dspace_id, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xd, &
          dims3, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)

     call h5screate_simple_f(rank, dims3, dspace_id, error)
     call h5dcreate_f(file_id, "Xt", H5T_NATIVE_DOUBLE, &
          dspace_id, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xt, &
          dims3, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)

     call h5screate_simple_f(rank, dims3, dspace_id, error)
     call h5dcreate_f(file_id, "X3he", H5T_NATIVE_DOUBLE, &
          dspace_id, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x3he, &
          dims3, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)

     call h5screate_simple_f(rank, dims3, dspace_id, error)
     call h5dcreate_f(file_id, "X4li", H5T_NATIVE_DOUBLE, &
          dspace_id, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x4li, &
          dims3, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     
  endif

  if (Hyperon) then
     call h5screate_simple_f(rank, dims3, dspace_id, error)
     call h5dcreate_f(file_id, "XL", H5T_NATIVE_DOUBLE, &
          dspace_id, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xL, &
          dims3, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
  endif
     
  call h5fclose_f(file_id,error)
  call h5close_f(error)


  write(*,*) "Done! :-)"


end subroutine write_table
