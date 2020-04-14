!-*-f90-*-
subroutine setup_final_table

  use eos_table_module
  implicit none

  integer i,j,k !local counters

  real*8 dlogrho,dlogtemp,dye

  !allocate arrays
  allocate(logrho(nlogrho))
  allocate(logtemp(nlogtemp))
  allocate(ye(nye))
  
  !initialize
  logrho(:) = 0.0d0
  logtemp(:) = 0.0d0
  ye(:) = 0.0d0

  dlogrho = (logrhomax-logrhomin)/(nlogrho-1)
  dlogtemp = (logtempmax-logtempmin)/(nlogtemp-1)
  dye = (yemax-yemin)/(nye-1)

  write(*,"(a,1P10E18.9)") "table spacings",dlogrho,dlogtemp,dye

  !now allocate final eos table
  !nuclear part
  allocate(eos_table(nlogrho,nlogtemp,nye,final_nvars))
     
  !electrons
  allocate(eos_electrons(nlogrho,nlogtemp,nye,n_electron_vars))
  allocate(eos_photons(nlogrho,nlogtemp,nye,n_electron_vars))

  !initialize
  eos_table(:,:,:,:) = 0.0d0
  eos_table(:,:,:,imasseff) = lowdensity_effectivemass
  eos_electrons = 0.0d0
  eos_photons = 0.0d0

  !now figure out which parts of the table fall into which category,
  !ie. set the transition indices indices
  
  irho_start_comps = -1
  irho_start_nuclear = -1
  do i=1,nlogrho
     !also use this opportunity to set final table densities
     logrho(i) = logrhomin + dlogrho*dble(i-1)
     if(logrho(i).ge.logrhomin_comps.and.irho_start_comps.eq.-1) then
        irho_start_comps = i
     endif
     if(logrho(i).ge.logrhomin_nuclear.and.irho_start_nuclear.eq.-1) then
        irho_start_nuclear = i
     endif
  enddo
  
  do j=1,nlogtemp

  enddo

  itemp_start_nuclear = -1
  itemp_stop_nuclear = -1
  do j=1,nlogtemp
     !set final table temperatures
     logtemp(j) = logtempmin + dlogtemp*dble(j-1)
     if(logtemp(j).ge.logtempmin_nuclear.and.itemp_start_nuclear.eq.-1) then
        itemp_start_nuclear = j
     endif
     if(logtemp(j).gt.logtempmax_nuclear.and.itemp_stop_nuclear.eq.-1) then
        itemp_stop_nuclear = j-1 !meaning the last index was the last
                                 !one that could be calculated from
                                 !the table
     endif
  enddo

  if (itemp_stop_nuclear.eq.-1) then
     !if we do not need to go to higher temps this will be the case,
     !then set itemp_stop_nuclear to nlogtemp, i.e. nlogtemp is the
     !last point to be calculated from the table (and the last point!)
     itemp_stop_nuclear = nlogtemp
  endif

  write(6,*) "irho_start_comps: ",irho_start_comps,logrho(irho_start_comps)
  write(6,*) "irho_start_nuclear: ",irho_start_nuclear,logrho(irho_start_nuclear)
  write(6,*) "itemp_start_nuclear: ",itemp_start_nuclear,logtemp(itemp_start_nuclear)
  write(6,*) "itemp_stop_nuclear: ",itemp_stop_nuclear,logtemp(itemp_stop_nuclear)

  do k=1,nye
     ye(k) = yemin + dye*dble(k-1)
  enddo
  
  !check settings for Hempel, Hyperon
  if (Hempel) then
     if (.not.Hempel_massfrac_reconstruction) then
        write(*,*) "Probably should be using Hempel Reconstruction"
        stop
     endif
     if (lump_into_np) then
        !this must be done in reading in, then set Hempel to
        !false to proceed as normal, hempel reconstruction will
        !proceed as normal
     endif
  else if (Hyperon) then
     if (.not.Hempel_massfrac_reconstruction) then
        write(*,*) "Probably should be using Hempel Reconstruction"
        stop
     endif
     if (lump_into_np) then
        !this will be done in reading in, then Hempel will be set to
        !false to proceed as normal, hempel reconstruction will
        !proceed as normal
     endif
  else
     if (lump_into_np) then
        write(*,*) "must be using Hempel EOS of Hyperons to lump lights into n & p", &
             " if another EOS has lights, things will need to be adjusted accordingly"
        stop
     endif
     if (Hempel_massfrac_reconstruction) then
        write(*,*) "Proceed with caution regarding mass fraction reconstruction", &
             " it should be good, but untested"
        stop
     endif
  endif



end subroutine setup_final_table
