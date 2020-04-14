!-*-f90-*-
program eos_table_maker

  use eos_table_module
  implicit none
  integer irho,itemp,iye !table counters
  integer i,j,k !other counters
  real*8 temp
  
  real*8 :: min_press,min_energy ! we take logs of pressure and energy, use these along the way...
  integer :: imin_energy(3) !location of minimum energy, once found
  character(256) :: initial_eos_table_name
  character(256) :: h5filename

  !extra pressure
  logical :: add_extra_pressure
  real*8 :: extrapressure_factor
  real*8 :: baryon_dens
  
  !versioning
  real*8 :: timestamp
  character(8) :: date
  integer :: values(8)
  character(100) :: base,vnum,srho,stemp,sye

  !!!!!!!!!!!!!!!!!!The following are EOS specific, update for your EOS
  !!!!!!!!!!!!!!!!!!It is currently formatted for H. Shen 2010 no leptons, no lambdas

  initial_eos_table_name = "./original_tables/FSU1.7eosb1.01.dat"

  nt_temp = 109 !number of temperature entries in initial table                                                               
  nt_yp = 53 !number of yp entries in initial table                                                                          
  nt_rho = 336 !number of density entries in initial table                                                                    
  nvars = 13 !columns in the initial table - 3 (don't care about ye
             !and rho etc., know those)

  !size you want to final table in
  !denstiy
  nlogrho = 280 !number of density bins
  logrhomin = 3.0d0 !global minimum of final table
  logrhomax = 15.39d0 !global maximum of final table

  !temperature
  nlogtemp = 180 !number of temperature bins
  logtempmin = -2.0d0
  logtempmax = 2.5d0
  
  !ye
  nye = 52 !number of ye bins
  yemin = 0.05d0
  yemax = 0.56d0

  m_ref = 939.0d0 !reference mass of free nucleons
  baryon_to_gram_conversion = m_ref*mev_to_gram !mass unit of EOS

  base="GShenFSU_1.7"
  vnum="_1.1"

  final_nvars = 14

  if (final_nvars.lt.nvars+1) stop "final_nvars must be at least nvars + mue variable"

  !specify nasty transitions........  at the following density we stop
  !using the thermodynamic results of the nuclear table and transition
  !to an ion eos.  The determine, try to find a density where the
  !pressures of the two EOS match up nicely and smoothly.  We require
  !so we also shift the energies at this density for every temperature
  !and ye to ensure continuity, this is described later on
  logrhomin_nuclear = 7.25d0
  
  ! We use the compositions of the nuclear table down to this log
  ! density...  Well out of NSE.... They do not do much
  logrhomin_comps = 7.25d0  
    
  !The nuclear table may not go low enough in temperature, here is the
  !transition temp of the table, below we extrapolate or use ion eos
  logtempmin_nuclear= -0.799d0

  !Or not high enough, we like black holes, so does Shen so we are
  !good here, GShen doesn't
  logtempmax_nuclear = 1.875d0

  !set column numbers of initial table.  We assume the first three
  !columns are not going to be inlcuded (i.e. they are density, ye,
  !and temp, or something similar.  Adjust read_iniital_data if
  !required
  ifree    = 1 !we don't use this
  ienergy  = 13 !specific internal energy code assumes units of
                !MeV/baryon, GShen is special, they don't give E so
                !I'll overwrite meff after reading and set it to F-TS
  ientropy = 3 !in units of /k_B
  iabar    = 7 !mass number of heavies, average
  izbar    = 8 !atomic number of heavies, average
  imasseff = 13 !we don't use this, see above
  ixn      = 9 !mass fraction of neutrons
  ixp      = 10 !mass fraction of protons
  ixa      = 11 !mass fraction of alphas
  ixh      = 12 !mass fraction of heavies
  ipress   = 2 !press assume units of MeV/fm^3

  !Note, must be careful regarding rest mass in/out of chemical
  !potential.  We assume incoming chemical potentials DO NOT include
  !rest mass difference.  This is true for H. Shen and G. Shen where
  !the mu's are defined relative to the nucleon mass M, which is the
  !same for both neutron and proton
  imun     = 4 !chemical potential of neutron, MeV
  imup     = 5 !chemical potential of proton, MeV
  imue     = 6 !chemical potential of electron, MeV, includes rest
               !mass, not part of H. Shen initial table but part of
               !G. Shen and part of final

  use_saved_electrons = .false.


  !first read in data
  call read_initial_table(initial_eos_table_name)

  add_extra_pressure = .false.
  extrapressure_factor = 2.0d-5
  
  if (add_extra_pressure) then
     do i=1,nt_temp
        do j=1,nt_yp
           do k=1,nt_rho
              baryon_dens = 10.0d0**table_logrho(k)/baryon_to_gram_conversion/1.0d39
              if (baryon_dens.gt.0.2d0) then
                 table_data(i,j,k,ipress) = table_data(i,j,k,ipress) + &
                      197.326938d0**3*extrapressure_factor*(baryon_dens**2-0.2d0**2)
                 table_data(i,j,k,ienergy) = table_data(i,j,k,ienergy) + &
                      197.326938d0**3*extrapressure_factor*(baryon_dens-2.0d0*0.2d0+0.2d0**2/baryon_dens)                      

              endif
           enddo
        enddo
     enddo
   
  endif

  !only thing left for you to do is read fix_table_units... important
  !stuff in there and adjust reading file for specific table format

  !!!!!!!!!!!!!!!!!!!!Done EOS specific stuff

  !setup up output table, then our routines interpolate the initial
  !table to this table with cubic hermite interpolation (one of the
  !points of this program is to take a sparse table and do a good job
  !at interpolating to a finer table where only linear interpolation
  !is needed during a computational simulation)
  write(*,*) "Setting up final table"
  call setup_final_table

  !take a first stab at cubicly-hermitely interpolating the initial
  !table to the final table.  This routine also sets up the
  !compositions everywhere (including the loweos part) and takes care
  !of the low-T high density region.
  write(*,*) "Interpolating table"
  call interpolate_initial_table_cubic_hermite

  ! tweak things to ensure mass and charge conservation of nuclear
  ! composition (neutrons,protons,alphas,heavies), this is need
  ! because the original table was good to only 1 part in ~1000 and we
  ! need better then that for our electron eos.  Also, the
  ! interpolations we did would not necessarily keep the sum of the
  ! mass fractions = 1 or the sum of the protons to ye
  write(*,*) "Fixing Compositions"
  call fix_compositions

  ! some units need converting (i.e. pressure need to be converted
  ! from MeV/fm^-3 to dynes/cm^2).  But also we mess with the overall
  ! energy scale to get correct gravitational mass, by adding and
  ! subtracting different rest masses.  We also add in the rest mass
  ! difference into the chemical potential so that \mu_n - \mu_p
  ! includes rest mass difference between neutron and proton
  !!!!!This routine will also be EOS specific, please read comments
  !!!!!regarding energy scaling and update with your normalization
  write(*,*) "Fixing table units"
  call fix_table_units

  !this routine calculates electrons and photons everywhere and adds
  !them to high density EOS component (starting at
  !irho_start_nuclear+1, all elese gets done in loweos
  write(*,*) "Electron & Photon EOS with saved?:", use_saved_electrons
  call electron_photon_eos
  
  !now call the loweos to fill in stuff below irho_start_nuclear.
  !Will add in nuclear (i.e. ion) part and also electron and photon
  !part.
  write(*,*) "Low Density EOS"
  call loweos

  !ok we are almost done.  The final table will be in log energy and
  !log pressure.  By now the pressure should be positive everywhere so
  !there is no longer a need to shift it, the energy will not be > 0
  !everywhere.  We will calculate the minimum energy and shift
  !everything by 1.01*|min energy|.  Then we store the energy shift in
  !the final table so it can be properly handled

  write(*,*) "Energy Shift and table logging"
  min_energy = minval(eos_table(:,:,:,ienergy))
  imin_energy = minloc(eos_table(:,:,:,ienergy))
  energy_shift = -1.01d0*min(min_energy,0.0d0)*massn_cgs/mev_to_erg 
  write(*,*) "min_energy (MeV/baryon): ", min_energy*massn_cgs/mev_to_erg
  write(*,*) "energy shift (MeV/baryon): ", energy_shift
  write(*,*) "location, indices",imin_energy(1:3)
  write(*,*) "location, values", 10.0d0**logrho(imin_energy(1)),10.0d0**logtemp(imin_energy(2)),ye(imin_energy(3))

  ! shift energy to be compatible with log
  do k=1,nye
     do j=1,nlogtemp
        do i=1,nlogrho
           !convert to MeV/baryon
           eos_table(i,j,k,ienergy) = &
                eos_table(i,j,k,ienergy)*massn_cgs/mev_to_erg
           !add shift
           eos_table(i,j,k,ienergy) = eos_table(i,j,k,ienergy) &
                + energy_shift
           !convert back to erg/g
           eos_table(i,j,k,ienergy) = &
                eos_table(i,j,k,ienergy)/massn_cgs*mev_to_erg

        enddo
     enddo
  enddo

  min_energy = minval(eos_table(:,:,:,ienergy))
  min_press = minval(eos_table(:,:,:,ipress))
  
  if(min_press.le.0.0d0.or.min_energy.le.0.0d0) then
     write(6,"(1P10E15.6)") min_press, min_energy
     stop "myeos: min_press or min_energy <= 0.0d0"
  endif
  eos_table(:,:,:,ipress) = log10(eos_table(:,:,:,ipress))
  eos_table(:,:,:,ienergy) = log10(eos_table(:,:,:,ienergy))

  !now calculate lots of things we need in the table but are not
  !provided (sound speed for example)
  write(*,*) "Producing table derivatives"
  call derivatives_production

  write(*,*) "Checking for NaNs"
  !check data for NaN's, Inf's
  call checkdata3d(eos_table(:,:,:,ipress), &
       nlogrho,nlogtemp,nye,"total pressure")
  call checkdata3d(eos_table(:,:,:,ienergy), &
       nlogrho,nlogtemp,nye,"total energy")
  call checkdata3d(eos_table(:,:,:,ientropy), &
       nlogrho,nlogtemp,nye,"total entropy")
  call checkdata3d(dedt,&
       nlogrho,nlogtemp,nye,"total dedt")
  call checkdata3d(cs2, &
       nlogrho,nlogtemp,nye,"total cs2")
  call checkdata3d(dpdrho, &
       nlogrho,nlogtemp,nye,"total dpdrho")
  call checkdata3d(dpde,&
       nlogrho,nlogtemp,nye,"total dpde")

  !write out table in H5 format
  call date_and_time(DATE=date,VALUES=values)
  write(srho,*) nlogrho
  write(stemp,*) nlogtemp
  write(sye,*) nye
  timestamp = dble(values(1))*10000.0d0+dble(values(2))*100.0+dble(values(3)) + &
       (dble(values(5))+dble(values(6))/60.0d0 + dble(values(7))/3600.0d0 )/24.0

  h5filename = trim(adjustl(base))//"EOS_rho"//trim(adjustl(srho))// &
       "_temp"//trim(adjustl(stemp))//"_ye"//trim(adjustl(sye))// &
       "_version"//trim(adjustl(vnum))//"_"//trim(adjustl(date))//".h5"

  call write_table(h5filename,timestamp)

  write(*,*) "This was fun! We are done!"

end program eos_table_maker

