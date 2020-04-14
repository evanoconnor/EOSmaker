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
  real*8 dm,dz
  integer c1,c2,c3,c4,c5

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
  !!!!!!!!!!!!!!!!!!It is currently formatted for Hempel Stellar Collapse format, no leptons, no lots of extra nuclei

  initial_eos_table_name = "./original_tables/sfho_frdm_sc_v1.03.bin"
  Hempel = .true.
  Hempel_massfrac_reconstruction = .true.
  lump_into_np = .false.

  nt_temp = 81 !number of temperature entries in initial table
  nt_yp = 60 !number of yp entries in initial table
  nt_rho = 326 !number of density entries in initial table
  nvars = 21 !columns in the initial table - 3 (don't care about ye
             !and rho etc., know those)

  !size you want to final table in
  !denstiy
  nlogrho = 222 !number of density bins (326+25-25+7), 326 original table, plus 25 for low density, -25+7 to remove zeros at top
  logrhomin = 2.2202492d0 !global minimum of final table
  logrhomax = 15.220249d0+7.0d0/25.0d0 !global maximum of final table

  !temperature
  nlogtemp = 180 !number of temperature bins
  logtempmin = -2.0d0
  logtempmax = 2.2d0
  
  !ye
  nye = 60 !number of ye bins
  yemin = 0.01d0
  yemax = 0.60d0

  m_ref = amu_mev !reference mass of free nucleons
  baryon_to_gram_conversion = amu_cgs !mass unit of EOS

  base="Hempel_SFHo"
  vnum="_1.1"

  final_nvars = 22

  if (final_nvars.lt.nvars+1) stop "final_nvars must be at least nvars + mue variable"

  !specify nasty transitions........  at the following density we stop
  !using the thermodynamic results of the nuclear table and transition
  !to an ion eos.  The determine, try to find a density where the
  !pressures of the two EOS match up nicely and smoothly.  We require
  !so we also shift the energies at this density for every temperature
  !and ye to ensure continuity, this is described later on
  logrhomin_nuclear = 7.1d0
  
  ! We use the compositions of the nuclear table down to this log
  ! density...  Well out of NSE.... They do not do much
  logrhomin_comps = 5.0d0
    
  !The nuclear table may not go low enough in temperature, here is the
  !transition temp of the table, below we extrapolate or use ion eos
  logtempmin_nuclear= -1.0d0

  !Or not high enough, we like black holes, so does Shen so we are good here.
  logtempmax_nuclear = 2.2d0

  !set column numbers of initial table.  We assume the first three
  !columns are not going to be inlcuded (i.e. they are density, ye,
  !and temp, or something similar.  Adjust read_iniital_data if
  !required
  !.....1         logarithm of baryon mass density: log(rho) [g/cm^3]
!.....2         baryon number density: nb [fm^-3]
!.....3         logarithm of total proton fraction: log(yp)
!.....4         total proton fraction: yp (=ye)
!.....5         free energy per baryon wrt 938 MeV: f/nb-938 [MeV]
!.....6         internal energy per baryon wrt amu=931.49432 MeV: e/nb-amu [MeV]
!.....7         entropy per baryon: s/nb [kB]
!.....8         averaged mass number of heavy nucleii with A>4: <A>
!.....9         averaged charge number of heavy nuclei with A>4: <Z>
!.....10        effective mass: m^* [MeV]
!.....11        unbound neutron mass fraction: X_n
!.....12        unbound proton mass fraction: X_p
!.....13        alpha particle mass fraction: X_alpha
!.....14        heavy nuclei mass fraction of nuclei with A>4: X_heavy
!.....15        pressure p: [MeV/fm^3]
!.....16        neutron chemical potential wrt m_n: mu_n-m_n [MeV]
!.....17        proton chemical potential wrt m_n: mu_p-m_p [MeV]
!.....18        deuteron mass fraction: X_d
!.....19        triton mass fraction X_t
!.....20        helion (3He) mass fraction X_h
!.....21        4Li mass fraction X_li4

  ifree    = 5 !we don't use this
  ienergy  = 6 !specific internal energy code assumes units of MeV/baryon
  ientropy = 7 !in units of /k_B
  iabar    = 8 !mass number of heavies, average
  izbar    = 9 !atomic number of heavies, average
  imasseff = 10 !we don't use this
  ixn      = 11 !mass fraction of neutrons
  ixp      = 12 !mass fraction of protons
  ixa      = 13 !mass fraction of alphas
  ixh      = 14 !mass fraction of heavies
  ipress   = 15 !press assume units of MeV/fm^3
  

  !Note, must be careful regarding rest mass in/out of chemical
  !potential.  We assume incoming chemical potentials DO NOT include
  !rest mass difference.  This is true for H. Shen where the mu's are
  !defined relative to the nucleon mass M, which is the same for both
  !neutron and proton
  imun     = 16 !chemical potential of neutron, MeV
  imup     = 17 !chemical potential of proton, MeV
  imue     = 22 !chemical potential of electron, MeV, includes rest
                !mass, not part of initial table but part of final

  !extra Hemple variables
  ixd = 18
  ixt = 19
  ix3he = 20
  ix4li = 21

  use_saved_electrons = .false.
  noelectrons = .false.

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
  !stuff in there

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
  if (Hempel) then
     write(*,*) "Just checking Hempel conservations, should be perfect"
     c1 = 0
     c2 = 0
     c3 = 0
     c4 = 0
     c5 = 0
     do i=irho_start_comps,nlogrho
        do j=1,nlogtemp
           do k=1,nye

              dm = eos_table(i,j,k,ixn) + eos_table(i,j,k,ixp) + &
                   eos_table(i,j,k,ixa) + eos_table(i,j,k,ixh) + eos_table(i,j,k,ixt) + &
                   eos_table(i,j,k,ixd) + eos_table(i,j,k,ix3he) + &
                   eos_table(i,j,k,ix4li) - 1.0d0
              
              dz = eos_table(i,j,k,ixp) + 0.5d0*eos_table(i,j,k,ixa) + &
                   eos_table(i,j,k,izbar)*eos_table(i,j,k,ixh)/eos_table(i,j,k,iabar) + &
                   eos_table(i,j,k,ixt)/3.0d0 + 0.5d0*eos_table(i,j,k,ixd) + &
                   2.0d0*eos_table(i,j,k,ix3he)/3.0d0 + eos_table(i,j,k,ix4li)*0.75d0 - ye(k)

              if(abs(dm).gt.1.0d-4 .or. abs(dz)/ye(k).gt.1.0d-4) then
                 c1 = c1 + 1
                 if(abs(dm).gt.0.001d0.or.abs(dz)/ye(k).gt.0.001d0) then
                    c2 = c2 + 1
                 endif
                 if(abs(dm).gt.0.01d0.or.abs(dz)/ye(k).gt.0.01d0) then
                    c3 = c3 + 1
                 endif
                 if(abs(dm).gt.0.05d0.or.abs(dz)/ye(k).gt.0.05d0) then
                    c4 = c4 + 1
                 endif
                 if(abs(dm).gt.1.0d-1.or.abs(dz)/ye(k).gt.1.0d-1) then
                    c5 = c5+1
                 endif
              endif
              
           enddo
        enddo
     enddo
     
     write(6,*) "Mass/charge conservation check"
     write(6,*) "total points: ", nlogrho*nlogtemp*nye
     write(6,"(A26,i9,f11.7,A1)") "|dm| or |dz|/ye > 0.0001: ", c1,100.0*c1/(nlogrho*nlogtemp*nye),"%"
     write(6,"(A26,i9,f11.7,A1)") "|dm| or |dz|/ye > 0.0010: ", c2,100.0*c2/(nlogrho*nlogtemp*nye),"%"
     write(6,"(A26,i9,f11.7,A1)") "|dm| or |dz|/ye > 0.0100: ", c3,100.0*c3/(nlogrho*nlogtemp*nye),"%"
     write(6,"(A26,i9,f11.7,A1)") "|dm| or |dz|/ye > 0.0500: ", c4,100.0*c4/(nlogrho*nlogtemp*nye),"%"
     write(6,"(A26,i9,f11.7,A1)") "|dm| or |dz|/ye > 0.1000: ", c5,100.0*c5/(nlogrho*nlogtemp*nye),"%"

     if (lump_into_np) then

        do i=irho_start_comps,nlogrho
           do j=1,nlogtemp
              do k=1,nye
                 eos_table(i,j,k,ixp) = eos_table(i,j,k,ixp) + 0.5d0*eos_table(i,j,k,ixd)+eos_table(i,j,k,ixt)/3.0d0+ &
                      2.0d0*eos_table(i,j,k,ix3he)/3.0d0+3.0d0*eos_table(i,j,k,ix4li)/4.0d0
                 eos_table(i,j,k,ixn) = eos_table(i,j,k,ixn) + 0.5d0*eos_table(i,j,k,ixd)+2.0d0*eos_table(i,j,k,ixt)/3.0d0+ &
                      eos_table(i,j,k,ix3he)/3.0d0+eos_table(i,j,k,ix4li)/4.0d0
              enddo
           enddo
        enddo

     endif

  else
     write(*,*) "Fixing Compositions"
     call fix_compositions
  endif
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
  min_energy = minval(eos_table(:,:,:,ipress))
  imin_energy = minloc(eos_table(:,:,:,ipress))
  write(*,*) "min_pressure: ", min_energy
  write(*,*) "location, indices",imin_energy(1:3)
  write(*,*) "location, values", 10.0d0**logrho(imin_energy(1)),10.0d0**logtemp(imin_energy(2)),ye(imin_energy(3))

  ! shift energy to be compatible with log
  do i=1,nlogrho
     do j=1,nlogtemp
        do k=1,nye

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

