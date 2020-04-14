!-*-f90-*-
module eos_table_module

  implicit none

  integer :: nt_temp !number of temperature entries in initial table
  integer :: nt_yp !number of yp entries in initial table
  integer :: nt_rho !number of density entries in initial table
  integer :: nvars !number of variables for initial table , in general
                   !this ignores rho,T,ye.  You may have to adjust the
                   !read_table.F90 file to accomodiate various formats

  integer :: final_nvars !number of variables we need to make the
                         !final table, note this must be greater or at
                         !leat equal to nvars as every nvar gets
                         !interpolated to the final table.

  integer :: nlogrho, nlogtemp, nye !number of rho/T/yp entries in
                                    !final table

  real*8  :: logrhomin,logrhomax !min and max of final table, rho
  real*8  :: logtempmin,logtempmax !min and max of final table, temp
  real*8  :: yemin,yemax !min and max of final table, ye

  real*8  :: logrhomin_comps,logrhomin_nuclear !transition densities
  real*8  :: logtempmin_nuclear,logtempmax_nuclear !transition densities

  integer :: irho_start_nuclear, irho_start_comps !transition density
                                                  !indexes, most
                                                  !tables do not go as
                                                  !low as you want
                                                  !them too, this gets
                                                  !set based on
                                                  !choosen transition
                                                  !densties
  integer :: itemp_start_nuclear,itemp_stop_nuclear !transition
                                                    !temperature
                                                    !indexes, again,
                                                    !some tables do
                                                    !not go low
                                                    !enough.  With T
                                                    !sometimes tables
                                                    !do not go high
                                                    !enough for the
                                                    !needs of people

  logical :: noelectrons = .false. !experimental
  logical :: use_saved_electrons = .false. ! electrons take a while to
                                           ! compute, once final table
                                           ! size is set use this set
                                           ! to .true. to avoid
                                           ! recalculating them all
                                           ! the time

  integer :: total_hacks  = 0 !This is to keep track of the number of
                              !times the mass fractions differ by more
                              !then 20%.  Normally such a situtition
                              !would stop the program but if charge
                              !conservation is off by more then 15% we
                              !accept this

  logical :: Hempel = .false. !used for adding light nuclei in the
                              !final .h5 file and adding them to the
                              !mass fraction reconstruction
  logical :: Hempel_massfrac_reconstruction = .false. ! set to true to
                                                      ! use
                                                      ! conservative
                                                      ! mass fraction
                                                      ! interpolation
  logical :: lump_into_np = .false. !this splits the mass fraction of
                                    !light nuclei into free neutrons
                                    !and protons right from the
                                    !beginning, this must have Hempel
                                    != .true., but then turns it off
                                    !after it lumps
  logical :: Hyperon = .false. !similar to Hempel, but for Hyperons,
                               !Hempel_massfrac_reconstruction must be
                               !true


  real*8,allocatable :: table_logtemp(:) !array for initial table temp
                                         !(stores the log)
  real*8,allocatable :: table_logrho(:) !array for initial table rho
                                        !(stores the log)
  real*8,allocatable :: table_lintemp(:) !array for initial table temp
  real*8,allocatable :: table_linrho(:) !array for initial table rho
  real*8,allocatable :: table_yp(:) !array for initial table yp

  real*8,allocatable :: table_data(:,:,:,:) !array to store all
                                            !initial table data,
                                            !indexes are
                                            !temp,yp,rho,variable (see
                                            !below)

  real*8,allocatable :: logrho(:),logtemp(:), ye(:) !final table coordinates

  real*8,allocatable :: eos_table(:,:,:,:) !final table indexes are
                                           !temp,yp,rho,variable (see
                                           !below)

  real*8,allocatable :: eos_electrons(:,:,:,:) ! electrons contribution to everything
  real*8,allocatable :: eos_photons(:,:,:,:) !photons contribution to everything

  !some arrays for along the way, these are things we want to
  !calculate that are not in the initial table but we want in the
  !final table
  real*8,allocatable :: dpdrho(:,:,:) 
  real*8,allocatable :: dpde(:,:,:)
  real*8,allocatable :: dedt(:,:,:)
  real*8,allocatable :: dpdt(:,:,:)
  real*8,allocatable :: gamma(:,:,:)
  real*8,allocatable :: cs2(:,:,:)

  real*8, allocatable :: electron_gamma(:,:,:)
  real*8, allocatable :: electron_dedT(:,:,:)
  real*8, allocatable :: electron_dpdT(:,:,:)
  real*8, allocatable :: electron_dpdrhoT(:,:,:)

  ! electron number densities and derivatives
  ! needed for Coulomb corrections
  real*8, allocatable :: electron_xne(:,:,:)
  real*8, allocatable :: electron_dxnedd(:,:,:)
  real*8, allocatable :: electron_dxnedt(:,:,:)
  real*8, allocatable :: electron_dxneda(:,:,:)
  real*8, allocatable :: electron_dxnedz(:,:,:)

  integer :: n_electron_vars = 5 !for array allocation
  integer :: n_photon_vars = 3 !for array allocation

  !these variables are where we set the numeric index for the various
  !EOS variables.  Instead of remembering what index, just use these
  !variables.  Will be unique for each EOS
  integer :: ifree
  integer :: ienergy
  integer :: ipress
  integer :: iabar
  integer :: izbar
  integer :: imasseff
  integer :: ientropy
  integer :: imun
  integer :: imup
  integer :: imue
  integer :: ixn
  integer :: ixp
  integer :: ixa
  integer :: ixh
  integer :: ixd
  integer :: ixt
  integer :: ix3he
  integer :: ix4li
  integer :: ixL
  integer :: iLambdamasseff

  !these are internal
  integer :: electron_ipress = 1
  integer :: electron_ieps = 2
  integer :: electron_ientropy = 3
  integer :: electron_imue = 4
  integer :: electron_imuep = 5

  integer :: photon_ipress = 1
  integer :: photon_ieps = 2
  integer :: photon_ientropy = 3


  !some unit conversions
  real*8, parameter :: m_n = 939.566d0
  real*8, parameter :: m_p = 938.272d0
  real*8, parameter :: mev_to_erg = 1.60217733d-6
  real*8, parameter :: erg_to_mev = 6.24150636d5
  real*8, parameter :: amu_cgs = 1.66053873d-24
  real*8, parameter :: massn_cgs = 1.674927211d-24
  real*8, parameter :: mev_to_gram = 1.7826627d-27
  real*8, parameter :: amu_mev = 931.494d0
  real*8, parameter :: kb_erg = 1.380658d-16
  real*8, parameter :: kb_mev = 8.61738568d-11
  real*8, parameter :: temp_mev_to_kelvin = 1.1604447522806d10
  real*8, parameter :: planck = 6.626176d-27
  real*8, parameter :: pi = 3.14159265358979d0
  real*8, parameter :: clight = 2.99792458d10
  real*8, parameter :: avo = 6.0221367d23

  !this will be for storing things as log
  real*8 :: energy_shift = 0.0d0
  real*8 :: m_ref
  real*8 :: baryon_to_gram_conversion
  real*8 :: lowdensity_effectivemass = 939.0d0
  logical :: anal_mu_lowden = .true.    

contains
  
  !allocate initial table
  subroutine allocate_eos
    
    allocate(table_logtemp(nt_temp))
    allocate(table_lintemp(nt_temp))
    allocate(table_yp(nt_yp))
    allocate(table_logrho(nt_rho))
    allocate(table_linrho(nt_rho))
    allocate(table_data(nt_temp,nt_yp,nt_rho,nvars))

   
    table_data(:,:,:,:) = 0.0d0


  end subroutine allocate_eos
  

end module eos_table_module
