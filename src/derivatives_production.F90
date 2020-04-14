
!-*-f90-*-
subroutine derivatives_production
  
  use eos_table_module
  implicit none

  integer i,j,k
  real*8 :: dx,x1,x2,f1,f2,z,zz

  real*8, allocatable :: dsdlnT(:,:,:)
  real*8, allocatable :: dsdlnrho(:,:,:)
  real*8, allocatable :: dpdrhoT(:,:,:)
  real*8, allocatable :: dedrhoT(:,:,:)

  integer :: igamma  = 1
  ! igamma: There are multiple ways to compute gamma1, via the entropy
  !         turns out to be pretty good.
  ! 1 -> Gamma1 = &
  !      dlnP/dlnrho|T,Y_e - ds/dlnrho|T,Y_e * (dlnP/dlnT / ds/dlnT)_rho,Y_e
  ! 2 -> Gamma1 = P/rho * ( dpdrho|e,Y_e + P/rho**2 dpde|rho,Y_e
  ! 3 -> Gamma1 = dlnP/dlnrho|T,Y_e + T * (dpdT|rho,Y_e)**2 / &
  !      (P*rho * dedT|rho,Y_e)

!###########################################################################  

  ! dedT, e is erg in log10, T is in MeV in log 10
  allocate(dedT(nlogrho,nlogtemp,nye))
  dedT = 0.0d0

  do k=1,nye
     do i=1,nlogrho
        do j=2,nlogtemp-1
           x1 = logtemp(j-1)
           f1 = eos_table(i,j-1,k,ienergy)
           x2 = logtemp(j+1)
           f2 = eos_table(i,j+1,k,ienergy)
           dedT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
                * 10.0d0**eos_table(i,j,k,ienergy)
        enddo

        ! boundaries: one-sided derivative
        x2 = logtemp(2)
        x1 = logtemp(1)
        f2 = eos_table(i,2,k,ienergy)
        f1 = eos_table(i,1,k,ienergy)
        dedT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1)  &
                * 10.0d0**eos_table(i,1,k,ienergy)
        x2 = logtemp(nlogtemp)
        x1 = logtemp(nlogtemp-1)
        f2 = eos_table(i,nlogtemp,k,ienergy)
        f1 = eos_table(i,nlogtemp-1,k,ienergy)
        dedT(i,nlogtemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(nlogtemp) &
             * 10.0d0**eos_table(i,nlogtemp,k,ienergy) 
     enddo
  enddo

!########################################################################### 
  
  !dpdT, p is in dyn/cm^2 in log10, T is in MeV in log 10
  !dsdlnT, s is in k_B/baryon not in log 10, T is in meV in log 10
  allocate(dpdT(nlogrho,nlogtemp,nye))
  allocate(dsdlnT(nlogrho,nlogtemp,nye))
  dpdT = 0.0d0
  dsdlnT = 0.0d0

  
  do k=1,nye
     do i=1,nlogrho
        do j=2,nlogtemp-1
           x1 = logtemp(j-1)
           f1 = eos_table(i,j-1,k,ipress)
           x2 = logtemp(j+1)
           f2 = eos_table(i,j+1,k,ipress)
           dpdT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
                 * 10.0d0**eos_table(i,j,k,ipress)

           x1 = logtemp(j-1)
           f1 = log10(eos_table(i,j-1,k,ientropy))
           x2 = logtemp(j+1)
           f2 = log10(eos_table(i,j+1,k,ientropy))
           dsdlnT(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)

        enddo

        ! boundaries: one-sided derivative
        x1 = logtemp(1)
        f1 = eos_table(i,1,k,ipress)
        x2 = logtemp(2)
        f2 = eos_table(i,2,k,ipress)
        dpdT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
             10.0d0**eos_table(i,1,k,ipress)

        x1 = logtemp(1)
        f1 = log10(eos_table(i,1,k,ientropy))
        x2 = logtemp(2)
        f2 = log10(eos_table(i,2,k,ientropy))
        dsdlnT(i,1,k) = (f2-f1)/(x2-x1) * eos_table(i,1,k,ientropy)


        x1 = logtemp(nlogtemp-1)
        f1 = eos_table(i,nlogtemp-1,k,ipress)
        x2 = logtemp(nlogtemp)
        f2 = eos_table(i,nlogtemp,k,ipress)
        dpdT(i,nlogtemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
             10.0d0**eos_table(i,1,k,ipress)

        x1 = logtemp(nlogtemp-1)
        f1 = log10(eos_table(i,nlogtemp-1,k,ientropy))
        x2 = logtemp(nlogtemp)
        f2 = log10(eos_table(i,nlogtemp,k,ientropy))
        dsdlnT(i,nlogtemp,k) = (f2-f1)/(x2-x1) & 
             * eos_table(i,nlogtemp,k,ientropy)

     enddo
  enddo

!########################################################################### 
  ! dp/drho|T, p is in dyn/cm^2 in log 10, T is in MeV in log 10
  ! ds/dlnrho|T, s is in k_b/baryon not in log 10, T is in MeV in log 10
  ! de/drho|T, e is in erg in log 10, T is in MeV in log 10
  ! dp/drho|e = dp/drho|T + dp/dT * (-de/drho|T) / de/dT
  ! dp/de|rho = dp/dT / de/dT
  
  allocate(dpdrho(nlogrho,nlogtemp,nye))
  allocate(dpde(nlogrho,nlogtemp,nye))
  allocate(dsdlnrho(nlogrho,nlogtemp,nye))
  allocate(dpdrhoT(nlogrho,nlogtemp,nye))
  allocate(dedrhoT(nlogrho,nlogtemp,nye))
  dpdrho   = 0.0d0
  dpde     = 0.0d0
  dsdlnrho = 0.0d0
  dpdrhoT  = 0.0d0
  dedrhoT  = 0.0d0

  do k=1,nye
     do j=1,nlogtemp
        do i=2,nlogrho-1
           x1 = logrho(i-1)
           x2 = logrho(i+1)
           f1 = eos_table(i-1,j,k,ipress)
           f2 = eos_table(i+1,j,k,ipress)
           dpdrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
                * 10.0d0**eos_table(i,j,k,ipress)

           x1 = logrho(i-1)
           x2 = logrho(i+1)
           f1 = eos_table(i-1,j,k,ienergy)
           f2 = eos_table(i+1,j,k,ienergy)
           dedrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
                * 10.0d0**eos_table(i,j,k,ienergy)

           x1 = logrho(i-1)
           x2 = logrho(i+1)
           f1 = log10(eos_table(i-1,j,k,ientropy))
           f2 = log10(eos_table(i+1,j,k,ientropy))
           dsdlnrho(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
          enddo
          
          ! boundaries: one-sided derivative
          x1 = logrho(1)
          x2 = logrho(2)
          f1 = eos_table(1,j,k,ipress)
          f2 = eos_table(2,j,k,ipress)
          dpdrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
               * 10.0d0**eos_table(1,j,k,ipress)

          f1 = eos_table(1,j,k,ienergy)
          f2 = eos_table(2,j,k,ienergy)
          dedrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
               * 10.0d0**eos_table(1,j,k,ienergy)

          x1 = logrho(nlogrho-1)
          x2 = logrho(nlogrho)
          f1 = eos_table(nlogrho-1,j,k,ipress)
          f2 = eos_table(nlogrho,j,k,ipress)
          dpdrhoT(nlogrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nlogrho) &
               * 10.0d0**eos_table(nlogrho,j,k,ipress)

          f1 = eos_table(nlogrho-1,j,k,ienergy)
          f2 = eos_table(nlogrho,j,k,ienergy)
          dedrhoT(nlogrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nlogrho) &
               * 10.0d0**eos_table(nlogrho,j,k,ienergy)

          x1 = logrho(1)
          x2 = logrho(2)
          f1 = log10(eos_table(1,j,k,ientropy))
          f2 = log10(eos_table(2,j,k,ientropy))
          dsdlnrho(1,j,k) = (f2-f1)/(x2-x1) * eos_table(1,j,k,ientropy)

          x1 = logrho(nlogrho-1)
          x2 = logrho(nlogrho)
          f1 = log10(eos_table(nlogrho-1,j,k,ientropy))
          f2 = log10(eos_table(nlogrho,j,k,ientropy))
          dsdlnrho(nlogrho,j,k) = (f2-f1)/(x2-x1) * eos_table(nlogrho,j,k,ientropy)

     enddo
  enddo

  do k=1,nye
     do j=1,nlogtemp
        do i=1,nlogrho
           dpdrho(i,j,k) = dpdrhoT(i,j,k) + dpdT(i,j,k) * &
                (-dedrhoT(i,j,k))/dedT(i,j,k)
           
           dpde(i,j,k) = dpdT(i,j,k) / dedT(i,j,k)
        enddo
     enddo
  enddo

  !########################################################################### 
  allocate(gamma(nlogrho,nlogtemp,nye))
  allocate(cs2(nlogrho,nlogtemp,nye))
  cs2 = 0.0d0
  gamma = 0.0d0

  do k=1,nye
     do j=1,nlogtemp
        do i=1,nlogrho
           z = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress)
           zz = 10.0d0**logtemp(j)/10.0d0**eos_table(i,j,k,ipress)

           if(igamma .eq. 1) then
              gamma(i,j,k) = z * dpdrhoT(i,j,k) & 
                   - dsdlnrho(i,j,k) * zz * dpdT(i,j,k) / dsdlnT(i,j,k)
           endif

           if(igamma .eq. 2) then
              gamma(i,j,k) =  dpdrho(i,j,k) + 10.0d0**eos_table(i,j,k,ipress) / &
                   (10.0d0**logrho(i))**2 * dpde(i,j,k)
              gamma(i,j,k) = z * gamma(i,j,k)
           endif


           if(igamma .eq. 3) then
              gamma(i,j,k) = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress) &
                   * dpdrhoT(i,j,k) &
                   + 10.0d0**logtemp(j)*dpdT(i,j,k)**2 / &
                   (10.0d0**eos_table(i,j,k,ipress)*10.0d0**logrho(i) * dedT(i,j,k))
           endif

           cs2(i,j,k) = gamma(i,j,k) / z

        enddo
     enddo
  enddo

end subroutine derivatives_production
