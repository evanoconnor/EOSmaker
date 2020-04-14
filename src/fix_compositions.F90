!-*-f90-*-
#define DEBUG_FRACTIONS 0
subroutine fix_compositions

  use eos_table_module
  implicit none

  real*8 norm
  real*8 :: tiny = 1.0d-30
  integer i,j,k

  real*8 K1,K2
  real*8 xp,xa,xh,xn
  real*8 xps,xas,xhs,xns
  real*8 xps0,xas0,xhs0,xns0
  real*8 zbar,abar,zbars,abars,zbars0,abars0
  
  real*8 mass_cons, charge_cons
  real*8 dm,dz
  real*8 dxh, dxa

  integer tf(2),rf(2),yf(2)
  real*8 min_x,max_x

  integer c1,c2,c3,c4,c5
  integer d1,d2,d3,d4,d5

  integer ordering(4)

  real*8 local_ye

  real*8 csort(4,2),xtemp(1,2)

  character(len=64) :: cstring = "main"
  logical done,swapped

  ! The problem is the following:
  ! In the original Shen table, all variables are
  ! given to 6/7 significant figures only, hence
  ! mass conservation and charge conservation of
  ! the NSE compositions are good to only 6 digits or
  ! so -- this leads to problem with the low density
  ! and electron EOS that need Abar/Zbar = Ye, where
  ! Abar and Zbar are constructed from all components
  ! as in electron_eos.F90
  !
  ! In this routine, we fix mass and charge conservation
  ! by fixing x_n and x_h, adjusting x_alpha and x_p
  ! and the average Z of the heavies if problems show
  ! up in computing x_alpha and x_p (such as negative
  ! values that are, of course, not allowed).
  
  ! ensure minimum of mass fractions is indeed 0.0d0

  do i=irho_start_comps,nlogrho
     do j=1,nlogtemp
        do k=1,nye
           eos_table(i,j,k,ixn) = max(eos_table(i,j,k,ixn),0.0d0)
           eos_table(i,j,k,ixp) = max(eos_table(i,j,k,ixp),0.0d0)
           eos_table(i,j,k,ixa) = max(eos_table(i,j,k,ixa),0.0d0)
           eos_table(i,j,k,ixh) = max(eos_table(i,j,k,ixh),0.0d0)
        enddo
     enddo
  enddo

  ! check how bad things are
  c1 = 0
  c2 = 0
  c3 = 0
  c4 = 0
  c5 = 0

  d1 = 0
  d2 = 0
  d3 = 0
  d4 = 0
  d5 = 0

  do i=irho_start_comps,nlogrho
     do j=1,nlogtemp
        do k=1,nye
           
           xn = eos_table(i,j,k,ixn)
           xp = eos_table(i,j,k,ixp)
           xa = eos_table(i,j,k,ixa)
           xh = eos_table(i,j,k,ixh)
           zbar = eos_table(i,j,k,izbar)
           abar = eos_table(i,j,k,iabar)

           xns = xn
           xps = xp
           xas = xa
           xhs = xh
           abars = abar
           zbars = zbar

           if(xh.lt.1.0d-14) then 
              eos_table(i,j,k,izbar) = tiny
              eos_table(i,j,k,iabar) = tiny
           endif
           
           ! make sure we have no pathological heavy nucleus (A <= 4)
           if(abar.le.4.0d0) then
              xh = 0.0d0
              xhs = 0.0d0
              zbar = tiny
              abar = tiny
           endif

           ! make sure we have no pathological heavy nucleus (Z <= 0)
           if(zbar.le.0.0d0) then
              xh = 0.0d0
              xhs = 0.0d0
              zbar = tiny
              abar = tiny
           endif


           dm = mass_cons(xn,xp,xa,xh)
           dz = charge_cons(xn,xp,xa,xh,abar,zbar,ye(k))

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

           norm = xn + xp + xa + xh
           xn = xn/norm
           xp = min(xp/norm,ye(k))
           xa = xa/norm
           xh = xh/norm

           dm = mass_cons(xn,xp,xa,xh)
           dz = charge_cons(xn,xp,xa,xh,abar,zbar,ye(k))

           if(abs(dm).gt.1.0d-4 .or. abs(dz)/ye(k).gt.1.0d-4) then
              d1 = d1 + 1
              if(abs(dm).gt.0.001d0.or.abs(dz)/ye(k).gt.0.001d0) then
                 d2 = d2 + 1
              endif
              if(abs(dm).gt.0.01d0.or.abs(dz)/ye(k).gt.0.01d0) then
                 d3 = d3 + 1
              endif
              if(abs(dm).gt.0.05d0.or.abs(dz)/ye(k).gt.0.05d0) then
                 d4 = d4 + 1
              endif
              if(abs(dm).gt.1.0d-1.or.abs(dz)/ye(k).gt.1.0d-1) then
                 d5 = d5+1
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

  write(6,*) " "
  write(6,*) "Mass/charge conservation check with mass normalization"
  write(6,*) "total points: ", nlogrho*nlogtemp*nye
  write(6,"(A19,i9,f11.7,A1)") "|dz|/ye > 0.0001: ", d1,100.0*d1/(nlogrho*nlogtemp*nye),"%"
  write(6,"(A19,i9,f11.7,A1)") "|dz|/ye > 0.0010: ", d2,100.0*d2/(nlogrho*nlogtemp*nye),"%"
  write(6,"(A19,i9,f11.7,A1)") "|dz|/ye > 0.0100: ", d3,100.0*d3/(nlogrho*nlogtemp*nye),"%"
  write(6,"(A19,i9,f11.7,A1)") "|dz|/ye > 0.0500: ", d4,100.0*d4/(nlogrho*nlogtemp*nye),"%"
  write(6,"(A19,i9,f11.7,A1)") "|dz|/ye > 0.1000: ", d5,100.0*d5/(nlogrho*nlogtemp*nye),"%"

  ! Fixing mass fractions while obeying
  ! mass conservation is a difficult task.
  ! We will only mess with cases in which |dm|
  ! or |dz| > 0.01.

  ! normalize fractions, but ensure 
  ! (a) mass conservation
  ! xn + xp + xa + xh = 1
  ! (b) charge conservation
  ! xp + 0 * xn + 2/4 * xa + Z/A * xh = Ye

  ! What are mass fractions used for?  1. They are needed to solved
  ! the ion eos at low densities 2. They are extensively needed in for
  ! neutrino interactions as each interaction depends on the
  ! composition.  This means we must take care when adjusting the mass
  ! fractions to not change a mass fraction by a relatively large
  ! amount.  For example, protons in the collapse phase contribute
  ! significantly to the neutrino emission while they're mass fraction
  ! is very small.


  do i=irho_start_comps,nlogrho
     do j=1,nlogtemp
        do k=nye,1,-1

           xn = eos_table(i,j,k,ixn)
           xp = eos_table(i,j,k,ixp)
           xa = eos_table(i,j,k,ixa)
           xh = eos_table(i,j,k,ixh)
           zbar = eos_table(i,j,k,izbar)
           abar = eos_table(i,j,k,iabar)

           local_ye = ye(k)

           call one_time_fix_comps(xn,xp,xa,xh,abar,zbar,local_ye)

           eos_table(i,j,k,ixn) = xn
           eos_table(i,j,k,ixp) = xp
           eos_table(i,j,k,ixa) = xa
           eos_table(i,j,k,ixh) = xh
           eos_table(i,j,k,izbar) = zbar
           eos_table(i,j,k,iabar) = abar
           
        enddo
     enddo
  enddo

end subroutine fix_compositions

subroutine one_time_fix_comps(xn,xp,xa,xh,abar,zbar,ye)

  use eos_table_module, only : total_hacks
  implicit none
  real*8, intent(inout) :: xn,xp,xa,xh
  real*8, intent(inout) :: abar,zbar
  
  real*8, intent(in) :: ye
  
  real*8 :: xn0,xp0,xa0,xh0,abar0,zbar0
  real*8 :: xns,xps,xas,xhs,abars,zbars
  real*8 :: tiny = 1.0d-14
  integer :: ordering(4)

  real*8 :: dz,dm,norm,dz_old
  real*8 :: mass_cons
  real*8 :: charge_cons
  integer :: smallest_error
  real*8 :: beta
  logical :: cont = .false.

  real*8 :: errors(5),neutron_error,proton_error,alpha_error,heavy_error,zbar_error,abar_error

  !sometimes negative...
  xn = max(xn,0.0d0)
  xp = max(xp,0.0d0)
  xa = max(xa,0.0d0)
  xh = max(xh,0.0d0)
  zbar = max(zbar,0.0d0)
  abar = max(abar,0.0d0)

  !takes all the compositions and fixes them
  !original values
  xn0 = xn
  xp0 = xp
  xa0 = xa
  xh0 = xh
  abar0 = abar
  zbar0 = zbar
  
  !to prevent NaNs
  if(xh.lt.1.0d-14) then
     zbar = tiny
     abar = tiny
  endif
           
  ! make sure we have no pathological heavy nucleus (A <= 4)
  if(abar.le.4.0d0) then
     xh = 0.0d0
     zbar = tiny
     abar = tiny
  endif

  ! make sure we have no pathological heavy nucleus (Z <= 0)
  if(zbar.le.0.0d0) then
     xh = 0.0d0
     xhs = 0.0d0
     zbar = tiny
     abar = tiny
  endif
  
  !normalize to 1
  norm = xn + xp + xa + xh
  xn = xn/norm
  xp = xp/norm
  xa = xa/norm
  xh = xh/norm

  !normalized original
  xns = xn
  xps = xp
  xas = xa
  xhs = xh
  abars = abar
  zbars = zbar
  
  dm = mass_cons(xn,xp,xa,xh)
  dz = charge_cons(xn,xp,xa,xh,abar,zbar,ye)
  dz_old = dz

  if(abs(dz)/ye.gt.0.0001d0) then

     !will set to true if we have a violating case but we approve the change
     cont = .false.

     !first get order of massfractions order returns array
     !ordering(4) first entry is most abundant fractions, 1
     !is neutron, 2 is proton, 3 is alpha, 4 is heavy.  This
     !means if ordering(1).eq.1 then neutron is most
     !abundant, if ordering(3).eq.1 then neutron is third
     !most abundant
     call sort_comps(xn,xp,xa,xh,ordering)

     if (ordering(1).eq.1) then
        !if neutrons are dominate
        !check proton error.
        errors(:) = 101.0d0
        if (xps.gt.0.0d0) then
           xp = ye-0.5d0*xas-zbars/abars*xhs
           xn = 1.0d0-xp-xas-xhs
           proton_error = (xp-xps)/xps
           neutron_error = (xn-xns)/xns
           errors(1) = sqrt(proton_error**2+neutron_error**2)
           if (xp.lt.0.0d0) errors(1) = 100.0d0
           if (xn.lt.0.0d0) errors(1) = 100.0d0
           xp = xps
           xn = xns
        endif
        !check alpha error
        if (xas.gt.0.0d0) then
           xa = (ye-xps-zbars/abars*xhs)*2.0d0
           xn = 1.0d0-xps-xa-xhs
           alpha_error = (xa-xas)/xas
           neutron_error = (xn-xns)/xns
           errors(2) = sqrt(alpha_error**2+neutron_error**2)
           if (xa.lt.0.0d0) errors(2) = 100.0d0
           if (xn.lt.0.0d0) errors(2) = 100.0d0
           xa = xas
           xn = xns
        endif
        !check heavy error
        if (xhs.gt.0.0d0) then
           xh = (ye-xps-0.5d0*xas)*abars/zbars
           xn = 1.0d0-xps-xas-xh
           heavy_error = (xh-xhs)/xhs
           neutron_error = (xn-xns)/xns
           errors(3) = sqrt(heavy_error**2+neutron_error**2)
           if (xh.lt.0.0d0) errors(3) = 100.0d0
           if (xn.lt.0.0d0) errors(3) = 100.0d0
           xh = xhs
           xn = xns
        endif        
        !check zbar error
        if (zbar.gt.0.0d0) then
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
           zbar_error = (zbar-zbars)/zbars
           errors(4) = abs(zbar_error)
           if (zbar.lt.0.0d0) errors(4) = 100.0d0
           if (zbar.gt.abars) errors(4) = 102.0d0
           zbar = zbars
        endif  
        !check abar error
        if (abar.gt.0.0d0) then
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
           abar_error = (abar-abars)/abars
           errors(5) = abs(abar_error)
           if (abar.lt.0.0d0) errors(5) = 100.0d0
           if (abar.lt.zbars) errors(5) = 102.0d0
           abar = abars
        endif   
        !now find and use smallest error
        smallest_error = minloc(errors,1)
        if (errors(smallest_error).gt.0.1d0) then
           !attempt to fix in some cases
#if DEBUG_FRACTIONS
           write(*,"(1P10E18.9)") xns,xps,xas,xhs,abars,zbars,abs(dz)/ye,ye
           write(*,"(1P10E18.9)") errors
           write(*,"(1P10E18.9)") xn,xp,xa,xh,zbar,abar
           write(*,"(1P10E18.9)") (ye-xps-0.5d0*xas)*abars/zbars,1.0d0-xps-xas-(ye-xps-0.5d0*xas)*abars/zbars
#endif
           if (errors(smallest_error).gt.0.2d0.and.abs(dz)/ye.gt.0.15d0) then
              write(*,"(a70,2P10E18.9)") "Conservation is horrible, something we just have to accept...", &
                   abs(dz)/ye,errors(smallest_error)
              total_hacks = total_hacks+1
              write(*,"(a20,i10,a40)") "This error occured ",total_hacks," times, this should be a low number..."
              cont = .true.
           else if (ordering(2)*ordering(3).eq.12.and.xa.gt.0.01.and.xh.gt.0.01) then
              !attempt to change xa and xh by the same realtive amount, equations follow.
              beta = (ye-xps-0.5d0*xas-zbars/abars*xhs)/(xas*0.5d0+xhs*zbars/abars)
              xa = (1.0d0+beta)*xas
              xh = (1.0d0+beta)*xhs
              xn = 1.0d0-xps-xa-xh
#if DEBUG_FRACTIONS
              write(*,"(1P10E18.9)") xn,xp,xa,xh,zbar,abar
              write(*,"(1P10E18.9)") ye-xps-0.5d0*xa-zbars/abars*xh,abs(xn-xns)/xns,abs(xa-xas)/xas,abs(xh-xhs)/xhs
#endif
              smallest_error = -1 !special case that means fractions set
           else
#if DEBUG_FRACTIONS
              write(*,*) "Hope error is < 25%..."
#endif
           endif
        endif
        if(smallest_error.eq.1) then
           !neutron and protons give the smallest error
           xp = ye-0.5d0*xas-zbars/abars*xhs
           xn = 1.0d0-xp-xas-xhs
        else if (smallest_error.eq.2) then
           !neutrons and alphas give smallest error
           xa = (ye-xps-zbars/abars*xhs)*2.0d0
           xn = 1.0d0-xps-xa-xhs
        else if (smallest_error.eq.3) then
           !neutrons and heavies give smallest error
           xh = (ye-xps-0.5d0*xas)*abars/zbars
           xn = 1.0d0-xps-xas-xh
        else if (smallest_error.eq.4) then
           !zbar gives smallest error
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
        else if (smallest_error.eq.5) then
           !abar gives smallest error
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
        else if (smallest_error.eq.-1) then
#if DEBUG_FRACTIONS
           write(*,*) "special case, fractions set"
#endif
        else
           stop "whats the smallest error?"
        endif
     else if (ordering(1).eq.2) then
        !if protons are dominate
        !check neutron error.
        errors(:) = 101.0d0
        if (xns.gt.0.0d0) then
           xp = ye-0.5d0*xas-zbars/abars*xhs
           xn = 1.0d0-xp-xas-xhs
           proton_error = (xp-xps)/xps
           neutron_error = (xn-xns)/xns
           errors(1) = sqrt(proton_error**2+neutron_error**2)
           if (xp.lt.0.0d0) errors(1) = 100.0d0
           if (xn.lt.0.0d0) errors(1) = 100.0d0
           xp = xps
           xn = xns
        endif
        !check alpha error
        if (xas.gt.0.0d0) then
           xa = (1.0d0-xns-ye+(zbars/abars-1.0d0)*xhs)*2.0d0
           xp = 1.0d0-xns-xa-xhs
           alpha_error = (xa-xas)/xas
           proton_error = (xp-xps)/xps
           errors(2) = sqrt(alpha_error**2+proton_error**2)
           if (xa.lt.0.0d0) errors(2) = 100.0d0
           if (xp.lt.0.0d0) errors(2) = 100.0d0
           xa = xas
           xp = xps
        endif
        !check heavy error
        if (xhs.gt.0.0d0) then
           xh = (ye-1.0d0+xns+0.5d0*xas)/(zbars/abars-1.0d0)
           xp = 1.0d0-xns-xas-xh
           heavy_error = (xh-xhs)/xhs
           proton_error = (xp-xps)/xps
           errors(3) = sqrt(heavy_error**2+proton_error**2)
           if (xh.lt.0.0d0) errors(3) = 100.0d0
           if (xp.lt.0.0d0) errors(3) = 100.0d0
           xh = xhs
           xp = xps
        endif        
        !check zbar error
        if (zbar.gt.0.0d0) then
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
           zbar_error = (zbar-zbars)/zbars
           errors(4) = abs(zbar_error)
           if (zbar.lt.0.0d0) errors(4) = 100.0d0
           if (zbar.gt.abars) errors(4) = 102.0d0
           zbar = zbars
        endif  
        !check abar error
        if (abar.gt.0.0d0) then
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
           abar_error = (abar-abars)/abars
           errors(5) = abs(abar_error)
           if (abar.lt.0.0d0) errors(5) = 100.0d0
           if (abar.lt.zbars) errors(5) = 102.0d0
           abar = abars
        endif   
        !now find and use smallest error
        smallest_error = minloc(errors,1)
        if (errors(smallest_error).gt.0.1) then
!           write(*,*) "Hope the error is smaller than 0.2, protons: smallest error > 10%"
        endif
        if(smallest_error.eq.1) then
           !protons and neutrons give the smallest error
           xp = ye-0.5d0*xas-zbars/abars*xhs
           xn = 1.0d0-xp-xas-xhs
        else if (smallest_error.eq.2) then
           !protons and alphas give smallest error
           xa = (1.0d0-xns-ye+(zbars/abars-1.0d0)*xhs)*2.0d0
           xp = 1.0d0-xns-xa-xhs
        else if (smallest_error.eq.3) then
           !protons and heavies give the smallest error
           xh = (ye-1.0d0+xns+0.5d0*xas)/(zbars/abars-1.0d0)
           xp = 1.0d0-xns-xas-xh
        else if (smallest_error.eq.4) then
           !zbar gives smallest error
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
        else if (smallest_error.eq.5) then
           !abar gives smallest error
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
        else
           stop "whats the smallest error?"
        endif        
     else if (ordering(1).eq.3) then
        !if alphas are dominate
        !check neutron error.
        errors(:) = 101.0d0
        if (xns.gt.0.0d0) then
           xa = (ye-xps-zbars/abars*xhs)*2.0d0
           xn = 1.0d0-xps-xa-xhs
           alpha_error = (xa-xas)/xas
           neutron_error = (xn-xns)/xns
           errors(1) = sqrt(alpha_error**2+neutron_error**2)
           if (xa.lt.0.0d0) errors(1) = 100.0d0
           if (xn.lt.0.0d0) errors(1) = 100.0d0
           xa = xas
           xn = xns
        endif
        !check proton error
        if (xps.gt.0.0d0) then
           xa = (1.0d0-xns-ye+(zbars/abars-1.0d0)*xhs)*2.0d0
           xp = 1.0d0-xns-xa-xhs
           alpha_error = (xa-xas)/xas
           proton_error = (xp-xps)/xps
           errors(2) = sqrt(alpha_error**2+proton_error**2)
           if (xa.lt.0.0d0) errors(2) = 100.0d0
           if (xp.lt.0.0d0) errors(2) = 100.0d0
           xa = xas
           xp = xps
        endif
        !check heavy error
        if (xhs.gt.0.0d0) then
           xh = (ye-0.5d0*xps-0.5d0+xns*0.5d0)/(zbars/abars-0.5d0)           
           xa = 1.0d0-xns-xps-xh
           heavy_error = (xh-xhs)/xhs
           alpha_error = (xa-xas)/xas
           errors(3) = sqrt(heavy_error**2+alpha_error**2)
           if (xh.lt.0.0d0) errors(3) = 100.0d0
           if (xa.lt.0.0d0) errors(3) = 100.0d0
           xh = xhs
           xa = xas
        endif        
        !check zbar error
        if (zbar.gt.0.0d0) then
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
           zbar_error = (zbar-zbars)/zbars
           errors(4) = abs(zbar_error)
           if (zbar.lt.0.0d0) errors(4) = 100.0d0
           if (zbar.gt.abars) errors(4) = 102.0d0
           zbar = zbars
        endif  
        !check abar error
        if (abar.gt.0.0d0) then
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
           abar_error = (abar-abars)/abars
           errors(5) = abs(abar_error)
           if (abar.lt.0.0d0) errors(5) = 100.0d0
           if (abar.lt.zbars) errors(5) = 102.0d0
           abar = abars
        endif   
        !now find and use smallest error
        smallest_error = minloc(errors,1)
        if (errors(smallest_error).gt.0.1d0) then
#if DEBUG_FRACTIONS
           write(*,"(1P10E18.9)") ye,xn0,xp0,xa0,xh0,abar0,zbar0
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           write(*,"(1P10E18.9)") errors
           write(*,"(1P10E18.9)") abs(dz)/ye,abs(xn-xn0)/xn0
#endif
           if (ordering(2)*ordering(3).eq.2.and.xa.gt.0.9d0) then           
              !instead of changing proton by a large amount
              !and potentially messing up neutrino stuff, change
              !neutron fractions, it is the next most abundant
              xa = (ye-xps-zbars/abars*xhs)*2.0d0
              xn = 1.0d0-xps-xa-xhs
              smallest_error = -1
#if DEBUG_FRACTIONS
              write(*,*) "changing xn instead"
#endif
              cont = .true.
           endif
           if (ordering(2)*ordering(3).eq.4.and.xa.gt.0.9d0) then           
              !instead of changing heavy by a large amount
              !and potentially messing up neutrino stuff, change
              !neutron fractions, it is the next most abundant
              xa = (ye-xps-zbars/abars*xhs)*2.0d0
              xn = 1.0d0-xps-xa-xhs
              smallest_error = -1
#if DEBUG_FRACTIONS
              write(*,*) "changing xn instead"
#endif
              cont = .true.
           endif
           if (ordering(2)*ordering(3).eq.8.and.xp.gt.0.01.and.xh.gt.0.01) then
              !attempt to change xp and xh by the same relative amount, equations follow.
              beta = (ye-xps/2.0d0-0.5d0+xns/2.0d0-(zbars/abars-0.5d0)*xhs)/(xps/2.0d0+xhs*(zbars/abars-0.5d0))
              xp = (1.0d0+beta)*xps
              xh = (1.0d0+beta)*xhs
              xa = 1.0d0-xns-xp-xh
#if DEBUG_FRACTIONS
              write(*,*) "attemping to change xp and xh at the smae time"
              write(*,"(1P10E18.9)") xn,xp,xa,xh,zbar,abar
              write(*,"(1P10E18.9)") ye-xp-0.5d0*xa-zbars/abars*xh,abs(xa-xas)/xas,abs(xp-xps)/xps,abs(xh-xhs)/xhs
#endif
              if (sqrt((abs(xa-xas)/xas)**2+(abs(xp-xps)/xps)**2+(abs(xh-xhs)/xhs)**2).lt.errors(smallest_error)) then
                 smallest_error = -1 !special case that means fractions set
              else
                 !reset to before this other attempt.
                 if (xas.gt.0.9d0) then
#if DEBUG_FRACTIONS
                    write(*,*) "pump the excess into neutrons"
#endif                   
                    smallest_error = 1
                    xp = xps
                    xh = xhs
                    xa = xas
                    cont = .true.
                 else
                    xp = xps
                    xh = xhs
                    xa = xas
                 endif
              endif
                 
           else if (errors(smallest_error).gt.0.2d0.and.abs(dz)/ye.gt.0.15d0) then
              write(*,"(a70,2P10E18.9)") "Conservation is horrible, something we just have to accept...", &
                   abs(dz)/ye,errors(smallest_error)
              cont = .true.
           else
#if DEBUG_FRACTIONS              
              write(*,*) "Hope the error is < 0.2..., alphas error > 10%"
#endif
           endif

        endif
        if(smallest_error.eq.1) then
           !alphas and neutrons give the smallest error
           xa = (ye-xps-zbars/abars*xhs)*2.0d0
           xn = 1.0d0-xps-xa-xhs
        else if (smallest_error.eq.2) then
           !protons and alphas give smallest error
           xa = (1.0d0-xns-ye+(zbars/abars-1.0d0)*xhs)*2.0d0
           xp = 1.0d0-xns-xa-xhs
        else if (smallest_error.eq.3) then
           !alphas and heavies give the smallest error
           xh = (ye-0.5d0*xps-0.5d0+xns*0.5d0)/(zbars/abars-0.5d0)
           xa = 1.0d0-xns-xps-xh
        else if (smallest_error.eq.4) then
           !zbar gives smallest error
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
        else if (smallest_error.eq.5) then
           !abar gives smallest error
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
        else if (smallest_error.eq.-1) then
#if DEBUG_FRACTIONS
           write(*,*) "special case, fractions set"
#endif
        else
           stop "whats the smallest error?"
        endif  
     else if(ordering(1).eq.4) then
        !if heavies are dominate
        !check neutron error.
        errors(:) = 101.0d0
        if (xns.gt.0.0d0) then
           xh = (ye-xps-0.5d0*xas)*abars/zbars
           xn = 1.0d0-xps-xas-xh
           heavy_error = (xh-xhs)/xhs
           neutron_error = (xn-xns)/xns
           errors(3) = sqrt(heavy_error**2+neutron_error**2)
           if (xh.lt.0.0d0) errors(3) = 100.0d0
           if (xn.lt.0.0d0) errors(3) = 100.0d0
           xh = xhs
           xn = xns
        endif
        !check proton error
        if (xps.gt.0.0d0) then
           xh = (ye-1.0d0+xns+0.5d0*xas)/(zbars/abars-1.0d0)
           xp = 1.0d0-xns-xas-xh
           heavy_error = (xh-xhs)/xhs
           proton_error = (xp-xps)/xps
           errors(3) = sqrt(heavy_error**2+proton_error**2)
           if (xh.lt.0.0d0) errors(3) = 100.0d0
           if (xp.lt.0.0d0) errors(3) = 100.0d0
           xh = xhs
           xp = xps
        endif
        !check alpha error
        if (xas.gt.0.0d0) then
           xh = (ye-0.5d0*xps-0.5d0+xns*0.5d0)/(zbars/abars-0.5d0)           
           xa = 1.0d0-xns-xps-xh
           heavy_error = (xh-xhs)/xhs
           alpha_error = (xa-xas)/xas
           errors(3) = sqrt(heavy_error**2+alpha_error**2)
           if (xh.lt.0.0d0) errors(3) = 100.0d0
           if (xa.lt.0.0d0) errors(3) = 100.0d0
           xh = xhs
           xa = xas
        endif        
        !check zbar error
        if (zbar.gt.0.0d0) then
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
           zbar_error = (zbar-zbars)/zbars
           errors(4) = abs(zbar_error)
           if (zbar.lt.0.0d0) errors(4) = 100.0d0
           if (zbar.gt.abars) errors(4) = 102.0d0
           zbar = zbars
        endif  
        !check abar error
        if (abar.gt.0.0d0) then
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
           abar_error = (abar-abars)/abars
           errors(5) = abs(abar_error)
           if (abar.lt.0.0d0) errors(5) = 100.0d0
           if (abar.lt.zbars) errors(5) = 102.0d0
           abar = abars
        endif   
        !now find and use smallest error
        smallest_error = minloc(errors,1)
        if (errors(smallest_error).gt.0.1) then
#if DEBUG_FRACTIONS
           write(*,*) "Hope error is less then 25% should check overall dz convergence here, might be bigger than error, heavies: smallest error > 10%"
#endif
        endif
        if(smallest_error.eq.1) then
           !heavies and neutrons give the smallest error
           xh = (ye-xps-0.5d0*xas)*abars/zbars
           xn = 1.0d0-xps-xas-xh
        else if (smallest_error.eq.2) then
           !heavies and protons give smallest error
           xh = (ye-1.0d0+xns+0.5d0*xas)/(zbars/abars-1.0d0)
           xp = 1.0d0-xns-xas-xh
        else if (smallest_error.eq.3) then
           !alphas and heavies give the smallest error
           xh = (ye-0.5d0*xps-0.5d0+xns*0.5d0)/(zbars/abars-0.5d0)
           xa = 1.0d0-xns-xps-xh
        else if (smallest_error.eq.4) then
           !zbar gives smallest error
           zbar = (ye-xps-0.5d0*xas)*abars/xhs
        else if (smallest_error.eq.5) then
           !abar gives smallest error
           abar = zbars*xhs/(ye-xps-0.5d0*xas)
        else
           stop "whats the smallest error?"
        endif 
     else
        stop "whats the ordering??"
     endif

     ! some checks:
     dm = mass_cons(xn,xp,xa,xh)
     dz = charge_cons(xn,xp,xa,xh,abar,zbar,ye)

     ! mass conservation
     if(abs(dm).gt.1.0d-13) then
        stop "mass conservation should not be violated at all"
     endif

     ! charge conservation
     if(abs(dz)/ye.gt.0.01) then
        write(*,*) xn,xp,xa,xh,abar,zbar
        stop "charge conservation should be better than |dz|/ye(k) = 0.01"
     endif

     !maximum change in quantities
     if (xns.gt.0.0d0) then
        if(abs(xn-xns)/xns.gt.0.25d0.and.(.not.cont)) then
           write(*,"(1P10E18.9)") ye,xns,xps,xas,xhs,abars,zbars
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           write(*,"(1P10E18.9)") abs(dz_old)/ye,abs(xn-xn0)/xn0
           stop "xn is > 25%"
        endif
     endif
     if (xps.gt.0.0d0) then
        if(abs(xp-xps)/xps.gt.0.25d0.and.(.not.cont)) then
           write(*,"(1P10E18.9)") ye,xns,xps,xas,xhs,abars,zbars
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           write(*,"(1P10E18.9)") abs(dz_old)/ye,abs(xp-xp0)/xp0
           stop "xp is > 25%"
        endif
     endif
     if (xas.gt.0.0d0) then
        if(abs(xa-xas)/xas.gt.0.25d0.and.(.not.cont)) then
           write(*,"(1P10E18.9)") ye,xns,xps,xas,xhs,abars,zbars
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           write(*,"(1P10E18.9)") abs(dz_old)/ye,abs(xa-xa0)/xa0
           stop "change in xa is > 25%"
        endif
     endif
     if (xhs.gt.0.0d0) then
        if(abs(xh-xhs)/xhs.gt.0.25d0.and.(.not.cont)) then
           write(*,"(1P10E18.9)") ye,xns,xps,xas,xhs,abars,zbars
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           stop "xh is > 25%"
        endif
     endif
     if  (zbars.gt.0.0d0) then
        if(abs(zbar-zbars)/zbars.gt.0.25d0.and.(.not.cont)) then
           write(*,"(1P10E18.9)") ye,xns,xps,xas,xhs,abars,zbars
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           write(*,"(1P10E18.9)") abs(dz_old)/ye,abs(zbar-zbar0)/zbar0
           stop "change in zbar is > 25%"
        endif
     endif
     if  (abars.gt.0.0d0) then
        if(abs(abar-abars)/abars.gt.0.25d0.and.(.not.cont)) then
           write(*,"(1P10E18.9)") ye,xn0,xp0,xa0,xh0,abar0,zbar0
           write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
           write(*,"(1P10E18.9)") abs(dz_old)/ye,abs(abar-abar0)/abar0
           stop "change in abar is > 25%"
        endif
     endif     
     
     ! positivity of mass fractions 
     if(xn < 0.0d0) then
        stop "xn < 0"
     endif
     if(xp < 0.0d0) then
        write(*,"(1P10E18.9)") ye,xn0,xp0,xa0,xh0,abar0,zbar0
        write(*,"(1P10E18.9)") ye,xn,xp,xa,xh,abar,zbar
        
        stop "xp < 0"
     endif
     if(xa < 0.0d0) then
        stop "xa < 0"
     endif
     if(xh < 0.0d0) then
        stop "xh < 0"
     endif
     
     !non NaNness
     if(xn.ne.xn) then
        stop "xn is NaN"
     endif
     if(xp.ne.xp) then
        stop "xn is NaN"
     endif
     if(xa.ne.xa) then
        stop "xn is NaN"
     endif
     if(xh.ne.xh) then
        stop "xn is NaN"
     endif
     if(zbar.ne.zbar) then
        stop "zbar is NaN"
     endif
     if(abar.ne.abar) then
        stop "abar is NaN"
     endif

  endif
     
     
end subroutine one_time_fix_comps
   
function charge_cons(xn,xp,xa,xh,abar,zbar,ye)

  implicit none
  real*8 :: charge_cons
  real*8 :: xn,xp,xa,xh,abar,zbar,ye

  charge_cons = ye - xp - 0.5d0*xa - zbar/abar*xh

end function charge_cons

function mass_cons(xn,xp,xa,xh)

  implicit none
  real*8 :: mass_cons
  real*8 :: xn,xp,xa,xh

  mass_cons = 1.0d0 - xn - xp - xa - xh

end function mass_cons

subroutine sort_comps(xn,xp,xa,xh,ordering)

  implicit none

  integer ordering(4)
  real*8 fordering(4)
  
  real*8 xn,xp,xa,xh

  ordering = 0
  fordering = 0.0d0

  fordering(1) = xn
  fordering(2) = xp
  fordering(3) = xa
  fordering(4) = xh

  ordering(1) = maxloc(fordering,1)
  fordering(ordering(1)) = -1.0d0

  ordering(2) = maxloc(fordering,1)
  fordering(ordering(2)) = -1.0d0

  ordering(3) = maxloc(fordering,1)
  fordering(ordering(3)) = -1.0d0

  ordering(4) = maxloc(fordering,1)
  fordering(ordering(4)) = -1.0d0

  fordering(1) = xn
  fordering(2) = xp
  fordering(3) = xa
  fordering(4) = xh

end subroutine sort_comps

