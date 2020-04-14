!-*-f90-*-
subroutine loweos

  use eos_table_module
  implicit none

  real*8 :: xabar, xzbar

  real*8 pion,eion,sion
  real*8 aion(4),zion(4),xion(4)

  real*8 xrho,xye,xtemp
  integer irho,iye,itemp

  real*8 xh,xa,xn,xp,mu_n,mu_p
  real*8 xd,xt,x3he,x4li,xL
  real*8 xxne,xdxnedd,xdxnedt,xdxnedz,xdxneda
  real*8 pcoul,ecoul,scoul

  real*8 pmatch,pep
  real*8 ematch,epep
  real*8 smatch,sep

  real*8, allocatable:: energy_offset(:,:)
  real*8, allocatable:: mun_offset(:,:)
  real*8, allocatable:: mup_offset(:,:)

  write(6,*) "##################################"
  write(6,*) "irho_start_nuclear:  ",irho_start_nuclear
  write(6,*) "irho_start_comps: ",irho_start_comps

  ! for irho < irho_start_comps, keep composition of irho_start_comps
  ! at irho_start_nuclear: 
  ! (i) get energy offset between loweos and nuclear EOS
  ! (ii) average things
  ! (iii) get chemical potential offset 

  allocate(energy_offset(nlogtemp,nye))
  allocate(mun_offset(nlogtemp,nye))
  allocate(mup_offset(nlogtemp,nye))

  energy_offset = 0.0d0
  mun_offset = 0.0d0
  mup_offset = 0.0d0
  ! handle transition first:
  do iye = 1,nye
     do itemp = 1,nlogtemp
        irho = irho_start_nuclear
        xrho = 10.0d0**logrho(irho)
        xtemp = 10.0d0**logtemp(itemp) * temp_mev_to_kelvin
        xye = ye(iye)
        
        xabar = eos_table(irho,itemp,iye,iabar)
        xzbar = eos_table(irho,itemp,iye,izbar)
        xh = eos_table(irho,itemp,iye,ixh)
        xa = eos_table(irho,itemp,iye,ixa)
        xn = eos_table(irho,itemp,iye,ixn)
        xp = eos_table(irho,itemp,iye,ixp)
        if (Hempel) then
           !for ion EOS, put lights in np, there shouldn't be any at these low densities
           xd = eos_table(irho,itemp,iye,ixd)
           xt = eos_table(irho,itemp,iye,ixt)
           x3he = eos_table(irho,itemp,iye,ix3he)
           x4li = eos_table(irho,itemp,iye,ix4li)
           xp = xp + 0.5d0*xd + xt/3.0d0 + 2.0d0*x3he/3.0d0 + 3.0d0*x4li/4.0d0
           xn = xn + 0.5d0*xd + 2.0d0*xt/3.0d0 + x3he/3.0d0 + x4li/4.0d0
        endif
        if (Hyperon) then
           !for ion EOS, put hyperons in n, there shouldn't be any at these low densities
           xL = eos_table(irho,itemp,iye,ixL)
           xn = xn + xL
        endif
        

        ! calculate the offset of mu_n and mu_p at transition density
        ! this mainly to take care of the low region T < 0.1 MeV
        if(anal_mu_lowden) then
           call mu_np(xrho,xtemp,xn,xp,xye,mu_n,mu_p)
           mun_offset(itemp,iye) =  eos_table(irho_start_nuclear,itemp,iye,imun) - mu_n
           mup_offset(itemp,iye) = eos_table(irho_start_nuclear,itemp,iye,imup) - mu_p
        endif
        
        xion(1) =  xn
        xion(2) =  xp
        xion(3) =  xa
        xion(4) =  xh
        
        aion(1) = 1.0
        zion(1) = 0.0
        
        aion(2) = 1.0
        zion(2) = 1.0
        
        aion(3) = 4.0
        zion(3) = 2.0
        
        aion(4) = xabar
        zion(4) = xzbar

        xxne = electron_xne(irho,itemp,iye)
        xdxnedd = electron_dxnedd(irho,itemp,iye)
        xdxnedt = electron_dxnedt(irho,itemp,iye)
        xdxnedz = electron_dxnedz(irho,itemp,iye)
        xdxneda = electron_dxneda(irho,itemp,iye)
        
        call ioneos(xrho,xtemp,xion,aion,zion,pion,eion,sion,&
             xxne,xdxnedd,xdxnedt,xdxnedz,xdxneda,pcoul,ecoul,scoul)
        
        sion = sion / kb_erg * amu_cgs
        scoul = scoul / kb_erg * amu_cgs

        eion = eion+ecoul
        pion = pion+pcoul
        sion = sion+scoul

        epep = eos_electrons(irho,itemp,iye,electron_ieps) &
                + eos_photons(irho,itemp,iye,photon_ieps) 

        ematch = eos_table(irho,itemp,iye,ienergy) + &
             eos_electrons(irho,itemp,iye,electron_ieps) &
             + eos_photons(irho,itemp,iye,photon_ieps)

        pep = eos_electrons(irho,itemp,iye,electron_ipress) &
                + eos_photons(irho,itemp,iye,photon_ipress) 

        pmatch = eos_table(irho,itemp,iye,ipress) + &
             eos_electrons(irho,itemp,iye,electron_ipress) &
             + eos_photons(irho,itemp,iye,photon_ipress)

        sep = eos_electrons(irho,itemp,iye,electron_ientropy) &
                + eos_photons(irho,itemp,iye,photon_ientropy) 

        smatch = eos_table(irho,itemp,iye,ientropy) + &
             eos_electrons(irho,itemp,iye,electron_ientropy) &
             + eos_photons(irho,itemp,iye,photon_ientropy)

        energy_offset(itemp,iye) = ematch - (epep+eion)

        eos_table(irho,itemp,iye,ienergy) =  &
             0.5d0 * (ematch + (epep+eion+energy_offset(itemp,iye) ))
        
        eos_table(irho,itemp,iye,ipress) = &
             0.5d0 * (pmatch + (pep+pion))

        eos_table(irho,itemp,iye,ientropy) = &
             0.5d0 * (smatch + (sep+sion))

        eos_table(irho,itemp,iye,imue) = eos_electrons(irho,itemp,iye,electron_imue)

     enddo
  enddo


  do iye = 1,nye
     do itemp = 1,nlogtemp
        do irho=1,irho_start_nuclear-1

           xrho = 10.0d0**logrho(irho)
           xtemp = 10.0d0**logtemp(itemp) * temp_mev_to_kelvin
           xye = ye(iye)
           
           if(irho.ge.irho_start_comps) then
              xabar = eos_table(irho,itemp,iye,iabar)
              xzbar = eos_table(irho,itemp,iye,izbar)
              xh = eos_table(irho,itemp,iye,ixh)
              xa = eos_table(irho,itemp,iye,ixa)
              xn = eos_table(irho,itemp,iye,ixn)
              xp = eos_table(irho,itemp,iye,ixp)
              if (Hempel) then
                 !for ion EOS, put lights in np, there shouldn't be any at these low densities
                 xd = eos_table(irho,itemp,iye,ixd)
                 xt = eos_table(irho,itemp,iye,ixt)
                 x3he = eos_table(irho,itemp,iye,ix3he)
                 x4li = eos_table(irho,itemp,iye,ix4li)
                 xp = xp + 0.5d0*xd + xt/3.0d0 + 2.0d0*x3he/3.0d0 + 3.0d0*x4li/4.0d0
                 xn = xn + 0.5d0*xd + 2.0d0*xt/3.0d0 + x3he/3.0d0 + x4li/4.0d0
              endif
              if (Hyperon) then
                 !for ion EOS, put hyperons in n, there shouldn't be any at these low densities
                 xL = eos_table(irho,itemp,iye,ixL)
                 xn = xn + xL
              endif
           else
              xabar = eos_table(irho_start_comps,itemp,iye,iabar)
              xzbar = eos_table(irho_start_comps,itemp,iye,izbar)
              xh = eos_table(irho_start_comps,itemp,iye,ixh)
              xa = eos_table(irho_start_comps,itemp,iye,ixa)
              xn = eos_table(irho_start_comps,itemp,iye,ixn)
              xp = eos_table(irho_start_comps,itemp,iye,ixp)
              if (Hempel) then
                 !for ion EOS, put lights in np, there shouldn't be any at these low densities
                 xd = eos_table(irho_start_comps,itemp,iye,ixd)
                 xt = eos_table(irho_start_comps,itemp,iye,ixt)
                 x3he = eos_table(irho_start_comps,itemp,iye,ix3he)
                 x4li = eos_table(irho_start_comps,itemp,iye,ix4li)
                 xp = xp + 0.5d0*xd + xt/3.0d0 + 2.0d0*x3he/3.0d0 + 3.0d0*x4li/4.0d0
                 xn = xn + 0.5d0*xd + 2.0d0*xt/3.0d0 + x3he/3.0d0 + x4li/4.0d0
              endif
              if (Hyperon) then
                 !for ion EOS, put hyperons in n, there shouldn't be any at these low densities
                 xL = eos_table(irho,itemp,iye,ixL)
                 xn = xn + xL
              endif
           endif

           xion(1) =  xn
           xion(2) =  xp
           xion(3) =  xa
           xion(4) =  xh
     
           aion(1) = 1.0
           zion(1) = 0.0
           
           aion(2) = 1.0
           zion(2) = 1.0
           
           aion(3) = 4.0
           zion(3) = 2.0
     
           aion(4) = xabar
           zion(4) = xzbar
           
           xxne = electron_xne(irho,itemp,iye)
           xdxnedd = electron_dxnedd(irho,itemp,iye)
           xdxnedt = electron_dxnedt(irho,itemp,iye)
           xdxnedz = electron_dxnedz(irho,itemp,iye)
           xdxneda = electron_dxneda(irho,itemp,iye)
     
           call ioneos(xrho,xtemp,xion,aion,zion,pion,eion,sion,&
                xxne,xdxnedd,xdxnedt,xdxnedz,xdxneda,pcoul,ecoul,scoul)
           
           sion = sion / kb_erg * amu_cgs
           scoul = scoul / kb_erg * amu_cgs

           eion = eion+ecoul
           pion = pion+pcoul
           sion = sion+scoul
           
           epep = eos_electrons(irho,itemp,iye,electron_ieps) &
                + eos_photons(irho,itemp,iye,photon_ieps) 
      
           pep = eos_electrons(irho,itemp,iye,electron_ipress) &
                + eos_photons(irho,itemp,iye,photon_ipress) 

           sep = eos_electrons(irho,itemp,iye,electron_ientropy) &
                + eos_photons(irho,itemp,iye,photon_ientropy) 

           ! get analytical mu_n and mu_p 
           if(anal_mu_lowden) then
             call mu_np(xrho,xtemp,xn,xp,xye,mu_n,mu_p)
             eos_table(irho,itemp,iye,imun) = mu_n + mun_offset(itemp,iye)
             eos_table(irho,itemp,iye,imup) = mu_p + mup_offset(itemp,iye)
           endif


           eos_table(irho,itemp,iye,ienergy) =  &
                epep+eion+energy_offset(itemp,iye)

           if (itemp.eq.itemp_start_nuclear) then
              write(*,*) xrho,xye,epep/1e18,eion/1e18,energy_offset(itemp,iye)/1e18
           endif

           eos_table(irho,itemp,iye,ipress) = &
                (pep+pion)

           eos_table(irho,itemp,iye,ientropy) = &
                (sep+sion)

           eos_table(irho,itemp,iye,imue) = eos_electrons(irho,itemp,iye,electron_imue)

           if(irho.ge.irho_start_comps) then
              if(.not.anal_mu_lowden) then
                 eos_table(irho,itemp,iye,imun) = eos_table(irho,itemp,iye,imun)
                 eos_table(irho,itemp,iye,imup) = eos_table(irho,itemp,iye,imup)
              endif
              eos_table(irho,itemp,iye,iabar) = eos_table(irho,itemp,iye,iabar)
              eos_table(irho,itemp,iye,izbar) = eos_table(irho,itemp,iye,izbar)
              eos_table(irho,itemp,iye,ixn) = eos_table(irho,itemp,iye,ixn)
              eos_table(irho,itemp,iye,ixp) = eos_table(irho,itemp,iye,ixp)
              eos_table(irho,itemp,iye,ixa) = eos_table(irho,itemp,iye,ixa)
              eos_table(irho,itemp,iye,ixh) = eos_table(irho,itemp,iye,ixh)
              if (Hempel) then
                 eos_table(irho,itemp,iye,ixd) = eos_table(irho,itemp,iye,ixd)
                 eos_table(irho,itemp,iye,ixt) = eos_table(irho,itemp,iye,ixt)
                 eos_table(irho,itemp,iye,ix3he) = eos_table(irho,itemp,iye,ix3he)
                 eos_table(irho,itemp,iye,ix4li) = eos_table(irho,itemp,iye,ix4li)
              endif
              if (Hyperon) then
                 eos_table(irho,itemp,iye,ixL) = eos_table(irho,itemp,iye,ixL)
              endif
           else
              if(.not.anal_mu_lowden) then
                 eos_table(irho,itemp,iye,imun) = eos_table(irho_start_comps,itemp,iye,imun)
                 eos_table(irho,itemp,iye,imup) = eos_table(irho_start_comps,itemp,iye,imup)
              endif
              eos_table(irho,itemp,iye,iabar) = eos_table(irho_start_comps,itemp,iye,iabar)
              eos_table(irho,itemp,iye,izbar) = eos_table(irho_start_comps,itemp,iye,izbar)
              eos_table(irho,itemp,iye,ixn) = eos_table(irho_start_comps,itemp,iye,ixn)
              eos_table(irho,itemp,iye,ixp) = eos_table(irho_start_comps,itemp,iye,ixp)
              eos_table(irho,itemp,iye,ixa) = eos_table(irho_start_comps,itemp,iye,ixa)
              eos_table(irho,itemp,iye,ixh) = eos_table(irho_start_comps,itemp,iye,ixh)
              if (Hempel) then
                 eos_table(irho,itemp,iye,ixd) = eos_table(irho_start_comps,itemp,iye,ixd)
                 eos_table(irho,itemp,iye,ixt) = eos_table(irho_start_comps,itemp,iye,ixt)
                 eos_table(irho,itemp,iye,ix3he) = eos_table(irho_start_comps,itemp,iye,ix3he)
                 eos_table(irho,itemp,iye,ix4li) = eos_table(irho_start_comps,itemp,iye,ix4li)
              endif
              if (Hyperon) then
                 eos_table(irho,itemp,iye,ixL) = eos_table(irho_start_comps,itemp,iye,ixL)
              endif
           endif
           
        enddo
     enddo
  enddo

contains

  subroutine ioneos(den,temp,xion,aion,zion,pion,eion,sion,&
       xne,dxnedd,dxnedt,dxnedz,dxneda,pcoul,ecoul,scoul)

    ! This outine is based on Frank Timmes'
    ! public_timmes.f

    use eos_table_module,only : kb_erg
    implicit none

    real*8 den,temp
    real*8 xion(4),aion(4)
    real*8 zion(4),pion,eion,sion

    real*8 abar,zbar,ymass,zbarxx

    ! internal vars
    real*8 kt,ktinv,deninv,tempinv
    real*8 xni,dxnidd,dxnidt,dxnida,dxnidz
    real*8 ytot1 
    real*8 dpiondd,dpiondt,dpionda,dpiondz
    real*8 deiondd,deiondt,deionda,deiondz
    real*8 y,yy,z,sifac,etaion

    ! for Coulomb corrections
    real*8 plasg,xne,dxnedd,dxnedt,dxneda,dxnedz,&
         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,&
         ecoul,decouldd,decouldt,decoulda,decouldz,&
         scoul,dscouldd,dscouldt,dscoulda,dscouldz  

    integer i

    sifac = 8.6322745944370191d-45

    pion = 0.0d0
    eion = 0.0d0
    sion = 0.0d0

    ! calculate abar and zbar
    ytot1 = 0.0d0
    zbarxx = 0.0d0
    do i=1,4
       ymass = xion(i)/(max(aion(i),1.0d0))
       ytot1 = ytot1 + ymass
       zbarxx = zbarxx + zion(i) * ymass
    enddo
    abar = 1.0d0/ytot1
    zbar = zbarxx * abar

    kt = kb_erg * temp
    ktinv = 1.0d0/kt
    deninv = 1.0d0/den
    tempinv = 1.0d0/temp

    ytot1 = 1.0d0/abar
    xni     = avo * ytot1 * den
    dxnidd  = avo * ytot1
    dxnidt  = 0.0d0
    dxnida  = -xni * ytot1
    dxnidz  = 0.0d0 
    
    pion    = xni * kt 
    dpiondd = dxnidd * kt 
    dpiondt = xni * kb_erg
    dpionda = -pion * ytot1 
    dpiondz = 0.0d0
    
    eion    = 1.5d0 * pion*deninv 
    deiondd = (1.5d0 * dpiondd - eion)*deninv 
    deiondt = 1.5d0 * dpiondt*deninv 
    deionda = 1.5d0 * dpionda*deninv 
    deiondz = 0.0d0
    
    y         = 1.0d0/(abar*kt)
    yy        = y * sqrt(y)
    z         = xni * sifac * yy
    etaion    = log(z)
    sion    = (eion + pion*deninv)*tempinv - etaion * kb_erg*avo*ytot1


    ! for the coulomb part, we need information about the electrons
    call coulomb(den,temp,abar,zbar,pion,dpiondd,dpiondt,dpionda,dpiondz,&
         xne,dxnedd,dxnedt,dxneda,dxnedz,&
         plasg, &
         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,&
         ecoul,decouldd,decouldt,decoulda,decouldz,&
         scoul,dscouldd,dscouldt,dscoulda,dscouldz)
  

  end subroutine ioneos

  subroutine coulomb(den,temp,abar,zbar,&
       pion,dpiondd,dpiondt,dpionda,dpiondz,&
       xne, dxnedd,dxnedt,dxneda,dxnedz,&
       plasg,&
       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,&
       ecoul,decouldd,decouldt,decoulda,decouldz,&
       scoul,dscouldd,dscouldt,dscoulda,dscouldz)
    ! This is a F90 adaptation of the coulomb corrections
    ! subroutine in public_timmes.f
    implicit none
!  plasg    = ratio of electrostatic to thermal energy
!  pcoul    = coulomb pressure correction
!  coulmult = coulomb component multiplier
!  dpcouldd = derivative of the coulomb pressure with density
!  dpcouldt = derivative of the coulomb pressure with temperature
!  dpcoulda = derivative of the coulomb pressure with abar
!  dpcouldz = derivative of the coulomb pressure with zbar

!  ecoul    = coulomb energy correction
!  decouldd = derivative of the coulomb energy with density
!  decouldt = derivative of the coulomb energy with temperature
!  decoulda = derivative of the coulomb energy with abar
!  decouldz = derivative of the coulomb energy with zbar

!  scoul    = coulomb entropy correction
!  dscouldd = derivative of the coulomb entropy with density
!  dscouldt = derivative of the coulomb entropy with temperature
!  dscoulda = derivative of the coulomb entropy with abar
!  dscouldz = derivative of the coulomb entropy with zbar

    real*8 den,temp,abar,zbar,&
         pion,dpiondd,dpiondt,dpionda,dpiondz,&
         plasg,xne,dxnedd,dxnedt,dxneda,dxnedz,&
         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,&
         ecoul,decouldd,decouldt,decoulda,decouldz,&
         scoul,dscouldd,dscouldt,dscoulda,dscouldz  

!  local variables
!  for the uniform background coulomb correction
     real*8 ytot1,kt,ktinv,&
          s,dsdd,dsdt,dsda,dsdz,sinv,&
          aele,daeledd,daeledt,daeleda,daeledz,aeleinv,&
          eplasg,eplasgdd,eplasgdt,eplasgda,eplasgdz,&
          aion,&
          lami,inv_lami,lamidd,lamida,lamidz,&
          plasgdd,plasgdt,plasgda,plasgdz,&
          x,y,z


      real*8 :: u0,rho6
      real*8,parameter :: a1 = -0.898004d0
      real*8,parameter :: b1 = 0.96786d0
      real*8,parameter :: c1 = 0.220703d0
      real*8,parameter :: d1 = -0.86097d0
      real*8,parameter :: e1 = 2.5269d0
      real*8,parameter :: a2 = 0.29561d0
      real*8,parameter :: b2 = 1.9885d0
      real*8,parameter :: c2 = 0.288675d0


!  various constants
      real*8, parameter :: third  = 1.0d0/3.0d0
      real*8, parameter :: forth  = 4.0d0/3.0d0
      real*8, parameter :: fiveth = 5.0d0/3.0d0
      real*8, parameter :: qe = 4.8032068d-10
      real*8, parameter :: esqu = qe*qe
      real*8, parameter :: pi = 3.1415926535897932384d0
      real*8, parameter :: kerg = 1.380658d-16
      real*8, parameter :: forthpi = forth*pi

!  common variables
       ytot1   = 1.0d0/abar
       kt      = kerg * temp
       ktinv   = 1.0d0/kt

!  yakovlev & shalybkov eqs 5, 9 and 10
       s        = forthpi * xne
       dsdd     = forthpi * dxnedd
       dsdt     = forthpi * dxnedt
       dsda     = forthpi * dxneda
       dsdz     = forthpi * dxnedz
       sinv     = 1.0d0/s

!  electron-sphere radius aele
       aele     = sinv**third
       z        = -third * aele * sinv
       daeledd  = z * dsdd
       daeledt  = z * dsdt
       daeleda  = z * dsda
       daeledz  = z * dsdz
       aeleinv  = 1.0d0/aele

!  electron coupling parameter eplasg
       eplasg   = esqu * ktinv * aeleinv
       z        = -eplasg * aeleinv
       eplasgdd = z * daeledd
       eplasgdt = z * daeledt - eplasg*ktinv*kerg
       eplasgda = z * daeleda
       eplasgdz = z * daeledz

!  ion-sphere radius aion
       x        = zbar**third
       aion     = x * aele

!  ion coupling parameter plasg
       z        = x*x*x*x*x
       plasg    = z * eplasg
       plasgdd  = z * eplasgdd
       plasgdt  = z * eplasgdt
       plasgda  = z * eplasgda
       plasgdz  = z * eplasgdz + fiveth*x*x * eplasg

!  yakovlev & shalybkov 1989 equations 82, 85, 86, 87
       if (plasg .ge. 1.0) then
        x        = plasg**(0.25d0) 
        u0       = a1*plasg + b1*x + c1/x + d1
        ecoul    = pion/den * u0
        pcoul    = third * ecoul * den
        scoul    = -avo*ytot1*kerg * &
             (3.0d0*b1*x - 5.0d0*c1/x &
             + d1 * (log(plasg) - 1.0d0) - e1)

        y        = avo/abar*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        decouldd = y * plasgdd 
        decouldt = y * plasgdt + ecoul/temp
        decoulda = y * plasgda - ecoul/abar
        decouldz = y * plasgdz

        y        = third * den
        dpcouldd = third * ecoul + y*decouldd
        dpcouldt = y * decouldt
        dpcoulda = y * decoulda
        dpcouldz = y * decouldz


        y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x +1.25d0*c1/x +d1)
        dscouldd = y * plasgdd
        dscouldt = y * plasgdt
        dscoulda = y * plasgda - scoul/abar
        dscouldz = y * plasgdz


!  yakovlev & shalybkov 1989 equations 102, 103, 104
       else if (plasg .lt. 1.0) then
        x        = plasg*sqrt(plasg)
        y        = plasg**b2
        z        = c2 * x - third * a2 * y
        pcoul    = -pion * z
        ecoul    = 3.0d0 * pcoul/den
        scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

        s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        dpcouldd = -dpiondd*z - pion*s*plasgdd
        dpcouldt = -dpiondt*z - pion*s*plasgdt
        dpcoulda = -dpionda*z - pion*s*plasgda
        dpcouldz = -dpiondz*z - pion*s*plasgdz

        s        = 3.0d0/den
        decouldd = s * dpcouldd - ecoul/den
        decouldt = s * dpcouldt
        decoulda = s * dpcoulda
        decouldz = s * dpcouldz

        s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x -a2*(b2-1.0d0)*y)
        dscouldd = s * plasgdd
        dscouldt = s * plasgdt
        dscoulda = s * plasgda - scoul/abar
        dscouldz = s * plasgdz
       end if



  end subroutine coulomb

  subroutine mu_np(xrho,xtemp,xn,xp,xye,mu_n,mu_p)
    ! calculate chemical potentials of neutrons and protons at low densities
    use eos_table_module, only: m_n,m_p,mev_to_erg,erg_to_mev,kb_erg,planck,pi,clight,avo
    implicit none
    
    real*8 :: xrho,xtemp,xn,xp,xye
    real*8 :: mu_n,mu_p
    real*8 :: n_n,n_p,a_e,Gamma_e,Gamma_p,mu_p_coul
    real*8 :: fac_1_erg,fac_1_mev,fac_2
    real*8,  parameter :: A_1 = -0.9052d0, &
         &                  A_2 = 0.6322d0,  &
         &                  A_3 = - sqrt(3.0d0)/2.0d0 - A_1/sqrt(A_2)
    real*8, parameter  :: m_n_cgs = m_n*mev_to_erg/clight**2,&
         &                  m_p_cgs = m_p*mev_to_erg/clight**2
    real*8, parameter  :: electron_charge = 4.803204184d-10
    real*8, parameter  :: tiny = 1.0d-20
    
    fac_1_erg = kb_erg*xtemp
    fac_1_mev = kb_erg*xtemp*erg_to_mev
    fac_2 = (planck**2/(2.0*pi*fac_1_erg))**1.5
    
    if(xn.le.0.0d0) then
       xn = tiny
    endif
    if(xp.le.0.0d0) then
       xp = tiny
    endif
    n_n = xn*xrho/m_n_cgs
    n_p = xp*xrho/m_p_cgs
    
    a_e = (4.0d0/3.0d0*pi*xrho*xye*avo)**(-1.0d0/3.0d0)
    Gamma_e = electron_charge**2/(fac_1_erg*a_e)
    Gamma_p = Gamma_e
    
    ! chemical potentials for Boltzmann gas (in unit of MeV)
    ! include rest mass but normalized by free neutron mass
    ! when n_n = 0 or n_p = 0, set mu_n = 0.0 and mu_p = 0, correspondingly
    ! 
    if(n_n.gt.0.0d0) then
       mu_n = fac_1_mev*log(n_n/2.0d0*fac_2/m_n_cgs**1.5) + m_n - m_n
    else 
       mu_n = 0.0d0
    endif
    if(n_p.gt.0.0d0) then
       mu_p = fac_1_mev*log(n_p/2.0d0*fac_2/m_p_cgs**1.5) + m_p - m_n
    else
       mu_p = 0.0d0
    endif
    ! coulomb correction for charged particles from Chabrier & Potkhin (1998)
    mu_p_coul = fac_1_mev * (A_1*(sqrt(Gamma_p*(A_2 + Gamma_p)) &
         - A_2 * log(sqrt(Gamma_p/A_2) + sqrt(1.0d0+Gamma_p/A_2))) &
         + 2.0d0*A_3*(sqrt(Gamma_p) - atan(sqrt(Gamma_p))))

    mu_p = mu_p + mu_p_coul
  end subroutine mu_np

end subroutine loweos
