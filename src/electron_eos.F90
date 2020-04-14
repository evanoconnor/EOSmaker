!-*-f90-*-
subroutine electron_photon_eos

  use eos_table_module
  implicit none

  integer iii,jjj,kkk,i,j,k

  integer ixnlogrho,jxnlogtemp,kxnye
  integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer OMP_GET_NUM_PROCS
  integer mythread
  real*8 abar,zbar,xrho,xtemp,pele,ppos,eele,epos,sele,&
       spos,prad,erad,srad,gamma1,dedT1,dpdT1,dpdrhoT1
  real*8 xxne,xdxnedd,xdxnedt,xdxnedz,xdxneda
  real*8 muep,mue
  real*8 ytot1,zbarxx

  if (noelectrons.eqv..true.) return

  write(6,*) "Adding Electrons, Positrons and Photons :-))))"
  
  allocate(electron_gamma(nlogrho,nlogtemp,nye))
  allocate(electron_dedT(nlogrho,nlogtemp,nye))
  allocate(electron_dpdT(nlogrho,nlogtemp,nye))
  allocate(electron_dpdrhoT(nlogrho,nlogtemp,nye))

  allocate(electron_xne(nlogrho,nlogtemp,nye))
  allocate(electron_dxnedd(nlogrho,nlogtemp,nye))
  allocate(electron_dxnedt(nlogrho,nlogtemp,nye))
  allocate(electron_dxnedz(nlogrho,nlogtemp,nye))
  allocate(electron_dxneda(nlogrho,nlogtemp,nye))

  if(use_saved_electrons) then
     write(6,*) "Using saved electron EOS data"
     open(unit=473,file="electrons_photons.bin.dat",status="unknown",form="unformatted")
     read(473) ixnlogrho,jxnlogtemp,kxnye
     if(ixnlogrho.ne.nlogrho .or. &
          jxnlogtemp.ne.nlogtemp .or. &
          kxnye.ne.nye) then
        write(6,*) ixnlogrho,jxnlogtemp,kxnye
        write(6,*) nlogrho,nlogtemp,nye
        call flush(6)
        stop "Supplied electron_photns.bin.dat not consistent with expected table dimensions!"
     endif

     read(473) eos_electrons
     read(473) eos_photons
     read(473) electron_gamma
     read(473) electron_dedT
     read(473) electron_dpdT
     read(473) electron_dpdrhoT
     read(473) electron_xne
     read(473) electron_dxnedd
     read(473) electron_dxnedt
     read(473) electron_dxnedz
     read(473) electron_dxneda
     close(473)

     write(*,*) "Read in Electron & photon table"

     if (noelectrons) then
        eos_electrons(i,j,k,electron_ipress) = 0.0d0
        eos_photons(i,j,k,photon_ipress) = 0.0d0

        eos_electrons(i,j,k,electron_ieps) = 0.0d0
        eos_photons(i,j,k,photon_ieps) = 0.0d0

        eos_electrons(i,j,k,electron_ientropy) = 0.0d0
        eos_photons(i,j,k,photon_ientropy) = 0.0d0

        eos_table(i,j,k,imue) = 0.0d0

     endif

     goto 474

  endif

  !if you don't have electrons saved then will have to calculate
  !it... will take most of the time.  This routine will write out a
  !table so next time you can use the save version
  do kkk=nye,1,-1
     do jjj=1,nlogtemp
        write(6,"(A20,1P1E15.6,A6,1P1E15.6)") "Working on Ye = ", &
             ye(kkk), " Temp = ",logtemp(jjj)
        do iii=1,nlogrho
           if(iii.ge.irho_start_comps) then
              
           ! compute abar and zbar
              if (Hempel.and..not.Hyperon) then
                 ytot1 = eos_table(iii,jjj,kkk,ixn)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixp)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar)) + &
                      eos_table(iii,jjj,kkk,ixd)/2.0d0 + &
                      eos_table(iii,jjj,kkk,ixt)/3.0d0 + &
                      eos_table(iii,jjj,kkk,ix3he)/3.0d0 + &
                      eos_table(iii,jjj,kkk,ix4li)/4.0d0
                 
                 zbarxx = 0.0d0 + eos_table(iii,jjj,kkk,ixp)*1.0d0 + &
                      2.0d0 * eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      1.0d0 * eos_table(iii,jjj,kkk,ixd)/2.0d0 + &
                      1.0d0 * eos_table(iii,jjj,kkk,ixt)/3.0d0 + &
                      2.0d0 * eos_table(iii,jjj,kkk,ix3he)/3.0d0 + &
                      3.0d0 * eos_table(iii,jjj,kkk,ix4li)/4.0d0 + &
                      eos_table(iii,jjj,kkk,izbar) &
                      * eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar))
              elseif (Hempel.and.Hyperon) then
                 ytot1 = eos_table(iii,jjj,kkk,ixn)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixp)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar)) + &
                      eos_table(iii,jjj,kkk,ixd)/2.0d0 + &
                      eos_table(iii,jjj,kkk,ixt)/3.0d0 + &
                      eos_table(iii,jjj,kkk,ix3he)/3.0d0 + &
                      eos_table(iii,jjj,kkk,ix4li)/4.0d0+ &
                      eos_table(iii,jjj,kkk,ixL)

                 zbarxx = 0.0d0 + eos_table(iii,jjj,kkk,ixp)*1.0d0 + &
                      2.0d0 * eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      1.0d0 * eos_table(iii,jjj,kkk,ixd)/2.0d0 + &
                      1.0d0 * eos_table(iii,jjj,kkk,ixt)/3.0d0 + &
                      2.0d0 * eos_table(iii,jjj,kkk,ix3he)/3.0d0 + &
                      3.0d0 * eos_table(iii,jjj,kkk,ix4li)/4.0d0 + &
                      eos_table(iii,jjj,kkk,izbar) & 
                      * eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar))
              else if (Hyperon) then
                 ytot1 = eos_table(iii,jjj,kkk,ixn)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixp)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar)) + &
                      eos_table(iii,jjj,kkk,ixL)/1.0d0
                 
                 zbarxx = 0.0d0 + eos_table(iii,jjj,kkk,ixp)*1.0d0 + &
                      2.0d0 * eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      eos_table(iii,jjj,kkk,izbar) & 
                      * eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar))

              else             
                 ytot1 = eos_table(iii,jjj,kkk,ixn)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixp)/1.0d0 + &
                      eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar))
                 
                 zbarxx = 0.0d0 + eos_table(iii,jjj,kkk,ixp)*1.0d0 + &
                      2.0d0 * eos_table(iii,jjj,kkk,ixa)/4.0d0 + &
                      eos_table(iii,jjj,kkk,izbar) & 
                      * eos_table(iii,jjj,kkk,ixh)/(eos_table(iii,jjj,kkk,iabar))
              endif
              abar = 1.0d0/ytot1
              zbar = zbarxx * abar

              if( abs(ye(kkk)-zbarxx).gt.1.0d-13) then
                 if(abs(ye(kkk)-zbarxx).le.0.01d0) then
                    zbarxx = ye(kkk)
                    zbar = zbarxx*abar
                 else
                    write(6,"(i4,i4,i4,1P10E15.6)") iii,jjj,kkk,ye(kkk),zbarxx,ye(kkk)-zbarxx
                    stop "Charge conservation violated!"
                 endif
              endif

           else
              !at these low density (irho < irho_start_comps) assume
              !ni56, the ion eos shouldn't make a big difference down
              !here...
              
              abar = 56.0d0
              zbar = 56.0d0*ye(kkk)
           endif

           xrho = 10.0d0**logrho(iii)
           xtemp = 10.0d0**logtemp(jjj) * temp_mev_to_kelvin

           call wrap_timmes(abar,zbar,xrho,xtemp,pele,ppos,eele,epos,sele,spos,&
                prad,erad,srad,mue,muep,gamma1,dedT1,dpdT1,dpdrhoT1,&
                xxne,xdxnedd,xdxnedt,xdxnedz,xdxneda)

           srad = srad / kb_erg * amu_cgs
           spos = spos / kb_erg * amu_cgs
           sele = sele / kb_erg * amu_cgs

           ! note that chemical potential of Timmes is: eta =
           ! mu_e/(k_B T) without restmass, we need: eta*k_B*T +
           ! 0.511d0

           mue = mue*10.0d0**logtemp(jjj) + 0.511d0

           eos_electrons(iii,jjj,kkk,electron_ipress) = &
                (pele + ppos)

           eos_electrons(iii,jjj,kkk,electron_ieps) = &
                (eele + epos)

           eos_electrons(iii,jjj,kkk,electron_ientropy) = &
                (sele + spos)

           eos_electrons(iii,jjj,kkk,electron_imue) = &
                (mue)

           eos_electrons(iii,jjj,kkk,electron_imuep) = &
                (-mue)

           eos_photons(iii,jjj,kkk,photon_ipress) = &
                (prad)

           eos_photons(iii,jjj,kkk,photon_ieps) = &
                (erad)

           eos_photons(iii,jjj,kkk,photon_ientropy) = &
                (srad)

           electron_gamma(iii,jjj,kkk) = gamma1
           electron_dedT(iii,jjj,kkk) = dedT1
           electron_dpdT(iii,jjj,kkk) = dpdT1
           electron_dpdrhoT(iii,jjj,kkk) = dpdrhoT1

           electron_xne(iii,jjj,kkk) = xxne
           electron_dxnedd(iii,jjj,kkk) = xdxnedd
           electron_dxnedt(iii,jjj,kkk) = xdxnedt
           electron_dxnedz(iii,jjj,kkk) = xdxnedz
           electron_dxneda(iii,jjj,kkk) = xdxneda

        enddo
     enddo
  enddo

  ! write electron and photon tables for next time
  write(*,*) "Writing electrons to table for next time"
  open(unit=473,file="electrons_photons.bin.dat",status="unknown",form="unformatted")
  write(473) nlogrho,nlogtemp,nye
  write(473) eos_electrons
  write(473) eos_photons
  write(473) electron_gamma
  write(473) electron_dedT
  write(473) electron_dpdT
  write(473) electron_dpdrhoT
  write(473) electron_xne
  write(473) electron_dxnedd
  write(473) electron_dxnedt
  write(473) electron_dxnedz
  write(473) electron_dxneda

  close(473)

474 continue

  write(*,*) "lets add these explicitly to the high density EOS table"
  ! add in electrons and photons to high density EOS component
  ! (irho_start_nuclear is the matching point, it is done in loweos)
  do k=1,nye
     do j=1,nlogtemp
        do i=irho_start_nuclear+1,nlogrho

           eos_table(i,j,k,ipress) = eos_table(i,j,k,ipress) &
                + eos_electrons(i,j,k,electron_ipress) &
                + eos_photons(i,j,k,photon_ipress)

           eos_table(i,j,k,ienergy) = eos_table(i,j,k,ienergy)  &
                + eos_electrons(i,j,k,electron_ieps) &
                + eos_photons(i,j,k,photon_ieps)

           eos_table(i,j,k,ientropy) = eos_table(i,j,k,ientropy) &
                + eos_electrons(i,j,k,electron_ientropy) & 
                + eos_photons(i,j,k,photon_ientropy)

           eos_table(i,j,k,imue) = eos_electrons(i,j,k,electron_imue) 

        enddo
     enddo
  enddo
  write(*,*) "Done that"


end subroutine electron_photon_eos

subroutine wrap_timmes(abar,zbar,rho,temp,pele,ppos,&
     eele,epos,sele,spos,prad,erad,srad,mue,muep,gamma,dedT,dpdT,dpdr,&
     xne,dxnedd,dxnedt,dxnedz,dxneda)

  implicit none
  include 'vector_eos.dekF90'

  integer ionmax
  parameter (ionmax=1)
  real*8 xmass(ionmax),ymass(ionmax),&
       aion(ionmax),zion(ionmax),temp,rho
  real*8 abar,zbar,pele,ppos,eele,epos,sele,spos,prad,erad,srad
  real*8 gamma,dedT,dpdT,dpdr,mue,muep
  real*8 xne,dxnedd,dxnedt,dxnedz,dxneda
  xmass(1) = 1.0d0

  aion(1) = abar
  zion(1) = zbar

  jlo_eos = 1
  jhi_eos = 1

  abar_row(1) = abar
  zbar_row(1) = zbar
  
  den_row(1) = rho
  temp_row(1) = temp

  !this calls a modified timmes eos where the ions are left out, we'll
  !add them later
  call eosfxt

  pele = pele_row(1)
  ppos = ppos_row(1)

  eele = eele_row(1)
  epos = epos_row(1)

  sele = sele_row(1)
  spos = spos_row(1)

  prad = prad_row(1)
  erad = erad_row(1)
  srad = srad_row(1)

  mue = etaele_row(1)
  muep = etapos_row(1)

  gamma = gam1_row(1)
  dedT = det_row(1)
  dpdT = dpt_row(1)
  dpdr = dpd_row(1)

  xne = xne_row(1)
  dxnedd = dxned_row(1)
  dxnedt = dxnet_row(1)
  dxnedz = dxnez_row(1)
  dxneda = dxnea_row(1)


end subroutine wrap_timmes
