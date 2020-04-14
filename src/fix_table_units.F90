!-*-f90-*-
subroutine fix_table_units

! * go to cgs for energy and pressure
! * fix energy normalization

  use eos_table_module
  implicit none

  integer i,j,k

  do k=1,nye
     do j=1,nlogtemp
        do i=1,nlogrho
! 1 Mev/fm = mev_to_erg * 1d13 dyn
! 1 fm*fm = 1d-26 cm^2
           eos_table(i,j,k,ipress) = &
                eos_table(i,j,k,ipress) * mev_to_erg*1.0d13/1.0d-26

! energy shift -- Shen et al. use the atomic mass unit as the
! normalization of their energy and density. We keep this to ensure
! proper adiabatic collapse.  We do not take into account the mass
! difference of the neutron and proton in the internal energy to
! remain consitant with the choice of 938MeV in the Shen EOS. The
! following is true for the H Shen. 2010 EOS: 
! 1. @ yp=0,T=0,n_B=0, the internal energy of the matter in the Shen
! EOS is 6.5MeV = 938MeV-931.494MeV.  
! 2. @ yp=0.5,T=0,n_B=0, the internal energy of the matter in the Shen
! EOS is -2.31MeV and the average A=53,Z=27.  This is roughly the
! excess binding energy, per nucleon of iron compared to carbon 12
! (i.e. the reference energy of c12, i.e. 931.494MeV)

           eos_table(i,j,k,ienergy) = eos_table(i,j,k,ienergy)

! energy: Mev/baryon --> erg/g
! Mev/baryon * n_b = MeV/cm^3
! Mev/cm^3 * 1.6022d-6 --> erg/cm^3
! erg/cm^3 / rho --> erg/g

           eos_table(i,j,k,ienergy) = &
                eos_table(i,j,k,ienergy) * &
                1.0d0/baryon_to_gram_conversion * mev_to_erg

! chemical potential: Shen et al. (both G. and H.) do not distinguish
! between nucleon masses for neutrons and protons, their values are
! relative to a nucleon mass of 938.0 MeV (939.0 for G. Shen), but we
! need to include the rest mass difference in muhat (to be consistent
! with the LS treatment).  We accomplish that by adding in the proton
! rest mass m_p into mu_p and the neutron rest mass m_n into mu_n,
! then subtract m_ref from both.  mu_n - mu_p has then a constant
! offset of ~1.294 MeV.

! Hemple does include the mass difference, but subtracts it out for
! the tables used in EOSmaker.  Therefore we have to add the
! difference back in

           eos_table(i,j,k,imup) = &
                eos_table(i,j,k,imup) + (m_p-m_ref)

           eos_table(i,j,k,imun) = &
                eos_table(i,j,k,imun) + (m_n-m_ref)

        enddo
     enddo
  enddo

end subroutine fix_table_units
