subroutine write_Burrows_table

  use eos_table_module
  implicit none

  ! Adam's units:
  ! pressure: non-log in Mev/fm (?)

  integer :: nt = 300
  integer :: nr = 350
  integer :: ny = 60
  integer :: numel = 16
  real*8,allocatable ::  adam_table(:,:,:,:)

  ! local stuff
  integer i,j,k,irecl,irec

  allocate(adam_table(numel,ny,nr,nt))

  if(nlogrho.ne.nr .or. &
     nlogtemp.ne.nt .or. &
     nye.ne.ny) then
     write(*,*) "Dimensions of Ott and Burrows tables do not match!"
     stop "This is bad!"
  endif

  do i=1,nlogrho
     do j=1,nlogtemp
        do k=1,nye

           ! energy per baryon with Adam's shift
           adam_table(1,k,i,j) = 10.d0**eos_table(i,j,k,ienergy) &
                * amu_cgs / mev_to_erg - energy_shift*amu_cgs/mev_to_erg &
                + 9.3d0

           if (adam_table(1,k,i,j)<0.0d0) write(*,*) k,i,j,adam_table(1,k,i,j)

#if 0
           if(i.eq.80.and.j.eq.55.and.k.eq.36) then
              write(6,*) 10.0d0**eos_table(i,j,k,ienergy) &
                   - energy_shift, &
                   10.0d0**eos_table(i,j,k,ienergy),&
                   energy_shift*amu_cgs/mev_to_erg
              stop "blah"
           endif
#endif
           ! pressure
           adam_table(2,k,i,j) = (10.0d0**eos_table(i,j,k,ipress)) / &
                1.60217733d33

           ! entropy per baryon
           adam_table(3,k,i,j) = eos_table(i,j,k,ientropy)

           ! C_v = dedt
           adam_table(4,k,i,j) = dedt(i,j,k) * &
                1.60217733d-6/1.381d-16

           ! xn
           adam_table(5,k,i,j) = eos_table(i,j,k,ixn)
           ! xp
           adam_table(6,k,i,j) = eos_table(i,j,k,ixp)
           ! xa
           adam_table(7,k,i,j) = eos_table(i,j,k,ixa)
           ! xh
           adam_table(8,k,i,j) = eos_table(i,j,k,ixh)

           ! za
           adam_table(9,k,i,j) = eos_table(i,j,k,izbar)/eos_table(i,j,k,iabar)

           ! aw
           adam_table(10,k,i,j) = eos_table(i,j,k,iabar)

           ! muhat
           adam_table(11,k,i,j) = eos_table(i,j,k,imun) - &
                eos_table(i,j,k,imup)

           ! gammac
           adam_table(12,k,i,j) = max(1.0,min(3.5,gamma(i,j,k)))

           ! dhy
           adam_table(13,k,i,j) = -1.0d55

           ! zht
           adam_table(14,k,i,j) = -1.0d55

           ! mue
           adam_table(15,k,i,j) = eos_table(i,j,k,imue)

           ! dpde
           adam_table(16,k,i,j) = dpde(i,j,k)

        enddo
     enddo
  enddo

  write(*,*) minval(adam_table(1,:,:,:))

  irecl = 8*ny*nr*nt
  open(10,file='SFHo_table_30035060_stellarcollapse_gammalimited_b.dat',&
       form='unformatted',access='direct',recl=irecl)
  do irec=1,numel
     write(10,rec=irec) (((adam_table(irec,k,i,j),&
          j=1,nt),i=1,nr),k=1,ny)
  enddo
  close(10)
  write(6,*) "Wrote EOS table in Burrows format!"

end subroutine write_Burrows_table

!The way Adam likes it:
!--------------------------------------------
!  table(1,y,r,t)  :: energy per baryon      |
!  table(2,y,r,t)  :: pressure               |
!  table(3,y,r,t)  :: entropy per baryon     |
!  table(4,y,r,t)  :: cv                     |
!  table(5,y,r,t)  :: xn                     |
!  table(6,y,r,t)  :: xp                     |
!  table(7,y,r,t)  :: xa                     |
!  table(8,y,r,t)  :: xh                     |
!  table(9,y,r,t)  :: za                     | 
!  table(10,y,r,t) :: aw                     |
!  table(11,y,r,t) :: muhat                  |
!  table(12,y,r,t) :: gammac                 |
!  table(13,y,r,t) :: dhy                    |
!  table(14,y,r,t) :: zht                    |
!  table(15,y,r,t) :: mue                    |
!  table(16,y,r,t) :: dpde                   |
!--------------------------------------------
