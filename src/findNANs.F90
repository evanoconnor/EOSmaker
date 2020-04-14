!-*-f90-*-
subroutine checkdata3d(array,nx,ny,nz,aname)

  implicit none
  integer nx,ny,nz
  integer i,j,k
   integer zerocounter
  integer nancounter
  integer infcounter
  integer totalcounter

  real*8 array(nx,ny,nz)
  character(*) aname
  logical myisnan
  
  zerocounter=0
  nancounter=0
  infcounter=0
  totalcounter=0

  do i=1,nx
     do j=1,ny
        do k=1,nz
           if(array(i,j,k).ne.array(i,j,k)) then
              nancounter=nancounter+1
              write(*,*) "NaN in ",aname
              write(*,"(i4,i4,i4,1P10E16.7)") i,j,k,array(i,j,k)
              stop
           endif
           if(abs(array(i,j,k)).ge.1.d99) then
              infcounter=infcounter+1
              write(*,"(i4,i4,i4,1P10E16.7)") i,j,k,array(i,j,k)
           endif
           if(array(i,j,k).eq.0.0d0) then
              zerocounter=zerocounter+1
           endif
        enddo
     enddo
  enddo

  write(*,*) zerocounter, "0s in ", aname," Total: ",totalcounter
  write(*,*) nancounter, "NANs in ", aname," Total: ",totalcounter
  write(*,*) infcounter, "INFs in ", aname," Total: ",totalcounter


end subroutine checkdata3d

