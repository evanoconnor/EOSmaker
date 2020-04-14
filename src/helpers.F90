!-*-f90-*-
subroutine extrapolate_lin(x,f,xx,ff)

  implicit none
  

  real*8 x(2), f(2), xx, ff

  ff = (f(2)-f(1))/(x(2)-x(1)) * (xx-x(1)) + f(1)

end subroutine extrapolate_lin

subroutine linterp_1D(n,data,x,xx,ff)
  
  implicit none
  integer :: n
  real*8  :: data(n)
  real*8  :: x(n)
  real*8  :: xx,ff

  ! functions
  real*8  :: psi0,psi1

  ! local vars
  integer :: nvars = 1
  integer :: i
  real*8  :: deriv(n)
  integer :: if(2)

  real*8  :: m,b
  
  ! computer derivatives
  deriv = 0.0d0
  if(n.le.3) stop "linterp_1D: need n > 3"
!  if(xx.lt.x(1)) stop "cubic_hermit_1D: need xx >= x(1)"
  if(xx.lt.x(1)) then
     write(6,"(1P10E18.9)") xx,x(1)
     stop "linterp_1D: need xx >= x(1)"
  endif

  ! find bracketing points
  call two_friends(xx,x,n,if)

  ! set up interpolation
  b = data(if(1))
  m = (data(if(2)) - data(if(1)))/(x(if(2)) - x(if(1)))

  ff = m * (xx-x(if(1))) + b


end subroutine linterp_1D

subroutine cubic_hermite_1D(n,data,x,xx,ff)
  
  implicit none
  integer :: n
  real*8  :: data(n)
  real*8  :: x(n)
  real*8  :: xx,ff

  ! functions
  real*8  :: psi0,psi1

  ! local vars
  integer :: nvars = 1
  integer :: i
  real*8  :: deriv(n)
  integer :: if(2)

  real*8  :: z,y,dx
  
  ! computer derivatives
  deriv = 0.0d0
  if(n.le.3) stop "cubic_hermite_1D: need n > 3"
!  if(xx.lt.x(1)) stop "cubic_hermit_1D: need xx >= x(1)"
  if(xx.lt.x(1)) then
     write(6,"(1P10E18.9)") xx,x(1)
     stop "cubic_hermit_1D: need xx >= x(1)"
  endif

  if(xx.ge.x(n)) stop "cubic_hermit_1D: need xx < x(n)"
  do i=2,n-1
     call deriv1ne(x(i-1),       &
          x(i),                  &
          x(i+1),                &
          data(i-1),             &
          data(i),               &
          data(i+1),             &
          deriv(i),              &
          nvars)
  enddo
  ! boundary: one-sided deriv
  deriv(1) = (data(2)-data(1))/(x(2)-x(1))
  deriv(n) = (data(n)-data(n-1))/(x(n)-x(n-1))

  ! find bracketing points
  call two_friends(xx,x,n,if)

  ! set up interpolation
  z = (xx-x(if(1)))/(x(if(2))-x(if(1)))
  y = 1.0d0 - z
  dx = x(if(2))-x(if(1))

  ff = data(if(1))*psi0(z)     &
     + data(if(2))*psi0(y)     &
     + deriv(if(1))*dx*psi1(z) &
     - deriv(if(2))*dx*psi1(y)

end subroutine cubic_hermite_1D


subroutine deriv1ne(x1,x2,x3,f1,f2,f3,fp,nvars)

  implicit none
  integer nvars, ivar
  real*8 x1,x2,x3,f1(nvars),f2(nvars),f3(nvars),fp(nvars)

  real*8 h1 ! h_i-1
  real*8 h2 ! h_i

  real*8 s1(nvars) ! s_i-1
  real*8 s2(nvars) ! s_i

  h1 = x2-x1
  h2 = x3-x2

  s1 = (f2-f1)/(x2-x1)
  s2 = (f3-f2)/(x3-x2)

  fp = (s1*h2 + s2*h1)/(h1+h2)
  
  ! limiter:
  do ivar = 1,nvars 
     if (s1(ivar)*s2(ivar) .le. 0.0d0) then
        fp(ivar) = 0.0d0
     else if ( (abs(fp(ivar)).gt.2.0d0*abs(s1(ivar))) .or. &
          (abs(fp(ivar)).gt.2.0d0*abs(s2(ivar))) ) then
        fp(ivar) = 2.0d0*sign(1.0d0,s1(ivar))*min(abs(s1(ivar)),abs(s2(ivar)))
     endif
  enddo


end subroutine deriv1ne

function psi0(z)

  implicit none
  real*8 psi0
  real*8 z

  psi0 = 2.0d0*z**3 - 3.0d0*z**2 + 1.0d0

end function psi0


function psi1(z)

  implicit none
  real*8 psi1
  real*8 z

  psi1 = z**3 - 2.0d0*z**2 + z

end function psi1



subroutine two_friends(val,tval,n,friends)

  implicit none
  integer n
  real*8 val,tval(n)
  integer friends(2)

  integer i

  !special case
  if (val.eq.tval(1)) then
     friends(1) = 1
     friends(2) = 2
     return
  else if (val.eq.tval(n)) then
     friends(1) = n-1
     friends(2) = n
     return
  endif

  i=1
  do while(tval(i).le.val.and.i.le.n)
     i=i+1
  enddo

  if(i.eq.1) then
     write(6,"(1P10E15.6)") val, tval(1) 
     stop "BAD! Out of table range!"
  endif

  if(i.gt.n) then
     write(*,*) "i .eq. n+1 in two_friends",val,tval(n)
     stop
  endif
  
  friends(1) = i-1
  friends(2) = i


end subroutine two_friends

subroutine findindex(n,array,value,iminus,iplus)
  implicit none
  
  integer n,iminus,iplus
  double precision array(n)
  double precision value
  double precision buffer
        
  buffer = (value-array(1))/(array(n)-array(1)) * (n-1)*1.0d0 
  iminus = int(buffer)+1
  iplus  = int(buffer)+2

  
end subroutine findindex
      
