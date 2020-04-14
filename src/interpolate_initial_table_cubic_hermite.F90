!-*-f90-*-
subroutine interpolate_initial_table_cubic_hermite

  use eos_table_module
  implicit none

  integer i,j,k,ii,jj,kk !counters
  real*8 min_energy, min_press !min energy and pressure variables
  
  !will store the d/drho, d/dtemp, and d/drhodtemp for interpolation
  real*8 ddrho(nt_temp,nt_yp,nt_rho,nvars)
  real*8 ddtemp(nt_temp,nt_yp,nt_rho,nvars)
  real*8 ddtempdrho(nt_temp,nt_yp,nt_rho,nvars)

  !interpolating stuff
  real*8, external :: psi0,psi1
  real*8 r2,r1,r 
  integer ir1,ir2

  real*8 :: d1(nvars), d2(nvars), z

  real*8 t,t1,t2,x,y,dt,dr
  real*8 w0x,w1x,w0mx,w1mx
  real*8 w0y,w1y,w0my,w1my
  real*8, allocatable :: temp_data(:,:,:,:)
  real*8 tdata11(nvars),ddtemp11(nvars),ddrho11(nvars),ddtempdrho11(nvars)
  real*8 tdata12(nvars),ddtemp12(nvars),ddrho12(nvars),ddtempdrho12(nvars)
  real*8 tdata21(nvars),ddtemp21(nvars),ddrho21(nvars),ddtempdrho21(nvars)
  real*8 tdata22(nvars),ddtemp22(nvars),ddrho22(nvars),ddtempdrho22(nvars)
  integer tf(2),rf(2)

  real*8 y1,y2,y3,my
  integer yf(2)

  real*8 y_lin,rho_lin,temp_lin
  real*8 alpha_ye,alpha_rho,alpha_temp
  real*8 c00,c01,c10,c11,c0,c1,c

  real*8 lr,lt,diff,res,d1a,d2a,d3a

  integer count_diff(16)
  integer irho,ik
  real*8 xx(2),yy(2)

  logical :: nopresswarning = .true.


  write(*,*) "Starting to interpolate..."
  ! make pressure and energy logarithmic
  
  ! Due to coulomb corrections the pressure and due to binding
  ! energies, the energy is negative at some places; hence can't do
  ! log without shifting things; find minima here:

  min_energy = minval(table_data(:,:,:,ienergy))
  min_press = minval(table_data(:,:,:,ipress))
  
  !shift table energy up by a constant, don't worry, we take it off
  !after interpolation
  if(min_energy.lt.0.0d0) then 
     table_data(:,:,:,ienergy) = table_data(:,:,:,ienergy) - min_energy*1.1d0
  endif
  
  !shift table pressure up by a constant, don't worry, we take it off
  !after interpolation
  if(min_press.lt.0.0d0) then 
     table_data(:,:,:,ipress) = table_data(:,:,:,ipress) - min_press*1.1d0
  endif

  !now log it
  table_data(:,:,:,ienergy) = log10(table_data(:,:,:,ienergy))
  table_data(:,:,:,ipress) = log10(table_data(:,:,:,ipress))

  !First, we interpolate in the rho-T plane, for this we need to
  ! compute derivatives with respect to rho, T and rho,T
  ddrho = 0.0d0
  ddtemp = 0.0d0
  ddtempdrho = 0.0d0
  
  ! 1st derivative in rho direction
  do k=1,nt_yp
     do j=1,nt_temp
        do i=2,nt_rho-1
           call deriv1ne(table_logrho(i-1),       &
                table_logrho(i),                  &
                table_logrho(i+1),                &
                table_data(j,k,i-1,:),            &
                table_data(j,k,i,:),              &
                table_data(j,k,i+1,:),            &
                ddrho(j,k,i,:),                   &
                nvars)
        enddo

        ! Boundary values for first derivative
        ddrho(j,k,1,:) = (table_data(j,k,2,:) - table_data(j,k,1,:)) &
             / (table_logrho(2) - table_logrho(1))

        ddrho(j,k,nt_rho,:) = (table_data(j,k,nt_rho,:) - table_data(j,k,nt_rho-1,:)) &
             / (table_logrho(nt_rho) - table_logrho(nt_rho-1))

     enddo
  enddo
  write(*,*) "Done 1st derivative in rho"

  ! 1st derivative in T direction
  do k=1,nt_yp
     do i=1,nt_rho
        do j=2,nt_temp-1
           call deriv1ne(table_logtemp(j-1),       &
                table_logtemp(j),                  &
                table_logtemp(j+1),                &
                table_data(j-1,k,i,:),             &
                table_data(j,k,i,:),               &
                table_data(j+1,k,i,:),             &
                ddtemp(j,k,i,:),                   &
                nvars)

        enddo
        ! Boundary values for first derivative
        ddtemp(1,k,i,:) = (table_data(2,k,i,:) - table_data(1,k,i,:)) &
             / (table_logtemp(2) - table_logtemp(1))
        
        ddtemp(nt_temp,k,i,:) = (table_data(nt_temp,k,i,:) - table_data(nt_temp-1,k,i,:)) &
             / (table_logtemp(nt_temp) - table_logtemp(nt_temp-1))
        
     enddo
  enddo
  write(*,*) "Done 1st derivative in temp"
  
  ! second derivative in T, rho direction
  ! standard centered differences
  do k=1,nt_yp
     do i=1,nt_rho
        do j=2,nt_temp-1
           call deriv1ne(table_logtemp(j-1),        &
                table_logtemp(j),                  &
                table_logtemp(j+1),                &
                ddrho(j-1,k,i,:),                  &
                ddrho(j,k,i,:),                    &
                ddrho(j+1,k,i,:),                  &
                ddtempdrho(j,k,i,:),               &
                nvars)
        enddo
        ! Boundary values
        ddtempdrho(1,k,i,:) =  ( ddrho(2,k,i,:) - ddrho(1,k,i,:))  &
             / (table_logtemp(2) - table_logtemp(1))
        
        ddtempdrho(nt_temp,k,i,:) =  ( ddrho(nt_temp,k,i,:) - ddrho(nt_temp-1,k,i,:))  &
             / (table_logtemp(nt_temp) - table_logtemp(nt_temp-1))
        
     enddo
  enddo
  write(*,*) "Done 2nd derivative in rhotemp"

  allocate(temp_data(nlogrho,nlogtemp,nt_yp,nvars))
  temp_data = 0.0d0

  do k=1,nt_yp
     do jj=itemp_start_nuclear,itemp_stop_nuclear
        do ii=irho_start_comps,nlogrho

           r = logrho(ii)
           t = logtemp(jj)
           ! find bracketing points, two_friends take a value (r) and
           ! finds the indices (tf(1) and tf(2)) in the supplied array
           ! (table_logtemp) of length (nt_temp)
           call two_friends(t,table_logtemp,nt_temp,tf)
           call two_friends(r,table_logrho,nt_rho,rf)

           !get those value of temp and rho of the initial table
           t1 = table_logtemp(tf(1))
           t2 = table_logtemp(tf(2))
           r1 = table_logrho(rf(1))
           r2 = table_logrho(rf(2))

           !get derivatives we calculated above
           ddrho11 = ddrho(tf(1),k,rf(1),:)
           ddrho12 = ddrho(tf(1),k,rf(2),:)
           ddrho21 = ddrho(tf(2),k,rf(1),:)
           ddrho22 = ddrho(tf(2),k,rf(2),:)

           ddtemp11 = ddtemp(tf(1),k,rf(1),:)
           ddtemp12 = ddtemp(tf(1),k,rf(2),:)
           ddtemp21 = ddtemp(tf(2),k,rf(1),:)
           ddtemp22 = ddtemp(tf(2),k,rf(2),:)

           ddtempdrho11 = ddtempdrho(tf(1),k,rf(1),:)
           ddtempdrho12 = ddtempdrho(tf(1),k,rf(2),:)
           ddtempdrho21 = ddtempdrho(tf(2),k,rf(1),:)
           ddtempdrho22 = ddtempdrho(tf(2),k,rf(2),:)

           !and get the table values themselves
           tdata11 = table_data(tf(1),k,rf(1),:)
           tdata12 = table_data(tf(1),k,rf(2),:)
           tdata21 = table_data(tf(2),k,rf(1),:)
           tdata22 = table_data(tf(2),k,rf(2),:)
           
           !do the interpolation
           x = (r - r1)/(r2-r1)
           y = (t - t1)/(t2-t1)
           
           dt = t2-t1
           dr = r2-r1

           w0x = psi0(x)
           w0mx = psi0(1.0d0-x)

           w1x = psi1(x)
           w1mx = psi1(1.0d0-x)

           w0y = psi0(y)
           w0my = psi0(1.0d0-y)

           w1y = psi1(y)
           w1my = psi1(1.0d0-y)


           temp_data(ii,jj,k,:) = &
!                
                tdata11(:)              *w0x          *w0y                 &
              + tdata12(:)              *w0mx         *w0y                 &
              + tdata21(:)              *w0x          *w0my                &
              + tdata22(:)              *w0mx         *w0my                &
!                                     
              + ddtemp11(:)*dt          *w0x          *w1y                 &
              + ddtemp12(:)*dt          *w0mx         *w1y                 &
              - ddtemp21(:)*dt          *w0x          *w1my                &
              - ddtemp22(:)*dt          *w0mx         *w1my                &
!                                     
              + ddrho11(:)*dr           *w1x          *w0y                 &
              - ddrho12(:)*dr           *w1mx         *w0y                 &
              + ddrho21(:)*dr           *w1x          *w0my                &
              - ddrho22(:)*dr           *w1mx         *w0my                &
!
              + ddtempdrho11(:)*dr*dt   *w1x          *w1y                 &
              - ddtempdrho12(:)*dr*dt   *w1mx         *w1y                 &
              - ddtempdrho21(:)*dr*dt   *w1x          *w1my                &
              + ddtempdrho22(:)*dr*dt   *w1mx         *w1my               

           ! sanity check:
           if(temp_data(ii,jj,k,ipress).ne.temp_data(ii,jj,k,ipress)) then
              stop "Failed NAN/Inf sanity check in interpolate"
           endif
        enddo
     enddo
  enddo

  write(*,*) "Done interpolating calculations in rho-T space..."

  ! now let's interpolate in Y_e :-)
  ! Again, let's use a cubic Hermite method

  do kk=1,nye
     do jj=itemp_start_nuclear,itemp_stop_nuclear
        do ii=irho_start_comps,nlogrho
           
           y = ye(kk)
           call two_friends(y,table_yp,nt_yp,yf)
           y = ye(kk)
           y1 = table_yp(yf(1))
           y2 = table_yp(yf(2))

           if(yf(2).eq.nt_yp) then
              d2(:) = (temp_data(ii,jj,yf(2),:) - temp_data(ii,jj,yf(2)-1,:)) &
                   / (table_yp(yf(2)) - table_yp(yf(2)-1))
           else
              call deriv1ne(table_yp(yf(2)-1),    &
                   table_yp(yf(2)),               &
                   table_yp(yf(2)+1),             &
                   temp_data(ii,jj,yf(2)-1,:),    &
                   temp_data(ii,jj,yf(2),:),      &
                   temp_data(ii,jj,yf(2)+1,:),    &
                   d2,                            &
                   nvars)
           endif
           if(yf(1).eq.1) then
              d1(:) = (temp_data(ii,jj,yf(1)+1,:) - temp_data(ii,jj,yf(1),:)) &
                   / (table_yp(yf(1)+1) - table_yp(yf(1)))
           else
              call deriv1ne(table_yp(yf(1)-1),    &
                   table_yp(yf(1)),            &
                   table_yp(yf(1)+1),          &
                   temp_data(ii,jj,yf(1)-1,:),    &
                   temp_data(ii,jj,yf(1),:),      &
                   temp_data(ii,jj,yf(1)+1,:),    &
                   d1,                            &
                   nvars)
           endif


           z = (y - y1) / (y2 - y1)
           x = 1.0d0-z

           !last step in interpolations
           eos_table(ii,jj,kk,1:nvars) = temp_data(ii,jj,yf(1),1:nvars)*psi0(z) &
                + temp_data(ii,jj,yf(2),1:nvars)*psi0(x) &
                + d1(1:nvars)*(y2-y1)*psi1(z) &
                - d2(1:nvars)*(y2-y1)*psi1(x)

        enddo
     enddo
  enddo

  write(*,*) "Done interpolating to final table.."


  if (Hempel_massfrac_reconstruction) then
     write(*,*) "Doing Hemple reconstruction"

     table_linrho = 10.0d0**table_logrho
     table_lintemp = 10.0d0**table_logtemp

     write(*,*) itemp_start_nuclear,itemp_stop_nuclear
     write(*,*) irho_start_comps,nlogrho

     do kk=1,nye
        do jj=itemp_start_nuclear,itemp_stop_nuclear
           do ii=irho_start_comps,nlogrho     

              y_lin = ye(kk)
              call two_friends(y_lin,table_yp,nt_yp,yf)
              y_lin = ye(kk)
              y1 = table_yp(yf(1))
              y2 = table_yp(yf(2))
              alpha_ye = (y_lin-y1)/(y2-y1)

              rho_lin = 10.0d0**logrho(ii)
              call two_friends(rho_lin,table_linrho,nt_rho,rf)
              rho_lin = 10.0d0**logrho(ii)
              r1 = table_linrho(rf(1))
              r2 = table_linrho(rf(2))
              alpha_rho = (rho_lin-r1)/(r2-r1)

              temp_lin = 10.0d0**logtemp(jj)
              call two_friends(temp_lin,table_lintemp,nt_temp,tf)
              temp_lin = 10.0d0**logtemp(jj)
              t1 = table_lintemp(tf(1))
              t2 = table_lintemp(tf(2))
              alpha_temp = (temp_lin-t1)/(t2-t1)

              c00 = table_data(tf(1),yf(1),rf(1),ixn)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixn)*alpha_ye
              c10 = table_data(tf(1),yf(1),rf(2),ixn)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixn)*alpha_ye
              c01 = table_data(tf(2),yf(1),rf(1),ixn)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixn)*alpha_ye
              c11 = table_data(tf(2),yf(1),rf(2),ixn)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixn)*alpha_ye

              c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
              c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

              c = c0*(1.0d0-alpha_temp) + c1*alpha_temp

              eos_table(ii,jj,kk,ixn) = c


              c00 = table_data(tf(1),yf(1),rf(1),ixp)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixp)*alpha_ye
              c10 = table_data(tf(1),yf(1),rf(2),ixp)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixp)*alpha_ye
              c01 = table_data(tf(2),yf(1),rf(1),ixp)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixp)*alpha_ye
              c11 = table_data(tf(2),yf(1),rf(2),ixp)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixp)*alpha_ye

              c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
              c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

              c = c0*(1.0d0-alpha_temp) + c1*alpha_temp

              eos_table(ii,jj,kk,ixp) = c

              c00 = table_data(tf(1),yf(1),rf(1),ixa)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixa)*alpha_ye
              c10 = table_data(tf(1),yf(1),rf(2),ixa)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixa)*alpha_ye
              c01 = table_data(tf(2),yf(1),rf(1),ixa)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixa)*alpha_ye
              c11 = table_data(tf(2),yf(1),rf(2),ixa)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixa)*alpha_ye

              c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
              c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

              c = c0*(1.0d0-alpha_temp) + c1*alpha_temp

              eos_table(ii,jj,kk,ixa) = c

              c00 = table_data(tf(1),yf(1),rf(1),ixh)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixh)*alpha_ye
              c10 = table_data(tf(1),yf(1),rf(2),ixh)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixh)*alpha_ye
              c01 = table_data(tf(2),yf(1),rf(1),ixh)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixh)*alpha_ye
              c11 = table_data(tf(2),yf(1),rf(2),ixh)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixh)*alpha_ye

              c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
              c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

              c = c0*(1.0d0-alpha_temp) + c1*alpha_temp

              eos_table(ii,jj,kk,ixh) = c

              if (Hempel) then
                 c00 = table_data(tf(1),yf(1),rf(1),ixd)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixd)*alpha_ye
                 c10 = table_data(tf(1),yf(1),rf(2),ixd)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixd)*alpha_ye
                 c01 = table_data(tf(2),yf(1),rf(1),ixd)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixd)*alpha_ye
                 c11 = table_data(tf(2),yf(1),rf(2),ixd)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixd)*alpha_ye

                 c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
                 c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

                 c = c0*(1.0d0-alpha_temp) + c1*alpha_temp
                 
                 eos_table(ii,jj,kk,ixd) = c

                 c00 = table_data(tf(1),yf(1),rf(1),ixt)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixt)*alpha_ye
                 c10 = table_data(tf(1),yf(1),rf(2),ixt)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixt)*alpha_ye
                 c01 = table_data(tf(2),yf(1),rf(1),ixt)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixt)*alpha_ye
                 c11 = table_data(tf(2),yf(1),rf(2),ixt)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixt)*alpha_ye
                 
                 c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
                 c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

                 c = c0*(1.0d0-alpha_temp) + c1*alpha_temp
                 
                 eos_table(ii,jj,kk,ixt) = c
                 
                 c00 = table_data(tf(1),yf(1),rf(1),ix3he)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ix3he)*alpha_ye
                 c10 = table_data(tf(1),yf(1),rf(2),ix3he)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ix3he)*alpha_ye
                 c01 = table_data(tf(2),yf(1),rf(1),ix3he)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ix3he)*alpha_ye
                 c11 = table_data(tf(2),yf(1),rf(2),ix3he)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ix3he)*alpha_ye
                 
                 c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
                 c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho
                 
                 c = c0*(1.0d0-alpha_temp) + c1*alpha_temp

                 eos_table(ii,jj,kk,ix3he) = c
                 
                 c00 = table_data(tf(1),yf(1),rf(1),ix4li)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ix4li)*alpha_ye
                 c10 = table_data(tf(1),yf(1),rf(2),ix4li)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ix4li)*alpha_ye
                 c01 = table_data(tf(2),yf(1),rf(1),ix4li)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ix4li)*alpha_ye
                 c11 = table_data(tf(2),yf(1),rf(2),ix4li)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ix4li)*alpha_ye
                 
                 c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
                 c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho
                 
                 c = c0*(1.0d0-alpha_temp) + c1*alpha_temp
                 
                 eos_table(ii,jj,kk,ix4li) = c
                 
              endif

              if (Hyperon) then

                 c00 = table_data(tf(1),yf(1),rf(1),ixL)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(1),ixL)*alpha_ye
                 c10 = table_data(tf(1),yf(1),rf(2),ixL)*(1.0d0-alpha_ye) + table_data(tf(1),yf(2),rf(2),ixL)*alpha_ye
                 c01 = table_data(tf(2),yf(1),rf(1),ixL)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(1),ixL)*alpha_ye
                 c11 = table_data(tf(2),yf(1),rf(2),ixL)*(1.0d0-alpha_ye) + table_data(tf(2),yf(2),rf(2),ixL)*alpha_ye
                 
                 c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
                 c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho
                 
                 c = c0*(1.0d0-alpha_temp) + c1*alpha_temp
              
                 eos_table(ii,jj,kk,ixL) = c
              endif

              c00 = table_data(tf(1),yf(1),rf(1),ixh)/table_data(tf(1),yf(1),rf(1),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(1),yf(2),rf(1),ixh)/table_data(tf(1),yf(2),rf(1),iabar)*alpha_ye
              c10 = table_data(tf(1),yf(1),rf(2),ixh)/table_data(tf(1),yf(1),rf(2),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(1),yf(2),rf(2),ixh)/table_data(tf(1),yf(2),rf(2),iabar)*alpha_ye
              c01 = table_data(tf(2),yf(1),rf(1),ixh)/table_data(tf(2),yf(1),rf(1),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(2),yf(2),rf(1),ixh)/table_data(tf(2),yf(2),rf(1),iabar)*alpha_ye
              c11 = table_data(tf(2),yf(1),rf(2),ixh)/table_data(tf(2),yf(1),rf(2),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(2),yf(2),rf(2),ixh)/table_data(tf(2),yf(2),rf(2),iabar)*alpha_ye
              
              c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
              c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho
              
              c = c0*(1.0d0-alpha_temp) + c1*alpha_temp
              
              eos_table(ii,jj,kk,iabar) = eos_table(ii,jj,kk,ixh)/c
                 
              c00 = (table_data(tf(1),yf(1),rf(1),ixh)*table_data(tf(1),yf(1),rf(1),izbar)/ &
                   table_data(tf(1),yf(1),rf(1),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(1),yf(2),rf(1),ixh)*table_data(tf(1),yf(2),rf(1),izbar)/ &
                   table_data(tf(1),yf(2),rf(1),iabar)*alpha_ye)
              c10 = (table_data(tf(1),yf(1),rf(2),ixh)*table_data(tf(1),yf(1),rf(2),izbar)/ &
                   table_data(tf(1),yf(1),rf(2),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(1),yf(2),rf(2),ixh)*table_data(tf(1),yf(2),rf(2),izbar)/ &
                   table_data(tf(1),yf(2),rf(2),iabar)*alpha_ye)
              c01 = (table_data(tf(2),yf(1),rf(1),ixh)*table_data(tf(2),yf(1),rf(1),izbar)/ &
                   table_data(tf(2),yf(1),rf(1),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(2),yf(2),rf(1),ixh)*table_data(tf(2),yf(2),rf(1),izbar)/ &
                   table_data(tf(2),yf(2),rf(1),iabar)*alpha_ye)
              c11 = (table_data(tf(2),yf(1),rf(2),ixh)*table_data(tf(2),yf(1),rf(2),izbar)/ &
                   table_data(tf(2),yf(1),rf(2),iabar)*(1.0d0-alpha_ye) + &
                   table_data(tf(2),yf(2),rf(2),ixh)*table_data(tf(2),yf(2),rf(2),izbar)/ &
                   table_data(tf(2),yf(2),rf(2),iabar)*alpha_ye)


              c0 = c00*(1.0d0-alpha_rho) + c10*alpha_rho
              c1 = c01*(1.0d0-alpha_rho) + c11*alpha_rho

              c = c0*(1.0d0-alpha_temp) + c1*alpha_temp

              eos_table(ii,jj,kk,izbar) = eos_table(ii,jj,kk,iabar)/eos_table(ii,jj,kk,ixh)*c

              if (eos_table(ii,jj,kk,ixh).eq.0.0d0) then
                 eos_table(ii,jj,kk,iabar) = 1.0d0
                 eos_table(ii,jj,kk,izbar) = 1.0d0
              endif


           enddo
        enddo
     enddo
  endif

  count_diff = 0
  ! quick & dirty sanity check, compare trilinear interpolation to
  ! what we just did.  If numbers are crazy then perhaps something
  ! wrong...

  do k=2,nt_yp-1
     do i=2,nt_rho-1
        do j=2,nt_temp-4

           lr = table_logrho(i)
           lt = table_logtemp(j)
           y =  table_yp(k)

           if (y.lt.yemin.or.y.gt.yemax) then

           else if (lr.gt.logrhomax) then

           else
              ! if there is error here, check if the routine is not trying to test the table outside its range
              call trilin_interp_table(lr,lt,y,&
                   res,eos_table(:,:,:,ipress),&
                   nlogrho,nlogtemp,nye,&
                   logrho,logtemp,ye,d1a,d2a,d3a)
              
              diff = abs(res - table_data(j,k,i,ipress)) &
                   / table_data(j,k,i,ipress)
              
              if(diff.lt.1.0d-2) then
                 count_diff(1) = count_diff(1)+1
              else if(diff.lt.2.0d-2) then
                 count_diff(2) = count_diff(2)+1
              else if (diff.lt.3.0d-2) then
                 count_diff(3) = count_diff(3)+1
              else if(diff.lt.4.0d-2) then
                 count_diff(4) = count_diff(4)+1
              else if(diff.lt.5.0d-2) then
                 count_diff(5) = count_diff(5)+1              
              else if(diff.lt.6.0d-2) then
                 count_diff(6) = count_diff(6)+1
              else if(diff.lt.7.0d-2) then
                 count_diff(7) = count_diff(7)+1              
              else if(diff.lt.8.0d-2) then
                 count_diff(8) = count_diff(8)+1
              else if(diff.lt.9.0d-2) then
                 count_diff(9) = count_diff(9)+1              
              else if(diff.lt.10.0d-2) then
                 count_diff(10) = count_diff(10)+1
              else if(diff.lt.20.0d-2) then
                 count_diff(11) = count_diff(11)+1
              else if(diff.lt.30.0d-2) then
                 count_diff(12) = count_diff(12)+1
              else if(diff.lt.40.0d-2) then
                 count_diff(13) = count_diff(13)+1
              else if(diff.lt.50.0d-2) then
                 count_diff(14) = count_diff(14)+1
              else if(diff.lt.60.0d-2) then
                 count_diff(15) = count_diff(15)+1
              else
                 count_diff(16) = count_diff(16)+1
              endif
           endif
    
        enddo
     enddo
  enddo
  
  write(6,*) "Diff counters: use your own judgement..."
  write(6,"(A7)") "1","2","3","4","5","6","7","8","9","10","20","30","40","50","60","rest"
  write(6,"(I7)") count_diff(1:16)

  ! unlog for now and re-adjust minima
  if(min_energy.lt.0.0d0) then 
    eos_table(:,:,:,ienergy) = 10.0d0**eos_table(:,:,:,ienergy) + min_energy*1.1d0
  endif

  if(min_press.lt.0.0d0) then 
     eos_table(:,:,:,ipress) = 10.0d0**eos_table(:,:,:,ipress) + min_press*1.10d0
  endif

  ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ! extrapolation to higher temperatures if needed

  if(itemp_stop_nuclear.le.nlogtemp) then
     ! linear extrapolation of energy, pressure, and entropy.
     ! linear extrapolation of the chemical potentials (is this good or bad?)
     ! composition fixe
     
     ! if you really want to do this right you need to understand
     ! nuclear physics a lot better then we (man kind that is) do.
     ! Perhaps a more thermodymaically consistant way would be to
     ! solve the uniform matter RMF fields (for RMF calculations) that
     ! is what most people do for high temps (>~ 15MeV)

     do i=nlogrho,irho_start_comps,-1
        do k=1,nye
           !itemp_stop_nuclear+1 is the first point that the table cannot handle
           do j=itemp_stop_nuclear+1,nlogtemp,1
              !just keep compositions at highest initial table temp
              eos_table(i,j,k,ixn) = eos_table(i,itemp_stop_nuclear,k,ixn)
              eos_table(i,j,k,ixp) = eos_table(i,itemp_stop_nuclear,k,ixp)
              eos_table(i,j,k,ixa) = eos_table(i,itemp_stop_nuclear,k,ixa)
              eos_table(i,j,k,ixh) = eos_table(i,itemp_stop_nuclear,k,ixh)
              eos_table(i,j,k,iabar) = eos_table(i,itemp_stop_nuclear,k,iabar)
              eos_table(i,j,k,izbar) = eos_table(i,itemp_stop_nuclear,k,izbar)
              if (Hempel) then
                 eos_table(i,j,k,ixt) = eos_table(i,itemp_stop_nuclear,k,ixt)
                 eos_table(i,j,k,ixd) = eos_table(i,itemp_stop_nuclear,k,ixd)
                 eos_table(i,j,k,ix3he) = eos_table(i,itemp_stop_nuclear,k,ix3he)
                 eos_table(i,j,k,ix4li) = eos_table(i,itemp_stop_nuclear,k,ix4li)
              endif
              if (Hyperon) then
                 eos_table(i,j,k,ixL) = eos_table(i,itemp_stop_nuclear,k,ixL)
              endif


           enddo
        enddo
     enddo

     do i=nlogrho,irho_start_nuclear,-1
        do k=1,nye
           !itemp_stop_nuclear+1 is the first point that the table cannot handle
           do j=itemp_stop_nuclear+1,nlogtemp,1

              xx(1) = 10.0d0**logtemp(itemp_stop_nuclear-1)
              xx(2) = 10.0d0**logtemp(itemp_stop_nuclear)
              x = 10.0d0**logtemp(j)

              ! pressure
              ! the extrapolation should not be done in log(quantity);
              yy(1) = eos_table(i,itemp_stop_nuclear-1,k,ipress)
              yy(2) = eos_table(i,itemp_stop_nuclear,k,ipress)
              call extrapolate_lin(xx,yy,x,y)
              y1 = max(y,0.0d0)
              eos_table(i,j,k,ipress) = y1

              ! energy
              ! the extrapolation should not be done in log(quantity);
              yy(1) = eos_table(i,itemp_stop_nuclear-1,k,ienergy)
              yy(2) = eos_table(i,itemp_stop_nuclear,k,ienergy)
              call extrapolate_lin(xx,yy,x,y)
              y1 = max(y,0.0d0)
              eos_table(i,j,k,ienergy) = y1

              ! entropy
              ! the extrapolation should not be done in log(quantity);
              yy(1) = eos_table(i,itemp_stop_nuclear-1,k,ientropy)
              yy(2) = eos_table(i,itemp_stop_nuclear,k,ientropy)
              call extrapolate_lin(xx,yy,x,y)
              y1 = max(y,0.0d0)
              eos_table(i,j,k,ientropy) = y1

              ! proton chemical potential
              ! the extrapolation should not be done in log(quantity);
              yy(1) = eos_table(i,itemp_stop_nuclear-1,k,imup)
              yy(2) = eos_table(i,itemp_stop_nuclear,k,imup)
              call extrapolate_lin(xx,yy,x,y)
              y1 = max(y,0.0d0)
              eos_table(i,j,k,imup) = y1

              ! neutron chemical potential
              ! the extrapolation should not be done in log(quantity);
              yy(1) = eos_table(i,itemp_stop_nuclear-1,k,imun)
              yy(2) = eos_table(i,itemp_stop_nuclear,k,imun)
              call extrapolate_lin(xx,yy,x,y)
              y1 = max(y,0.0d0)
              eos_table(i,j,k,imun) = y1

           enddo
        enddo
     enddo
  endif
  ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  ! in high-density, low temperature region, fix composition to value
  ! at itemp_start_nuclear linear extrapolate the pressure, energy,
  ! entropy, mun and mup
  if(irho_start_nuclear.ne.1) then
     eos_table(1:irho_start_nuclear-1,:,:,ienergy) = 0.0d0
     eos_table(1:irho_start_nuclear-1,:,:,ipress) = 0.0d0
     eos_table(1:irho_start_nuclear-1,:,:,ientropy) = 0.0d0
  endif

  do i=nlogrho,irho_start_comps,-1
     do k=1,nye
        do j=itemp_start_nuclear-1,1,-1

           xx(1) = 10.0d0**logtemp(itemp_start_nuclear)
           xx(2) = 10.0d0**logtemp(itemp_start_nuclear+1)
           x = 10.0d0**logtemp(j)

           ! pressure
           yy(1) = eos_table(i,itemp_start_nuclear,k,ipress)
           yy(2) = eos_table(i,itemp_start_nuclear+1,k,ipress)
           call extrapolate_lin(xx,yy,x,y)
           if (y.lt.0.0d0.and.nopresswarning) then
              write(*,*) "check here if you have problems with negative pressure"
              nopresswarning = .false.
           endif
           y1 = y
           eos_table(i,j,k,ipress) = y1

           ! energy
           yy(1) = eos_table(i,itemp_start_nuclear,k,ienergy)
           yy(2) = eos_table(i,itemp_start_nuclear+1,k,ienergy)
           call extrapolate_lin(xx,yy,x,y)
           y2 = y
           eos_table(i,j,k,ienergy) = y2
           
           ! entropy
           yy(1) = eos_table(i,itemp_start_nuclear,k,ientropy)
           yy(2) = eos_table(i,itemp_start_nuclear+1,k,ientropy)
           call extrapolate_lin(xx,yy,x,y)
           y3 = max(y,0.0d0)
           if (y3.eq.0.0d0) then
              y3 = 0.0d0
           endif
           eos_table(i,j,k,ientropy) = y3

           
           ! copy over compositions at low T -- just keep em fixed at
           ! the value that they have a the lowest temperature in the
           ! initial table
           eos_table(i,j,k,ixn) = eos_table(i,itemp_start_nuclear,k,ixn)
           eos_table(i,j,k,ixp) = eos_table(i,itemp_start_nuclear,k,ixp)
           eos_table(i,j,k,ixa) = eos_table(i,itemp_start_nuclear,k,ixa)
           eos_table(i,j,k,ixh) = eos_table(i,itemp_start_nuclear,k,ixh)
           eos_table(i,j,k,iabar) = eos_table(i,itemp_start_nuclear,k,iabar)
           eos_table(i,j,k,izbar) = eos_table(i,itemp_start_nuclear,k,izbar)
           if (Hempel) then
              eos_table(i,j,k,ixt) = eos_table(i,itemp_start_nuclear,k,ixt)
              eos_table(i,j,k,ixd) = eos_table(i,itemp_start_nuclear,k,ixd)
              eos_table(i,j,k,ix3he) = eos_table(i,itemp_start_nuclear,k,ix3he)
              eos_table(i,j,k,ix4li) = eos_table(i,itemp_start_nuclear,k,ix4li)
           endif
           if (Hyperon) then
              eos_table(i,j,k,ixL) = eos_table(i,itemp_start_nuclear,k,ixL)
           endif

           ! mup
           yy(1) = eos_table(i,itemp_start_nuclear,k,imup)
           yy(2) = eos_table(i,itemp_start_nuclear+1,k,imup)
           call extrapolate_lin(xx,yy,x,y)
           y3 = y
           eos_table(i,j,k,imup) = y3

           ! mun
           yy(1) = eos_table(i,itemp_start_nuclear,k,imun)
           yy(2) = eos_table(i,itemp_start_nuclear+1,k,imun)
           call extrapolate_lin(xx,yy,x,y)
           y3 = y
           eos_table(i,j,k,imun) = y3
        enddo
     enddo
  enddo

  ! reset pressure/energy/entropy etc. to zero where not needed
  ! (between irho_start_comps and irho_start_nuclear, keep
  ! compositional information. Note the difference between
  ! irho_start_comps and irho_start_nuclear, see setup_my_table.F90
  
  write(*,*) "interpolation done"

end subroutine interpolate_initial_table_cubic_hermite

