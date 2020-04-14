!-*-f90-*-
subroutine read_initial_table(eosfilename)

  use eos_table_module
  implicit none

  integer i,j,k,l !counters
  real*8 dbuffer1,dbuffer2 !buffer

  character(len=256) eosfilename
  character(len=1024) string_buffer1,string_buffer2

  real*8 my_xn,my_xp,my_xa,my_xh,my_abar,my_zbar
  real*8 my_one,my_ye

  call allocate_eos

  open(unit=473,file=trim(adjustl(eosfilename)),form='unformatted', status='old')

  write(*,*) "Reading in table data"
  read(473) table_data
  write(*,*) "Done reading in table data"

  !just a check for unknownness
  do j=1,nt_temp
     table_logtemp(j) = -1.0d0 + (j-1)*0.04d0
     do k=1,nt_yp
        table_yp(k) =  0.01d0 + (k-1)*0.01d0
        do i=1,nt_rho
           if (table_data(j,k,i,ixh).eq.0.0d0) then
              table_data(j,k,i,iabar) = 1.0d0
              table_data(j,k,i,izbar) = 1.0d0
           endif
           if (j.eq.1.and.k.eq.1) then
              table_logrho(i) = log10(10.0d0**(-12.0d0 + &
                   (i-1)*0.04d0)*1.66053928d-24*1.0d39)
           endif
           if (lump_into_np) then
              !neutrons
              table_data(j,k,i,11) = table_data(j,k,i,11) + &
                   table_data(j,k,i,18)*0.5d0 + & !deuterons
                   table_data(j,k,i,19)*2.0d0/3.0d0 + & !tritons
                   table_data(j,k,i,20)*1.0d0/3.0d0 + & !helions
                   table_data(j,k,i,21)*1.0d0/4.0d0 + & !lithions
                   table_data(j,k,i,22) !lambdas
              !protons
              table_data(j,k,i,12) = table_data(j,k,i,12) + &
                   table_data(j,k,i,18)*0.5d0 + & !deuterons
                   table_data(j,k,i,19)*1.0d0/3.0d0 + & !tritons
                   table_data(j,k,i,20)*2.0d0/3.0d0 + & !helions
                   table_data(j,k,i,21)*3.0d0/4.0d0  !lithions
                   
              table_data(j,k,i,18) = 0.0d0
              table_data(j,k,i,19) = 0.0d0
              table_data(j,k,i,20) = 0.0d0
              table_data(j,k,i,21) = 0.0d0
              table_data(j,k,i,22) = 0.0d0
           endif
           do l=1,nvars
              if(table_data(j,k,i,l).ne.table_data(j,k,i,l)) then
                 write(*,*) table_data(j,k,i,l)
                 write(*,*) j,k,i,l," something's wrong here!!!"
                 stop "read_table.F90: BAD!"
              endif
           enddo
        enddo
     enddo
  enddo

  if (lump_into_np) then
     Hempel = .false.
  endif

  close(473)

end subroutine read_initial_table
