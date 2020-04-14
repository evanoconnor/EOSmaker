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

  open(unit=473,file=trim(adjustl(eosfilename)))

  !this routine is EOS format specific... adjust as necessary

  !loop over T
  do j=1,nt_temp
     read(473,*) string_buffer1
     read(473,*) string_buffer1,string_buffer2
     read(473,*) dbuffer1,dbuffer2
     read(473,*)
     !dbuffer is temperature of the first chunk
     table_logtemp(j) = log10(dbuffer2)

     write(6,*) "Reading in temperature: ",table_logtemp(j)

     !loop over yp next
     do k=1,nt_yp
        !loop over rho last
        do i=1,nt_rho

           read(473,*) table_logrho(i),dbuffer1, &
                table_yp(k),(table_data(j,k,i,l), l=1,nvars)

           my_xn = table_data(j,k,i,ixn)
           my_xp = table_data(j,k,i,ixp)
           my_xa = table_data(j,k,i,ixa)
           my_xh = table_data(j,k,i,ixh)
           my_abar = table_data(j,k,i,iabar)
           my_zbar = table_data(j,k,i,izbar)

           call one_time_fix_comps(my_xn,my_xp,my_xa,my_xh,my_abar,my_zbar,table_yp(k))
           
           table_data(j,k,i,ixn) = my_xn
           table_data(j,k,i,ixp) = my_xp
           table_data(j,k,i,ixa) = my_xa
           table_data(j,k,i,ixh) = my_xh
           table_data(j,k,i,iabar) = my_abar
           table_data(j,k,i,izbar) = my_zbar

           my_ye = table_data(j,k,i,ixp)+0.5d0*table_data(j,k,i,ixa)+table_data(j,k,i,izbar)*table_data(j,k,i,ixh)/table_data(j,k,i,iabar)
           my_one =  table_data(j,k,i,ixp)+table_data(j,k,i,ixn)+table_data(j,k,i,ixa)+table_data(j,k,i,ixh)

           !just a check for unknownness
           do l=1,nvars
              if(table_data(j,k,i,l).ne.table_data(j,k,i,l)) then
                 write(*,*) j,k,i,l," something's wrong here!!!"
                 stop "read_table.F90: BAD!"
              endif
           enddo

        enddo
     enddo
  enddo
  

  close(473)

end subroutine read_initial_table
