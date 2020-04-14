!-*-f90-*-
      subroutine trilin_interp_table                &
          (lr,lt,y,value,array,nrho,ntemp,          &
           nye,logrho,logtemp,ye,d1,d2,d3)

      implicit none

      integer nrho,ntemp,nye
      double precision lr,lt,y,value
      double precision logrho(nrho)
      double precision logtemp(ntemp)
      double precision ye(nye)
      double precision d1,d2,d3
      double precision array(nrho,ntemp,nye)

      integer ixplus,ixminus
      integer iyplus,iyminus
      integer izplus,izminus

      double precision ddxr,ddyr,ddzr
      double precision q1,q2,q3,q4
      double precision q5,q6,q7,q8

      call trilin_findindex(nrho,logrho,lr,ixminus,ixplus)
      call trilin_findindex(ntemp,logtemp,lt,iyminus,iyplus)
      call trilin_findindex(nye,ye,y,izminus,izplus)

      ddxr = (lr - logrho(ixminus)) / &
          (logrho(ixplus) - logrho(ixminus))

      ddyr = (lt - logtemp(iyminus)) / &
          (logtemp(iyplus) - logtemp(iyminus))

      ddzr = (y - ye(izminus)) / &
          (ye(izplus) - ye(izminus))


      !set up base points

      q1 = array(ixminus  ,iyminus  ,izminus  )
      q2 = array(ixminus+1,iyminus  ,izminus  )
      q3 = array(ixminus  ,iyminus+1,izminus  )
      q4 = array(ixminus  ,iyminus  ,izminus+1)
      q5 = array(ixminus+1,iyminus+1,izminus  )
      q6 = array(ixminus+1,iyminus  ,izminus+1)
      q7 = array(ixminus  ,iyminus+1,izminus+1)
      q8 = array(ixminus+1,iyminus+1,izminus+1)


      d1 = (q2-q1)/(logrho(ixplus)-logrho(ixminus))
      d2 = (q3-q1)/(logtemp(iyplus)-logtemp(iyminus))
      d3 = (q4-q1)/(ye(izplus)-ye(izminus))

      call trilinterp(ddxr,ddyr,ddzr,q1,q2,q3,q4,&
          q5,q6,q7,q8,value)


      end subroutine

      subroutine trilin_findindex(n,array,value,iminus,iplus)
        
        implicit none 
        integer n,iminus,iplus
        double precision array(n)
        double precision value
        double precision buffer
        
        buffer = (value-array(1))/(array(n)-array(1)) * (n-1)*1.0d0 
        iminus = int(buffer)+1
        iplus  = int(buffer)+2
        
        if(iminus.lt.1.or.iplus.lt.1) then
           write(6,*) n,value,array(n),array(1),iminus,iplus
           call flush(6)
           stop "this is very bad"
        endif
        
      end subroutine trilin_findindex


      subroutine trilinterp(ddxr,ddyr,ddzr, &
          q1,q2,q3,q4,q5,q6,q7,q8,res)

      implicit none

      double precision ddxr,ddyr,ddzr
      double precision q1,q2,q3,q4
      double precision q5,q6,q7,q8
      double precision res

      
      res = q1 +    &                                   
     
          ddxr *      &
          (q2 - q1) + &
     
          ddyr *      &
          (q3 - q1) + &
     
          ddzr *      &
          (q4 - q1) + &
     
          ddxr * ddyr * &
          (q5 - q3 - q2 + q1) + &
     
          ddxr * ddzr * &
          (q6 - q4 - q2 + q1) + &
     
          ddyr * ddzr * &
          (q7 - q4 - q3 + q1) + &
     
          ddxr * ddyr * ddzr *  &                   
          (q8 - q7 - q6 + q4 - q5 + q3 + q2 - q1)



      end subroutine


