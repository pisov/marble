        program npv1
c     Simulates a binary system with vacancies, they play the only (main) role
c     In this version are not used pointers to the coordinates of the vacancies
c     but they are placed on a lattice.
c     Thus it is specially suited for the cases with large concentrations of 
c     vacancies
c     Well branched system of monitoring is to be developed.
c
c     (c) Vesselin Tonchev, Marcin Minkowski, 2016, Sofia - Warsaw
c     (c) Vesselin Tonchev, 2017, Sofia   
        integer LS                 ! Size of the lattice
        parameter (LS = 200)
        integer CA(LS,LS)          ! The lattice itself
        integer ip(LS), im(LS)     ! The vectors with the PBC
        integer it, nt             ! Time steps
        integer ans                ! writing switch
        parameter (nt = 10000)     ! Time steps
c        integer exru(2,-1:1,-1:1)  ! NO Exclusion rule in this version !!!
        real cv, cm                ! Concentration of vacancies and minuses
        real r                     ! random number
        integer idum               ! current seed of the RNG
        common /ixx/ idum
        parameter (cv = 0.01, cm = (1.0-cv)/2.0) ! conc
c   The coding of the lattice sites is as follows:
c      -1 : one of the particles kind
c       0 : vacancy
c       1 : the other of the kinds of particles
c end of the declarations block  
c    Initialize the lattice AND the random number generator
        call ofi                     ! Open files
        call initCA1(LS, CA, cv, cm) ! this is initialization
                                     ! with totally random distribution
        ans = 0
        call wrCA(ans, LS, CA)       ! write the lattice    
c   NO!     call definerule(exru)        ! define the exclusion rule
        call prepBC(LS, im, ip)      ! prepare the boundary conditions
        
c    Time loop
        do it = 1, nt
         call Latex(LS,CA,im,ip) !,exru) ! Make one lattice exchange
         ans = 1
         call wrCA(ans, LS, CA)       ! .. and save the lattice
        enddo
c
c end of executions block
        end
c-----------------------------------------
        subroutine prepBC(LS, im, ip) ! prepare the boundary conditions
        integer LS, is
        integer ip(LS), im(LS)
        im(1) = LS
        do is = 2, LS
          im(is) = is -1
        enddo
        do is = 1, LS-1
          ip(is) = is +1
        enddo
          ip(IS) = 1
        return
        end
c-----------------------------------------
        subroutine ofi  ! open the files , not tweaked well ..
         open (11, file = 'redini.dat')
         open (12, file = 'bluini.dat')
         open (13, file = 'vacini.dat')
         open (21, file = 'red.dat')
         open (22, file = 'blu.dat')
         open (23, file = 'vac.dat')
        return
        end
c-----------------------------------------
       subroutine wrCA(ans, LS, CA)  ! writes the lattice
       integer ans
       integer LS
       integer CA(LS,LS)
       close(21)
       close(22)
       close(23)
        open (21, file = 'red.dat')
         open (22, file = 'blu.dat')
         open (23, file = 'vac.dat')
       if (ans.eq.0) then
        do ii = 1, LS
          do jj=1,LS
            Cat=CA(ii,jj)
            if (Cat.eq.-1) then
             write (11,*) ii,jj
            elseif(Cat.eq.1) then
             write (12,*) ii,jj
            else
             write (13,*) ii,jj
            endif
          enddo
        enddo
       else
       do ii = 1, LS
          do jj=1,LS
            Cat=CA(ii,jj)
            if (Cat.eq.-1) then
             write (21,*) ii,jj
            elseif(Cat.eq.1) then
             write (22,*) ii,jj
            else
             write (23,*) ii,jj
            endif
          enddo
        enddo
       endif
       return
       end
c-----------------------------------------
       subroutine Latex(LS,CA,im,ip) !,exru) ! Performs the vacancy exchanges
       integer LS, LS2, is, ix, iy, il, ir
       integer CA(LS,LS)
       integer ip(LS),im(LS)
c       integer exru(2,-1:1,-1:1)
       integer idum, iaaa
       common /ixx/ idum
       LS2 = LS*LS                 ! Total number of attempts = lattice area
       do is = 1, LS2
c       generate random coordinates of the current site to be checked
       r = ran2r(idum)
       ix = int(r*LS)+1
       r = ran2r(idum)
       iy = int(r*LS)+1
       if (CA(ix,iy).eq.0) then    ! .. whether there is a vacancy there
          r = ran2r(idum)          ! then choose the direction of exchange
          if (r.lt.0.25) then ! exclude right
             il = CA(im(ix),iy)
             ir = CA(ip(ix),iy)
             iaaa = il*ir
             if (iaaa.eq.1) then
                CA(ix,iy)=ir
                CA(ip(ix),iy)=0
             endif
          elseif(r.lt.0.5) then ! exclude down
             il = CA(ix,im(iy))
             ir = CA(ix,ip(iy))
             iaaa = il*ir
             if (iaaa.eq.1) then
               CA(ix,iy)=il
               CA(ix,im(iy))=0
             endif
          elseif (r.lt.0.75) then !exclude left
             il = CA(ip(ix),iy)
             ir = CA(im(ix),iy)
             iaaa = il*ir
             if (iaaa.eq.1) then
               CA(ix,iy)=ir
               CA(im(ix),iy)=0 
             endif
          else                     ! exlcude up
             il = CA(ix,ip(iy))
             ir = CA(ix,im(iy))
             iaaa = il*ir
             if (iaaa.eq.1) then
              CA(ix,iy)=ir
              CA(ix,ip(iy))=0
             endif          
          endif
       endif
       enddo
       return
       end
c .. IF SOMETHING is wrong it is here in this procedure !!
c-----------------------------------------
       subroutine definerule(exru)   ! Defines the exchange rule
       integer exru(2,-1:1,-1:1)
c   Define the exclusions
        do ir = -1,1
            do il = -1,1
              if (ir.eq.il) then ! Exchange !
                exru(1,il,ir) = ir
                exru(2,il,ir) = 0
              else
                exru(1,il,ir) = 0 ! No exchange !
                exru(2,il,ir) = ir  ! Should be also defined    
              endif                 ! because the code is UNCONDITIONAL
            enddo
        enddo
       return
       end
c-----------------------------------------
        subroutine initCA1(LS, CA, cv, cm)
c      Totally random spread of the three components -1, 0 (vacancy), +1
        integer LS
        integer CA(LS,LS)
        real cv, cm, acm
        integer idum
        common /ixx/ idum
        call randomize
        acm = cm + cv
        do ii = 1, LS
           do jj = 1, LS
           r = ran2r(idum)
           if (r.lt.cm) then
              CA(ii,jj) = -1
           elseif (r.lt.acm) then
              CA(ii,jj) =  0
           else
              CA(ii,jj) = 1
           endif
           enddo
        enddo
        return
        end
c----------------------------------------- 
        subroutine randomize  ! randomization of the RNG
        integer ittt(3)
        integer idum
        common /ixx/ idum
        real r
        call itime (ittt)
        idum=(ittt(1)+1)*(ittt(2)+2)*(ittt(3)+3)+(ittt(3)+1)*(ittt(1)+4)
        idum = - idum
        r = ran2init(idum) ! note the version of the RNG used
        return
        end
c----------------------------------------------------------------------------------        
c   Random Number Generator by Numerical Recipes:
c   (to be called for first time with a negative integer that is kept further in the
c    common area gnrsvw)
      FUNCTION ran2r(idum) ! from this code is removed the if block
c                            producing the initial seed
c
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2r,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2r=min(AM*iy,RNMX)
      return
      END
c---------------------------------------------------      
C  (C) Copr. 1986-92 Numerical Recipes Software 
      FUNCTION ran2init(idum) ! only used in the randomization code
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2init,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2init=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .
