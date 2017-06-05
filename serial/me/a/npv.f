        program npv1
c     Simulates a binary system with vacancies, they play the only (main) role
c     In this version are not used pointers to the coordinates of the vacancies
c     but they are placed on a lattice.
c     Thus it is specially suited for the cases with large concentrations of 
c     vacancies
c     Well branched system of monitoring is to be developed.
c
c     (c) Vesselin Tonchev, Marcin Minkowski, 2016, Sofia - Warsaw
        integer LS                 ! Size of the lattice
        parameter (LS = 128)
        integer CA(LS,LS)          ! The lattice itself
        integer ip(LS), im(LS)     ! The vectors with the PBC
        integer it, nt             ! Time steps
        integer ans                ! writing switch
        parameter (nt = 10000)     ! Time steps
        integer exru(2,-1:1,-1:1)  ! Exclusion rule
        real cv, cm                ! Concentration of vacancies and minuses
        real r                     ! random number
        integer idum               ! current seed of the RNG
        common /ixx/ idum
        parameter (cv = 0.01, cm = (1.0-cv)/2.0) ! conc
        integer aa,bb,cc,dd,remain
        real enpyx,enpyy
        integer nuclu
        character*32 :: filename

         write(*,*)'LS=',LS,'nt=',nt,'mod='
         write(*,*)' it    ',' nuclu   ',' enpyx   ',' enpyy  ','remain'
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
!        call wrCA(ans, LS, CA)       ! write the lattice    
        call prepBC(LS, im, ip)      ! prepare the boundary conditions
        
       
c    Time loop
        do it = 0, nt+1
         call Latex(LS,CA,im,ip,exru) ! Make one lattice exchange
         ans = 1
        
!        call wrCA(ans, LS, CA)       ! .. and save the lattice
        if (MOD(it,10000).eq.0)then
     

        call enpy(LS,CA,enpyx,enpyy)

              write(*,*)'iter # ', it
!        call draw(LS,CA,it)

        remain=0
        do, i=1,LS
           do, j=1,LS
           aa=CA(ip(i),j)
           bb=CA(im(i),j)
           cc=CA(i,ip(j))
           dd=CA(i,im(j))
           if((CA(i,j).eq.0).and.(((aa.eq.bb).and.(aa.ne.0))
     &     .or.(cc.eq.dd).and.(cc.ne.0)) )then
              remain=remain+1
           endif
           enddo
        enddo
        call clus(LS,CA,nuclu) 
         write(*,*)it,nuclu,enpyx,enpyy,remain
!        write(*,*)'remain',remain
!        write(*,*)'enpyx', enpyx,'enpyy', enpyy,'sum',enpyx+enpyy
           end if
        enddo
        
c
c end of executions block
        end
c-----------------------------------------
        subroutine clus(LS,CA,nuclu) 
       implicit none
       integer LS,   ix, iy,i,j
       integer CA(LS,LS), label(LS,LS,2), clusters(LS*LS)
       integer max_label,nuclu,newlabel,oldlabel
       integer, dimension(:,:),allocatable :: clusterdata, shortclusdata
       integer :: counter,counter2,nusizes,typeclus
       integer,dimension(:), allocatable :: sizes
       logical, dimension(:), allocatable :: mask
        
!label(*,*,1) -> labels
!label(*,*,2) -> type of cluster


       label(1:LS,1:LS,1:2)=0

       max_label = 1
       do, ix=1,LS
          do, iy=1,LS
          ! if bouth up and left are the same
            if((CA(ix-1,iy).eq.CA(ix,iy)).and.(CA(ix,iy).eq.
     &        CA(ix,iy-1)).and.(ix.ne.1).and.(iy.ne.1))then
                   oldlabel=label(ix-1,iy,1)
                   newlabel=label(ix,iy-1,1)
                   label(ix,iy,1)=newlabel

                   do i=1,ix
                      do j=1, iy
                         if(label(i,j,1).eq.oldlabel)then
                            label(i,j,1)=newlabel
                         endif
                      enddo
                   enddo

          !if up is the same
            elseif((CA(ix-1,iy).eq.CA(ix,iy)).and.(ix.ne.1))then

                   label(ix,iy,1)=label(ix-1,iy,1)
                   label(ix,iy,2)=label(ix-1,iy,2)


           !if left is the same
             elseif((CA(ix,iy-1).eq.CA(ix,iy)).and.(iy.ne.1))then
                   label(ix,iy,1)=label(ix,iy-1,1)
                   label(ix,iy,2)=label(ix,iy-1,2)
             else
                  label(ix,iy,1)= max_label
                   max_label=max_label+1
             endif
                   if(CA(ix,iy).eq.-1)then
                      label(ix,iy,2)=1
                   elseif(CA(ix,iy).eq.0)then
                      label(ix,iy,2)=2
                   else
                      label(ix,iy,2)=3
                   endif
                  
          enddo
       enddo
        
        !boundary conditions
        do ix=1,LS
           !down
           if(CA(LS,ix).eq.CA(1,ix))then
                   oldlabel=label(LS,ix,1)
                   newlabel=label(1,ix,1)
                   label(LS,ix,1)=newlabel

                   do i=1,LS
                      do j=1, LS
                         if(label(i,j,1).eq.oldlabel)then
                            label(i,j,1)=newlabel
                         endif
                      enddo
                   enddo
           endif
           !right
           if(CA(ix,LS).eq.CA(ix,LS))then
                   oldlabel=label(ix,LS,1)
                   newlabel=label(ix,LS,1)
                   label(LS,ix,1)=newlabel

                   do i=1,LS
                      do j=1, LS
                         if(label(i,j,1).eq.oldlabel)then
                            label(i,j,1)=newlabel
                         endif
                      enddo
                   enddo
           endif     

        enddo
!make a list of the cluster names
        do i=1,LS
           do,j=1, LS
              clusters((i-1)*LS+j)=label(i,j,1)
           enddo
        enddo
!sort it        
        call buble(LS*LS,clusters)
        oldlabel=-1
!count it
        nuclu=0
        do i=1,LS*LS
              newlabel=clusters(i)
              if(oldlabel.ne.newlabel)then
                 nuclu=nuclu+1
              endif

              oldlabel=clusters(i)
        enddo

!fulling clusterdata with info about the clusters
!clusterdata(*,1) -> label
!clusterdata(*,2) -> size
!clusterdata(*,3) -> type

        !labels
        allocate(clusterdata(nuclu,3))
        counter=1
        oldlabel=-1
        do i=1,LS*LS
           newlabel=clusters(i)
           if(oldlabel.ne.newlabel)then
              clusterdata(counter,1)=newlabel
              counter=counter+1
           endif
           oldlabel=clusters(i)
        enddo
        

        !size and type of clusters
        do i=1, nuclu
           counter=0
           typeclus=3
           newlabel=clusterdata(i,1)

           do ix=1, LS
              do iy=1, LS
                 if(newlabel.eq.label(ix,iy,1))then
                    counter = counter+1
                    typeclus = label(ix,iy,2)!not very smart :/
                 endif
              enddo
           enddo
           
          clusterdata(i,2)=counter
          clusterdata(i,3)=typeclus
        enddo

        open (unit = 111, file='clusterdata.txt')
        do i=1,nuclu
           write(111,*) clusterdata(i,:)
        end do
        close(unit = 111)
c--  --  -- --  --  --  --

        allocate(sizes(nuclu))       
        allocate(mask(nuclu))       
!       shortclusdata(typesize,cout,type)
!       sizes(nuclu)
!       buble(sizesort)
!       
!make a list of the cluster sizes
       


       sizes=clusterdata(:,2)




!sort it        
        call buble(nuclu,sizes)

         do i=1, size(sizes)
           write(*,*)sizes(i)
        end do 

       
        oldlabel=clusters(1)
!count it
        nusizes=0
!reusing the variables new/oldlabel as new/oldsize of clusters        
        oldlabel=-1
        do i=1,nuclu
              newlabel=sizes(i)
              if(oldlabel.ne.newlabel)then
                 nusizes=nusizes+1
              endif
              oldlabel=sizes(i)
        enddo
        
        
       allocate(shortclusdata(nusizes,3))

!shortclusterdata(*,1) -> size of the clusters
!shortclusterdata(*,2) -> count of the clusters of this size
!shortclusterdata(*,3) -> type of the clusters of this size 

!fulling shortclusterdata with data

       !sizes
       oldlabel = -1
       counter = 1
       

       do i=1,nuclu
          newlabel=sizes(i)
       mask=newlabel.eq.sizes
          if(newlabel.ne.oldlabel)then
             shortclusdata(counter,1)=newlabel
             shortclusdata(counter,2)=count(mask)
             counter=counter+1
          elseif(newlabel.eq.oldlabel)then
          endif

          oldlabel=sizes(i)
       enddo

       open (unit = 112, file='shortclusdata.txt')
        do i=1,nusizes
           write(112,*) shortclusdata(i,1:2 )
        end do
        close(unit = 112)


        deallocate(shortclusdata)
        deallocate(clusterdata)
        deallocate(sizes)       
        deallocate(mask)       
        return
        end
 
c---------------------
      subroutine buble(LS,clusters)
      implicit none
      integer LS
      integer,dimension(:) :: clusters(LS)
      integer :: temp
      INTEGER :: i, j
      LOGICAL :: swapped

      do j = LS-1, 1, -1
      swapped = .FALSE.
      do i = 1, j
      if (clusters(i) .gt. clusters(i+1)) then
        temp = clusters(i)
        clusters(i) = clusters(i+1)
        clusters(i+1) = temp
        swapped = .TRUE.
      endif
      end do
      if (.not. swapped)exit
      end do
      
      
      return
      end


c-----------------------------------------
        subroutine enpy(LS,CA,enpyx,enpyy) 
       integer LS,   ix, iy,i
       integer CA(LS,LS)
       integer enpy1(LS)
       integer enpy2(LS)
       real enpyx,enpyy
       integer br1,br2,aa,bb,cc,dd

       bb=-2
       dd=-2
       do ix=1,LS
          br1=0
          br2=0

          do iy=1,LS
             aa=CA(ix,iy)
             cc=CA(iy,ix)
           

             if(aa.eq.0)then
                i=iy
                do while (aa.eq.0.and.i.ne.1)
                   i=i-1
                   aa=CA(ix,i)
                enddo
             endif
             if(aa.eq.bb.and.aa.ne.0)then
                br1=br1+1
             endif


             if(cc.eq.0)then
                i=ix
                do while (cc.eq.0.and.i.ne.1)
                   i=i-1
                   cc=CA(i,iy)
                enddo
             endif
             if(cc.eq.dd.and.cc.ne.0)then
                br2=br2+1
             endif


             dd= CA(iy,ix)
             bb= CA(ix,iy)
          enddo
          enpy1(ix)=LS-1-br1
          enpy2(ix)=LS-1-br2
       enddo
        enpyx=SUM(enpy1)/real(LS)
        enpyy=SUM(enpy2)/real(LS)

        return
        end

c-----------------------------------------
        subroutine draw(LS,CA,it)
        use utils
        integer LS,it
        integer CA(LS,LS)
        double precision buf(LS,LS)
        integer i,j,row,column
        character*32 :: filename


        buf(:,:)=0
        do, i=1,LS
           do, j=1,LS
              do, row=0, 0
                 do, column=0,0
                    if (CA(i,j).eq.-1)then
                       buf(i*4-row,j*4-column)=255.0
                   elseif(CA(i,j).eq.0)then
                       buf(i*4-row,j*4-column)=127.0
                   elseif(CA(i,j).eq.1)then
                      buf(1*4-row,j*4-column)=0.0
                   end if
                 end do
              end do
            end do
         end do


!      write(filename,'(A8,I,A4)')'img/step',it,'.ppm'
      write(filename,'(A4I0.8A4)')'img/step',it,'.ppm'
      call ppmwrite(filename,buf(:,:),LS,LS)
        return
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
       subroutine Latex(LS,CA,im,ip,exru) ! Performs the vacancy exchanges
       integer LS, LS2, is, ix, iy, il, ir
       integer CA(LS,LS)
       integer ip(LS),im(LS)
       integer exru(2,-1:1,-1:1)
       integer idum
       integer el,er,eu,ed
       common /ixx/ idum

       el=0
       er=0
       eu=0
       ed=0
       LS2 = LS*LS                 ! Total number of attempts = lattice area
       r=ran2r(idum)
        if (r.lt.0.857)then
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
             if(il.eq.ir.and.il.ne.0)then
                CA(ix,iy)=il
                CA(ip(ix),iy)=0
                er=er+1
             endif
          elseif(r.lt.0.5) then ! exclude down
             il = CA(ix,im(iy))
             ir = CA(ix,ip(iy))
             if(il.eq.ir.and.il.ne.0)then
                CA(ix,iy)=il
                CA(ix,ip(iy))=0
                ed=ed+1
             endif
          elseif (r.lt.0.75) then !exclude left
             il = CA(ip(ix),iy)
             ir = CA(im(ix),iy)
             if(il.eq.ir.and.il.ne.0)then
                CA(ix,iy)=il
                CA(im(ix),iy)=0 
                el=el+1
             endif
          else                     ! exlcude up
             il = CA(ix,ip(iy))
             ir = CA(ix,im(iy))
             if(il.eq.ir.and.il.ne.0)then
                CA(ix,iy)=il
                CA(ix,im(iy))=0
                eu=eu+1
             endif
          endif
       endif
       enddo

       else
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
                CA(ix,iy)=ir
                CA(ip(ix),iy)=0
          elseif(r.lt.0.5) then ! exclude down
             il = CA(ix,im(iy))
             ir = CA(ix,ip(iy))
                CA(ix,iy)=ir
                CA(ix,ip(iy))=0
          elseif (r.lt.0.75) then !exclude left
             il = CA(ip(ix),iy)
             ir = CA(im(ix),iy)
                CA(ix,iy)=ir
                CA(im(ix),iy)=0 
          else                     ! exlcude up
             il = CA(ix,ip(iy))
             ir = CA(ix,im(iy))
                CA(ix,iy)=ir
                CA(ix,im(iy))=0
          endif
       endif
      
       enddo
       endif
!
!       write(*,*)el,er,eu,ed
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
