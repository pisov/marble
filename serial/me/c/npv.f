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
        parameter (LS = 256)
        integer CA(LS,LS)          ! The lattice itself
        integer nuvac              ! The number of vacancies
        integer, dimension(:,:), allocatable :: vac !vacancy list
        integer ip(LS), im(LS)     ! The vectors with the PBC
        integer it, nt             ! Time steps
        parameter (nt = 100000)     ! Time steps
        real cv, cm                ! Concentration of vacancies and minuses
        real r                     ! random number
        parameter (cv = 0.01, cm = (1.0-cv)/2.0) ! conc
        character*32 :: filename

c   The coding of the lattice sites is as follows:
c      -1 : one of the particles kind
c       0 : vacancy
c       1 : the other of the kinds of particles
c end of the declarations block  
c    Initialize the lattice AND the random number generator
        call initCA1(LS, CA, cv, cm, nuvac) ! this is initialization
                                     ! with totally random distribution
        allocate(vac(nuvac,2))       ! allocate the vacancy list

        call updatevaclist(CA,LS,vac,nuvac)  ! make list of the vacancies
        
        call prepBC(LS, im, ip)      ! prepare the boundary conditions
        
        open (unit = 213, file='siminfo.txt')
              write(213,*) LS
        close(unit = 213)

      
c    Time loop
        do it = 0, nt+1
           call Latex(LS,CA,im,ip,exru,vac,nuvac) ! Make one lattice exchange
           call updatevaclist(CA,LS,vac,nuvac)  ! update the list of the vacancies
           if (MOD(it,10000).eq.0)then
     
              call infoparticle(LS,CA) !print info about the particles
             !dump CA to file
             write(filename,'(A5I0.8A4)')'data/step',it,'.dat'
              open (unit = 212, file=filename)
                 do i=1,LS
                    write(212,*) CA(i,:)
                 end do
              close(unit = 212)
          endif

        enddo

      deallocate(vac) 
c end of executions block
        end
c----------------------------------------        
!prints the number of each particles (usefull for monitoring)        
        subroutine  infoparticle(LS,CA)
        
        integer LS
        integer CA(LS,LS)
        integer a,b,c,i,j

        a=0
        b=0
        c=0

        !count them
        do i=1,LS
           do j=1,LS
              if(CA(i,j).eq.-1)then
              a=a+1
              elseif(CA(i,j).eq.0)then
              b=b+1
              elseif(CA(i,j).eq.1)then
              c=c+1
              endif
           enddo
        enddo
        !print them
        write(*,*)a,b,c


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
       subroutine Latex(LS,CA,im,ip,exru, vac,nuvac) ! Performs the vacancy exchanges
       integer LS, LS2,nuvac, is, ix, iy, il, ir
       integer CA(LS,LS)
       integer, dimension(:,:) :: vac(nuvac,2)
       integer ip(LS),im(LS)

       call random_number(r)
       if (r.lt.0.857)then
       do is = 1, nuvac
c       generate random coordinates of the current site to be checked
       ix = vac(is,1)
       iy = vac(is,2)
       if(CA(ix,iy).ne.0)then
!         write(*,*)'arrrrr'
       endif
                    ! then choose the direction of exchange
          call random_number(r)
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
       enddo

       else
       do is = 1, nuvac
c       generate random coordinates of the current site to be checked
       ix = vac(is,1)
       iy = vac(is,2)
                    ! then choose the direction of exchange
          call random_number(r)
          if (r.lt.0.25) then ! move right
             il = CA(im(ix),iy)
             ir = CA(ip(ix),iy)
                CA(ix,iy)=ir
                CA(ip(ix),iy)=0
          elseif(r.lt.0.5) then ! move down
             il = CA(ix,im(iy))
             ir = CA(ix,ip(iy))
                CA(ix,iy)=ir
                CA(ix,ip(iy))=0

          elseif (r.lt.0.75) then !move left
             il = CA(ip(ix),iy)
             ir = CA(im(ix),iy)
                CA(ix,iy)=ir
                CA(im(ix),iy)=0 
          else                     ! move up
             il = CA(ix,ip(iy))
             ir = CA(ix,im(iy))
                CA(ix,iy)=ir
                CA(ix,im(iy))=0
          endif
      
       enddo
       endif

       return
       end
c .. IF SOMETHING is wrong it is here in this procedure !!
c-----------------------------------------
        subroutine initCA1(LS, CA, cv, cm, nuvac)
        implicit none 
c      Totally random spread of the three components -1, 0 (vacancy), +1
        integer LS,  nuvac
        integer CA(LS,LS)
        real cv, cm, acm
        integer :: ii,jj
        real :: r
        
        acm = cm + cv

        nuvac=0
        do ii = 1, LS
           do jj = 1, LS
           call random_number(r)
           if (r.lt.cm) then
              CA(ii,jj) = -1
           elseif (r.lt.acm) then
              CA(ii,jj) =  0
              nuvac = nuvac + 1
           else
              CA(ii,jj) = 1
           endif
           enddo
        enddo

        return
        end
c----------------------------------------- 
        subroutine updatevaclist(CA,LS,vac,nuvac)  ! make list of the vacancies
        integer LS,counter
        integer,dimension(:,:) :: vac(nuvac,2)
        integer CA(LS,LS)
        
        counter=1
        do ii=1, LS
           do jj=1, LS
              if(CA(ii,jj) .eq. 0)then
                 vac(counter,1)=ii
                 vac(counter,2)=jj
                 counter=counter+1
              endif
          enddo
       enddo
        

        return
        end
