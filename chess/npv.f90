        program npv1
!c     Simulates a binary system with vacancies, they play the only (main) role
!c     In this version are not used pointers to the coordinates of the vacancies
!c     but they are placed on a lattice.
!c     Thus it is specially suited for the cases with large concentrations of 
!c     vacancies
!c     Well branched system of monitoring is to be developed.
!c
!c     (c) Vesselin Tonchev, Marcin Minkowski, 2016, Sofia - Warsaw
        
        use utils
        use mpi
        implicit none

        integer LS                 ! Size of the lattice
        parameter (LS = 128)
        integer ip(LS), im(LS)     ! The vectors with the PBC
        integer it, nt             ! Time steps
        parameter (nt = 1000000)     ! Time steps
        integer exru(2,-1:1,-1:1)  ! Exclusion rule
        real cv, cm                ! Concentration of vacancies and minuses
        real r                     ! random number
        integer idum               ! current seed of the RNG
        common /ixx/ idum
        parameter (cv = 0.01, cm = (1.0-cv)/2.0) ! conc
!c   The coding of the lattice sites is as follows:
!c      -1 : one of the particles kind
!c       0 : vacancy
!c       1 : the other of the kinds of particles
!c end of the declarations block  
!c    Initialize the lattice AND the random number generator

        integer seed, len, junk,streamnum, nstreams
        integer gtype,stream
!mpi

        integer :: i, j, n,  ncount , mcount ,m , sendcnts , int_size,a
       integer, dimension(:,:),allocatable :: CA,sCA
       integer, dimension(:),allocatable :: displs, sendcounts
       integer(kind=mpi_address_kind) :: lb, extent
       
       character*32 :: filename
       double precision, dimension(:,:), allocatable :: buf
       integer :: row, column
       integer :: ii,jj, nuvac
       integer,dimension(:),allocatable :: black,wite

       integer MPI_COMM_CART,MPI_COMM_TwoD,MPI_ONE_ROW, MPI_ONE_COL, MPI_BLOCK, MPI_BLOCKA
       integer :: mpi_rowo,mpi_colo,mpi_rowon
        integer,dimension(2) :: crd
        integer, dimension(2) :: dims
        logical, dimension(2) :: periods
        integer :: ierror, my_rank, comm_size, up, down, left, right
        double precision :: time
 

        call MPI_Init(ierror)
        call MPI_Comm_size(MPI_COMM_WORLD,comm_size,ierror)


        allocate(displs(comm_size))
        allocate(sendcounts(comm_size))


        call  MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierror)
        call  MPI_Dims_create(comm_size,2,dims,ierror)
        periods(2) = .true.
        periods(1) = .true.

        call  MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,.true.,MPI_COMM_TwoD, ierror)
        call  MPI_Comm_rank(MPI_COMM_TwoD,my_rank,ierror)
      


         !Calculate stripes
        mcount = LS / dims(1)
        ncount = LS / dims(2)
        m=LS
         
        do n=1,comm_size
        call MPI_Cart_coords(MPI_COMM_TwoD,n-1,2,crd,ierror)
        i=crd(1)
        j=crd(2)
        if(my_rank.eq.0)then
        endif
        displs(n)=(j*m*ncount +i*mcount)
        sendcounts(n)=1
        end do
        
       
        ! Get directions in up/down
        call MPI_Cart_shift(MPI_COMM_TwoD,0,1,up,down,ierror)
        ! Get directions in left/right
        call MPI_Cart_shift(MPI_COMM_TwoD,1,1,left,right,ierror)
!        write(*,*)my_rank,"side lr",left,right,"ud", up,down 

        ! Define new type COLUMN
        call MPI_Type_contiguous(mcount+2,MPI_INTEGER,MPI_ONE_COL,ierror)
        call  MPI_Type_commit(MPI_ONE_COL, ierror)
        
        ! Define new type ROW
        call  MPI_Type_vector(ncount+2,1,mcount+2,MPI_INTEGER,MPI_ONE_ROW,ierror)
        call  MPI_Type_commit(MPI_ONE_ROW,ierror)

       call MPI_Type_size(MPI_INTEGER,a,ierror)
        extent=a
        lb=0
        ! Define BLOCK type
       call  MPI_Type_vector(ncount,mcount,LS,MPI_INTEGER,MPI_BLOCKA, ierror)
        call  MPI_Type_create_resized(MPI_BLOCKA,lb,extent,MPI_BLOCK,ierror )
        call  MPI_Type_commit(MPI_BLOCK,ierror)
!define the sending blocks
        call mpi_type_contiguous(2*(mcount+2),mpi_integer,mpi_colo,ierror)
        call mpi_type_commit(mpi_colo,ierror)

        call mpi_type_vector(2*(ncount+2),2,mcount+2,mpi_integer,mpi_rowo,ierror)
!        call mpi_type_create_resized(mpi_rowon,lb,extent,mpi_rowo,ierror)
        call mpi_type_commit(mpi_rowo,ierror)

        write(*,*)'comm_size', comm_size

        allocate(black(comm_size/2))
        allocate(wite(comm_size/2))
!define the chess board        
        i=1
        do,j=0,comm_size-1,2
           if(mod(j/(nint(sqrt(real(comm_size)))),2).eq.0)then
              black(i)=j
              wite(i)=j+1
              i=i+1
           else
              black(i)=j+1
              wite(i)=j
              i=i+1
           endif
        enddo
    
        write(*,*)'black',black
        write(*,*)'wite',wite

        allocate(CA(LS,LS))
        allocate(sCA(0:mcount+1,0:ncount+1))
        allocate(buf(LS*4,LS*4))
        

        time = MPI_Wtime()


        if(my_rank.eq.0)then
        call initCA1(LS, CA, cv, cm) ! this is initialization
        end if
        ca(ls,ls)=0
        ca(1,1)=0

        call  MPI_Scatterv(CA, sendcounts, displs, MPI_BLOCK,&
                CA, 1, MPI_BLOCK, 0, MPI_COMM_TwoD,ierror)
!       
       
        sCA(1:mcount,1:ncount)=CA(1:mcount,1:ncount)

        call  MPI_Sendrecv(sCA(1,0)      ,1,MPI_ONE_ROW, up ,0,&
                sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,&
                MPI_COMM_TwoD,MPI_STATUS_IGNORE,ierror)

        call  MPI_Sendrecv(sCA(mcount,0),1,MPI_ONE_ROW,down,0,&
                sCA(0,0),        1,MPI_ONE_ROW,up  ,0,&
                MPI_COMM_TwoD,MPI_STATUS_IGNORE,ierror)

        call  MPI_Sendrecv(sCA(0,1)      ,1,MPI_ONE_COL,left ,0,&
                sCA(0,ncount+1),1,MPI_ONE_COL,right,0,&
                MPI_COMM_TwoD,MPI_STATUS_IGNORE,ierror)

        call MPI_Sendrecv(sCA(0,ncount),1,MPI_ONE_COL,right,0,&
                sCA(0,0)        ,1,MPI_ONE_COL,left ,0,&
                MPI_COMM_TwoD,MPI_STATUS_IGNORE,ierror)

!c    Time loop
        do it = 0, nt

        CA(1:mcount,1:ncount)=sCA(1:mcount,1:ncount)


           if (MOD(it,10000).eq.0)then

        call  MPI_Gatherv(CA          , 1, MPI_BLOCK, CA, sendcounts,&
                displs, MPI_BLOCK, 0, MPI_COMM_TwoD,ierror)


        if(my_rank.eq.0)then
         
                  nuvac=0
                  i=0
                  j=0
                  a=0
                  do,ii=1,LS
                  do,jj=1,LS
                  if(CA(ii,jj).eq.0)then
                          nuvac=nuvac+1
                  elseif(ca(ii,jj).eq.1)then
                          i=i+1
                  elseif(ca(ii,jj).eq.-1)then
                          j=j+1
                  else
                          a=a+1
                          write(*,*)'kor',ii,jj,ca(ii,jj)
                  endif
                  enddo
                  enddo
              write(*,*)'iter # ', it,nuvac,i,j,a



       buf(:,:)=0
       do, i=1,LS
          do, j=1,LS
             do, row=0, 3
                do, column=0,3
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
       

     ! write(filename,'(A4I0.8A4)')'img/step',it,'.ppm'
     write(filename,'(A4I0.8A4)')'img/step',it,'.ppm'
     call ppmwrite(filename,buf(:,:),4*LS,4*LS)


           end if
        end if


         if (mod(it,2).eq.0)then 

        if(any(black==my_rank))then 

            call yLatex(mcount,ncount,sCA) ! Make one lattice exchange

      call MPI_send(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
!---
      call MPI_send(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)

!------------------

     call MPI_send(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)

!------------------

     call MPI_recv(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)

      else

      call MPI_recv(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
!---
      call MPI_recv(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)

!--------------------
     call MPI_recv(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)

!--------------------

     call MPI_send(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
      endif

         else

        if(any(black==my_rank))then 
            call xLatex(mcount,ncount,sCA) ! Make one lattice exchange

     call MPI_send(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)
!---
     call MPI_send(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
!---------------


      call MPI_send(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
!---------------


      call MPI_recv(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      else

      call MPI_recv(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
!---
     call MPI_recv(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
!------------------


      call MPI_recv(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)


!------------------
      call MPI_send(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)
         endif
      endif

         if (mod(it,2).eq.0)then 

        if(any(wite==my_rank))then 

            call yLatex(mcount,ncount,sCA) ! Make one lattice exchange

      call MPI_send(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
!---
      call MPI_send(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)

!------------------

     call MPI_send(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)

!------------------

     call MPI_recv(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)

      else

      call MPI_recv(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
!---
      call MPI_recv(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)

!--------------------
     call MPI_recv(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)

!--------------------

     call MPI_send(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
      endif

         else

        if(any(wite==my_rank))then 
            call xLatex(mcount,ncount,sCA) ! Make one lattice exchange

     call MPI_send(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,ierror)
!---
     call MPI_send(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
     call MPI_send(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,ierror)
!---------------


      call MPI_send(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
!---------------


      call MPI_recv(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      else

      call MPI_recv(sCA(mcount,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(mcount+1,0),1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
!---
     call MPI_recv(sCA(0,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(sCA(1,0),1,MPI_ONE_ROW,up,0,MPI_COMM_TwoD,&
             MPI_STATUS_IGNORE,ierror)
!------------------


      call MPI_recv(sCA(0,0),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(sCA(0,ncount+1),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,&
              MPI_STATUS_IGNORE,ierror)


!------------------
      call MPI_send(sCA(0,1),1,MPI_ONE_COL,left,0,MPI_COMM_TwoD,ierror)
      call MPI_send(sCA(0,ncount),1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,ierror)
         endif
      endif

        enddo

        
        if(my_rank.eq.0)then
                time = MPI_Wtime() - time

                write(*,*)time 
       
  

         
        endif

        deallocate(CA)
        deallocate(sCA)
        deallocate(displs)
        deallocate(sendcounts)
        deallocate(buf)
        deallocate(black)
        deallocate(wite)


        call MPI_Finalize(ierror)

!c
!c end of executions block
        end
!c-----------------------------------------
        subroutine xLatex(mcount,ncount,sCA) ! Performs the vacancy exchanges
       integer mcount,ncount, LS2, is, ix, iy, il, ir
       integer sCA(0:mcount+1,0:ncount+1)
!       integer ip(LS),im(LS)
      real r 
        integer seed, len, junk,streamnum, nstreams
        integer gtype,stream
       LS2 = mcount * ncount                 ! Total number of attempts = lattice area

!        seed = make_sprng_seed()
!        stream = init_sprng(seed,SPRNG_DEFAULT,1)
!

       !0.857
       call random_number(r)
       if(r.lt.0.857)then
       do is = 1, LS2
!c       generate random coordinates of the current site to be checked
       call random_number(r)
       ix = int(r*mcount)+1
       call random_number(r)
       iy = int(r*ncount)+1
       if (sCA(ix,iy).eq.0) then    ! .. whether there is a vacancy there
                   ! then choose the direction of exchange

       call random_number(r)
          if (r.lt.0.50) then ! exclude right
             il = sCA(ix-1,iy)
             ir = sCA(ix+1,iy)
             if(il.eq.ir.and.il.ne.0)then
                sCA(ix,iy)=il
                sCA(ix+1,iy)=0
             endif
          elseif(r.lt.0.0) then ! exclude down
             il = sCA(ix,iy-1)
             ir = sCA(ix,iy+1)
             if(il.eq.ir.and.il.ne.0)then
             sCA(ix,iy)=il
             sCA(ix,iy+1)=0
             endif
          elseif (r.lt.1.00) then !exclude left
             il = sCA(ix+1,iy)
             ir = sCA(ix-1,iy)
             if(il.eq.ir.and.il.ne.0)then
             sCA(ix,iy)=il
             sCA(ix-1,iy)=0
             endif
          else                     ! exlcude up
             il = sCA(ix,iy+1)
             ir = sCA(ix,iy-1)
             if(il.eq.ir.and.il.ne.0)then
             sCA(ix,iy)=il
             sCA(ix,iy-1)=0        
             endif
          endif
       endif
       enddo

       else
!diffusion
       do is = 1, LS2
!c       generate random coordinates of the current site to be checked
       call random_number(r)
       ix = int(r*mcount)+1
       call random_number(r)
       iy = int(r*ncount)+1
       if (sCA(ix,iy).eq.0) then    ! .. whether there is a vacancy there
                   ! then choose the direction of exchange

       call random_number(r)

          if (r.lt.0.5) then ! exclude right
             il = sCA(ix-1,iy)
             ir = sCA(ix+1,iy)
             sCA(ix,iy)=ir
             sCA(ix+1,iy)=0
          elseif(r.lt.0.00) then ! exclude down
             il = sCA(ix,iy-1)
             ir = sCA(ix,iy+1)
             sCA(ix,iy)=ir
             sCA(ix,iy+1)=0
          elseif (r.lt.1.00) then !exclude left
             il = sCA(ix+1,iy)
             ir = sCA(ix-1,iy)
             sCA(ix,iy)=ir
             sCA(ix-1,iy)=0 
           elseif(r.lt.0.00) then                    ! exlcude up
             il = sCA(ix,iy+1)
             ir = sCA(ix,iy-1)
             sCA(ix,iy)=ir
             sCA(ix,iy-1)=0          
          endif
       endif
       enddo
       endif
       return
       end
!c-----------------------------------------
        subroutine yLatex(mcount,ncount,sCA) ! Performs the vacancy exchanges
       integer mcount,ncount, LS2, is, ix, iy, il, ir
       integer sCA(0:mcount+1,0:ncount+1)
!       integer ip(LS),im(LS)
      real r 
       integer exru(2,-1:1,-1:1)
!#include "sprng_f.h"
!#define SIMPLE_SPRNG	
        integer seed, len, junk,streamnum, nstreams
        integer gtype,stream
       LS2 = mcount * ncount                 ! Total number of attempts = lattice area

!        seed = make_sprng_seed()
!        stream = init_sprng(seed,SPRNG_DEFAULT,1)
!

       !0.857
       call random_number(r)
       if(r.lt.0.857)then
       do is = 1, LS2
!c       generate random coordinates of the current site to be checked
       call random_number(r)
       ix = int(r*mcount)+1
       call random_number(r)
       iy = int(r*ncount)+1
       if (sCA(ix,iy).eq.0) then    ! .. whether there is a vacancy there
                   ! then choose the direction of exchange

       call random_number(r)
          if (r.lt.0.00) then ! exclude right
             il = sCA(ix-1,iy)
             ir = sCA(ix+1,iy)
             if(il.eq.ir.and.il.ne.0)then
                sCA(ix,iy)=il
                sCA(ix+1,iy)=0
             endif
          elseif(r.lt.0.5) then ! exclude down
             il = sCA(ix,iy-1)
             ir = sCA(ix,iy+1)
             if(il.eq.ir.and.il.ne.0)then
             sCA(ix,iy)=il
             sCA(ix,iy+1)=0
             endif
          elseif (r.lt.0.00) then !exclude left
             il = sCA(ix+1,iy)
             ir = sCA(ix-1,iy)
             if(il.eq.ir.and.il.ne.0)then
             sCA(ix,iy)=il
             sCA(ix-1,iy)=0
             endif
          else                     ! exlcude up
             il = sCA(ix,iy+1)
             ir = sCA(ix,iy-1)
             if(il.eq.ir.and.il.ne.0)then
             sCA(ix,iy)=il
             sCA(ix,iy-1)=0        
             endif
          endif
       endif
       enddo

       else
!diffusion
       do is = 1, LS2
!c       generate random coordinates of the current site to be checked
       call random_number(r)
       ix = int(r*mcount)+1
       call random_number(r)
       iy = int(r*ncount)+1
       if (sCA(ix,iy).eq.0) then    ! .. whether there is a vacancy there
                   ! then choose the direction of exchange

       call random_number(r)

          if (r.lt.0.00) then ! exclude right
             il = sCA(ix-1,iy)
             ir = sCA(ix+1,iy)
             sCA(ix,iy)=ir
             sCA(ix+1,iy)=0
          elseif(r.lt.0.50) then ! exclude down
             il = sCA(ix,iy-1)
             ir = sCA(ix,iy+1)
             sCA(ix,iy)=ir
             sCA(ix,iy+1)=0
          elseif (r.lt.0.00) then !exclude left
             il = sCA(ix+1,iy)
             ir = sCA(ix-1,iy)
             sCA(ix,iy)=ir
             sCA(ix-1,iy)=0 
           elseif(r.lt.1.0) then                    ! exlcude up
             il = sCA(ix,iy+1)
             ir = sCA(ix,iy-1)
             sCA(ix,iy)=ir
             sCA(ix,iy-1)=0          
          endif
       endif
       enddo
       endif
       return
       end
!c-----------------------------------------
        subroutine initCA1(LS, CA, cv, cm)
!c      Totally random spread of the three components -1, 0 (vacancy), +1
        integer LS
        integer CA(LS,LS)
        real cv, cm, acm
        real *8 r
        acm = cm + cv
        do ii = 1, LS
           do jj = 1, LS
           call random_number(r)
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
