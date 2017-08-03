program marble
  use globvars
  use mc
  use utils
  implicit none
  !use mpi
  !include 'mpif.h'


  !Define variables
  !Interaction matrix where energy between each sample is stored
  !blue-blue, blue-red, red-red, red-vacancy, blue-cavancy, vacancy-vacancy
  double precision, dimension(3) :: Jmtx
  !Grid variables
  integer, allocatable, dimension(:,:) :: mesh, cbuf, rbuf, wbuf
  !Vacancy list
  integer, allocatable, dimension(:,:) :: vac
  
  !chess grid
  integer,dimension(:),allocatable :: black,wite

  integer :: i, j,jj, k, ierror, req, msg_status, nvac, nvacob, rbuf_size
  integer :: totnvac
  integer(kind=8) :: cnt
  character(len=32) :: write_fmt, filename

  double precision :: startTime, endTime

  !Initialize the MPI environment
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, commsize, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

  !Read the config file
  if (rank.eq.0) then
    call read_input_file(inpt_filename)
  end if

  !Set write format
  write(write_fmt,'(A2,I0,A9)')'(A',len_trim(chkp_filename),',A1,I0.6)'

  !Distribute the initial variables
  call MPI_Bcast(nit,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(nout,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(nsize,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(Temp,    1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(pbvac,   1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(pbratio, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(pbfrmv,  1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)

  !Create Cartesian 2D topology
  dims(:) = 0
  call MPI_Dims_create(commsize, 2, dims, ierror)
  periods(:) = .true.
  call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., MPI_COMM_2D, ierror)
  call MPI_Comm_rank(MPI_COMM_2D, rank, ierror)
  !Get directions in up/down
  call MPI_Cart_shift(MPI_COMM_2D, 0, 1, up, down, ierror)
  !Get directions in left/right
  call MPI_Cart_shift(MPI_COMM_2D, 1, 1,left,right, ierror)
  ! write(*,*)rank,'ud',up,down,'l/r',left,right
  !Check grid decomposition
!  if ((mod(nsize,dims(1))+mod(nsize,dims(2))) > 0) then
!    if (rank.eq.0) then
!      write(0,*)'Grid decomposition is not possible with grid size: ',nsize,' and decomposition: ',nrow,'x',ncol
!      write(0,*)'Please try to run your code with different number ot processes'
!    end if
!    call MPI_Finalize(ierror)
!    call exit
!  end if
  !Calculate the local mesh size
  nrow = nsize / dims(1)
  ncol = nsize / dims(2)

  if (rank.eq.0) then
!    write(0,*)'Decomposition: ',nrow,'x',ncol,'[',dims(1),',',dims(2),']'
  end if

  !Create MPI types
  !Define new type COLUMN
   call MPI_Type_contiguous(nrow+2, MPI_INTEGER, MPI_TYPE_COL, ierror)
   call MPI_Type_commit(MPI_TYPE_COL, ierror)
  !Define new type ROW
   call MPI_Type_vector(ncol+2, 1, nrow+2, MPI_INTEGER, MPI_TYPE_ROW, ierror)
   call MPI_Type_commit(MPI_TYPE_ROW, ierror)
   !Define BLOCK type
   call MPI_Type_vector(ncol, nrow, nsize, MPI_INTEGER, MPI_BLOCK2, ierror)
   call MPI_Type_commit(MPI_BLOCK2, ierror)
   call MPI_Type_create_resized(MPI_BLOCK2, 0, 4, MPI_BLOCK, ierror)
   call MPI_Type_commit(MPI_BLOCK, ierror)

  !Allocate displs and scounts arrays give dims
  allocate(scounts(dims(1)*dims(2)))
  allocate(displs(dims(1)*dims(2)))
  allocate(black(commsize/2))
  allocate(wite(commsize/2))
  !Calculate possible maximal vacancies and allocate vaclist array


  scounts(:) = 1
  do k = 1, commsize
    call MPI_Cart_coords(MPI_COMM_2D, k - 1, 2, crd, ierror)
    i = crd(1) ! nrow
    j = crd(2) ! ncol
    !displs(k) = (i* nsize * nrow + j* ncol)*4
    displs(k) = (i * nrow + j * ncol * nsize)*4
    !if (rank.eq.0) then
    !  write(0,*)k-1, i, j, displs(k)
    !end if
  end do
  call MPI_Cart_coords(MPI_COMM_2D, rank, 2, crd, ierror)
 jj=1 
       do i=0, dims(1)-1
          do j=0, dims(2)-1,2
           
           if(mod(i,2).eq.0)then
                   black(jj)=j+dims(2)*i
                   wite(jj)=j+dims(2)*i+1
                   jj=jj+1
           else
                   black(jj)=j+dims(2)*i+1
                   wite(jj)=j+dims(2)*i
                   jj=jj+1
           endif

          enddo
       enddo

!       write(*,*)'black',black
!       write(*,*)'wite',wite
        
        if(mod(dims(1),2)+mod(dims(2),2).ne.0.or.dims(1)*dims(2).ne.commsize)then

                if(rank.eq.0)then
                   write(*,*)"please choose another number of processes"
                endif

                call mpi_finalize(ierror)
               call exit
        endif


  !Allocate memory
  allocate(mesh(0:nrow+1, 0:ncol+1))
  allocate(wbuf(nrow, ncol))

  !seed the random number generator
  call MC_init_random_seed(rank)
  !Initialize the grid
  call MC_init_mesh(mesh,  nvac, pbvac, pbratio,rank)
  
  call MPI_Allreduce(nvac, totnvac, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_2D, ierror)
  allocate(vac(totnvac,2))
  !Main iteration loop
  startTime = MPI_Wtime()

  !File writes counter initialization
!!!!!!!!!!!!!
!!!!!!x`

   call  MPI_Sendrecv(mesh(1,0)      ,1,MPI_TYPE_ROW, up ,0,&
        mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,&
        MPI_COMM_2D,MPI_STATUS_IGNORE,ierror)

   call  MPI_Sendrecv(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,&
        mesh(0,0),        1,MPI_TYPE_ROW,up  ,0,&
        MPI_COMM_2D,MPI_STATUS_IGNORE,ierror)

   call  MPI_Sendrecv(mesh(0,1),1,MPI_TYPE_COL,left ,0,&
        mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,&
        MPI_COMM_2D,MPI_STATUS_IGNORE,ierror)

   call MPI_Sendrecv(mesh(0,ncol),1,MPI_TYPE_COL,right,0,&
        mesh(0,0),1,MPI_TYPE_COL,left,0,&
        MPI_COMM_2D,MPI_STATUS_IGNORE,ierror)

  cnt = 0
  do it = 0, nit

         if (mod(it,2).eq.0)then 

        if(any(black==rank))then 

      call updatevaclist (nrow,ncol, mesh,totnvac, vac )
            call yLatex(nrow,ncol,mesh, vac,totnvac ) ! Make one lattice exchange
!     write(*,*)nvac 
      call MPI_send(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
!---
      call MPI_send(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)

!------------------

     call MPI_send(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)

!------------------

     call MPI_recv(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)

      else

      call MPI_recv(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
!---
      call MPI_recv(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)

!--------------------
     call MPI_recv(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)

!--------------------

     call MPI_send(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
      endif

         else

        if(any(black==rank))then 
     call updatevaclist (nrow,ncol, mesh,totnvac, vac )
            call xLatex(nrow,ncol,mesh, vac, totnvac ) ! Make one lattice exchange

     call MPI_send(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)
!---
     call MPI_send(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
!---------------


      call MPI_send(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
!---------------


      call MPI_recv(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      else

      call MPI_recv(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
!---
     call MPI_recv(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
!------------------


      call MPI_recv(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)


!------------------
      call MPI_send(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)
         endif
      endif

         if (mod(it,2).eq.0)then 

        if(any(wite==rank))then 

      call updatevaclist (nrow,ncol, mesh,totnvac, vac )
            call yLatex(nrow,ncol,mesh, vac, totnvac ) ! Make one lattice exchange

      call MPI_send(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
!---
      call MPI_send(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)

!------------------

     call MPI_send(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)

!------------------

     call MPI_recv(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)

      else

      call MPI_recv(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
!---
      call MPI_recv(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)

!--------------------
     call MPI_recv(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)

!--------------------

     call MPI_send(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
      endif

         else

        if(any(wite==rank))then 

     call updatevaclist (nrow,ncol, mesh,totnvac, vac )
     call xLatex(nrow,ncol,mesh, vac, totnvac ) ! Make one lattice exchange

     call MPI_send(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,ierror)
!---
     call MPI_send(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
     call MPI_send(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,ierror)
!---------------


      call MPI_send(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
!---------------


      call MPI_recv(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      else

      call MPI_recv(mesh(nrow,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
!---
     call MPI_recv(mesh(0,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
     call MPI_recv(mesh(1,0),1,MPI_TYPE_ROW,up,0,MPI_COMM_2D,&
             MPI_STATUS_IGNORE,ierror)
!------------------


      call MPI_recv(mesh(0,0),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)
      call MPI_recv(mesh(0,ncol+1),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,&
              MPI_STATUS_IGNORE,ierror)


!------------------
      call MPI_send(mesh(0,1),1,MPI_TYPE_COL,left,0,MPI_COMM_2D,ierror)
      call MPI_send(mesh(0,ncol),1,MPI_TYPE_COL,right,0,MPI_COMM_2D,ierror)
         endif
      endif


!!!!!!!
    if (mod(it, nout).eq.0) then
      call MPI_Reduce(nvac, totnvac, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_2D, ierror)
      if (rank.eq.0) then
!       write(6,'(A7,I10,A8,I10)')' step = ',it,'nvac = ',totnvac
      end if
      !wbuf(1:nrow,1:ncol) = abs(mod(crd(1), 2)-mod(crd(2), 2))+1
      write(filename, write_fmt)chkp_filename,'-',cnt
      cnt = cnt + 1
      wbuf(1:nrow,1:ncol) = mesh(1:nrow, 1:ncol)
      call write_out_mesh(wbuf, filename)
    end if
  end do
  endTime = MPI_Wtime()
  if (rank.eq.0) then
    !write(0,'(A16,F10.2)'),'Execution time: ',endTime-startTime
    write(*,*)commsize,endTime-startTime
  end if
  !Deallocate the memory
  deallocate(mesh)
  deallocate(wbuf)
  deallocate(scounts)
  deallocate(displs)
  deallocate(vac)
  deallocate(black)
  deallocate(wite)

  !Finish the job
  call MPI_Type_free(MPI_TYPE_COL, ierror)
  call MPI_Type_free(MPI_TYPE_ROW, ierror)
  call MPI_Type_free(MPI_BLOCK2, ierror)
  call MPI_Type_free(MPI_BLOCK, ierror)
  call MPI_Finalize(ierror)
end program marble
