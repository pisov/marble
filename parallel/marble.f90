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
  integer, allocatable, dimension(:,:) :: mesh, nmesh, cbuf, rbuf
  !Vacancy list
  integer, allocatable, dimension(:,:) :: vaclist, vcl_lft, vcl_rgt, vcl_up, vcl_dwn

  integer :: i, j, k, ierror, req, msg_status, maxvacnum, nvac, nvacob, rbuf_size
  integer :: nvac_lft, nvac_rgt, nvac_up, nvac_dwn

  !Initialize the MPI environment
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, commsize, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

  !Read the config file
  if (rank.eq.0) then
    call read_input_file(inpt_filename)
  end if

  !Distribute the initial variables
  call MPI_Bcast(nit,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(nout,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(nsize,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(Temp,    1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(pdvac,   1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(pdratio, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierror)

  !Create Cartesian 2D topology
  dims(:) = 0
  call MPI_Dims_create(commsize, 2, dims, ierror)
  periods(:) = .true.
  call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, MPI_COMM_2D, ierror)
  call MPI_Comm_rank(MPI_COMM_2D, rank, ierror)
  !Get directions in up/down
  call MPI_Cart_shift(MPI_COMM_2D, 0, 1, up, down, ierror)
  !Get directions in left/right
  call MPI_Cart_shift(MPI_COMM_2D, 1, 1,left,right, ierror)

  !Calculate the local mesh size
  nrow = nsize / dims(1)
  ncol = nsize / dims(2)

  !Create MPI types
  !Define new type COLUMN
   call MPI_Type_contiguous(ncol+2, MPI_INTEGER, MPI_TYPE_COL, ierror)
   call MPI_Type_commit(MPI_TYPE_COL, ierror)
  !Define new type ROW
   call MPI_Type_vector(nrow+2, 1, ncol+2, MPI_INTEGER, MPI_TYPE_ROW, ierror)
   call MPI_Type_commit(MPI_TYPE_ROW, ierror)
   !Define BLOCK type
   call MPI_Type_vector(nrow, ncol, nsize, MPI_INTEGER, MPI_BLOCK2, ierror)
   call MPI_Type_commit(MPI_BLOCK2, ierror)
   call MPI_Type_create_resized(MPI_BLOCK2, 0, 4, MPI_BLOCK, ierror)
   call MPI_Type_commit(MPI_BLOCK, ierror)

  !Allocate displs and scounts arrays give dims
  allocate(scounts(dims(1)*dims(2)))
  allocate(displs(dims(1)*dims(2)))

  !Calculate possible maximal vacancies and allocate vaclist array
  maxvacnum = nint(nsize * nsize * pdvac * 1.2e0)

  allocate(vaclist(4,maxvacnum))
  allocate(vcl_lft(4,maxvacnum))
  allocate(vcl_rgt(4,maxvacnum))
  allocate(vcl_up(4,maxvacnum))
  allocate(vcl_dwn(4,maxvacnum))
  allocate(cbuf(4,maxvacnum))
  allocate(rbuf(4,maxvacnum))

  scounts(:) = 1
  do k = 1, commsize
    call MPI_Cart_coords(MPI_COMM_2D, k - 1, 2, crd, ierror)
    i = crd(1)
    j = crd(2)
    displs(k) = (i* nsize * nrow + j* ncol)*4
  end do

  !Allocate memory
  allocate(mesh(0:nrow+1, 0:ncol+1))
  allocate(nmesh(0:nrow+1, 0:ncol+1))

  !Initialize the grid
  call MC_init_mesh(mesh(1:nrow, 1:ncol), vaclist, nvac, pdvac, pdratio)

  !Main iteration loop
  do it = 1, nit
    !Generate velocities
    call MC_genvel(vaclist, nvac, vcl_lft, nvac_lft, vcl_rgt, nvac_rgt,&
     vcl_up, nvac_up, vcl_dwn, nvac_dwn)

    !copy boundaries
    call select_boundary_vac(vcl_up, nvac_up, 'up', cbuf, nvacob)
    call MPI_Sendrecv(mesh(1, 0), 1, MPI_TYPE_ROW, up, 0,&
      mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Isend(cbuf, 4 * nvacob, MPI_INTEGER, up, 0, MPI_COMM_2D, req, ierror)
    call MPI_Probe(down, MPI_ANY_TAG, MPI_COMM_2D, msg_status, ierror)
    call MPI_Get_count(msg_status, MPI_INTEGER, rbuf_size, ierror)
    call MPI_Recv(rbuf, rbuf_size, MPI_INTEGER, down, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    if (rbuf_size > 0) then
      call update_vaclist(vaclist, nvac, rbuf, rbuf_size/4)
    end if
    !make one MC step up
    call MC_Step(mesh(1:nrow, 1:ncol), vaclist, nvac, 'up')

    !call select_boundary_vac(vcl_dwn, nvac_dwn, 'down', cbuf, nvacob)
    !call MPI_Sendrecv(mesh(nrow, 0), 1, MPI_TYPE_ROW, down, 0,&
    !  mesh(0, 0), 1, MPI_TYPE_ROW, up, 0, MPI_COMM_2D,&
    !  MPI_STATUS_IGNORE, ierror)
    !!make one MC step up
    !call MC_Step(mesh(1:nrow, 1:ncol), vaclist, nvac, 'up')


!    call MC_Step(mesh(1:nrow, 1:ncol), vaclist, nvac, 'left')

!    call MPI_Sendrecv(mesh(0, 1), 1, MPI_TYPE_COL, left , 0,&
!      mesh(0, ncol+1), 1, MPI_TYPE_COL, right, 0, MPI_COMM_2D,&
!      MPI_STATUS_IGNORE, ierror)
!    call MPI_Sendrecv(mesh(0, ncol), 1, MPI_TYPE_COL, right, 0,&
!      mesh(0, 0), 1, MPI_TYPE_COL, left , 0, MPI_COMM_2D,&
!      MPI_STATUS_IGNORE, ierror)


    if (mod(it, nout).eq.0) then
      if (rank.eq.0) then
        write(6,'(A7,I10)')'step = ',it
      end if
      call write_out_mesh(mesh(1:nrow, 1:ncol))
    end if
  end do
  !Deallocate the memory
  deallocate(mesh)
  deallocate(nmesh)
  deallocate(cbuf)
  deallocate(rbuf)
  deallocate(scounts)
  deallocate(displs)
  deallocate(vaclist)
  deallocate(vcl_lft)
  deallocate(vcl_rgt)
  deallocate(vcl_up)
  deallocate(vcl_dwn)

  !Finish the job
  call MPI_Type_free(MPI_TYPE_COL, ierror)
  call MPI_Type_free(MPI_TYPE_ROW, ierror)
  call MPI_Type_free(MPI_BLOCK2, ierror)
  call MPI_Type_free(MPI_BLOCK, ierror)
  call MPI_Finalize(ierror)
end program marble
