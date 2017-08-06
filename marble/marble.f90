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
  integer, allocatable, dimension(:,:) :: vaclist, vcl_lft, vcl_rgt, vcl_up, vcl_dwn, vac2move, vcl_still

  integer :: i, j, k, ierror, req, msg_status, nvac, nvacob, rbuf_size, nvac2move, nvac_still
  integer :: nvac_lft, nvac_rgt, nvac_up, nvac_dwn, totnvac
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
  call MPI_Bcast(nequib,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(neout,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
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

  !Check grid decomposition
  if ((mod(nsize,dims(1))+mod(nsize,dims(2))) > 0) then
    if (rank.eq.0) then
      write(0,*)'Grid decomposition is not possible with grid size: ',nsize,' and decomposition: ',nrow,'x',ncol
      write(0,*)'Please try to run your code with different number ot processes'
    end if
    call MPI_Finalize(ierror)
  end if
  !Calculate the local mesh size
  nrow = nsize / dims(1)
  ncol = nsize / dims(2)

  if (rank.eq.0) then
    write(0,*)'Decomposition: ',nrow,'x',ncol,'[',dims(1),',',dims(2),']'
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

  !Calculate possible maximal vacancies and allocate vaclist array
  maxvacnum = nint(nsize * nsize * pbvac * 4e0)
  if (rank.eq.0) then
    write(0,*)'vac size list: ', maxvacnum
  end if

  allocate(vaclist(4,maxvacnum))
  allocate(vcl_lft(4,maxvacnum))
  allocate(vcl_rgt(4,maxvacnum))
  allocate(vcl_up(4,maxvacnum))
  allocate(vcl_dwn(4,maxvacnum))
  allocate(vac2move(4,maxvacnum))
  allocate(vcl_still(4,maxvacnum))
  allocate(cbuf(4,maxvacnum))
  allocate(rbuf(4,maxvacnum))

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

  !Allocate memory
  allocate(mesh(0:nrow+1, 0:ncol+1))
  allocate(wbuf(nrow, ncol))

  !seed the random number generator
  call MC_init_random_seed(rank)
  !Initialize the grid
  call MC_init_mesh(mesh, vaclist, nvac, pbvac, pbratio)

  !Main iteration loop
  startTime = MPI_Wtime()

  !File writes counter initialization
  cnt = 0
  do it = 0, nequib+nit
    !Generate velocities
    call MC_genvel(mesh, vaclist, nvac, vcl_lft, nvac_lft, vcl_rgt, nvac_rgt,&
     vcl_up, nvac_up, vcl_dwn, nvac_dwn, vcl_still, nvac_still)

    call print_vac_list(vaclist, nvac,         'begin', mesh)
    call print_vac_list(vcl_up, nvac_up,       'up',    mesh)
    call print_vac_list(vcl_dwn, nvac_dwn,     'down',  mesh)
    call print_vac_list(vcl_lft, nvac_lft,     'left',  mesh)
    call print_vac_list(vcl_rgt, nvac_rgt,     'right', mesh)
    call print_vac_list(vcl_still, nvac_still, 'still', mesh)

    !copy boundaries
    write(0,'(A12)')'-----bb-----'
    do i = 0, nrow+1
      write(0,'(I4,A1,I2,A1,I3,A1,<ncol+2>I1)')cnt,':',rank,':',i,':',(mesh(i,j),j=0,ncol+1)
    end do
    write(0,'(A12)')'-----be-----'

    call MC_select_boundary_vac(vcl_up, nvac_up, 'up', cbuf, nvacob)
    !write(0, *)'it:',it,rank,':[up]',up,'->[down]',down
    call MPI_Sendrecv(mesh(1, 0), 1, MPI_TYPE_ROW, up, 0,&
      mesh(nrow+1,0),1,MPI_TYPE_ROW,down,0,MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Isend(cbuf, 4 * nvacob, MPI_INTEGER, up, 0, MPI_COMM_2D, req, ierror)
    call MPI_Probe(down, MPI_ANY_TAG, MPI_COMM_2D, msg_status, ierror)
    call MPI_Get_count(msg_status, MPI_INTEGER, rbuf_size, ierror)
    call MPI_Recv(rbuf, rbuf_size, MPI_INTEGER, down, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Wait(req, MPI_STATUS_IGNORE, ierror)
    if (rbuf_size > 0) then
!      call print_vac_list(   vcl_up, nvac_up, 'v_d_bef', mesh)
      call MC_update_vaclist(vcl_up, nvac_up, rbuf, rbuf_size/4)
!      call print_vac_list(   vcl_up, nvac_up, 'v_d_aft', mesh)
    end if
    nvac2move = 0
    call MC_join_vac(vac2move, nvac2move, vcl_up,  nvac_up)
    !Sort vacancy list
    call MC_sort_vclist(vac2move, nvac2move, 1)
    nvac = 0
    call MC_Step(mesh, vac2move, nvac2move, vaclist, nvac)


    call MC_select_boundary_vac(vcl_dwn, nvac_dwn, 'down', cbuf, nvacob)
    !write(0, *)'it:',it,rank,':[down]',down,'->[up]',up
    call MPI_Sendrecv(mesh(nrow, 0), 1, MPI_TYPE_ROW, down, 0,&
      mesh(0, 0), 1, MPI_TYPE_ROW, up, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Isend(cbuf, 4 * nvacob, MPI_INTEGER, down, 0, MPI_COMM_2D, req, ierror)
    call MPI_Probe(up, MPI_ANY_TAG, MPI_COMM_2D, msg_status, ierror)
    call MPI_Get_count(msg_status, MPI_INTEGER, rbuf_size, ierror)
    call MPI_Recv(rbuf, rbuf_size, MPI_INTEGER, up, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Wait(req, MPI_STATUS_IGNORE, ierror)
    if (rbuf_size > 0) then
      call print_vac_list(   vcl_dwn, nvac_dwn, 'v_d_bef', mesh)
      call MC_update_vaclist(vcl_dwn, nvac_dwn, rbuf, rbuf_size/4)
      call print_vac_list(   vcl_dwn, nvac_dwn, 'v_d_aft', mesh)
    end if
    nvac2move = 0
    call MC_join_vac(vac2move, nvac2move, vcl_dwn,  nvac_dwn)
    !Sort vacancy list
    call MC_sort_vclist(vac2move, nvac2move, 1)
    call MC_Step(mesh, vac2move, nvac2move, vaclist, nvac)
    write(0,'(A12)')'-----ab-----'
    do i = 0, nrow+1
      write(0,'(I4,A1,I2,A1,I3,A1,<ncol+2>I1)')cnt,':',rank,':',i,':',(mesh(i,j),j=0,ncol+1)
    end do
    write(0,'(A12)')'-----ae-----'



    call MC_select_boundary_vac(vcl_lft, nvac_lft, 'left', cbuf, nvacob)
    !write(0, *)'it:',it,rank,':[left]',left,'->[right]',right
    call MPI_Sendrecv(mesh(0, 1), 1, MPI_TYPE_COL, left , 0,&
      mesh(0, ncol+1), 1, MPI_TYPE_COL, right, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Isend(cbuf, 4 * nvacob, MPI_INTEGER, left, 0, MPI_COMM_2D, req, ierror)
    call MPI_Probe(right, MPI_ANY_TAG, MPI_COMM_2D, msg_status, ierror)
    call MPI_Get_count(msg_status, MPI_INTEGER, rbuf_size, ierror)
    call MPI_Recv(rbuf, rbuf_size, MPI_INTEGER, right, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Wait(req, MPI_STATUS_IGNORE, ierror)
    if (rbuf_size > 0) then
      !call print_vac_list(   rbuf, rbuf_size/4, 'b lft', mesh)
      !call print_vac_list(   vcl_lft, nvac_lft, 'vcl_lft', mesh)
      call MC_update_vaclist(vcl_lft, nvac_lft, rbuf, rbuf_size/4)
      !call print_vac_list(   vcl_lft, nvac_lft, 'vcl_lft', mesh)
    end if
    nvac2move = 0
    call MC_join_vac(vac2move, nvac2move, vcl_lft,  nvac_lft)
    !Sort vacancy list
    call MC_sort_vclist(vac2move, nvac2move, 1)
    call MC_Step(mesh, vac2move, nvac2move, vaclist, nvac)

    call MC_select_boundary_vac(vcl_rgt, nvac_rgt, 'right', cbuf, nvacob)
    !write(0, *)'it:',it,rank,':[right]',right,'->[left]',left
    call MPI_Sendrecv(mesh(0, ncol), 1, MPI_TYPE_COL, right, 0,&
      mesh(0, 0), 1, MPI_TYPE_COL, left , 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Isend(cbuf, 4 * nvacob, MPI_INTEGER, right, 0, MPI_COMM_2D, req, ierror)
    call MPI_Probe(left, MPI_ANY_TAG, MPI_COMM_2D, msg_status, ierror)
    call MPI_Get_count(msg_status, MPI_INTEGER, rbuf_size, ierror)
    call MPI_Recv(rbuf, rbuf_size, MPI_INTEGER, left, 0, MPI_COMM_2D,&
      MPI_STATUS_IGNORE, ierror)
    call MPI_Wait(req, MPI_STATUS_IGNORE, ierror)
    if (rbuf_size > 0) then
      !call print_vac_list(   rbuf, rbuf_size/4, 'b rgt', mesh)
      call MC_update_vaclist(vcl_rgt, nvac_rgt, rbuf, rbuf_size/4)
    end if
    nvac2move = 0
    call MC_join_vac(vac2move, nvac2move, vcl_rgt,  nvac_rgt)
    !Sort vacancy list
    call MC_sort_vclist(vac2move, nvac2move, 1)
    call MC_Step(mesh, vac2move, nvac2move, vaclist, nvac)


!    write(0,'(A12)')'------------'
!    do i = 0, nrow+1
!      write(0,'(I4,A1,I2,A1,I3,A1,<ncol+2>I1)')cnt,':',rank,':',i,':',(mesh(i,j),j=0,ncol+1)
!    end do
!    write(0,'(A12)')'------------'

!    nvac2move = 0
!    call MC_join_vac(vac2move, nvac2move, vcl_up,  nvac_up)
!    call MC_join_vac(vac2move, nvac2move, vcl_dwn, nvac_dwn)
!    !Sort vacancy list
!    call MC_sort_vclist(vac2move, nvac2move, 1)
!    !call print_vac_list(vac2move, nvac2move, 'u<->d', mesh)
!    nvac = 0
!    call MC_Step(mesh, vac2move, nvac2move, vaclist, nvac)


!    nvac2move = 0
!    call MC_join_vac(vac2move, nvac2move, vcl_lft, nvac_lft)
!    call MC_join_vac(vac2move, nvac2move, vcl_rgt, nvac_rgt)
!    !Sort vacancy list
!    call MC_sort_vclist(vac2move, nvac2move, 1)
!
!    call MC_Step(mesh, vac2move, nvac2move, vaclist, nvac)

    !write(0,'(A12)')'------------'
    !do i = 0, nrow+1
    !  write(0,'(I4,A1,I2,A1,I3,A1,<ncol+2>I1)')cnt,':',rank,':',i,':',(mesh(i,j),j=0,ncol+1)
    !end do
    !write(0,'(A12)')'------------'

    call MC_join_vac(vaclist, nvac, vcl_still, nvac_still)

    call print_vac_list(vaclist, nvac, 'end  ', mesh)

    if (((mod(it-nequib, nout).eq.0).and.(it.ge.nequib)).or.(((mod(it, neout).eq.0).and.(it.lt.nequib)))) then
      call MPI_Reduce(nvac, totnvac, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_2D, ierror)
      if (rank.eq.0) then
        if (it.lt.nequib) then
          write(6,'(A20,I10,A8,I10)')'equlibration step = ',it,'nvac = ',totnvac
        else
          write(6,'(A20,I10,A8,I10)')'production   step = ',it - nequib,'nvac = ',totnvac
        end if
      end if
      if (it.gt.nequib) then
        write(filename, write_fmt)chkp_filename,'-',cnt
        cnt = cnt + 1
        wbuf(1:nrow,1:ncol) = mesh(1:nrow, 1:ncol)
        call write_out_mesh(wbuf, filename)

        !call print_vac_list(vaclist, nvac, 'end  ', mesh)
      end if
    end if
  end do
  endTime = MPI_Wtime()
  if (rank.eq.0) then
    write(0,'(A16,F10.2)'),'Execution time: ',endTime-startTime
  end if
  !Deallocate the memory
  deallocate(mesh)
  deallocate(wbuf)
  deallocate(cbuf)
  deallocate(rbuf)
  deallocate(scounts)
  deallocate(displs)
  deallocate(vaclist)
  deallocate(vcl_lft)
  deallocate(vcl_rgt)
  deallocate(vcl_up)
  deallocate(vcl_dwn)
  deallocate(vac2move)
  deallocate(vcl_still)

  !Finish the job
  call MPI_Type_free(MPI_TYPE_COL, ierror)
  call MPI_Type_free(MPI_TYPE_ROW, ierror)
  call MPI_Type_free(MPI_BLOCK2, ierror)
  call MPI_Type_free(MPI_BLOCK, ierror)
  call MPI_Finalize(ierror)
end program marble
