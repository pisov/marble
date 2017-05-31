module utils
  use globvars
  implicit none


  contains
  subroutine write_out_mesh(mesh)
    integer, dimension(:,:), intent(inout) :: mesh

    integer :: ierror
    integer :: file_ptr
    integer :: bufsize

    !write(0,*)'--->',mesh(5, 11)
    bufsize = ncol * nrow

    call MPI_File_open(MPI_COMM_2D, chkp_filename, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, file_ptr, ierror)

    call MPI_File_set_view(file_ptr, displs(rank+1), MPI_INTEGER, &
                           MPI_BLOCK2, 'native', &
                           MPI_INFO_NULL, ierror)
    call MPI_File_write(file_ptr, mesh, bufsize, MPI_INTEGER, &
                        MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(file_ptr, ierror)

  end subroutine write_out_mesh

  subroutine read_input_file(filename)
    character(len=64), intent(in) :: filename

    !Control file read variables
    character(len=128) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios = 0
    integer :: line = 0

    open(fh, file=filename)
    do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '=')
        label = trim(buffer(1:pos-1))
        buffer = trim(buffer(pos+1:))

        select case (trim(label))
        case ('nit')
           read(buffer, *, iostat=ios) nit
        case ('nout')
           read(buffer, *, iostat=ios) nout
        case ('Temp')
            read(buffer, *, iostat=ios) Temp
        case ('pdvac')
            read(buffer, *, iostat=ios) pdvac
        case ('pdratio')
            read(buffer, *, iostat=ios) pdratio
        case ('nsize')
            read(buffer, *, iostat=ios) nsize
        case default
          !write(0,*)'Skip line:',buffer
        end select
     end if
  end do
  end subroutine read_input_file

  subroutine sort_vclist(vclist, nvac, row_id)
    implicit none
    integer, dimension(:,:), intent(inout) :: vclist
    integer, intent(in) :: nvac, row_id

    integer :: idx, i
    integer, dimension(4) :: buf

    do i = 1,  nvac
      idx = minloc(vclist(row_id,i:nvac),dim=1) + i - 1
      buf(:) = vclist(:, idx)
      vclist(:, idx) = vclist(:, i)
      vclist(:, i) = buf(:)
    end do
  end subroutine sort_vclist

  subroutine select_boundary_vac(vclist, nvac, direction, buf, nvacob)
    implicit none
    integer, dimension(:,:), intent(in) :: vclist
    integer, intent(in) :: nvac
    character(len=*), intent(in) :: direction
    integer, dimension(:,:), intent(out) :: buf
    integer, intent(out) :: nvacob

    integer :: i, id, side_id

    if (trim(direction).eq.'left') then
      id = 1
      side_id = 1
    else if (trim(direction).eq.'right') then
      id = 1
      side_id = ncol
    else if (trim(direction).eq.'up') then
      id = 2
      side_id = 1
    else if (trim(direction).eq.'down') then
      id = 2
      side_id = nrow
    else
      id = 1
      side_id = 1
    end if

    nvacob = 0
    do i = 1, nvac
      if (vclist(id,i).eq.side_id) then
        nvacob = nvacob + 1
        buf(i:, nvacob) = vclist(:,i)
      end if
    end do
    
  end subroutine select_boundary_vac
end module utils
