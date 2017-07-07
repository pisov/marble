module utils
  use globvars
  implicit none


  contains
  subroutine write_out_mesh(mesh, filename)
    integer, dimension(:,:), intent(inout) :: mesh
    character(len=*), intent(in) :: filename

    integer :: ierror
    integer :: file_ptr
    integer :: bufsize


    !write(0,*)'--->',mesh(5, 11)
    bufsize = ncol * nrow

    call MPI_File_open(MPI_COMM_2D, trim(filename), &
                        MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                        MPI_INFO_NULL, file_ptr, ierror)

    call MPI_File_set_view(file_ptr, displs(rank+1), MPI_INTEGER, &
                           MPI_BLOCK2, 'native', &
                           MPI_INFO_NULL, ierror)
    call MPI_File_write(file_ptr, mesh, bufsize, MPI_INTEGER, &
                        MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(file_ptr, ierror)

  end subroutine write_out_mesh

  subroutine check_out_file()
    implicit none
    logical :: res

    inquire( file=trim(chkp_filename), exist=res )

    if (res.eq..true.) then
      !
    end if

  end subroutine check_out_file

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
        case ('pbvac')
            read(buffer, *, iostat=ios) pbvac
        case ('pbratio')
            read(buffer, *, iostat=ios) pbratio
        case ('pbfrmv')
            read(buffer, *, iostat=ios) pbfrmv
        case ('nsize')
            read(buffer, *, iostat=ios) nsize
        case default
          !write(0,*)'Skip line:',buffer
        end select
     end if
  end do
  end subroutine read_input_file

  subroutine print_vac_list(vaclist, nvac, listname, mesh)
    implicit none
    integer, dimension(4,maxvacnum), intent(in) :: vaclist
    integer, intent(in) :: nvac
    character(len=*) :: listname
    integer, dimension(0:nrow+1, 0:ncol+1) :: mesh

    integer :: i

    do i = 1, nvac
        write(0,'(I3,A3,A1,A5,A10,I3,A1,I3,A2,I3,A2,I3,A4,I1)')rank,' : ','[',listname,'] step: ',it,':',i,' (',vaclist(1, i),', ',vaclist(2, i),') : ',mesh(vaclist(1, i), vaclist(2, i))
    end do

  end subroutine print_vac_list

end module utils
