program dump
use globvars
use utils
implicit none

  integer :: i, j
  integer :: fh, cval,thresh
  integer :: cr, cg, cb


  thresh = 255
  fh = 15
  call read_input_file(inpt_filename)
  open(fh, file=chkp_filename, form='binary')
  write(6,'(A2)')'P3'
  write(6,'(2I9)')nsize,nsize
  write(6,'(I3)')thresh
  do j = 1, nsize
    do i = 1, nsize
      read(fh)cval
      if (cval.eq.type_vac) then
        cr = 255
        cg = 255
        cb = 255
      elseif (cval.eq.type_red) then
        cr = 255
        cg = 0
        cb = 0
      elseif (cval.eq.type_blue) then
        cr = 0
        cg = 0
        cb = 255
      else
        cr = 0
        cg = 0
        cb = 0
      end if
      write(6,'(3I4)')cr,cg,cb
    end do
  end do
  close(fh)
end program dump
