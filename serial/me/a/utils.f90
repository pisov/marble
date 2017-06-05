module utils
implicit none
contains
subroutine pgmwrite(filename, x, nx, ny)

  implicit none

  character*(*) :: filename
  integer :: nx, ny

  real,    dimension(nx, ny) :: x

  real,    dimension(nx, ny) :: tmp
  integer, dimension(nx, ny) :: grey

  real :: tmin, tmax
  real, parameter :: thresh = 255.0

  integer, parameter :: iounit = 10

  integer :: i, j

  tmp(:,:) = x(:,:)

!  Find the max and min absolute values of the array

  tmin = minval(abs(tmp(:,:)))
  tmax = maxval(abs(tmp(:,:)))

!  Scale the values appropriately so the lies between 0 and thresh

  if (tmin .lt. 0 .or. tmax .gt. thresh) then

    tmp(:,:) = int((thresh*((abs(tmp(:,:)-tmin))/(tmax-tmin))) + 0.5)

  else

    tmp(:,:) = int(abs(tmp(:,:)) + 0.5)

  end if

!  Increase the contrast by boosting the lower values

  grey(:,:) = thresh * sqrt(tmp(:,:)/thresh)

  open(unit=iounit, file=filename)

  write(iounit,fmt='(''P2''/''# Written by pgmwrite'')')
  write(iounit,*) nx, ny
  write(iounit,*) int(thresh)
  write(iounit,fmt='(16(i3,'' ''))') ((grey(i,j), i=1,nx), j=1,ny)

  close(unit=iounit)

end subroutine pgmwrite

subroutine ppmwrite(filename, buf, nx, ny)
  implicit none
  integer, intent(in) :: nx, ny
  double precision, dimension(nx, ny), intent(in) :: buf
!  character, dimension(:) :: filename
  character (len=32) :: filename

  integer, dimension(nx, ny) :: red, blue, green
  double precision :: vmax, vmin, deltaV
  double precision :: thresh
  integer :: i, j, cr, cb, cg
  integer, parameter :: iounit = 10

  vmax = maxval(buf)
  vmin = minval(buf)
  thresh = 255.d0
  deltaV = vmax - vmin

  !open file for writing
  open(unit=iounit, file=filename)
  !write header info
  write(iounit,fmt='(''P3''/''# Written by ppmwrite'')')
  write(iounit,*) nx, ny
  write(iounit,*) int(thresh)

  red(:,:) = int((buf(:,:) - vmin) * thresh / (vmax - vmin))
  blue(:,:) = red(:,:)
  green(:,:) = red(:,:)

  do j = 1, ny
    do i = 1, nx
      cr = 255
      cg = 255
      cb = 255
      
      if (buf(i,j) .lt. (vmin + 0.25d0 * deltaV)) then
        cr = 0
        cg = 0
!        cb = 0
!        cg = int(thresh*(4 * (buf(i,j) - vmin) / deltaV))
      else if (buf(i,j) .lt. (vmin + 0.5d0 * deltaV)) then
!         cr = 0
!        cb=0
!        cg=0
!         cb = int(thresh*(1 + 4 * (vmin + 0.25d0*deltaV - buf(i,j)) / deltaV))
      else if (buf(i,j) .lt. (vmin + 0.75d0 * deltaV)) then
!        cr = int(thresh*(4 * (buf(i,j) - vmin - 0.5d0*deltaV) / deltaV))
!          cb = 0
!        cr=0
!        cg=0
      else
!        cg = int(thresh*(1+4*(vmin  + 0.75d0*deltaV-buf(i,j)) / deltaV))
!        cr = 0
        cb=0
        cg=0
      end if
      red(i,j)   = cr
      green(i,j) = cg
      blue(i,j)  = cb
    end do
  end do

  write(iounit,fmt='(5(3I4,'' ''))') ((red(i,j),green(i,j),blue(i,j), i=1,nx), j=1,ny)

  !close file
  close(unit=iounit)
end subroutine ppmwrite

end module utils
