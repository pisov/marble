module mc
  use globvars
  implicit none
  contains
  subroutine MC_init_mesh(mesh, vaclist, nvac, pdvac, pdratio)
    integer, dimension(:, :), intent(inout) :: mesh
    integer, dimension(:, :), intent(out) :: vaclist
    integer, intent(out) :: nvac
    double precision, intent(in) :: pdvac, pdratio

    integer :: nrow, ncol, i, j, vidx
    double precision :: prob

    nrow = size(mesh, 1)
    ncol = size(mesh, 2)
    nvac = size(vaclist, 2)

    vidx = 1
    vaclist(1, :) = -1
    vaclist(2, :) = -1

    do j = 1, nrow
      do i = 1, ncol
        call random_number(prob)
        if ((prob.lt.pdvac).and.(vidx.le.nvac)) then
          mesh(i, j) = type_vac
          vaclist(1, vidx) = i
          vaclist(2, vidx) = j
          !write(0,'(4I3)')vidx,i,j,mesh(i, j)
          vidx = vidx + 1
        else
          call random_number(prob)
          if(prob.lt.pdratio) then
            mesh(i, j) = type_red
          else
            mesh(i, j) = type_blue
          end if
        end if
      end do
    end do
    nvac = vidx - 1
  end subroutine MC_init_mesh

  subroutine MC_genvel(vaclist, nvac,vcl_lft, nvac_lft, vcl_rgt, nvac_rgt,&
      vcl_up, nvac_up, vcl_dwn, nvac_dwn)
    integer, dimension(:, :), intent(inout) ::  vaclist, vcl_lft, vcl_rgt, vcl_up, vcl_dwn
    integer, intent(out) :: nvac_lft, nvac_rgt, nvac_up, nvac_dwn
    integer, intent(in) :: nvac

    integer :: vidx, vx, vy, ii, jj
    double precision :: p

    nvac_lft = 0
    nvac_rgt = 0
    nvac_up  = 0
    nvac_dwn = 0

    do vidx = 1, nvac
      ii = vaclist(1, vidx)
      jj = vaclist(2, vidx)
      call random_number(p)
      if (p.lt.0.2) then
        vx = -1
        vy =  0
        nvac_lft = nvac_lft + 1
        vcl_lft(1, vidx) = ii
        vcl_lft(2, vidx) = jj
        vcl_lft(3, vidx) = vx
        vcl_lft(4, vidx) = vy
      else if (p.lt.0.4) then
        vx =  1
        vy =  0
        nvac_rgt = nvac_rgt + 1
        vcl_rgt(1, vidx) = ii
        vcl_rgt(2, vidx) = jj
        vcl_rgt(3, vidx) = vx
        vcl_rgt(4, vidx) = vy
      else if (p.lt.0.6) then
        vx =  0
        vy = -1
        nvac_dwn = nvac_dwn + 1
        vcl_dwn(1, vidx) = ii
        vcl_dwn(2, vidx) = jj
        vcl_dwn(3, vidx) = vx
        vcl_dwn(4, vidx) = vy
      else if (p.lt.0.8) then
        vx =  0
        vy =  1
        nvac_up = nvac_up + 1
        vcl_up(1, vidx) = ii
        vcl_up(2, vidx) = jj
        vcl_dwn(3, vidx) = vx
        vcl_dwn(4, vidx) = vy
      else
        vx =  0
        vy =  0
      end if
      !write(0,*)vidx,vx,vy
      vaclist(3, vidx) = vx
      vaclist(4, vidx) = vy
    end do

  end subroutine MC_genvel

  subroutine MC_Step(mesh, vaclist, nvac)
    integer, dimension(:, :), intent(inout) :: mesh
    integer, dimension(:, :), intent(inout) :: vaclist
    integer, intent(in) :: nvac
    character(len=*), intent(in) :: direction

    integer :: i, j, ii, jj, vidx, swap, k
    double precision :: p

    !Move vacancies    

    do vidx = 1, nvac
      i = vaclist(1, vidx)
      j = vaclist(2, vidx)

    end do

  end subroutine MC_Step

  subroutine update_vaclist(vaclist, nvac, buf, buf_size)
    integer, dimension(:, :), intent(inout) :: vaclist
    integer, intent(inout) :: nvac
    integer, dimension(:, :), intent(in) :: buf
    integer, intent(in) :: buf_size

    integer :: i

    do i = 1, buf_size
      nvac = nvac + 1
      vaclist(:, nvac) = buf(:, i)
    end do
  end subroutine update_vaclist
end module mc
