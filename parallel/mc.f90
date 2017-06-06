module mc
  use globvars
  implicit none
  contains
  subroutine MC_init_mesh(mesh, vaclist, nvac, pdvac, pdratio)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, allocatable, dimension(:, :), intent(inout) :: vaclist
    integer, intent(out) :: nvac
    double precision, intent(in) :: pdvac, pdratio

    integer :: i, j, vidx
    double precision :: prob

    nvac = size(vaclist, 2)

    vidx = 1
    vaclist(1, :) = -1
    vaclist(2, :) = -1

    do j = 1, ncol
      do i = 1, nrow
        call random_number(prob)
        if ((prob.lt.pdvac).and.(vidx.le.nvac)) then
          mesh(i, j) = type_vac
          vaclist(1, vidx) = i
          vaclist(2, vidx) = j
          !write(0,'(A6,I3,A3,4I3)')'Init: ', rank, ' : ', vidx, i, j, mesh(i, j)
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
      vcl_up, nvac_up, vcl_dwn, nvac_dwn, vcl_still, nvac_still)
    integer, allocatable, dimension(:, :), intent(inout) ::  vaclist, vcl_lft, vcl_rgt, vcl_up, vcl_dwn, vcl_still
    integer, intent(out) :: nvac_lft, nvac_rgt, nvac_up, nvac_dwn, nvac_still
    integer, intent(in) :: nvac

    integer :: vidx, vx, vy, ii, jj
    double precision :: p

    nvac_lft = 0
    nvac_rgt = 0
    nvac_up  = 0
    nvac_dwn = 0
    nvac_still = 0

    do vidx = 1, nvac
      ii = vaclist(1, vidx)
      jj = vaclist(2, vidx)

      call random_number(p)
      !p = p*0.4d0 + 0.4d0
      !write(0,*)'process',ii,jj,p
      if (p.lt.0.2) then
        vx = 0
        vy =  -1
        nvac_lft = nvac_lft + 1
        vcl_lft(1, nvac_lft) = ii
        vcl_lft(2, nvac_lft) = jj
        vcl_lft(3, nvac_lft) = vx
        vcl_lft(4, nvac_lft) = vy
      else if (p.lt.0.4) then
        vx =  0
        vy =  1
        nvac_rgt = nvac_rgt + 1
        vcl_rgt(1, nvac_rgt) = ii
        vcl_rgt(2, nvac_rgt) = jj
        vcl_rgt(3, nvac_rgt) = vx
        vcl_rgt(4, nvac_rgt) = vy
      else if (p.lt.0.6) then
        vx = 1
        vy = 0
        nvac_dwn = nvac_dwn + 1
        vcl_dwn(1, nvac_dwn) = ii
        vcl_dwn(2, nvac_dwn) = jj
        vcl_dwn(3, nvac_dwn) = vx
        vcl_dwn(4, nvac_dwn) = vy
      else if (p.lt.0.8) then
        vx =  -1
        vy =  0
        nvac_up = nvac_up + 1
        vcl_up(1, nvac_up) = ii
        vcl_up(2, nvac_up) = jj
        vcl_up(3, nvac_up) = vx
        vcl_up(4, nvac_up) = vy
      else
        vx =  0
        vy =  0
        nvac_still = nvac_still + 1
        vcl_still(1, nvac_still) = ii
        vcl_still(2, nvac_still) = jj
        vcl_still(3, nvac_still) = vx
        vcl_still(4, nvac_still) = vy
      end if
    end do

  end subroutine MC_genvel

  subroutine MC_Step(mesh, vaclist, nvac)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, allocatable, dimension(:, :), intent(inout) :: vaclist    
    integer, intent(in) :: nvac


    integer :: i, j, ii, jj, vidx, swap, k
    double precision :: p

    !write(0,*)'nvac = ',nvac
    !Move vacancies
    do vidx = 1, nvac
      !Get the vacancy coordinate
      i = vaclist(1, vidx)
      j = vaclist(2, vidx)
      !Calculate the position
      ii = i + vaclist(3, vidx)
      jj = j + vaclist(4, vidx)
      !Move the vacancy only if destination is NOT vacancy
      if ((mesh(ii, jj).ne.type_vac).and.((ii.ge.1).and.(ii.le.nrow)).and.((jj.ge.1).and.(jj.le.ncol))) then
        vaclist(1, vidx) = ii
        vaclist(2, vidx) = jj

        !write(0,'(I4,A3,A2 ,I3   ,A2 ,I1        ,A1 ,I3,A2  ,I3,A5     ,I1          ,A1 ,I3,A2  ,I3,A1)')&
        !          rank,' : ',' [',vidx,'] ',mesh(i, j),'(',i ,', ',j ,') -> ',mesh(ii, jj),'(',ii,', ',jj,')'
        !swap the elements

        swap = mesh(ii, jj)
        mesh(ii, jj) = mesh(i, j)
        mesh(i, j) = swap
      else
        !write(0,'(I4,A3,A2 ,I3   ,A2 ,I1        ,A1 ,I3,A2  ,I3,A5     ,I1          ,A1 ,I3,A2  ,I3,A1)')&
        !          rank,' : ','s[',vidx,'] ',mesh(i, j),'(',i ,', ',j ,') -> ',mesh(ii, jj),'(',ii,', ',jj,')'
      end if
    end do

  end subroutine MC_Step

  subroutine MC_update_vaclist(vaclist, nvac, buf, buf_size)
    integer, allocatable, dimension(:, :), intent(inout) :: vaclist
    integer, intent(inout) :: nvac
    integer, allocatable, dimension(:, :), intent(in) :: buf
    integer, intent(in) :: buf_size

    integer :: i

    do i = 1, buf_size
      nvac = nvac + 1
      !write(0,'(A8,I3,A2,I3,A2,I3,A1)')'update: ',nvac,' (',buf(1, i),', ',buf(2, i),')'
      vaclist(:, nvac) = buf(:, i)
    end do
  end subroutine MC_update_vaclist

  subroutine MC_sort_vclist(vclist, nvac, row_id)
    implicit none
    integer, allocatable, dimension(:,:), intent(inout) :: vclist
    integer, intent(in) :: nvac, row_id

    integer :: idx, i
    integer, dimension(4) :: buf

    do i = 1,  nvac
      idx = minloc(vclist(row_id,i:nvac),dim=1) + i - 1
      buf(:) = vclist(:, idx)
      vclist(:, idx) = vclist(:, i)
      vclist(:, i) = buf(:)
    end do
  end subroutine MC_sort_vclist

  subroutine MC_select_boundary_vac(vclist, nvac, direction, buf, nvacob)
    implicit none
    integer, allocatable, dimension(:,:), intent(in) :: vclist
    integer, intent(in) :: nvac
    character(len=*), intent(in) :: direction
    integer, allocatable, dimension(:,:), intent(inout) :: buf
    integer, intent(out) :: nvacob

    integer :: i, id, side_id, offset

    if (trim(direction).eq.'left') then
      id = 2
      side_id = 1
      offset = ncol+ 1
    else if (trim(direction).eq.'right') then
      id = 2
      side_id = ncol
      offset = 0
    else if (trim(direction).eq.'up') then
      id = 1
      side_id = 1
      offset = nrow + 1
    else if (trim(direction).eq.'down') then
      id = 1
      side_id = nrow
      offset = 0
    else
      write(0,*)'Unknown direction: ',trim(direction)
      stop
    end if

    nvacob = 0
    do i = 1, nvac
      if (vclist(id,i).eq.side_id) then
        !write(0,*)'[',nvacob,'] (',vclist(1,i),', ',vclist(2,i),')'
        nvacob = nvacob + 1
        buf(:, nvacob) = vclist(:,i)
        buf(id, nvacob) = offset
      end if
    end do
  end subroutine MC_select_boundary_vac

  subroutine MC_join_vac(vaclist, nvac, vac_sub_list, nvac_subc, direction)
    implicit none
    integer, allocatable, dimension(:,:), intent(inout) :: vaclist
    integer, intent(inout) :: nvac
    integer, allocatable, dimension(:,:), intent(inout) :: vac_sub_list
    integer, intent(in) :: nvac_subc
    character(len=*), intent(in) :: direction

    integer :: i, idx, ii, jj

    idx = 1
    do while (idx.le.nvac_subc)
      ii = vac_sub_list(1, idx)
      jj = vac_sub_list(2, idx)

      if (trim(direction).eq.'forward') then
       ii = ii + vac_sub_list(3, idx)
       jj = jj + vac_sub_list(4, idx)
     end if

      if (((ii.ge.1).and.(ii.le.nrow)).and.((jj.ge.1).and.(jj.le.ncol))) then
        nvac = nvac + 1
        vaclist(:, nvac) = vac_sub_list(:, idx)
      end if
      idx = idx + 1
    end do

  end subroutine MC_join_vac

  subroutine MC_init_random_seed(rank)
    implicit none
    integer, intent(in) :: rank
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(COUNT=clock)

    !seed = clock + rank * 37 * (/ (i - 1, i = 1, n) /)
    seed = rank * 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)

    deallocate(seed)
    end subroutine MC_init_random_seed
end module mc
