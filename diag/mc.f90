module mc
  use globvars
  implicit none
  contains
  subroutine MC_init_mesh(mesh, vaclist, nvac, pbvac, pbratio)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, allocatable, dimension(:, :), intent(inout) :: vaclist
    integer, intent(out) :: nvac
    double precision, intent(in) :: pbvac, pbratio

    integer :: i, j, vidx
    double precision :: prob

    nvac = size(vaclist, 2)

    vidx = 1
    vaclist(1, :) = -1
    vaclist(2, :) = -1

    do j = 1, ncol
      do i = 1, nrow
        call random_number(prob)
        if ((prob.lt.pbvac).and.(vidx.le.nvac)) then
          mesh(i, j) = type_vac
          vaclist(1, vidx) = i
          vaclist(2, vidx) = j
          !write(0,'(A6,I3,A3,4I3)')'Init: ', rank, ' : ', vidx, i, j, mesh(i, j)
          vidx = vidx + 1
        else
          call random_number(prob)
          if(prob.lt.pbratio) then
            mesh(i, j) = type_red
          else
            mesh(i, j) = type_blue
          end if
        end if
      end do
    end do
    nvac = vidx - 1
  end subroutine MC_init_mesh

  subroutine MC_genvel(mesh, vaclist, nvac)
    integer, allocatable, dimension(:, :), intent(in) :: mesh
    integer, allocatable, dimension(:, :), intent(inout) ::  vaclist
    integer, intent(in) :: nvac

    integer :: vidx, vx, vy, ii, jj
    double precision :: p, q

    do vidx = 1, nvac
      ii = vaclist(1, vidx)
      jj = vaclist(2, vidx)

      call random_number(q)
      if (q.lt.pbfrmv) then
        call random_number(p)
        if ((p.lt.0.25).and.(mesh(ii,jj-1).ne.type_vac)) then
          vx = 0
          vy =  -1
        else if ((p.lt.0.50).and.(mesh(ii,jj+1).ne.type_vac)) then
          vx =  0
          vy =  1
        else if ((p.lt.0.75).and.(mesh(ii+1,jj).ne.type_vac)) then
          vx = 1
          vy = 0
        else if (mesh(ii-1,jj).ne.type_vac) then
          vx =  -1
          vy =  0
        else
          vx =  0
          vy =  0
        end if
        vaclist(3, vidx) = vx
        vaclist(4, vidx) = vy
      else
        ! Move vacancies toward sametype
        vx = 0
        vy = 0
        if ((mesh(ii-1, jj).eq.mesh(ii+1, jj)).and.(mesh(ii-1,jj).ne.type_vac)) then
          vx = 1
        else
          vx = 0
        end if

        if ((mesh(ii, jj-1).eq.mesh(ii, jj+1)).and.(mesh(ii, jj-1).ne.type_vac)) then
          vy = 1
        else
          vy = 0
        end if
      
        if ( (vx*vy).eq.1 ) then
          call random_number(q)
          if (q.le.0.5) then
            vy = 0
          else
            vx = 0
          end if
        end if

        if(vx.eq.1) then
          call random_number(q)
          if (q.le.0.5) then
            vx = 1
            vy = 0
          else
            vx = -1
            vy = 0
          end if
        else if (vy.eq.1) then 
          call random_number(q)
          if (q.le.0.5) then
            vx =  0
            vy =  1
          else
            vx = 0
            vy =  -1
          end if
        else
          vx =  0
          vy =  0
        end if
        vaclist(3, vidx) = vx
        vaclist(4, vidx) = vy
      end if
    end do

  end subroutine MC_genvel

  subroutine MC_Step(mesh, vaclist, nvac, side)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, allocatable, dimension(:, :), intent(inout) :: vaclist
    integer, intent(in) :: nvac
    character(len=*), intent(in) :: side 

    integer :: i, j, ii, jj, vidx, swap 
    logical :: move

    !Move vacancies
    do vidx = 1, nvac
      !Get the vacancy coordinate
      i = vaclist(1, vidx)
      j = vaclist(2, vidx)
      !Decide to move based on side
      if (trim(side).eq.'down' ) then
        if (i.gt.j) then
          move = .true.
        else
          move = .false.
        end if 
      else
        if (i.le.j) then
          move = .true.
        else
          move = .false.
        end if  
      endif
      if (move) then
        !Calculate the new position
        ii = i + vaclist(3, vidx)
        jj = j + vaclist(4, vidx)

        !Update new vacancy list
        !Move the vacancy only if destination is NOT vacancy
        if ((mesh(ii, jj).ne.type_vac)) then
          vaclist(1, vidx) = ii
          vaclist(2, vidx) = jj

          swap = mesh(ii, jj)
          mesh(ii, jj) = mesh(i, j)
          mesh(i, j) = swap
        else
          ii = i
          jj = j
        end if

        vaclist(3, vidx) = 0
        vaclist(4, vidx) = 0
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

  subroutine MC_join_vac(vaclist, nvac, vac_sub_list, nvac_subc)
    implicit none
    integer, allocatable, dimension(:,:), intent(inout) :: vaclist
    integer, intent(inout) :: nvac
    integer, allocatable, dimension(:,:), intent(inout) :: vac_sub_list
    integer, intent(in) :: nvac_subc

    integer :: i, idx, ii, jj

    idx = 1
    do while (idx.le.nvac_subc)
      nvac = nvac + 1
      vaclist(:, nvac) = vac_sub_list(:, idx)
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
