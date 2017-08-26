module mc
  use globvars
  implicit none
  contains
  subroutine MC_init_mesh(mesh, nvac, pbvac, pbratio,rank)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, intent(out) :: nvac
    double precision, intent(in) :: pbvac, pbratio
    integer :: rank

    integer :: i, j, vidx
    double precision :: prob

    nvac=0
   do i = 1, nrow
   do j = 1, ncol
   call random_number(prob)
   if (prob.lt.pbvac) then
           mesh(i,j) = type_vac
        nvac = nvac + 1
   elseif (prob.lt.(1d0-pbvac)*pbratio) then
           mesh(i,j) =  type_red
   else
           mesh(i,j) =  type_blue
   endif
   enddo
   enddo

!  if(rank.eq.3)then
!          mesh=type_red
!  else
!          mesh=type_vac
!  endif

!  write(*,*)'ratio',(1d0-pbvac)*pbratio
  end subroutine MC_init_mesh
!!!! -  bez genvel + mc_step
!--------------------------------------
      subroutine xLatex(mcount,ncount,mesh, vac, nvac ) ! Performs the vacancy exchanges
       integer mcount,ncount, LS2, nvac,is,il,ir,ix,iy
       integer mesh(0:mcount+1,0:ncount+1), vac(nvac,2)
!       integer ip(LS),im(LS)
      real r 
        integer seed, len, junk,streamnum, nstreams
        integer gtype,stream
       LS2 = mcount * ncount                 ! Total number of attempts = lattice area

!        seed = make_sprng_seed()
!        stream = init_sprng(seed,SPRNG_DEFAULT,1)
!

       !0.857
       call random_number(r)
       if(r.lt.pbfrmv)then
       do is = 1, nvac

       ix = vac(is,1)
       iy = vac(is,2)
                   ! then choose the direction of exchange

       if(ix.eq.-1.or.iy.eq.-1)then
               return
       endif
       call random_number(r)
          if (r.lt.0.50) then ! exclude right
             il = mesh(ix-1,iy)
             ir = mesh(ix+1,iy)
             if(il.eq.ir.and.il.ne.0)then
                mesh(ix,iy)=il
                mesh(ix+1,iy)=0
             endif
          elseif (r.lt.1.00) then !exclude left
             il = mesh(ix+1,iy)
             ir = mesh(ix-1,iy)
             if(il.eq.ir.and.il.ne.0)then
                mesh(ix,iy)=il
                mesh(ix-1,iy)=0
             endif
          endif
       enddo

       else
!diffusion
       do is = 1, nvac
!c       generate random coordinates of the current site to be checked
       
       
       ix = vac(is,1)
       iy = vac(is,2)
                   ! then choose the direction of exchange

       if(ix.eq.-1.or.iy.eq.-1)then
               return
       endif
       call random_number(r)

          if (r.lt.0.5) then ! exclude right
             il = mesh(ix-1,iy)
             ir = mesh(ix+1,iy)
             mesh(ix,iy)=ir
             mesh(ix+1,iy)=0
          elseif (r.lt.1.00) then !exclude left
             il = mesh(ix+1,iy)
             ir = mesh(ix-1,iy)
             mesh(ix,iy)=ir
             mesh(ix-1,iy)=0 
          endif
       enddo
       endif
       return
       end



!c-----------------------------------------
        subroutine yLatex(mcount,ncount,mesh,vac,nvac) ! Performs the vacancy exchanges
       integer mcount,ncount, LS2,nvac, is, ix, iy, il, ir
       integer mesh(0:mcount+1,0:ncount+1),vac(nvac,2)
!       integer ip(LS),im(LS)
      real r 
!#include "sprng_f.h"
!#define SIMPLE_SPRNG	
        integer seed, len, junk,streamnum, nstreams
        integer gtype,stream
       LS2 = mcount * ncount                 ! Total number of attempts = lattice area

!        seed = make_sprng_seed()
!        stream = init_sprng(seed,SPRNG_DEFAULT,1)
!

       !0.857
       call random_number(r)
       if(r.lt.pbfrmv)then
       do is = 1, nvac
       ix = vac(is,1)
       iy = vac(is,2)

       if(ix.eq.-1.or.iy.eq.-1)then
               return
       endif
                   ! then choose the direction of exchange

       call random_number(r)
          if(r.lt.0.5) then ! exclude down
             il = mesh(ix,iy-1)
             ir = mesh(ix,iy+1)
             if(il.eq.ir.and.il.ne.0)then
                mesh(ix,iy)=il
                mesh(ix,iy+1)=0
             endif
          else                     ! exlcude up
             il = mesh(ix,iy+1)
             ir = mesh(ix,iy-1)
             if(il.eq.ir.and.il.ne.0)then
                mesh(ix,iy)=il
                mesh(ix,iy-1)=0        
             endif
          endif
       enddo

       else
!diffusion
       do is = 1, nvac
       ix = vac(is,1)
       iy = vac(is,2)
                   ! then choose the direction of exchange
       if(ix.eq.-1.or.iy.eq.-1)then
            return
       endif

       call random_number(r)

          if(r.lt.0.50) then ! exclude down
             il = mesh(ix,iy-1)
             ir = mesh(ix,iy+1)
             mesh(ix,iy)=ir
             mesh(ix,iy+1)=0
           elseif(r.lt.1.0) then                    ! exlcude up
             il = mesh(ix,iy+1)
             ir = mesh(ix,iy-1)
             mesh(ix,iy)=ir
             mesh(ix,iy-1)=0          
          endif
       enddo
       endif
       return
       end
!c--------------------------------------------


        subroutine updatevaclist (nrow,ncol, mesh,nvac, vac )

                integer :: nrow, ncol 
                integer :: nvac
                integer :: mesh(0:nrow+1,0:ncol+1), vac(nvac,2)

                integer :: i, j, counter 

                counter = 1 
                vac(:,:) = -1

                do i=1, ncol
                   do j=1, nrow
                   
                   if(mesh(j,i).eq.0)then
                           vac(counter,1)=j
                           vac(counter,2)=i
                           counter=counter+1
                   endif

                   enddo
                enddo

        return
        end


!c-----------------------------------------


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
