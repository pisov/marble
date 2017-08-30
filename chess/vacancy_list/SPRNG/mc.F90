module mc
  use globvars
  implicit none
        
#define SIMPLE_SPRNG
  contains

  subroutine MC_init_mesh(mesh, nvac, pbvac, pbratio,rank,stream)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, intent(out) :: nvac
    double precision, intent(in) :: pbvac, pbratio
    integer :: rank,commsize
#include "sprng_f.h"


    integer :: i, j, vidx
    real*8 :: prob
! SPRNG
    SPRNG_POINTER stream

    nvac=0
   do i = 1, nrow
   do j = 1, ncol
    prob=sprng()
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

  end subroutine MC_init_mesh
!--------------------------------------
  subroutine MC_add_vac(mesh, nvac, pbvac, pbratio,rank,stream)
    integer, allocatable, dimension(:, :), intent(inout) :: mesh
    integer, intent(out) :: nvac
    double precision, intent(in) :: pbvac, pbratio
    integer :: rank

    integer :: i, j, vidx
    real*8 :: prob

    SPRNG_POINTER stream
#include "sprng_f.h"

    nvac=0
   do i = 1, nrow
   do j = 1, ncol
   prob=sprng()
   if (prob.lt.pbvac) then
           mesh(i,j) = type_vac
        nvac = nvac + 1
   endif
   enddo
   enddo


  end subroutine MC_add_vac
!--------------------------------------
      subroutine xLatex(mcount,ncount,mesh, vac, nvac,stream ) ! Performs the vacancy exchanges
       integer mcount,ncount, LS2, nvac,is,il,ir,ix,iy
       integer mesh(0:mcount+1,0:ncount+1), vac(nvac,2)
       real*8 r 
#include "sprng_f.h"
! SPRNG
       SPRNG_POINTER stream
       LS2 = mcount * ncount                 ! Total number of attempts = lattice area
       
!

       !0.857
       r=sprng()
       if(r.lt.pbfrmv)then
       do is = 1, nvac

       ix = vac(is,1)
       iy = vac(is,2)
                   ! then choose the direction of exchange

       if(ix.eq.-1.or.iy.eq.-1)then
               return
       endif
       r=sprng()
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
       r=sprng()

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
        subroutine yLatex(mcount,ncount,mesh,vac,nvac,stream) ! Performs the vacancy exchanges
       integer mcount,ncount, LS2,nvac, is, ix, iy, il, ir
       integer mesh(0:mcount+1,0:ncount+1),vac(nvac,2)
       real*8 r 
#include "sprng_f.h"
       SPRNG_POINTER stream

       LS2 = mcount * ncount                 ! Total number of attempts = lattice area


       !0.857
       r=sprng()
       if(r.lt.pbfrmv )then
       do is = 1, nvac
       ix = vac(is,1)
       iy = vac(is,2)

       if(ix.eq.-1.or.iy.eq.-1)then
               return
       endif
                   ! then choose the direction of exchange

       r=sprng()
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

       r=sprng()

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


end module mc
