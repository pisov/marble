program dump
use globvars
use utils
implicit none
        real enpyx,enpyy      ! 
!        character*128 :: filename
        character*16 :: iter
!        character*128 :: path


  integer :: i, j,nuiter
  integer :: fh, cval,thresh
  integer :: cr, cg, cb
  character(len=32) :: filename
  integer, allocatable, dimension(:,:) :: CA
  integer :: nuclu


  fh = 15
  call read_input_file(inpt_filename)
  allocate(CA(nsize,nsize))

  if (iargc().eq.0) then
    open(fh, file=chkp_filename, form='binary')
  else
    call getarg(1, filename)
    open(fh, file=trim(filename), form='binary')
  end if
  
  do j=1, nsize
  do i=1, nsize
  read(fh)cval
  CA(j,i)=cval
  enddo
  enddo

       iter= trim(filename(len_trim(filename)-5:len_trim(filename)))
       read(iter,*) nuiter
        
        call clus(nsize,CA,iter,nuclu) 
        call draw(nsize,CA,iter)
        call enpy(nsize,CA,enpyx,enpyy) 
        
! write (*,'(f6.5,3X,I9,3X,I9,3X,f10.6,3X,f10.6)') real(pbvac) ,nuiter*nout,nuclu,enpyx,enpyy
  close(fh)
  deallocate(CA)
!c end of executions block
        end
!c-----------------------------------------
        subroutine clus(LS,CA,iter,nuclu) 
       use globvars
       implicit none

       integer LS,   ix, iy,i,j
       integer CA(LS,LS), label(LS,LS,2), clusters(LS*LS)
       integer max_label,nuclu,newlabel,oldlabel
       integer, dimension(:,:),allocatable :: clusterdata, shortclusdata
       integer :: counter,counter2,nusizes,typeclus
       integer,dimension(:), allocatable :: sizes
       logical, dimension(:), allocatable :: mask
       character*16 :: iter
       character*128 :: path       

!label(*,*,1) -> labels
!label(*,*,2) -> type of cluster

       label(1:LS,1:LS,1:2)=0

       max_label = 1
       do, ix=1,LS
          do, iy=1,LS
          ! if bouth up and left are the same
            if((CA(ix-1,iy).eq.CA(ix,iy)).and.(CA(ix,iy).eq. &
                    CA(ix,iy-1)).and.(ix.ne.1).and.(iy.ne.1))then
                   oldlabel=label(ix-1,iy,1)
                   newlabel=label(ix,iy-1,1)
                   label(ix,iy,1)=newlabel

                   do i=1,ix
                      do j=1, iy
                         if(label(i,j,1).eq.oldlabel)then
                            label(i,j,1)=newlabel
                         endif
                      enddo
                   enddo

          !if up is the same
            elseif((CA(ix-1,iy).eq.CA(ix,iy)).and.(ix.ne.1))then

                   label(ix,iy,1)=label(ix-1,iy,1)
                   label(ix,iy,2)=label(ix-1,iy,2)
!
!
           !if left is the same
             elseif((CA(ix,iy-1).eq.CA(ix,iy)).and.(iy.ne.1))then
                   label(ix,iy,1)=label(ix,iy-1,1)
                   label(ix,iy,2)=label(ix,iy-1,2)
             else
                  label(ix,iy,1)= max_label
                   max_label=max_label+1
             endif
                   if(CA(ix,iy).eq.1)then
                      label(ix,iy,2)=1
                   elseif(CA(ix,iy).eq.0)then
                      label(ix,iy,2)=0
                   else
                      label(ix,iy,2)=3
                   endif
!                   if(mod(ix,100).eq.0.and.mod(iy,100).eq.0)then
!                      write(*,*)'aaa',ix,iy
!                   endif
          enddo
       enddo
        
        !boundary conditions
        do ix=1,LS
           !down
           if(CA(LS,ix).eq.CA(1,ix))then
                   oldlabel=label(LS,ix,1)
                   newlabel=label(1,ix,1)
                   label(LS,ix,1)=newlabel

                   do i=1,LS
                      do j=1, LS
                         if(label(i,j,1).eq.oldlabel)then
                            label(i,j,1)=newlabel
                         endif
                      enddo
                   enddo
           endif
           !right
           if(CA(ix,LS).eq.CA(ix,LS))then
                   oldlabel=label(ix,LS,1)
                   newlabel=label(ix,LS,1)
                   label(LS,ix,1)=newlabel

                   do i=1,LS
                      do j=1, LS
!                   if(mod(i,100).eq.0.and.mod(j,100).eq.0)then
!                      write(*,*)'bbb',i,j
!                   endif
                         if(label(i,j,1).eq.oldlabel)then
                            label(i,j,1)=newlabel
                         endif
                      enddo
                   enddo
           endif     

        enddo

        
!make a list of the cluster names
        do i=1,LS
           do,j=1, LS
              clusters((i-1)*LS+j)=label(i,j,1)
           enddo
        enddo
!sort it        
        call buble(LS*LS,clusters)
        oldlabel=-1
!count it
        nuclu=0
        do i=1,LS*LS
              newlabel=clusters(i)
              if(oldlabel.ne.newlabel)then
                 nuclu=nuclu+1
              endif

              oldlabel=clusters(i)
        enddo

!fulling clusterdata with info about the clusters
!clusterdata(*,1) -> label
!clusterdata(*,2) -> size
!clusterdata(*,3) -> type

        !labels
        allocate(clusterdata(nuclu,3))
        counter=1
        oldlabel=-1
        do i=1,LS*LS
           newlabel=clusters(i)
           if(oldlabel.ne.newlabel)then
              clusterdata(counter,1)=newlabel
              counter=counter+1
           endif
           oldlabel=clusters(i)
        enddo
        

        !size and type of clusters
        do i=1, nuclu
           counter=0
           typeclus=5
           newlabel=clusterdata(i,1)

           do ix=1, LS
              do iy=1, LS
                 if(newlabel.eq.label(ix,iy,1))then
                    counter = counter+1
                    typeclus = label(ix,iy,2)!not very smart :/
                 endif
              enddo
           enddo
           
          clusterdata(i,2)=counter
          clusterdata(i,3)=typeclus
        enddo

        open (unit = 111, file='postdata/' // trim(iter) //'clusterdata.&
     &txt')
        do i=1,nuclu
           write(111,*) clusterdata(i,:)
        end do
        close(unit = 111)
!c--  --  -- --  --  --  --

        allocate(sizes(nuclu))       
        allocate(mask(nuclu))       
!       
!make a list of the cluster sizes
       
       sizes=clusterdata(:,2)

!sort it        
        call buble(nuclu,sizes)

        oldlabel=clusters(1)
!count it
        nusizes=0
!reusing the variables new/oldlabel as new/oldsize of clusters        
        oldlabel=-1
        do i=1,nuclu
              newlabel=sizes(i)
              if(oldlabel.ne.newlabel)then
                 nusizes=nusizes+1
              endif
              oldlabel=sizes(i)
        enddo
        
        
       allocate(shortclusdata(nusizes,2))

!shortclusterdata(*,1) -> size of the clusters
!shortclusterdata(*,2) -> count of the clusters of this size

!fulling shortclusterdata with data

       !sizes
       oldlabel = -1
       counter = 1
       

       do i=1,nuclu
          newlabel=sizes(i)
       mask=newlabel.eq.sizes
          if(newlabel.ne.oldlabel)then
             shortclusdata(counter,1)=newlabel
             shortclusdata(counter,2)=count(mask)
             counter=counter+1
          elseif(newlabel.eq.oldlabel)then
          endif

          oldlabel=sizes(i)
       enddo
        open (unit = 112,file='./postdata/'// trim(iter) //'shortclusdat&
     &a.txt')
        do i=1, nusizes
           write(112,*) shortclusdata(i,:)
        end do
        close(unit = 112)



        deallocate(shortclusdata)
        deallocate(clusterdata)
        deallocate(sizes)       
        deallocate(mask)       
        return
        end
 
!c---------------------
      subroutine buble(LS,clusters)
      implicit none
      integer LS
      integer,dimension(:) :: clusters(LS)
      integer :: temp
      INTEGER :: i, j
      LOGICAL :: swapped

      do j = LS-1, 1, -1
      swapped = .FALSE.
      do i = 1, j
      if (clusters(i) .gt. clusters(i+1)) then
        temp = clusters(i)
        clusters(i) = clusters(i+1)
        clusters(i+1) = temp
        swapped = .TRUE.
      endif
      end do
      if (.not. swapped)exit
      end do
      
      
      return
      end


!c-----------------------------------------
        subroutine enpy(LS,CA,enpyx,enpyy) 
       integer LS,   ix, iy,i
       integer CA(LS,LS)
       integer enpy1(LS)
       integer enpy2(LS)
       real enpyx,enpyy
       integer br1,br2,aa,bb,cc,dd

       bb=-2
       dd=-2
       do ix=1,LS
          br1=0
          br2=0

          do iy=1,LS
             aa=CA(ix,iy)
             cc=CA(iy,ix)
           

             if(aa.eq.0)then
                i=iy
                do while (aa.eq.0.and.i.ne.1)
                   i=i-1
                   aa=CA(ix,i)
                enddo
             endif
             if(aa.eq.bb.and.aa.ne.0)then
                br1=br1+1
             endif


             if(cc.eq.0)then
                i=ix
                do while (cc.eq.0.and.i.ne.1)
                   i=i-1
                   cc=CA(i,iy)
                enddo
             endif
             if(cc.eq.dd.and.cc.ne.0)then
                br2=br2+1
             endif


             dd= CA(iy,ix)
             bb= CA(ix,iy)
          enddo
          enpy1(ix)=LS-1-br1
          enpy2(ix)=LS-1-br2
       enddo
        enpyx=SUM(enpy1)/real(LS)
        enpyy=SUM(enpy2)/real(LS)

        return
        end

!c-----------------------------------------
        subroutine draw(LS,CA,iter)
        use utils
        implicit none
        integer LS,ir
        integer CA(LS,LS)
        double precision,dimension(:,:),allocatable :: buf
        integer i,j,row,column,scale_img
        character*16 :: iter
        character*128 :: path
        character*256 :: filename

        scale_img=1
        allocate(buf(scale_img*LS,scale_img*LS))
        buf(:,:)=0
        do, i=1,LS
           do, j=1,LS
              do, row=0, scale_img-1
                 do, column=0,scale_img-1
                    if (CA(i,j).eq.1)then
                      ! buf(i,j)=255.0
                      buf(i*scale_img-row,j*scale_img-column)=255.0
                   elseif(CA(i,j).eq.0)then
                     !  buf(i,j)=127.0
                      buf(i*scale_img-row,j*scale_img-column)=127.0
                   elseif(CA(i,j).eq.2)then
                     ! buf(i,j)=0.0
                     buf(i*scale_img-row,j*scale_img-column)=0.0
                   end if
                 end do
              end do
            end do
         end do

!      write(filename,'(A20)')'img/' //trim(iter) //'.ppm'
      write(filename,'(A16)')'img/' //trim(iter) //'.ppm'
      call ppmwrite(filename,buf(:,:),scale_img*LS,scale_img*LS)
      deallocate(buf)
        return
        end
!c-----------------------------------------
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
!        cg = 0
!        cb = 0
!        cr = 0
!         cb = int(thresh*(1 + 4 * (vmin + 0.25d0*deltaV - buf(i,j)) / deltaV))
      else if (buf(i,j) .lt. (vmin + 0.75d0 * deltaV)) then
!        cr = int(thresh*(4 * (buf(i,j) - vmin - 0.5d0*deltaV) / deltaV))
!        cb = 0
!        cr = 0
!        cg = 0
      else
!        cg = int(thresh*(1+4*(vmin  + 0.75d0*deltaV-buf(i,j)) / deltaV))
!        cb = 0

        cb = 0
!        cr = 0
        cg = 0
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

