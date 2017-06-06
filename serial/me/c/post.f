        program post
        implicit none
c     Simulates a binary system with vacancies, they play the only (main) role
c     In this version are not used pointers to the coordinates of the vacancies
c     but they are placed on a lattice.
c     Thus it is specially suited for the cases with large concentrations of 
c     vacancies
c     Well branched system of monitoring is to be developed.
c
c     (c) Vesselin Tonchev, Marcin Minkowski, 2016, Sofia - Warsaw
        integer LS                 ! Size of the lattice
        integer,dimension(:,:),allocatable :: CA          ! The lattice itself
        integer i             ! counter
        real enpyx,enpyy      ! 
        character*128 :: filename
        character*16 :: iter
        character*128 :: path

        !get the size of the CA
        open (unit = 213, file='siminfo.txt')
           read(213,*) LS
        close(unit = 213)
        
        allocate(CA(LS,LS))
        
        call getarg(1,filename)
        iter= trim(filename(len_trim(filename)-11:len_trim(filename)-4))
        path= trim(filename(1:(len_trim(filename)-12)))


        open (unit = 211, file=filename)
        do i=1,LS
           read(211,*) CA(i,:)
        end do
        close(unit = 211)
        
        call clus(LS,CA,path,iter) 
        call draw(LS,CA,path,iter)

        deallocate(CA)
        
c end of executions block
        end
c-----------------------------------------
        subroutine clus(LS,CA,path,iter) 
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
            if((CA(ix-1,iy).eq.CA(ix,iy)).and.(CA(ix,iy).eq.
     &        CA(ix,iy-1)).and.(ix.ne.1).and.(iy.ne.1))then
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
                   if(CA(ix,iy).eq.-1)then
                      label(ix,iy,2)=1
                   elseif(CA(ix,iy).eq.0)then
                      label(ix,iy,2)=2
                   else
                      label(ix,iy,2)=3
                   endif
                  
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
c--  --  -- --  --  --  --

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
 
c---------------------
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


c-----------------------------------------
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

c-----------------------------------------
        subroutine draw(LS,CA,path,iter)
        use utils
        implicit none
        integer LS,it
        integer CA(LS,LS)
        double precision buf(4*LS,4*LS)
        integer i,j,row,column
        character*16 :: iter
        character*128 :: path
        character*256 :: filename

        buf(:,:)=0
        do, i=1,LS
           do, j=1,LS
              do, row=0, 3
                 do, column=0,3
                    if (CA(i,j).eq.-1)then
                      ! buf(i,j)=255.0
                      buf(i*4-row,j*4-column)=255.0
                   elseif(CA(i,j).eq.0)then
                     !  buf(i,j)=127.0
                      buf(i*4-row,j*4-column)=127.0
                   elseif(CA(i,j).eq.1)then
                     ! buf(i,j)=0.0
                     buf(i*4-row,j*4-column)=0.0
                   end if
                 end do
              end do
            end do
         end do

!      write(filename,'(A20)')'img/' //trim(iter) //'.ppm'
      write(filename,'(A16)')'img/' //trim(iter) //'.ppm'
      call ppmwrite(filename,buf(:,:),4*LS,4*LS)
        return
        end
c-----------------------------------------
