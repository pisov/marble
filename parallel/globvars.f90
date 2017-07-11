module globvars
  implicit none
  include 'mpif.h'

  integer, parameter :: type_vac = 0
  integer, parameter :: type_red = 1
  integer, parameter :: type_blue = 2

  character(len=64) :: chkp_filename = 'marble.chkp'
  character(len=64) :: inpt_filename = 'marble.in'

  !MC variables
  !nit - number of productive iteration steps
  !it - iteration step
  !nout - number of step to writeout data
  !nquib - number of MC equilibrium steps before writeout
  !nsize - grid size
  !ncol - number of columns of local mesh
  !nrow - number of row of local mesh
  integer(kind=8) :: nit = 1, it, nout = 1, nequib = 0, neout = 1
  integer :: nsize = 128, ncol, nrow, maxvacnum
  !The Temperature of the sample
  !Total energy of the sample
  !pbvac - probability density of vacancies
  !pbratio - probability dencity ratio between type I and type II particles
  double precision :: Temp = 0.d0, Energy, pbvac = 0.01d0, pbratio = 0.5d0
  double precision :: pbfrmv = 0.1d0

  !MPI specific variables
  !commsize - size of MPI comunivator (number of processes)
  !rank - process rank
  integer :: commsize, rank

  !dims - 1d array with dimensions of Cartesian topology
  !period - array controling period boundary conditions of Cartesian topology

  !integer(KIND=MPI_OFFSET_KIND), dimension(2) :: dims
  integer, dimension(2) :: dims, crd
  logical, dimension(2) :: periods
  integer :: MPI_COMM_2D, MPI_TYPE_ROW, MPI_TYPE_COL, MPI_BLOCK, MPI_BLOCK2

  !scounts - array with number of cunks to distribute among processors
  !displs - array with calculate chunk displacement

  integer, allocatable, dimension(:) :: scounts
  integer(kind=MPI_OFFSET_KIND), allocatable, dimension(:) :: displs
  integer :: up, down, left, right

end module globvars
