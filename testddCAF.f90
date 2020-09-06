! program to test the domain decomposition 
! compile: caf testddCAF.f90 -Wall -Wextra -Wuninitialized -fbounds-check -o testddCAF
! run: echo cafrun -np 2 testddCAF
! Date: Sun Aug 23 21:41:12 EDT 2020
! vim:fdm=marker
module thingy

integer :: myImg, myNpts, ourNpts, Nim
logical :: PBC

real(8), dimension(3) :: BoxSize
real(8), allocatable :: pos(:,:)[:] ! coarray definition
real(8), allocatable :: pene(:)
real(8) :: PotEnergy 

type communicator
  integer :: NNcells ! number of neighbouring image cels
  integer, dimension(:), allocatable :: NeighCellID ! list of neighbour cell ID
  integer :: maxOurNpts ! for allocation
  integer, dimension(2) :: bulkrange
  integer :: Nsend, Nget
  integer, dimension(:), allocatable :: sourceCellID
  integer, dimension(:,:), allocatable :: sendrange, getrange, sourcerange
end type communicator
type(communicator) :: com

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readCOMfile( fil, Nim ) !{{{
implicit none
integer, intent(in) :: Nim
character(*) :: fil
integer :: i, j
character(40) :: ctmp

open(10,file=trim(fil))

read(10,*) com%maxOurNpts ! max image ourNpts for allocation of coarray
read(10,*) BoxSize ! global boxsize
do i = 1, Nim ! images are writen in order, so no need to check
  if (i.eq.myImg) then ! read the section
    read(10,*)  j, myNpts, ourNpts, com%NNcells
      allocate( com%NeighCellID( com%NNcells ) )
    read(10,*) com%NeighCellID(:)
    read(10,*) ctmp, com%bulkrange
    read(10,*) com%Nsend
      allocate( com%sendrange( 2, com%Nsend ) )
    read(10,*) ctmp, com%sendrange( 1, : ) ! min
    read(10,*) ctmp, com%sendrange( 2, : ) ! max
    read(10,*) com%Nget
      allocate( com%getrange( 2, com%Nget ) )
      allocate( com%sourceCellID( com%Nget ) )
      allocate( com%sourcerange( 2, com%Nget ) )
    read(10,*) ctmp, com%getrange( 1, : ) ! min
    read(10,*) ctmp, com%getrange( 2, : ) ! max
    read(10,*) ctmp, com%sourceCellID(:)
    read(10,*) ctmp, com%sourcerange( 1, : ) ! min
    read(10,*) ctmp, com%sourcerange( 2, : ) ! max
  else ! skip it
    do j = 1, 12
      read(10,*) ! skip
    enddo
  endif
enddo
close(10)

allocate( pos(3,com%maxOurNpts)[*] ) !allocate 
allocate( pene(OurNpts) )
end subroutine readCOMfile !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readMyConfFile(fil) !{{{
implicit none
character(4) :: cf
integer :: i
character(80) :: fil

write(cf,'(i4.4)') myImg
open(9,file=trim(fil)//'_'//cf//'.dd')
do i = 1, myNpts
  read(9,*) cf, pos(:,i)
enddo
close(9)

end subroutine readMyConfFile !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! communicate the position on the neighbouring image to this one as a ghost layer
subroutine getGhostPositions !{{{
implicit none
integer :: i, my(2), ur(2), urID

! loop over each block of ranges to get
do i = 1, com%Nget
  my = com%getrange(:,i)
  ur = com%sourcerange(:,i)
  urID = com%sourceCellID(i)
  pos(:, my(1):my(2) ) = pos(:, ur(1):ur(2) )[urID]
!  write(6,'(i2,A,2i5,A,2i5,A,i2,A)') MyImg, "pos(:, ",my," ) = pos(:, ",ur," )[",urID,"]"
enddo

end subroutine getGhostPositions !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readXYZfile( fil ) !{{{
implicit none
integer :: i
character(*) :: fil
character(2) :: ch

open(7,file=trim(fil))
read(7,*) myNpts
ourNpts = myNpts ! for Nim == 1

allocate( pos(3,ourNpts)[*] )
allocate( pene(OurNpts) )

read(7,*) BoxSize
do i = 1, myNpts
  read(7,*) ch, pos(:,i)
enddo
close(7)
end subroutine readXYZfile !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeMyConf !{{{
implicit none
character(4) :: ctmp
integer :: i
write(ctmp,'(i4.4)') myImg
open(10+myImg, file='dmp_'//ctmp//'.conf')

do i = 1, myNpts
  write(10+myImg,*) "dd", pos(:,i), pene(i), 1
enddo
do i = myNpts+1, ourNpts
  write(10+myImg,*) "dd", pos(:,i), 0, 0
enddo

close(10+myImg)
end subroutine writeMyConf !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate Fake forces based on relative particle positions
subroutine calculateForces( rc ) !{{{
implicit none
real(8), intent(in) :: rc
integer :: i, j, nb
real(8) :: dist, rij(3), ene

PotEnergy = 0.0d0

pene(:) = 0.0d0

if (NUM_IMAGES().gt.1) then
   nb = com%bulkrange(2)
else ! single image
   nb = myNpts
endif

!do i = 1, myNpts
do i = 1, nb
  ! loop over neighs for Newton Pair calculations
  do j = i+1, myNpts 
    rij = pos(:,j) - pos(:,i)
    if (PBC) then
       rij = rij/BoxSize ! normalize
       rij = rij -anint(rij) !fold
       rij = rij*BoxSize
    endif
    dist = sqrt(dot_product(rij, rij))
    if (dist.lt.rc) then
      ene = dist**2 ! some kind of function of distance, LOL
      pene(i) = pene(i) +ene ! some kind of function of distance, LOL
      pene(j) = pene(j) +ene ! some kind of function of distance, LOL
      PotEnergy = PotEnergy +2.0*ene ! some kind of function of distance, LOL
    endif
  enddo
enddo

SYNC ALL ! ensure all ghost neighbours got received

do i = nb+1, myNpts
  ! loop over neighs for Newton Pair calculations
  do j = i+1, myNpts 
    rij = pos(:,j) - pos(:,i)
    if (PBC) then
       rij = rij/BoxSize ! normalize
       rij = rij -anint(rij) !fold
       rij = rij*BoxSize
    endif
    dist = sqrt(dot_product(rij, rij))
    if (dist.lt.rc) then
      ene = dist**2 ! some kind of function of distance, LOL
      pene(i) = pene(i) +ene ! some kind of function of distance, LOL
      pene(j) = pene(j) +ene ! some kind of function of distance, LOL
      PotEnergy = PotEnergy +2.d0*ene ! some kind of function of distance, LOL
    endif
  enddo
  ! loop over Ghost neighbours without Newton Pair
  do j = myNpts+1, ourNpts
    rij = pos(:,j) - pos(:,i)
    if (PBC) then
       rij = rij/BoxSize ! normalize
       rij = rij -anint(rij) !fold
       rij = rij*BoxSize
    endif
    dist = sqrt(dot_product(rij, rij))
    if (dist.lt.rc) then
      ene = dist**2 ! some kind of function of distance, LOL
      pene(i) = pene(i) +ene ! some kind of function of distance, LOL
      PotEnergy = PotEnergy +ene ! some kind of function of distance, LOL
    endif
  enddo
enddo


end subroutine calculateForces !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module thingy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program testddCAF
use thingy
implicit none
real(8) :: rc ! cutoff radius
character(80) :: ctmp, fil

if (iargc().ne.3) call help

call getarg(1,ctmp); read(ctmp,*) rc
call getarg(2,ctmp); read(ctmp,*) PBC
call getarg(iargc(),fil) ! last argument is file name base.
!rc = 0.867 ! interparticle cutoff radius
!PBC = .true.

myImg = this_image()
Nim = NUM_IMAGES()

!write(0,*) myImg, "myImg of ", Nim; call flush(0)
if (Nim.gt.1) call readCOMfile( "dd.com", Nim )
!write(0,*) myImg, "read the communication file"; call flush(0)

if (Nim.gt.1) then
   call readMyConfFile( fil )
else
   !call readXYZfile( "disk.xyz" )
   !call readXYZfile( "L.xyz" )
   !call readXYZfile( "Tri.xyz" )
   !call readXYZfile( "tube2.xyz" )
   call readXYZfile( fil )
endif
!write(0,*) myImg, "read the local configuration file"; call flush(0)

SYNC ALL ! must sync to be sure other images are done reading their data.
! need to share the positions across the images in the ghost layers first
if (Nim.gt.1) call getGhostPositions
!write(0,*) myImg, "got the ghost positions"; call flush(0)

call calculateForces( rc )
!write(0,*) myImg, "calculated forces"; call flush(0)

call writeMyConf

! total energy
write(6,*) MyImg, PotEnergy

SYNC ALL ! are they all done calculating PotEnergy?
call co_sum( PotEnergy, result_image=1 ) ! send result to img=1
if (myImg.eq.1) write(6,*) "total pot energy:", PotEnergy

end program testddCAF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine help
write(0,*) "Test program for the domdec program. This uses fortran coarrays."
write(0,*) "Usage: cafrun -np N testddCAF rc pbc FILE"
write(0,*) "  N    is the total number of images to run this program with."
write(0,*) "  rc   is the cutoff radius for the particle interactions, used in domdec"
write(0,*) "        to create the proper ghost domains for communication to neighbouring"
write(0,*) "        images."
write(0,*) "  pbc  Set to T for use of periodic boundary conditions, domdec must have -p"
write(0,*) "  FILE input file basename to read with _XXX.dd extension for each image."
write(0,*) "        for single image, FILE is just the file."
write(0,*) 
write(0,*) "Outputs"
write(0,*) " this program produces a set of files from each image: dmp_XXX.conf"
write(0,*) "  it contains the elemnt, position, energy, flag; for each particle,"
write(0,*) "  where flag is 1 if it is a local particle and 0 if a ghost particle."
STOP
end subroutine help
