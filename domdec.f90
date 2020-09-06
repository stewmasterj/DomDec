! Fri Aug 14 20:29:23 EDT 2020
! Sun Sep  6 13:43:01 EDT 2020
! vim:fdm=marker
! compile: gfortran heapSort.f90 -Wall -Wextra -Wuninitialized -fbounds-check ../ArgsMod/ArgsMod.f90 domdec.f90 -o domdec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module filMod

type cell
 integer :: n  !number of particles
 integer :: lvl ! recursion level
 character(8) :: nam !cell name
 real(4), dimension(3) :: BoxSize, psize ! cell size and material size
 real(4), dimension(2,3) :: extent
 real(4), dimension(:,:), allocatable :: pos, dat
 integer, dimension(:), allocatable :: id
 character(2), dimension(:), allocatable :: ele
 integer :: nnc ! number of neighbouring cells
 integer, dimension(50) :: ncl ! neighbour cell list
 !
 !integer, dimension(50) :: np2send ! number of points in cell to send to neigh cell
 !integer, dimension(50) :: np2get ! number of ghost points for cell to receive
 integer :: npsc, npgc ! number of permutations in psc & pgc (max 64)
 !  psc(0,:) collects how many particles are in a unique communication location
 !    and the 1:8 array holds what neighbour cells they are to be sent to.
 !  psc(-1,:) hold how many cells to send to.
 integer, dimension(-1:16,64) :: psc, pgc, pgcN ! various permutations of send and get
 integer, dimension(2) :: bulkrange  ! start and end particle IDs
 integer, dimension(:,:), allocatable :: sendrange, getrange, sourcerange
 integer, dimension(:), allocatable :: sourceCellID 
 ! some extra stuff
 integer :: Nim ! number of images to use
 logical :: PBC !either full PBC or none, so rare to have a mix
end type cell

type(cell) :: gcl  ! global cell
type(cell), dimension(:), allocatable :: cl ! local cells

contains
!!!!!!!!!!!!!!!!!!!!!! 
subroutine readXYZ( fil ) !{{{
implicit none
integer :: i
character(*) :: fil
character(2) :: ch
real(4), dimension(2,3) :: te

open(10,file=trim(fil))
read(10,*) gcl%n
allocate( gcl%pos(3,gcl%n), gcl%id(gcl%n) )
if (abs(product(gcl%BoxSize)).lt.1.e-12) then
   read(10,*) gcl%BoxSize
else
   read(10,*)
endif
! assuming the origin in the center of the BoxSize
gcl%extent(1,:) = -gcl%BoxSize*0.5
gcl%extent(2,:) =  gcl%BoxSize*0.5

do i = 1, gcl%n
  read(10,*) ch, gcl%pos(:,i)
  gcl%id(i) = i
enddo
close(10)

te(1,1) = minval(gcl%pos(1,:))
te(1,2) = minval(gcl%pos(2,:))
te(1,3) = minval(gcl%pos(3,:))
te(2,1) = maxval(gcl%pos(1,:))
te(2,2) = maxval(gcl%pos(2,:))
te(2,3) = maxval(gcl%pos(3,:))
gcl%psize = te(2,:) -te(1,:)

gcl%lvl = 0 ! recursion level 0
gcl%nam = "" !name

end subroutine readXYZ !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine cellsplit( c, d, o1, o2 ) !{{{
implicit none
integer, intent(in) :: d ! direction to sort and split across
type(cell) :: c, o1, o2
real(4) :: sp
real(4), allocatable, dimension(:) :: pp
integer :: i
real(4), dimension(2,3) :: tex

allocate( pp(c%n) )
! output: sorted positions and 2Boxsizes of sub domains
!
!   heap-sort the positions along the d direction
!     add partition at the halfway point of the sorted position array
! record the position of this partition. and two boxsize dimensions
pp = c%pos(d,:)
call hpsort_eps_epw( c%n, pp, c%id, 1.e-10 )

o1%n = c%n/2 ! intentional integer division
o2%n = c%n-o1%n
o1%lvl = c%lvl+1 ! increment recursion level
o2%lvl = c%lvl+1
o1%nam = trim(c%nam)//"1"
o2%nam = trim(c%nam)//"2"
if (allocated(o1%pos)) then;  deallocate( o1%pos, o1%id ); endif
if (allocated(o2%pos)) then;  deallocate( o2%pos, o2%id ); endif
allocate( o1%pos(3,o1%n), o1%id(o1%n) )
allocate( o2%pos(3,o2%n), o2%id(o2%n) )
do i = 1, o1%n
   o1%pos(:,i) = gcl%pos(:,c%id(i))
   o1%id(i) = c%id(i)
enddo
do i = 1, o2%n
   o2%pos(:,i) = gcl%pos(:,c%id(o1%n+i))
   o2%id(i) = c%id(o1%n+i)
enddo

! split location along d
sp = 0.5*( pp(o1%n) +pp(o1%n+1) )

! record the extents and boxsize of the two subdomains
o1%extent = c%extent
o1%extent(2,d) = sp
o1%BoxSize = c%BoxSize
o1%BoxSize(d) = o1%extent(2,d) -o1%extent(1,d)
tex(1,1) = minval(o1%pos(1,:))
tex(1,2) = minval(o1%pos(2,:))
tex(1,3) = minval(o1%pos(3,:))
tex(2,1) = maxval(o1%pos(1,:))
tex(2,2) = maxval(o1%pos(2,:))
tex(2,3) = maxval(o1%pos(3,:))
o1%psize = tex(2,:) -tex(1,:)

o2%extent = c%extent
o2%extent(1,d) = sp
o2%BoxSize = c%BoxSize
o2%BoxSize(d) = o2%extent(2,d) -o2%extent(1,d)
tex(1,1) = minval(o2%pos(1,:))
tex(1,2) = minval(o2%pos(2,:))
tex(1,3) = minval(o2%pos(3,:))
tex(2,1) = maxval(o2%pos(1,:))
tex(2,2) = maxval(o2%pos(2,:))
tex(2,3) = maxval(o2%pos(3,:))
o2%psize = tex(2,:) -tex(1,:)

end subroutine cellsplit !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine cellMsplit( c, d, m, o ) !{{{
implicit none
integer, intent(in) :: d, m ! direction to sort and number of splits
type(cell) :: c, o(m)
real(4) :: sp(m)
real(4), allocatable, dimension(:) :: pp
integer :: i, j, n
real(4), dimension(2,3) :: tex
character(4) :: ctmp

allocate( pp(c%n) )
! output: sorted positions and m-Boxsizes of sub domains
!
!   heap-sort the positions along the d direction
!     add partition at the halfway point of the sorted position array
! record the position of this partition. and two boxsize dimensions
pp = c%pos(d,:)
call hpsort_eps_epw( c%n, pp, c%id, 1.e-10 )

i = c%n/m
o(:)%n = i ! intentional integer division
i = c%n-m*i
if (i.lt.m) then ! ok
 do j = 1, i
  o(j)%n = o(j)%n +1 ! distribute remainder
 enddo
else
 write(0,*) "weird amount of particles left over after array split:",i,"of",m,"images"
endif

o(:)%lvl = c%lvl+1 ! increment recursion level
! allocate the split domains
do j = 1, m
 write(ctmp,'(i4)') j
 o(j)%nam = trim(adjustl(c%nam))//trim(adjustl(ctmp))
 if (allocated(o(j)%pos)) then;  deallocate( o(j)%pos, o(j)%id ); endif
 allocate( o(j)%pos(3,o(j)%n), o(j)%id(o(j)%n) )
enddo

n=0
do j = 1, m
 do i = 1, o(j)%n
   o(j)%pos(:,i) = gcl%pos(:,c%id(n+i))
   o(j)%id(i) = c%id(n+i)
 enddo
 n = n +o(j)%n
! split location along d
 if (j<m) sp(j) = 0.5*( pp(n) +pp(n+1) )
enddo

!do j = 1, m-1
!  sp(j) = 0.5*( pp(o(j)%n) +pp(o(j)%n+1) )
!enddo

! record the extents and boxsize of the two subdomains
do j = 1, m
  o(j)%extent = c%extent
  if(j>1) o(j)%extent(1,d) = sp(j-1)
  if(j<m) o(j)%extent(2,d) = sp(j)
  o(j)%BoxSize = c%BoxSize
  o(j)%BoxSize(d) = o(j)%extent(2,d) -o(j)%extent(1,d)
  tex(1,1) = minval(o(j)%pos(1,:))
  tex(1,2) = minval(o(j)%pos(2,:))
  tex(1,3) = minval(o(j)%pos(3,:))
  tex(2,1) = maxval(o(j)%pos(1,:))
  tex(2,2) = maxval(o(j)%pos(2,:))
  tex(2,3) = maxval(o(j)%pos(3,:))
  o(j)%psize = tex(2,:) -tex(1,:)
enddo

end subroutine cellMsplit !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine writeTopFile( fil ) !{{{
implicit none
character(*) :: fil
integer :: i, j
integer :: Nim
character(20) :: frm, cn, cll

Nim = gcl%Nim

write(cn,*) Nim
write(cll,*) len(trim(adjustl(cn)))+1

open(9,file=trim(fil))
! write global data
write(9,'(A)') "# [Nim, BoxSize(:)] \n [dir, min, max]"
write(9,*) Nim, gcl%BoxSize
do j = 1, 3
   write(9,*) j, gcl%extent(:,j)
enddo

! loop over each cell's data
write(9,'(A)') "# imID \n [dir, min, max] \n [Nneighcell, neighcellIDs(:)]"
do i = 1, Nim
 write(9,*) i ! image number (ID)
 do j = 1, 3
   write(9,*) j, cl(i)%extent(:,j)
 enddo
 frm = '(i'//trim(adjustl(cll))//')'
 write(9,frm) cl(i)%nnc
 write(cn,*) cl(i)%nnc
 frm = '('//trim(adjustl(cn))//'i'//trim(adjustl(cll))//')'
 write(9,frm) cl(i)%ncl(1:cl(i)%nnc)
enddo

close(9)
end subroutine writeTopFile !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeComFile( fil ) !{{{
implicit none
character(*) :: fil
integer :: i, Nim
integer, dimension(gcl%Nim) :: op
character(80) :: frm, word

Nim = gcl%Nim

open(8,file=trim(fil))

do i = 1, Nim
  op(i) = cl(i)%getrange(2,cl(i)%npgc)
enddo

write(8,*) maxval(op), "# max ourNpts across all images"
write(8,*) gcl%BoxSize, "# global BoxSize"
do i = 1, Nim
!!!!!!!!!!!!!!!
  write(8,*) i, cl(i)%n, op(i),  cl(i)%nnc, "# imID, myNpts, ourNpts, NNcells"
  write(word,'(i4)') cl(i)%nnc
  frm='('//trim(adjustl(word))//'i5,A)'
  write(8,frm) cl(i)%ncl(1:cl(i)%nnc), " # NeighCellIDs"
  !write(8,*) "nPts2send:",cl(i)%np2send(1:cl(i)%nnc)
  !write(8,*) cl(i)%np2get(1:cl(i)%nnc), "# NPts2get"

  write(8,*) "bulkminmax:",cl(i)%bulkrange
  write(8,*) cl(i)%npsc-1, "# NumPermutationsSendCells"
  write(word,'(i4)') cl(i)%npsc-1
  frm='(A,'//trim(adjustl(word))//'i10)'
  write(8,frm) "sendmin:",cl(i)%sendrange(1,:)
  write(8,frm) "sendmax:",cl(i)%sendrange(2,:)
  write(8,*) cl(i)%npgc, "# NumPermutationsGetCells"
  write(word,'(i4)') cl(i)%npgc
  frm='(A,'//trim(adjustl(word))//'i10)'
  write(8,frm) "getmin:",cl(i)%getrange(1,:)
  write(8,frm) "getmax:",cl(i)%getrange(2,:)
  frm='(A,'//trim(adjustl(word))//'i5)'
  write(8,frm) "sourceCellIDs:",cl(i)%sourceCellID(:)
  frm='(A,'//trim(adjustl(word))//'i10)'
  write(8,frm) "sourcemin:",cl(i)%sourcerange(1,:)
  write(8,frm) "sourcemax:",cl(i)%sourcerange(2,:)
enddo

call flush(8) ! why is this necessary?
close(8)


end subroutine writeComFile !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeXYZ( fil ) !{{{
! output the full set of points with each identified for their cell.
implicit none
character(*) :: fil
integer :: i, j
open(7,file=trim(fil))
write(7,*) gcl%n
write(7,*) gcl%BoxSize
do i = 1, gcl%Nim
   do j = 1, cl(i)%n
     write(7,*) "Cu",cl(i)%pos(:,j), cl(i)%id(j), i
   enddo
enddo 
close(7)
end subroutine writeXYZ !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findNeighcells( rc ) !{{{
implicit none
integer :: i, j, k, l
real(4), intent(in) :: rc
real(4), dimension(2,3) :: te, we ! temporary extent for jcell
logical :: match
!!!! find the neignbouring cells !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! not accounting for Periodic Boundary Conditions
! there's probably a faster way to do this with linked lists, but most
!  applications aren't going to be decomposing the domain beyond 1000.
! this is 0.5N^2, beware if trying to make OpenMP due to race conditions.
do i = 1, gcl%Nim
   cl(i)%nnc = 0 ! number of neighbour cells
enddo
!write(6,*) "# Nim, NumOfNeighCells, NeighCellID"
do i = 1, gcl%Nim
   ! extend jcell's extents by rc and see if its boundary is now inside icell.
   do j = i+1, gcl%Nim
     te(1,:) = cl(j)%extent(1,:) -rc ! move lower
     te(2,:) = cl(j)%extent(2,:) +rc ! move higher
     we = te
     if (gcl%PBC) then ! wrap extents?
       do l = 1, 3
         if (te(1,l).le.gcl%extent(1,l)) we(1,l) = te(1,l)+gcl%BoxSize(l)
         if (te(2,l).ge.gcl%extent(2,l)) we(2,l) = te(2,l)-gcl%BoxSize(l)
       enddo
     endif
     match = .false.
     do k = 1, 3
       match =            cl(i)%extent(1,k) < we(1,k) .and. cl(i)%extent(2,k) > we(1,k)
       match = match.or. (cl(i)%extent(1,k) < we(2,k) .and. cl(i)%extent(2,k) > we(2,k))
      if (match) then ! there is an extent that is within this cell, now it's 2D.
        l = 1+mod(k,3) ! next direction
        match=match.and.((cl(i)%extent(2,l) > we(1,l) .and. cl(i)%extent(1,l) < we(1,l)) &
        &  .or. (cl(i)%extent(2,l) > we(2,l) .and. cl(i)%extent(1,l) < we(2,l)) &
        &  .or. (cl(i)%extent(2,l) < te(2,l) .and. cl(i)%extent(1,l) > te(1,l)) )
        l = 1+mod(k+1,3) ! next direction
        match=match.and.((cl(i)%extent(2,l) > we(1,l) .and. cl(i)%extent(1,l) < we(1,l)) &
        &  .or. (cl(i)%extent(2,l) > we(2,l) .and. cl(i)%extent(1,l) < we(2,l)) &
        &  .or. (cl(i)%extent(2,l) < te(2,l) .and. cl(i)%extent(1,l) > te(1,l)) )
        if (match) exit ! end loop
      endif
     enddo
     if (match) then
      cl(i)%nnc = cl(i)%nnc + 1
      cl(j)%nnc = cl(j)%nnc + 1
      cl(i)%ncl(cl(i)%nnc) = j ! save this cell ID to i list
      cl(j)%ncl(cl(j)%nnc) = i ! save this cell ID to j list
     endif
   enddo
   !write(6,*) i, cl(i)%nnc, cl(i)%ncl(1:cl(i)%nnc)
enddo 
end subroutine findNeighcells !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
function inCellID( x ) !{{{
implicit none
real(4), dimension(3), intent(in) :: x
integer :: inCellID
integer :: i
inCellID = 0
do i = 1, gcl%Nim
  if (x(1) .lt. cl(i)%extent(1,1)) cycle
  if (x(1) .ge. cl(i)%extent(2,1)) cycle
  if (x(2) .lt. cl(i)%extent(1,2)) cycle
  if (x(2) .ge. cl(i)%extent(2,2)) cycle
  if (x(3) .lt. cl(i)%extent(1,3)) cycle
  if (x(3) .ge. cl(i)%extent(2,3)) cycle
  inCellID = i
  exit
enddo
end function inCellID !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
function inGhstCell( x, g, rc ) !{{{
! will need modification for periodic boundary conditions
implicit none
real(4), intent(in) :: rc
real(4), dimension(3), intent(in) :: x
integer, intent(in) :: g
logical :: inGhstCell, l(3)
real(4), dimension(2,3) :: te ! temporary extent for jcell
integer :: k

inGhstCell = .false.
te(1,:) = cl(g)%extent(1,:) -rc ! move lower
te(2,:) = cl(g)%extent(2,:) +rc ! move higher

if (gcl%PBC) then ! wrap extents?
  do k = 1, 3
    if (te(2,k).gt.gcl%extent(2,k).and.te(1,k).lt.gcl%extent(1,k)) cycle 
    if (te(2,k).gt.gcl%extent(2,k)) te(2,k) = te(2,k)-gcl%BoxSize(k)
    if (te(1,k).lt.gcl%extent(1,k)) te(1,k) = te(1,k)+gcl%BoxSize(k)
  enddo
endif

l = .false.
do k = 1, 3
 if (te(1,k).gt.te(2,k)) then ! it's wrapped
  !l(k) = x(k) < te(1,k) .or. x(k) >= te(2,k) ! too inclusive, full cells...
  l(k) = x(k) < te(2,k) .or. x(k) >= te(1,k) ! this works for many cells, not Nim=2
  !if (.not.l(k)) l(k) = x(k) < te(1,k) .and. x(k) >= te(2,k) ! try reverse?
 else
  !if (x(dir).lt.te(1,dir) .or. x(dir).ge.te(2,dir)) return
  l(k) = x(k) < te(2,k) .and. x(k) >= te(1,k)
 endif
enddo
if (l(1).and.l(2).and.l(3)) inGhstcell = .true.

end function inGhstCell !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
end module filMod 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
program domdec
use filMod
use ArgsMod
implicit none
integer :: Nim, i, j, rl, ns, nc, NNim(3), ndcol
integer, dimension(1) :: bd
character(80) :: fil, topfil, oxyz, finefil, comfil
type(cell) :: tc ! temporary cell domain
type(cell), dimension(:), allocatable :: tca ! temporary cell domain
logical :: found, PBC, verbose
real(4) :: rc ! cutoff radius
real(4) :: r, ra
! for comand line arguments
integer, dimension(1) :: intv1
integer, dimension(3) :: intv
logical, dimension(1) :: logv
real(4), dimension(1) :: realv1
real(4), dimension(3) :: realv
character(80), dimension(1) :: charv
! for phase II
character(2) :: ch
real(4), dimension(:), allocatable :: x
integer :: ic, jc, k, l, ORtype
character(80) :: ctmp, cerr
integer, dimension(16) :: sc ! send cell local IDs for a particle
integer, dimension(:), allocatable :: cnt

!!!! OPTIONS and DEFAULT PARAMETERS !!!!{{{
ORtype = 1 ! 1=ORB, 2=ORM
!call getarg(1,ch); read(ch,*) Nim
Nim = 2 ! number of final boxes (ideally a power of 2, esp for ORB)
!call getarg(2,ctmp); read(ctmp,*) rc
rc = 0.867 ! cutoff radius for interparticle interactions, determins ghost size
ndcol = 4 ! columns in file: XX, pos(3), dat(:)
!fil    = "../Ndyn_OMP/examples/test_larger/out.xyz" ! 64000
!fil    = "../Ndyn_OMP/examples/test/out.xyz" ! 1000
!fil    = "tube2.xyz" ! 1564 
!fil    = "disk.xyz" ! 1328
!fil    = "L.xyz" ! 768
!fil    = "Tri.xyz" ! 768
oxyz   = "dd-out.xyz" ! output xyz file with cell info [optional]
topfil = "dd.top" ! output topology file name
PBC = .false. ! periodic boundary conditions? 
verbose = .false. ! verbose output
! Phase II options
!finefil = trim(fil) ! high res model, needed for phase II only
comfil = "dd.com" ! output communication file name, for phase II
gcl%BoxSize = 0.0
!!!! END USER INPUTS !!!!
if (iargc().eq.0) call help

call getarg( 1, cerr )
if (trim(cerr).eq."-h".or.trim(cerr).eq."--help") call help

! last argument is the file name
call getarg(iargc(),finefil)
fil = trim(finefil) ! set file to coarse fil if not overwritten below.

! get number of images for ORB
call getOpt( "-n", 1, intv1, cerr )
if (cerr(1:1).eq.'0') then
  write(6,*) "Performing Orthogonal Recursive Bisection"
  Nim = intv1(1); ORtype=1
endif
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! get number of images for ORM
call getOpt( "-N", 3, intv, cerr )
if (cerr(1:1).eq.'0') then
  write(6,*) "Performing Orthogonal Recursive Multisection"
  NNim = intv
  Nim = PRODUCT(intv); ORtype=2
endif
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! get cutoff radius
call getOpt( "-r", 1, realv1, cerr )
if (cerr(1:1).eq.'0') rc = realv1(1)
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! set simulation model boxsize
call getOpt( "-b", 3, realv, cerr )
if (cerr(1:1).eq.'0') gcl%Boxsize = realv
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! set PBC=.true.
call getOpt( "-p", 0, logv, cerr )
if (cerr(1:1).eq.'0') PBC = .true.
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! set verbose=.true.
call getOpt( "-v", 0, logv, cerr )
if (cerr(1:1).eq.'0') verbose = .true.
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! set coarse file
call getOpt( "-f", 1, charv, cerr )
if (cerr(1:1).eq.'0') fil = trim(charv(1))
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

! set columns in input file to read
call getOpt( "-c", 1, intv1, cerr )
if (cerr(1:1).eq.'0') ndcol = intv1(1)
if (cerr(1:1).ne.'0'.and.cerr(1:1).ne.'4') write(0,*) cerr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70 

! these are just shared with the module filMod
gcl%Nim = Nim
gcl%PBC = PBC 
!}}}

!!! first read the positions of a skeleton model or coarse grid model
call readXYZ( fil ) ! populates a "global" cell

allocate( cl(Nim) ) ! allocate the total number of final cells
cl(:)%lvl = -1 ! set all cells recursion level
cl(:)%nam = "0" ! set all cells name

!!!! PHASE I !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
!!!! Orthogonal Recursive Bisection (ORB) !!!!!!!!!!!!!!!!!!!!!!!!!!70
if (ORtype == 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! {{{
! split the global cell first
bd = maxloc(gcl%psize)! mod(rl3)+1
call cellsplit( gcl, bd(1), cl(1), cl(2) )
nc = 2 ! number of currently occupied cells
rl = 1 ! set recursion level counter

! this is the population scheme being used
! C: 1, 2, 0, 0, 0, 0, 0, 0 	split(0, 1,2) 
! C: 11, 2, 12, 0, 0, 0, 0, 0 	split(1, 1,3)
! C: 11, 21, 12, 22, 0, 0, 0, 0	split(2, 2,4)
! C: 111, 21, 12, 22, 112, 0, 0, 0	split(1, 1,5)
! C: 111, 211, 12, 22, 112, 212, 0, 0	split(2, 2,6)
! C: 111, 211, 121, 22, 112, 212, 122, 0	split(3, 3,7)
! C: 111, 211, 121, 221, 112, 212, 122, 222	split(4, 4,8)
if (verbose) then
 write(6,*) "# Nim names, recursion level, bisection direction"
 do i = 1, Nim
   write(6,'(A)',advance="no") cl(i)%nam
 enddo
endif
! you could cycle through the bisection direction with mod(rl,3)+1
!  or based on the cell's longest edge
bd = maxloc(gcl%psize)! mod(rl,3)+1
if (verbose) write(6,*) rl, bd(1) !mod(rl,3)+1

! need to perform Nim-1 cellsplits
do ns = 2, Nim-1 ! loop cell splits
 ! what cell should we split?
100 continue ! if we need to redo this search for the given split
   found = .false. ! initialize this to false
   do i = 1, nc ! search for the next cell to split
     if (cl(i)%lvl.eq.rl) then ! if the cell's recursion_level==current_level
       found = .true. ! found a cell to split
       tc = cl(i) ! save this as temporary cell to split
       bd = maxloc(cl(i)%psize)! mod(rl3)+1
       call cellsplit( tc, bd(1), cl(i), cl(nc+1) )
       nc = nc+1 ! increment occupied cell index
       exit ! search loop
     endif
   enddo
   if (.not.found) then ! done with this recursion level, move to next.
     rl = rl+1 ! increment recursion level
     goto 100  ! search again at this new recursion level
   endif
   if (verbose) then
     !write(6,*) ns, i, nc-1
     do i = 1, Nim
       write(6,'(A)',advance="no") cl(i)%nam
     enddo
     write(6,*) rl, bd(1) !mod(rl,3)+1
   endif
enddo
!!!! END ORB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! }}}
!!!! Orthogonal Recursive Multisection (ORM) !!!!!!!!!!!!!!!!!!!!!!!70
elseif (ORtype == 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! {{{
! split the global cell first
bd = 1 !maxloc(gcl%psize)! mod(rl,3)+1
call cellMsplit( gcl, bd(1), NNim(bd(1)), cl(1:NNim(bd(1))) )
nc = NNim(bd(1)) ! number of currently occupied cells
rl = 1 ! set recursion level counter

allocate( tca(maxval(NNim)) )
! this is the population scheme being used
! C: 1, 2, 0, 0, 0, 0, 0, 0 	split(0, 1,2) 
! C: 11, 2, 12, 0, 0, 0, 0, 0 	split(1, 1,3)
! C: 11, 21, 12, 22, 0, 0, 0, 0	split(2, 2,4)
! C: 111, 21, 12, 22, 112, 0, 0, 0	split(1, 1,5)
! C: 111, 211, 12, 22, 112, 212, 0, 0	split(2, 2,6)
! C: 111, 211, 121, 22, 112, 212, 122, 0	split(3, 3,7)
! C: 111, 211, 121, 221, 112, 212, 122, 222	split(4, 4,8)
if (verbose) then
 write(6,*) "# Nim names, recursion level, bisection direction"
 do i = 1, Nim
  write(6,'(A)',advance="no") cl(i)%nam
 enddo
endif
! you could cycle through the bisection direction with mod(rl,3)+1
!  or based on the cell's longest edge
bd = maxloc(gcl%psize)! mod(rl,3)+1
if (verbose) write(6,*) rl, bd(1),nc !mod(rl,3)+1

! need to perform Nim-1 cellsplits
do ns = 2, NNim(1)*NNim(2) +NNim(1)+1 ! loop cell splits
 ! what cell should we split?
200 continue ! if we need to redo this search for the given split
   found = .false. ! initialize this to false
   do i = 1, nc ! search for the next cell to split
     if (cl(i)%lvl.eq.rl) then ! if the cell's recursion_level==current_level
       found = .true. ! found a cell to split
       tc = cl(i) ! save this as temporary cell to split
       bd = mod(rl,3)+1
       call cellMsplit( tc, bd(1), NNim(bd(1)), tca(1:NNim(bd(1))) )
       cl(i) = tca(1)
       if (nc+NNim(bd(1))-1.gt.Nim) then
         write(0,*) "rl:",rl,"bd:",bd(1),"Nim:",Nim,"v",nc+NNim(bd(1))-1
         write(0,*) "cname",cl(i)%nam
         STOP "ERROR: split cell overflow"
       endif
       cl(nc+1:nc+NNim(bd(1))-1) = tca(2:NNim(bd(1)))
       nc = nc+NNim(bd(1))-1 ! increment occupied cell index
       exit ! search loop
     endif
   enddo
   if (.not.found) then ! done with this recursion level, move to next.
     rl = rl+1 ! increment recursion level
     if (nc.lt.Nim)      goto 200  ! search again at this new recursion level
   endif
   if (verbose) then
     !write(6,*) ns, i, nc-1
     do i = 1, Nim
       write(6,'(A)',advance="no") cl(i)%nam
     enddo
     write(6,*) rl, bd(1),nc,ns !mod(rl,3)+1
   endif
enddo
endif ! end ORtype
!!!! END ORM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! }}}

if (verbose) then !{{{
  ! write some topology stuff to screen
  write(6,*) "# Nim, Npts, Lx, Ly, Lz"
  do i = 1, Nim
     write(6,*) i, cl(i)%n, cl(i)%BoxSize
  enddo
  
  write(6,*) "# Extents: Nim, dir, min, max"
  do i = 1, Nim
   do j = 1, 3
     write(6,*) i, j, cl(i)%extent(:,j)
   enddo
  enddo
endif !}}}

! write an XYZ file that includes particle ID and cell ID for debugging
if (oxyz.ne."") call writeXYZ( oxyz )

!!!! find the neignbouring cells !!!!
call findNeighCells( rc )
if (verbose) then !{{{
  write(6,*) "# Nim, NumOfNeighCells, NeighCellID"
  do i = 1, Nim
    write(6,*) i, cl(i)%nnc, cl(i)%ncl(1:cl(i)%nnc)
  enddo
endif !}}}

!!!! Write the topology file: Nim, cell extents, cell neighbours !!!!!!!!!!!!!80
call writeTopFile( topfil )
!!!! END OF PHASE I !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
write(6,*) "Wrote Topology file: dd.top"
write(6,*) "Commencing Phase II for communication Data."

!!!! PHASE II !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
!!!!  use sorted points in cells to identify ghost points !!!!!!!!!!70
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
! If you're using a high resolution model, read it in here and
!   bin the particles into each cell, so that they can be identified.
! This whole routine reads the file input stream.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70 {{{
write(6,*) "Reading fine resolution file: "//trim(finefil)
open(10,file=trim(finefil))
read(10,*) gcl%n
read(10,*) !boxsize?

! reset all cell point numbers
do i = 1, Nim
   cl(i)%n = 0
   !cl(i)%np2send(:) = 0
   !cl(i)%np2get(:) = 0
   cl(i)%npsc = 1 ! cardinality of psc (first is always a bulk pattern)
   cl(i)%psc = 0 ! initialize the purmutations of send cells (sc)
   cl(i)%npgc = 0 ! cardinality of pgc
   cl(i)%pgc = 0 ! initialize the purmutations of get cells (gc)
   write(ctmp,'(i4.4)') i
   open(10+i,file=trim(ctmp)//'.dd') !open temp files for each image
enddo

allocate( x(ndcol-1) ) ! data columns in data file
! just counting points in ghost domains for later allocation, etc.
do i = 1, gcl%n
   read(10,*) ch, x
   ! what cell is x in?
   ic = inCellID(x(1:3))  !binning points into cells
   cl(ic)%n = cl(ic)%n + 1 ! add to cell's total points
   ns = 0
   sc = 0

   do j = 1, cl(ic)%nnc ! loop over neighbouring cell IDs
     jc = cl(ic)%ncl(j)
     ! for this cell pair, we extend the extent of cell j by +rc
     !  the points in cell i that are inside that domain are saved.
     ! function inGhstCell needs to be modified for periodic boundaries
     if (inGhstCell( x(1:3), jc, rc )) then ! will need to be sent to cell j
       !cl(ic)%np2send(j) = cl(ic)%np2send(j) + 1 ! increment points to send cell j
       ns = ns + 1
       sc(ns) = j ! save local cell neighbour ID
     endif
   enddo
   ! sort the sc array for convenience (it's short so this might work ok)
   if (ns.gt.1) then
     ! sort algorithm by John Mahaffy March 10, 1995
     do k = 1, ns
       bd = MAXLOC( sc(k:ns) )
       if (bd(1)+k-1.ne.k) then
         jc = sc(k)
         sc(k) = sc(bd(1)+k-1)
         sc(bd(1)+k-1) = jc
       endif
     enddo
   endif
   ! Need to collect all unique purmutations of sc(:)
   found = .false.
   outer: do j = 1, cl(ic)%npsc
     do k = 1, 16
       if (sc(k).ne.cl(ic)%psc(k,j)) cycle outer
       if (sc(k).eq.0) then ! completed the comparison without failure
         found=.true.
         nc = j ! save the purmutation ID
       endif
     enddo
   enddo outer
   if (found) then
     cl(ic)%psc(0,nc) = cl(ic)%psc(0,nc) + 1 ! increment particle count for this pattern
   else ! make new pattern
     cl(ic)%npsc = cl(ic)%npsc + 1 ! increment number of purmutations of sc
     nc = cl(ic)%npsc ! the purmutation ID
     cl(ic)%psc(1:ns,nc) = sc(1:ns) ! save this pattern
     cl(ic)%psc(0,nc) = cl(ic)%psc(0,nc) + 1 ! increment particle count for this pattern
     cl(ic)%psc(-1,nc) = ns
   endif
   ! print the point Name, position(1:3), and purmutationID
   write(10+ic,*) ch, x, nc
enddo !}}}

!!!! counting points in subdomain blocks !!!! 
do i = 1, Nim !{{{
  !close(10+i) ! close the temp files
  call flush(10+i) ! flush them for debugging?
  REWIND(10+i) ! Rewind the temp files for rereading
  ! populate the arrays for points from cells to GET
  do j = 1, cl(i)%nnc ! loop over neighbour cells
    jc = cl(i)%ncl(j) ! neighbour cell ID
    do k = 1, cl(jc)%npsc ! loop over neighbour cell's send permutations
      do l = 1, cl(jc)%psc(-1,k) ! loop over send cell IDs in permutation
        if (i.eq.cl(jc)%ncl(cl(jc)%psc(l,k))) then ! if this perm sends to our Icell
          cl(i)%npgc = cl(i)%npgc + 1 ! increment total perm count
          cl(i)%pgc(-1,j) = cl(i)%pgc(-1,j) + 1 ! increment perm count
          cl(i)%pgcN(-1,j) = cl(i)%pgcN(-1,j) + 1 ! increment perm count
          nc = cl(i)%pgc(-1,j)
          cl(i)%pgc(0,j) = cl(i)%pgc(0,j) +cl(jc)%psc(0,k)  ! add the particle count 
          cl(i)%pgc(nc,j) = k ! save the local neigh cell's perm ID
          cl(i)%pgcN(nc,j) = cl(jc)%psc(0,k)  ! copy the particle count 
        endif
      enddo
    enddo
  enddo
  ! allocate the convinient arrays
  allocate( cl(i)%sendrange(2,cl(i)%npsc-1) )
  allocate( cl(i)%getrange(2,cl(i)%npgc) )
  allocate( cl(i)%sourceCellID(cl(i)%npgc) )
  allocate( cl(i)%sourcerange(2,cl(i)%npgc) )

  cl(i)%bulkrange(:) = (/ 1, cl(i)%psc(0,1) /) ! particle count of first permutation
  ns = 0
  do j = 2, cl(i)%npsc
    ns = ns + cl(i)%psc(0,j-1) ! cumulative IDs
    cl(i)%sendrange(:,j-1) = (/ ns+1, ns+cl(i)%psc(0,j) /) ! min and max
  enddo
  ns = ns + cl(i)%psc(0,cl(i)%npsc) ! cumulative IDs
  nc = 0
  do j = 1, cl(i)%nnc
    do k = 1, cl(i)%pgc(-1,j) ! loop perm count
      nc = nc + 1 ! should sum to cl(i)%npgc
      cl(i)%getrange(:,nc) = (/ ns+1, ns+cl(i)%pgcN(k,j) /) ! min and max
      cl(i)%sourceCellID(nc) = cl(i)%ncl(j) ! global cell ID
      ns = ns + cl(i)%pgcN(k,j) ! cumulative IDs
    enddo
  enddo
enddo

do i = 1, Nim
  nc = 0
  do j = 1, cl(i)%nnc ! local neighbour cell index
    jc = cl(i)%ncl(j) ! neighbour cell ID
    do k = 1, cl(i)%pgc(-1,j) ! loop perm count
      nc = nc + 1 ! should sum to cl(i)%npgc
      l = cl(i)%pgc(k,j) ! returns the local neigh cell's perm ID.. I think...
      cl(i)%sourcerange(:,nc) = cl(jc)%sendrange(:,l-1) ! min and max
    enddo
  enddo
enddo !}}}

!write(6,*) "Done counting particles in Cells and Ghost domains"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
call writeComFile( comfil )
write(6,*) "Wrote Communication file: dd.com"

if (verbose) then !{{{
  do i = 1, Nim
    write(6,*) cl(i)%npsc, "# NumPermutationsSendCells"
    write(6,*) "#PermID, Nsendcells, NptsInPerm, Permutation Pattern cells to be sent to"
    do j = 1, cl(i)%npsc
      write(6,'(i2,x,i2,x,i6,16i4)') j, cl(i)%psc(:,j)
    enddo
    write(6,*) "#NeiCellID, Ngetcells, NptsInPerm, Permutation Pattern IDs\nNpts to get"
    do j = 1, cl(i)%nnc
      write(6,'(i2,x,i2,x,i6,16i4)') j, cl(i)%pgc(:,j)
      write(6,'(i2,x,i2,x,i6,16i4)') j, cl(i)%pgcN(:,j)
    enddo
  enddo
  
  write(6,*) "# particle ratio: get/local"
endif

ra = 0.0; ic=0; jc=0
do i = 1, Nim
  j = cl(i)%npsc-1
  j = cl(i)%sendrange(2,j) ! num o particles local to i
  ic = ic + j ! collect for total sum
  r = real(j) ! num o particles local to i
  j = cl(i)%npgc
  j = cl(i)%getrange(2,j) -cl(i)%getrange(1,1) ! num o particles needed by i
  jc = jc + j ! collect for total sum
  r = real(j)/r ! local/get ratio
  write(6,*) i, r
enddo !}}}
ra = real(jc)/real(ic) ! total get/local ratio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!70
! Now read in each cell's temporary file, and allocate all of their particle positions
! make sure to deallocate the memory from the previously read cell data.
! Sort the points acording to the cells that need the information.
!  i.e. Lists of particle IDs to be sent to which image...
! Ideally you want a position array to be like:
!    [ bulk particles,  particles to be sent,  particles to be received ]
! Particle ID ranges for each block: [ bulkrange; sendrange(:); getrange(:) ]
!  the 'get' blocks also have source cell IDs and source particle ID ranges

! reread the temp files into cell types ! {{{
do i = 1, Nim
  if (allocated(cnt))        deallocate( cnt )
  allocate( cnt( cl(i)%npsc ) )
  cnt = 0
  if (allocated(cl(i)%pos))  deallocate( cl(i)%pos, cl(i)%id )
  ns = cl(i)%npsc -1 
  nc = cl(i)%sendrange(2, ns ) ! particle ID of last in local particle set.
  allocate( cl(i)%pos(3, nc ), cl(i)%dat(ndcol-4, nc), cl(i)%ele(nc) )
  ! read the file lines
  do j = 1, cl(i)%n ! loop over all points in temp cell file
    read(10+i,*) ch, x, nc
    ! where to put this point a the appropriate offset?
    cnt(nc) = cnt(nc) + 1 ! increment this block count
    if (nc.eq.1) then ! bulk block
      l = cl(i)%bulkrange(1)+cnt(nc)-1
      cl(i)%ele(l) = ch
      cl(i)%pos(:, l ) = x(1:3) ! save position
      if (ndcol.gt.4) cl(i)%dat(:, l ) = x(4:ndcol-1)
    else ! in a send block
      l = cl(i)%sendrange(1,nc-1)+cnt(nc)-1
      cl(i)%ele(l) = ch
      cl(i)%pos(:, l ) = x(1:3) ! save position
      if (ndcol.gt.4) cl(i)%dat(:, l ) = x(4:ndcol-1)
    endif
  enddo ! end file line loop
  close(10+i, STATUS="DELETE") ! close and delete temp file

  ! open new file to write the sorted positions
  write(ctmp,'(i4.4)') i
  ! try to remove the temporary cell file here...
  open(10+i,file=trim(finefil)//'_'//trim(ctmp)//'.dd') !open model files for each image
  do j = 1, cl(i)%n
    write(10+i,*) cl(i)%ele(j), cl(i)%pos(:,j), cl(i)%dat(:,j)
  enddo
  close(10+i)
  deallocate( cl(i)%ele, cl(i)%pos, cl(i)%dat ) ! free some memory for the next cells
enddo !}}}
write(6,*) "Wrote Model files: "//trim(finefil)//"_XXXX.dd"

! this is the last line, as it is probably the most important metric.
write(6,*) ra," # total get/local ratio (you want to minimize this)"

end program domdec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine help !{{{
  write(0,*) "Domain Decomposition for non-local particle interactions"
  write(0,*) "Usage: domdec [OPTIONS] FILE"
  write(0,*) "Options:"
  write(0,*) "  -n i   ORB: Number of images or subdomains to create. [2] Default"
  write(0,*) "  -N i3  ORM: Number of images along each direction. [2]"
  write(0,*) "  -r r   Non-local cutoff radius. [0.867]"
  write(0,*) "  -p     Use Periodic Boundary conditions. [F]"
  write(0,*) "  -v     Verbose output data to stdout. [F]"
  write(0,*) "  -f s   Use coarse grid file s to decompose. [FILE]"
  write(0,*) "  -c i   Number of data columns in FILE. [4]"
  write(0,*) "  -b r3  Override model BoxSize or supply if not in FILE"
  write(0,*)
  write(0,*) "Output files include:"
  write(0,*) " dd.top        topology file, data for each subdomain."
  write(0,*) " dd.com        particle IDs to send and get from each subdomain."
  write(0,*) " dd-out.xyz    config showing what subdomain each particle is in."
!  write(0,*) " XXXX.dd       temporary file of points in each subdomain XXXX."
  write(0,*) " FILE_XXXX.dd final output files to be read by simulation program."
  write(0,*) "                this is not an XYZ file, just column data."
  write(0,*) "                the metadata is in dd.com"
  write(0,*) 
  write(0,*) "Input file format (XYZ)"
  write(0,*) "  Number_of_points"
  write(0,*) "  LX LY LZ [optional]"
  write(0,*) "  XX x y z data"
  write(0,*) "  ..."
  write(0,*) "The first two lines are for the number of points and boxsize, "
  write(0,*) " then is a list of the particles, such as element name e.g. Cu"
  write(0,*) " with it's (x,y,z) coordinate and any other data about it, such"
  write(0,*) " charge or coordination, etc."
  write(0,*)
  write(0,*) "Notes:"
  write(0,*) " 1. ORB = Orthogonal Recursive Bisection"
  write(0,*) " 2. ORM = Orthogonal Recursive Multisection"
  write(0,*) " 3. for ORB (using option -n) best to use powers of 2 (e.g.: 1,2,4,8,16..)"
  write(0,*) " 4. if you can't use a power of 2, try using ORM, PRODUCT(i(3))=NumOImages"
  write(0,*) " 5. September 5th, 2020 by Ross J. Stewart"
  STOP
end subroutine help !}}}
