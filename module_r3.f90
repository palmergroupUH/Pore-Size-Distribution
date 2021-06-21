module pore_size_calculation
implicit none

! Total number of atoms in simulation box and number of frames in trajectory
integer, public :: natom
integer, public :: nframe 

!Framewise particles positions and simulation box size
real(8), allocatable, protected, save, public :: rxi(:,:)         !framewise x-coordinates
real(8), allocatable, protected, save, public :: ryi(:,:)         !framewise y-cordinates
real(8), allocatable, protected, save, public :: rzi(:,:)         !framewise z-coordinates
real(8), allocatable, protected, save, public :: box_length(:,:)  !framewise box size

real(8), public :: box_size(3) !simulation box size in a frame
real(8), public :: ibox_size(3) !inverse of simulation box size in a frame
real(8), public :: rp(3)       !pore-size calculation location point
real(8), public :: r_hs        !hard-sphere distance to avoid overlap between probe particle and monomer
real(8), allocatable, public :: rpos(:,:) !monomers positions in a frame

! Variables for linked-list method
real(8), public :: dcell_init                   !initial cell-size
integer, protected, save, public :: tot_cell    !total number of cells
integer, protected, save, public :: ncell(3)    !number of cells in each direction
real(8), protected, save, public:: dcell(3)     !cell size in each direction
integer, allocatable, protected, save, public :: head(:,:,:)
integer, allocatable, protected, save, public :: list(:)
integer, allocatable, protected, save, public :: atom_id(:,:)

contains
!########################################################################

!subroutine for initial random seed

SUBROUTINE init_random_seed()
implicit none
integer :: i, n, clock
integer, allocatable:: seed(:)
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!##########################################################################

!Subroutine to read input text-file

subroutine read_inputfile(sigma_p, sigma_np, dcell_init, nmax, frame_interval, traj_file, output_file)
implicit none
real(8), intent(out) :: sigma_p, sigma_np
real(8), intent(out) :: dcell_init
integer, intent(out) :: nmax, frame_interval
character*20, intent(out) :: traj_file
character*20, intent(out) :: output_file

open(unit=14,file='input.txt',action='read')
read(14,*)
read(14,*) sigma_p, sigma_np
read(14,*)
read(14,*) dcell_init
read(14,*)
read(14,*) nmax
read(14,*)
read(14,*) frame_interval
read(14,*)
read(14,*) traj_file
read(14,*)
read(14,*) output_file
close(14)
return
end subroutine read_inputfile
!##################################################################

! subroutine to read trajectory file in dcd format

subroutine dcd_reader(filename)
implicit none
character(len=*), intent(in) :: filename
character*4 :: header
integer :: control(20), i
real(8) :: box(6)
real(4),allocatable :: r_temp(:)
control(:) = 0
open(14, file=trim(filename), status='old',form='unformatted')
read(14) header, control
read(14)
read(14) natom
nframe = control(1)

allocate( rxi(nframe,natom) )
allocate( ryi(nframe,natom) )
allocate( rzi(nframe,natom) )
allocate( box_length(nframe,3) )

allocate( r_temp(natom) )
do i = 1, nframe
        read(14) box
        box_length(i,1) = box(1)
        box_length(i,2) = box(3)
        box_length(i,3) = box(6)
        read(14) r_temp
        rxi(i,:) = dble(r_temp)
        read(14) r_temp
        ryi(i,:) = dble(r_temp)
        read(14) r_temp
        rzi(i,:) = dble(r_temp)
enddo
close(14)
return
end subroutine dcd_reader

!#########################################################################

!Linked-list initialization - allocation of head and list arrays

subroutine init_list()
implicit none
ncell(:) = floor(box_size(:) / dcell_init)
if (any(ncell(:) < 3)) then
        print *, 'system is too small to use cell links'
        stop
endif
tot_cell = ncell(1) * ncell(2) * ncell(3)
dcell(:) = box_size(:) / dble(ncell(:))
allocate( head(0 : ncell(1)-1, 0 : ncell(2)-1, 0 : ncell(3)-1) )
allocate( list(natom) )
allocate( atom_id(3,natom) )
head(:,:,:) = 0
list(:) = 0
atom_id(:,:) = 0
return
end subroutine init_list
!##########################################################################

!Compute cell numbers of monomers

subroutine list_icell(ri,icell)
implicit none
real(8), intent(in) :: ri(3)
integer, intent(out) :: icell(3)

if (any(dabs(ri(:)/box_size(:)) > 0.5d0 )) then
        print *, 'atom not in the main-box'
        stop
end if

icell(:) = 0
icell(:) = floor( (ri(:)/box_size(:)+0.5d0) * dble(ncell(:)) )
icell(:) = modulo( icell(:), ncell(:) )
return
end subroutine list_icell
!###########################################################################

!Formation of head and list arrays

subroutine  link_list()
implicit none
real(8) :: ri(3)
integer :: i, cell_num, icell(3)

do i=1,natom
        ri(:) = rpos(:,i)
        icell(:) = 0
        call list_icell(ri,icell)
        list(i) = head(icell(1),icell(2),icell(3))
        head(icell(1),icell(2),icell(3)) = i
        atom_id(:,i) = icell(:)
enddo
return
end subroutine link_list
!######################################################################

!Check formation of head and list arrays

subroutine check_list()
implicit none
integer :: i, j, k, c, icell(3)
real(8) :: ri(3)
do i=1,natom
        ri(:) = rpos(:,i)
        icell(:) = 0
        call list_icell(ri,icell)
        if(any(icell(:) .ne. atom_id(:,i))) then
                print *, 'inconsistency1 found', i, icell, atom_id(:,i)
                stop
        endif
enddo
do i=0, ncell(1)-1
        do j=0, ncell(2)-1
                do k=0, ncell(3)-1
                        icell(:) = (/i,j,k/)
                        c = head(i,j,k)
                        do while ( c .ne. 0)
                                if(any( icell(:) .ne. atom_id(:,c))) then
                                        print *, 'inconsistency2 found', c, icell(:), atom_id(:,c)
                                        stop
                                endif
                                c = list(c)
                        enddo
                enddo
        enddo
enddo
return
end subroutine check_list
!#######################################################################

!Deallocation of arrays

subroutine finalize_list()
implicit none
deallocate(list)
deallocate(head)
deallocate(atom_id)
end subroutine finalize_list
!########################################################################

!Subroutine to compute pore size
subroutine fun(rin,f)
implicit none
real(8),intent(in) :: rin(3) 
real(8), intent(out) :: f
integer :: c, i, j, k, new_i, new_j, new_k, icell(3)
real(8) :: dist, rij(3), ri(3), pore_size, max_size

f=0.d0
max_size = dmax1(box_size(1),box_size(2),box_size(3))
pore_size = max_size

ri = rin
ri(:) = ri(:) - box_size(:) * dnint(ri(:)*ibox_size(:))

icell(:) = 0
call list_icell(ri,icell)

do i = icell(1)-1, icell(1)+1
        do j = icell(2)-1, icell(2)+1
                do k = icell(3)-1, icell(3)+1
                        new_i = modulo(i,ncell(1))
                        new_j = modulo(j,ncell(2))
                        new_k = modulo(k,ncell(3))
                        c = head(new_i,new_j,new_k)
                        do while(c .ne. 0)
                                rij = rpos(:,c) - ri
                                rij(:) = rij(:) - box_size(:)* dnint(rij(:)*ibox_size(:))
                                dist = dsqrt(sum(rij*rij))
                                if(dist .lt. r_hs) then
                                        return
                                else if (dist .lt. pore_size) then
                                        pore_size=dist     
                                end if
                                c = list(c)
                        enddo
                enddo
        enddo
enddo
if(pore_size .ge. dcell_init) then
        do i = 1,natom
                rij = rpos(:,i) - ri
                rij(:) = rij(:) - box_size(:) * dnint(rij(:)*ibox_size(:))
                dist = dsqrt(sum(rij*rij))
                if (dist .lt. pore_size) then
                        pore_size = dist
                endif
        enddo
endif
if(pore_size .gt. r_hs) f = -(pore_size-r_hs) * (pore_size-r_hs)
return
end subroutine fun
!############################################################################

!Subroutine to compute constraint value
subroutine func(rin,f)
implicit none
real(8),intent(in) :: rin(3)
real(8), intent(out) :: f
real(8) :: f_local, dr2
call fun(rin,f_local)
dr2 = sum((rp-rin)*(rp-rin)) + f_local
f = max(0.d0,dr2)
return
end subroutine func
!############################################################################

end module pore_size_calculation
