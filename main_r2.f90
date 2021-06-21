program main
use pore_size_calculation
implicit none
integer :: i, j
integer :: nmax  !Maximum number of points for pore size calculations in a frame
integer :: frame_interval !Frame-interval for pore size calculations
real(8) :: f1, f2 !local variables
real(8) :: poresize_1, poresize_2 !Pore sizes from Torquato and Gubbins method, respectively
real(8) :: ri(3)   !temporary local pore size calculation point
real(8) :: sigma_p, sigma_np   ! Monomer and prober particle diameter
external :: null
character* 20 :: traj_file, output_file !Monomers trajectory and output filenames 

call read_inputfile(sigma_p, sigma_np, dcell_init, nmax, frame_interval, traj_file, output_file) 

r_hs = 0.5d0*(sigma_p+sigma_np) !Hard sphere distance to avoid overlap between monomer and probe particles

call init_random_seed()         ! Random seed initialization
call dcd_reader(traj_file)      ! Read trajectory file

open(unit=14,file=output_file,action='write')

do i=1, nframe, frame_interval
        box_size(:) = box_length(i,:) !Read box_size in a given frame
        ibox_size(:) = 1.d0/box_size(:)
        allocate(rpos(3,natom))
        rpos(:,:)=0.d0
        do j=1,natom  !Read monomer positions in a given frame
                rpos(1,j)=rxi(i,j)-box_size(1)*dnint(rxi(i,j)*ibox_size(1))
                rpos(2,j)=ryi(i,j)-box_size(2)*dnint(ryi(i,j)*ibox_size(2))
                rpos(3,j)=rzi(i,j)-box_size(3)*dnint(rzi(i,j)*ibox_size(3))
        enddo
        ! Linked list initialization, formation, and check
        call init_list()
        call link_list()
        call check_list()
        ! Pore size calculations on random points in the box
        do j = 1, nmax
                f1=0.d0
                ! select random point in simulation box and check overlap with monomers
                do while( f1 == 0.d0)
                        call random_number(rp(1))
                        call random_number(rp(2))
                        call random_number(rp(3))
                        rp(1)=0.5d0*box_size(1)*(2.d0*rp(1)-1.d0)
                        rp(2)=0.5d0*box_size(2)*(2.d0*rp(2)-1.d0)
                        rp(3)=0.5d0*box_size(3)*(2.d0*rp(3)-1.d0)
                        ri(:) = rp(:)
                        call fun(ri,f1)
                enddo
                poresize_1 = 2.d0*dsqrt(-f1)
                f2=0.d0
                call solvopt(3, ri, f2, fun, .false., null, .true., func, .false., null)
                poresize_2 = 2.d0*dsqrt(-f2)
                write(14,'(8f16.5)')  rp, ri, poresize_1, poresize_2
        enddo
        deallocate(rpos)
        call finalize_list()
enddo
close(14)
end program main
