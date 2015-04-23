!============
! Author: Elmo Tempel
! Date: 22.04.2015
!============
module fil_cat
	use constants
	use parameters
	use cylinders
	use galaxies
	use cubesort
	implicit none
	!
	type fil_spine_points
		real(rk),dimension(1:3):: p ! coordinates
		real(rk),dimension(1:3):: dvis,dfil ! orientation of visitmap and fil points
		real(rk):: vis,wvis ! visit and weighted visit map values
		real(rk):: dstr ! orientation strength
		real(rk):: rad ! filament approximate radius
	end type fil_spine_points
	type fil_spine
		integer:: npts ! nr of points in a spine:: fpt(1:npts)
		type(fil_spine_points),dimension(:),allocatable:: fpt
		integer,dimension(1:2):: ecode ! exit code for spine detection
	end type fil_spine
	!
	integer:: nfil ! nr of filament spines:: fil(1:nfil)
	type(fil_spine),dimension(:),allocatable:: fil ! filament spines
	!
contains
	!
	subroutine write_basic_catalogue_to_file(fname)
		implicit none
		character(len=*),intent(in):: fname
		integer:: i,n,iunit
		real(rk):: wmax
		open(newunit=iunit, file=fname)
		write(iunit,fmt='(A)') "# idfil npts idpt x y z dx dy dz vis wvis dstr rad ecode1 ecode2"
		wmax=0.0
		do n=1,nfil
			wmax=max(wmax, maxval(fil(n)%fpt(:)%wvis) )
		end do
		do n=1,nfil
			do i=1,fil(n)%npts
				write(iunit,fmt='(3I8,3F13.5,6F10.6,F12.6,2I4)') n,fil(n)%npts,i,fil(n)%fpt(i)%p,&
				fil(n)%fpt(i)%dvis,fil(n)%fpt(i)%vis,fil(n)%fpt(i)%wvis/wmax,fil(n)%fpt(i)%dstr,&
				fil(n)%fpt(i)%rad,fil(n)%ecode(1:2)
			end do
		end do
		close(iunit)
	end subroutine write_basic_catalogue_to_file
end module fil_cat