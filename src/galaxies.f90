!============
! Author: Elmo Tempel
! Date (last changed): 21.01.2015
!============
!
module galaxies
    use constants
    use cubesort
    use parameters
	use utilities
    implicit none
    type galaxy_sample
        integer:: ngal ! number of galaxies
		integer:: ii_cubesort ! cubesort index
        real(rk),dimension(:,:),allocatable:: p ! coordinates(1:3 , 1:n)
        real(rk):: xmin,ymin,zmin,xmax,ymax,zmax ! min and max values in galaxy array
		! additional parameters... not actually used
		real(rk),dimension(:),allocatable:: lum ! heledused
    end type galaxy_sample
    !
    ! define a global galaxy type
    type(galaxy_sample):: gal
    !
contains
	!
	!==================================================
	! read and initialise galaxy data
	!==================================================
	!
	! read just coordinates: general subroutine
	! .. generate binary file... delete binary file manually for update
    subroutine read_points_xyzlum(filename,distlim)
        implicit none
		character(len=*),intent(in):: filename
		real(rk),intent(in),optional:: distlim ! maximum distance limit
        real(rk),dimension(1:4):: xdum
        integer:: n,iunit,n2,n1,ivar
		real(rk):: dist
		logical:: from_file,rlum
        print*, "Reading xyz(lum) galaxy sample... ", filename
		if (present(distlim)) print*, "Distance limit: ", distlim
		!
		! try to read galaxy file from binary file...automatically generated
		from_file=.false.
		open(newunit=iunit,file=filename//'.binary',status='old', form='unformatted',iostat=ivar)
		if (ivar==0) from_file=.true.
		close(iunit)
		!
		if (from_file) then
			open(newunit=iunit,file=filename//'.binary',status='old', form='unformatted',iostat=ivar)
			if (ivar/=0) stop "Error reading binary galaxy file!"
			read(iunit) gal%ngal
			print*, " ..valid points in file: ", gal%ngal
			call init_gal_sample(gal%ngal)
			read(iunit) gal%p
			read(iunit) rlum
			if (rlum) then
				read(iunit) gal%lum
			else
				gal%lum=1.0
			end if
			close(iunit)
		else
			! read three first columns (xyz)
			n1=0
			open(newunit=iunit,file=filename)
			do
				102 read(iunit,fmt=*,err=102,end=101) xdum(1:3)
				dist=sqrt(sum(xdum(1:3)**2))
				if ( .not.present(distlim) .or. (present(distlim).and.dist<=distlim) ) n1=n1+1
			end do
			101 continue
			close(iunit)
			! try to read four first columns (xyz lum)
			n2=0
			open(newunit=iunit,file=filename)
			do
				122 read(iunit,fmt=*,err=122,end=121) xdum(1:4)
				dist=sqrt(sum(xdum(1:3)**2))
				if ( .not.present(distlim) .or. (present(distlim).and.dist<=distlim) ) n2=n2+1
			end do
			121 continue
			close(iunit)
			!
			gal%ngal=n1
			print*, " ..valid points in file: ", gal%ngal
			call init_gal_sample(gal%ngal)
			!
			n=0
			open(newunit=iunit,file=filename)
			do
				if (n2==n1) then
					112 read(iunit,fmt=*,err=112,end=111) xdum(1:4)
				else
					113 read(iunit,fmt=*,err=113,end=111) xdum(1:3)
				end if
				dist=sqrt(sum(xdum(1:3)**2))
				if ( .not.present(distlim) .or. (present(distlim).and.dist<=distlim) ) then
					n=n+1
					gal%p(1:3,n)=xdum(1:3)*c_scale_factor
					if (n2==n1) then
						gal%lum(n)=xdum(4)
						rlum=.true.
					else
						rlum=.false.
						gal%lum(n)=1.0
					end if
				end if
			end do
			111 continue
			close(iunit)
			!
			open(newunit=iunit,file=filename//'.binary', form='unformatted')
			write(iunit) gal%ngal
			write(iunit) gal%p
			write(iunit) rlum
			if (rlum) then
				write(iunit) gal%lum
			end if
			close(iunit)
		end if
		!
		call init_gal_sample_after_reading()
    end subroutine read_points_xyzlum
	!
	!-----------------------------------------------------------
    ! initialise galaxy sample
    subroutine init_gal_sample(ngal)
        implicit none
        integer,intent(in):: ngal
        gal%ngal=ngal
        allocate(gal%p(1:3,1:ngal))
        allocate(gal%lum(1:ngal))
    end subroutine init_gal_sample
    !
    subroutine init_gal_sample_after_reading()
        implicit none
        gal%xmin=minval(gal%p(1,:)); gal%xmax=maxval(gal%p(1,:))
        gal%ymin=minval(gal%p(2,:)); gal%ymax=maxval(gal%p(2,:))
        gal%zmin=minval(gal%p(3,:)); gal%zmax=maxval(gal%p(3,:))
		!
		! initialise cubesort
		call init_cubesort(gal%p(1,:),gal%p(2,:),gal%p(3,:),c_cubesort_rsearch,gal%ii_cubesort,return_data=.true.)
    end subroutine init_gal_sample_after_reading
	!
	!
end module galaxies