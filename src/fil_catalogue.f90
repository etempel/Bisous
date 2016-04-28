!============
! Author: Elmo Tempel
! Date: 28.04.2016
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
	type fil_spine_cubesort
		integer:: ii
		real(rk),dimension(:),allocatable:: x,y,z
		integer,dimension(:),allocatable:: ifil,ipt
	end type fil_spine_cubesort
	!
	integer:: nfil ! nr of filament spines:: fil(1:nfil)
	type(fil_spine),dimension(:),allocatable:: fil ! filament spines
	type(fil_spine_cubesort):: fidx ! filamend indexes
	!
contains
	!
	subroutine calc_spine_statistics_for_input_points(file_in,file_out)
		implicit none
		character(len=*),intent(in):: file_in,file_out
		integer:: uin,uout ! unit in and out
		integer:: i,n,nin,ncnt,ierr,j,i0,n0
		real(rk),dimension(1:3):: crd,dori
		real(rk):: vm,wm,dst,rad
		integer,dimension(:),allocatable:: ivec
		print*
		print*, "Calculate spine statistics for points..."
		print*, "   Filename: ", file_in
		open(newunit=uin,file=file_in, status='old', iostat=ierr, action='read')
		if (ierr/=0) then
			print*, "ERR: calc_filament_statistics_for_input_points"
			print*, "Input file does not exist or is corrupted: ", file_in
			stop "Stopping the program..."
		end if
		open(newunit=uout,file=file_out)
		write(uout,fmt='(A)') '# x y z vmap wmap d_str  dx dy dz   dist_fil idfil idpt npts fil_rad'
		nin=0; ncnt=1
		do
			301 continue
			if (nin>0) ncnt=ncnt+1 ! jatame valja faili paise
			read(uin,fmt=*,err=301,end=302) crd
			nin=nin+1 ! valid xyz coordinates
			!
			do j=1,10 ! ten times the maximum radius
				call extract_points_index(fidx%ii,crd,rad=cdata_cyl_rad_max*j,ivec=ivec,sphere=.true.)
				if (size(ivec)>0) exit
			end do
			!print*, "done ", size(ivec), nin; read*
			n=0; i=0; rad=huge(1.0); n0=0; i0=0
			do j=1,size(ivec)
				n=fidx%ifil(ivec(j))
				i=fidx%ipt(ivec(j))
				! fil(n)%fpt(i) -- filamendi punkt
				dst=sum( (crd-fil(n)%fpt(i)%p)**2 )
				if (dst<rad) then
					n0=n; i0=i; rad=dst
				end if
			end do
			n=n0; i=i0 ! see on lahim filamendi punkt
			dst=sqrt(rad) ! see on kaugus sellest punktist
			!
			if (n>0 .and. i>0) then
				write(uout,fmt='(3F13.5, 3F8.5, 3F12.8, 1F10.5,3I7,F10.5)') crd(1:3),fil(n)%fpt(i)%vis,fil(n)%fpt(i)%wvis,&
				fil(n)%fpt(i)%dstr,fil(n)%fpt(i)%dvis(1:3), dst,n,i,fil(n)%npts,fil(n)%fpt(i)%rad
			else
				write(uout,fmt='(3F13.5, 3F8.5, 3F12.8, 1F10.5,3I7,F10.5)') crd(1:3), 0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0,0,0,0.0
			end if
		end do
		302 continue; ncnt=ncnt-1
		!
		if (nin/=ncnt) then
			print*, "WARNING: number of input and output points do not match (2)!"
			print*, "in, out: ", nin, ncnt
		end if
		close(uin)
		close(uout)
	end subroutine calc_spine_statistics_for_input_points
	!
	subroutine read_filament_catalogue_from_file(fname)
		implicit none
		character(len=*),intent(in):: fname
		integer:: iunit,idfil,npts,idpt,n,i,n2,i2
		integer:: ierr,cnt,cnt0
		!
		! allocate array for filaments
		if (allocated(fil)) stop "ERR: fil array already allocated!!!"
		allocate(fil(1:craw_max_nr_of_fil_spines)); nfil=0
		!
		! read data and allocate arrays
		open(newunit=iunit, file=fname, status='old', iostat=ierr, action='read')
		if (ierr/=0) then
			print*, "ERR: read_filament_catalogue_from_file"
			print*, "Input file does not exist or is corrupted: ", fname
			stop "Stopping the program..."
		end if
		read(iunit,fmt=*)
		cnt0=0
		do
			read(iunit,fmt=*,end=101) idfil,npts,idpt
			if (idpt==1) allocate(fil(idfil)%fpt(1:npts))
			fil(idfil)%npts=npts
			nfil=idfil
			if (idpt==1) cnt0=cnt0+npts
		end do
		101 continue
		close(iunit)
		!
		print*
		print*, "After first reading, npts: ", cnt0
		allocate(fidx%x(1:cnt0),fidx%y(1:cnt0),fidx%z(1:cnt0),fidx%ifil(1:cnt0),fidx%ipt(1:cnt0))
		cnt=0
		! read data to array
		open(newunit=iunit, file=fname, status='old', iostat=ierr, action='read')
		if (ierr/=0) stop "ERR: big problem!!!"
		read(iunit,fmt=*)
		do n=1,nfil
			if (fil(n)%npts/=size(fil(n)%fpt)) stop "ERR: read filament catalogue (1)!"
			do i=1,fil(n)%npts
				read(iunit,fmt=*) n2,npts,i2,fil(n)%fpt(i)%p,&
				fil(n)%fpt(i)%dvis,fil(n)%fpt(i)%vis,fil(n)%fpt(i)%wvis,fil(n)%fpt(i)%dstr,&
				fil(n)%fpt(i)%rad,fil(n)%ecode(1:2)
				if (n2/=n) stop "ERR: read filament catalogue (2)!"
				if (npts/=fil(n)%npts .or. i2>fil(n)%npts) stop "ERR: read filament catalogue (3)!"
				cnt=cnt+1
				fidx%x(cnt)=fil(n)%fpt(i)%p(1)
				fidx%y(cnt)=fil(n)%fpt(i)%p(2)
				fidx%z(cnt)=fil(n)%fpt(i)%p(3)
				fidx%ifil(cnt)=n
				fidx%ipt(cnt)=i
			end do
		end do
		close(iunit)
		print*, "Nr of spines, points: ", nfil,cnt
		if (cnt/=cnt0) stop "ERR: read filament catalogue (4)!"
		call init_cubesort(fidx%x,fidx%y,fidx%z,cdata_cyl_rad_max*3,fidx%ii)
		print*, "done init cubesort for filament spines!"
	end subroutine read_filament_catalogue_from_file
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