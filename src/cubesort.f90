!============
! Author: Elmo Tempel
! Date (last changed): 17.03.2013
!============
!
! some subroutines in this module were used only for testing
!
module cubesort
	use constants
	use quicksort
	implicit none
	private
	type datacube
        integer:: n=-1 ! number of points
        real(rk),dimension(:),allocatable:: x,y,z
		integer,dimension(:),allocatable:: idx ! index vector for data array
        real(rk):: xmin,ymin,zmin,xmax,ymax,zmax
		real(rk):: rsearch
        ! sorting index vectors
        integer,dimension(:),allocatable:: ivec_x
        integer,dimension(:,:),allocatable:: ivec_y
        integer,dimension(:,:,:),allocatable:: ivec_z
	end type datacube
	type(datacube),dimension(1:100):: dcube ! maximum nr of datasets is given here
	!
	public:: init_cubesort, free_cubesort, extract_points_index, extract_points_raw_index
	public:: extract_points_raw_index2
contains
	!
	subroutine init_cubesort(vecx,vecy,vecz,rsearch,ii,return_data)
		implicit none
		real(rk),dimension(:),intent(inout):: vecx,vecy,vecz
		real(rk),intent(in):: rsearch ! search radius for sorting
		integer,intent(out):: ii ! output index
		logical,intent(in),optional:: return_data ! return sorted data (default=.false.)
		integer:: i,n
		! find new index for datacube
		ii=0
		do i=1,size(dcube)
			if (dcube(i)%n==-1) then
				ii=i; exit
			end if
		end do
		if (ii==0) stop "Error: cubesort: no place for new dataset!"
		! allocate memory and store data
		n=size(vecx)
		dcube(ii)%n=n
		allocate(dcube(ii)%x(1:n),dcube(ii)%y(1:n),dcube(ii)%z(1:n))
		allocate(dcube(ii)%idx(1:n))
		forall(i=1:n) dcube(ii)%idx(i)=i
		dcube(ii)%x=vecx; dcube(ii)%y=vecy; dcube(ii)%z=vecz
		dcube(ii)%rsearch=rsearch
		! store the min and max values
        dcube(ii)%xmin=minval(vecx)
        dcube(ii)%ymin=minval(vecy)
        dcube(ii)%zmin=minval(vecz)
        dcube(ii)%xmax=maxval(vecx)
        dcube(ii)%ymax=maxval(vecy)
        dcube(ii)%zmax=maxval(vecz)
		! sort the data
		call sort_dcube(dcube(ii))
		!
		if (present(return_data) .and. return_data) then
			vecx=dcube(ii)%x(1:n)
			vecy=dcube(ii)%y(1:n)
			vecz=dcube(ii)%z(1:n)
		end if
	end subroutine init_cubesort
	!
	subroutine free_cubesort(ii)
		implicit none
		integer,intent(in):: ii
		dcube(ii)%n=-1
		deallocate(dcube(ii)%x,dcube(ii)%y,dcube(ii)%z,dcube(ii)%idx)
		deallocate(dcube(ii)%ivec_x,dcube(ii)%ivec_y,dcube(ii)%ivec_z)
	end subroutine free_cubesort
	!
	! return raw initial index vector.. useful when sorted data returned
	subroutine extract_points_raw_index(ii,p,rad,ivec,sphere)
		implicit none
		integer,intent(in):: ii
		real(rk),dimension(1:3),intent(in):: p
		logical,intent(in),optional:: sphere ! return data in a sphere (default=.true.)
        real(rk),intent(in):: rad ! centre and searching radius
        integer,dimension(:),allocatable,intent(out):: ivec ! output vector of galaxy indexes
		logical:: sph
		sph=.true.; if (present(sphere)) sph=sphere
		call extract_points_internal(dcube(ii),p,rad,ivec,sph)
	end subroutine extract_points_raw_index
	! return initial index vector
	subroutine extract_points_index(ii,p,rad,ivec,sphere)
		implicit none
		integer,intent(in):: ii
		real(rk),dimension(1:3),intent(in):: p
		logical,intent(in),optional:: sphere ! return data in a sphere (default=.true.)
        real(rk),intent(in):: rad ! centre and searching radius
        integer,dimension(:),allocatable,intent(out):: ivec ! output vector of galaxy indexes
		logical:: sph
		sph=.true.; if (present(sphere)) sph=sphere
		call extract_points_internal(dcube(ii),p,rad,ivec,sph)
		ivec=dcube(ii)%idx(ivec)
	end subroutine extract_points_index
	!
	subroutine extract_points_raw_index2(ii,p,rad,ivec,n,sphere)
		implicit none
		integer,intent(in):: ii
		real(rk),dimension(1:3),intent(in):: p
		logical,intent(in),optional:: sphere ! return data in a sphere (default=.true.)
        real(rk),intent(in):: rad ! centre and searching radius
        integer,dimension(:),intent(out):: ivec ! output vector of galaxy indexes
		integer,intent(out):: n
		integer,dimension(:),allocatable:: ivec0
		logical:: sph
		sph=.true.; if (present(sphere)) sph=sphere
		call extract_points_internal(dcube(ii),p,rad,ivec0,sph)
		n=size(ivec0)
		if (n>size(ivec)) stop "Error: please increase ivec size for galaxy extraction!"
		ivec(1:n)=ivec0
	end subroutine extract_points_raw_index2
	!
    ! return galaxies in a certain radius: return internal index (can be used if sorted data returned)
    subroutine extract_points_internal(dat,p,rad,ivec,sphere)
        implicit none
		type(datacube),intent(in):: dat
		real(rk),dimension(1:3),intent(in):: p
		logical,intent(in):: sphere ! return sphere, otherwise return cube
        real(rk),intent(in):: rad ! centre and searching radius
        integer,dimension(:),allocatable,intent(out):: ivec ! output vector of galaxy indexes
        integer,dimension(:),allocatable:: dumvec ! temporary vector
        integer:: ix,iy,iz
        integer:: i,j,k,di
        integer:: k1,k2,ii,kk
        real(rk):: dist2,rad2,x0,y0,z0
        if (allocated(ivec)) deallocate(ivec)
        ! some tests
		x0=p(1); y0=p(2); z0=p(3)
		if (p(1)<dat%xmin) x0=dat%xmin; if (p(1)>dat%xmax) x0=dat%xmax
		if (p(2)<dat%ymin) y0=dat%ymin; if (p(2)>dat%ymax) y0=dat%ymax
		if (p(3)<dat%zmin) z0=dat%zmin; if (p(3)>dat%zmax) z0=dat%zmax
        ! how many neighbouring volumes are needed
        di=floor( rad/dat%rsearch ) + 1
        rad2=rad**2
        ! indexes of xyz coordinates
        ix=floor( (x0-dat%xmin)/dat%rsearch ) + 1
        iy=floor( (y0-dat%ymin)/dat%rsearch ) + 1
        iz=floor( (z0-dat%zmin)/dat%rsearch ) + 1
        ! allocate temporary array
        kk=dat%ivec_x( min(size(dat%ivec_x)-1,ix+di)+1 ) - dat%ivec_x(max(1,ix-di))
        allocate(dumvec(1:kk))
        ! cycle over xyz coordinates
        kk=0
        do i=max(1,ix-di),min(size(dat%ivec_x)-1,ix+di)
            do j=max(1,iy-di),min(size(dat%ivec_y,dim=2)-1,iy+di)
                do k=max(1,iz-di),min(size(dat%ivec_z,dim=3)-1,iz+di)
                    k1=dat%ivec_z(i,j,k)
                    k2=dat%ivec_z(i,j,k+1)-1
                    if (k2-k1<0) cycle
                    do ii=k1,k2
                        if (sphere) then ! calculate the distances
                            dist2=(p(1)-dat%x(ii))**2 + (p(2)-dat%y(ii))**2 + (p(3)-dat%z(ii))**2
                            if (dist2<=rad2) then
                                kk=kk+1
                                dumvec(kk)=ii
                            end if
                        else ! return all indexes
                            kk=kk+1
                            dumvec(kk)=ii
                        end if
                    end do
                end do
            end do
        end do
        allocate(ivec(1:kk))
        ivec(1:kk)=dumvec(1:kk)
        deallocate(dumvec)
    end subroutine extract_points_internal
	!
    ! return galaxies in a certain radius: return internal index (can be used if sorted data returned)
    subroutine extract_points_crd_internal(dat,p,rad,ipts,sphere)
        implicit none
		type(datacube),intent(in):: dat
		real(rk),dimension(1:3),intent(in):: p
		logical,intent(in):: sphere ! return sphere, otherwise return cube
        real(rk),intent(in):: rad ! centre and searching radius
        real(rk),dimension(:,:),allocatable,intent(out):: ipts ! output vector of galaxy indexes
        real(rk),dimension(:,:),allocatable:: dumvec ! temporary vector
        integer:: ix,iy,iz
        integer:: i,j,k,di
        integer:: k1,k2,ii,kk
        real(rk):: dist2,rad2,x0,y0,z0
        if (allocated(ipts)) deallocate(ipts)
        ! some tests
		x0=p(1); y0=p(2); z0=p(3)
		if (p(1)<dat%xmin) x0=dat%xmin; if (p(1)>dat%xmax) x0=dat%xmax
		if (p(2)<dat%ymin) y0=dat%ymin; if (p(2)>dat%ymax) y0=dat%ymax
		if (p(3)<dat%zmin) z0=dat%zmin; if (p(3)>dat%zmax) z0=dat%zmax
        ! how many neighbouring volumes are needed
        di=floor( rad/dat%rsearch ) + 1
        rad2=rad**2
        ! indexes of xyz coordinates
        ix=floor( (x0-dat%xmin)/dat%rsearch ) + 1
        iy=floor( (y0-dat%ymin)/dat%rsearch ) + 1
        iz=floor( (z0-dat%zmin)/dat%rsearch ) + 1
        ! allocate temporary array
        kk=dat%ivec_x( min(size(dat%ivec_x)-1,ix+di)+1 ) - dat%ivec_x(max(1,ix-di))
        allocate(dumvec(1:3,1:kk))
        ! cycle over xyz coordinates
        kk=0
        do i=max(1,ix-di),min(size(dat%ivec_x)-1,ix+di)
            do j=max(1,iy-di),min(size(dat%ivec_y,dim=2)-1,iy+di)
                do k=max(1,iz-di),min(size(dat%ivec_z,dim=3)-1,iz+di)
                    k1=dat%ivec_z(i,j,k)
                    k2=dat%ivec_z(i,j,k+1)-1
                    if (k2-k1<0) cycle
                    do ii=k1,k2
                        if (sphere) then ! calculate the distances
                            dist2=(p(1)-dat%x(ii))**2 + (p(2)-dat%y(ii))**2 + (p(3)-dat%z(ii))**2
                            if (dist2<=rad2) then
                                kk=kk+1
								dumvec(1,kk)=dat%x(ii)
								dumvec(2,kk)=dat%y(ii)
								dumvec(3,kk)=dat%z(ii)
                            end if
                        else ! return all indexes
                            kk=kk+1
							dumvec(1,kk)=dat%x(ii)
							dumvec(2,kk)=dat%y(ii)
							dumvec(3,kk)=dat%z(ii)
                        end if
                    end do
                end do
            end do
        end do
        allocate(ipts(1:3,1:kk))
        ipts(1:3,1:kk)=dumvec(1:3,1:kk)
        deallocate(dumvec)
    end subroutine extract_points_crd_internal
	!
    ! sort datacube for fast selection
    subroutine sort_dcube(dat)
        implicit none
		type(datacube),intent(inout):: dat
		integer:: nx,ny,nz ! number of cube boxes in every direction
        integer:: i,k, k1,k2,j,l
        integer:: ii
        real(rk):: dum
        !
        ! find how many cubes are in every coordinate
        nx=floor( (dat%xmax-dat%xmin)/dat%rsearch ) + 1
        allocate(dat%ivec_x(1:nx+1))
        ny=floor( (dat%ymax-dat%ymin)/dat%rsearch ) + 1
        allocate(dat%ivec_y(1:nx+1,1:ny+1))
        nz=floor( (dat%zmax-dat%zmin)/dat%rsearch ) + 1
        allocate(dat%ivec_z(1:nx+1,1:ny+1,1:nz+1))
        !
        ! sort by x coordinate
        call qsort(dat%x,dat%y,dat%z,dat%idx)
        ! find the index vectors for x coordinate
        k=0
        dum=dat%xmin
        do i=1,dat%n
            if (dat%x(i)>=dum) then
                do ii=k+1,nx
                    dum=dat%xmin+ii*dat%rsearch
                    if (dat%x(i)<dum) exit
                    if (ii==nx) exit
                end do
                dat%ivec_x(k+1:ii)=i
                k=ii
                dum=dat%xmin+k*dat%rsearch
            end if
        end do
        k2=dat%n
        dat%ivec_x(k+1:nx+1)=k2+1
        !
        ! sort by y coordinate
        do i=1,nx
            k1=dat%ivec_x(i)
            k2=dat%ivec_x(i+1)-1
            if (k2-k1<0) then
                dat%ivec_y(i,:)=k1
                cycle
            end if
            call qsort(dat%y(k1:k2),dat%x(k1:k2),dat%z(k1:k2),dat%idx(k1:k2))
            dum=dat%ymin
            k=0
            do j=k1,k2
                if (dat%y(j)>=dum) then
                    do ii=k+1,ny
                        dum=dat%ymin+ii*dat%rsearch
                        if (dat%y(j)<dum) exit
                        if (ii==ny) exit
                    end do
                    dat%ivec_y(i,k+1:ii)=j
                    k=ii
                    dum=dat%ymin+k*dat%rsearch
                end if
            end do
            dat%ivec_y(i,k+1:ny+1)=k2+1
        end do
        !
        ! sort by z coordinate
        do i=1,nx
            do j=1,ny
                k1=dat%ivec_y(i,j)
                k2=dat%ivec_y(i,j+1)-1
                if (k2-k1<0) then
                    dat%ivec_z(i,j,:)=k1
                    cycle
                end if
                call qsort(dat%z(k1:k2),dat%y(k1:k2),dat%x(k1:k2),dat%idx(k1:k2))
                dum=dat%zmin
                k=0
                do l=k1,k2
                    if (dat%z(l)>=dum) then
                        do ii=k+1,nz
                            dum=dat%zmin+ii*dat%rsearch
                            if (dat%z(l)<dum) exit
                            if (ii==nz) exit
                        end do
                        dat%ivec_z(i,j,k+1:ii)=l
                        k=ii
                        dum=dat%zmin+k*dat%rsearch
                    end if
                end do
                dat%ivec_z(i,j,k+1:nz+1)=k2+1
            end do
        end do
    end subroutine sort_dcube
end module cubesort