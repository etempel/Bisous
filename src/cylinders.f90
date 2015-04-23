!============
! Author: Elmo Tempel
! Date (last changed): 21.01.2015
!============
!
! Module, defining the cylinders and procedures for them
module cylinders
    use constants
	use parameters
    implicit none
	! internal parameters for this module
	integer,parameter:: nr_max_con=6 ! max three connections at both ends. Usually only one
	! ... it cannot be variable, so, I'm storing it in this module.
    ! define the parameters for every cylinder
    ! ... only permanent variables, no allocatable arrays (there were memory problems? (in omp?))
    type cylinder
		logical:: fixed=.false. ! fixed boundary cylinder, cannot be removed or changed
        logical:: inuse=.true. ! active cylinder? -- internal parameter, overwrites inactive cylinders
        ! basic parameters to define the cylinder and connections
        real(rk),dimension(1:3):: p ! cylinder position, centre
        real(rk):: u,t  ! u:0..1   t:0..2pi -- cylinder orientation
        real(rk):: h,r ! length and radius of a cylinder
		! data potential rmin, rmax
		real(rk):: rmin,rmax
        ! parameters calculated from xyz,ut,hr. This params are changed automatically
        real(rk),dimension(1:3):: eu ! unit vector definition - determined by u,t
        real(rk),dimension(1:3):: endp1,endp2 ! cylinder endpoints (coordinates)
        ! connection info
        integer:: nr_con ! number of connected cylinders
        integer, dimension(1:nr_max_con):: id_con ! ids of connected cylinders
        integer, dimension(1:nr_max_con):: id_con_end ! endpoints of connected cylinders (1 or 2)
        integer, dimension(1:nr_max_con):: id_con_here ! endpoint of this cylinder (1 or 2)
		! parameters for fast cylinder finding
		integer:: ix,iy,iz ! cylinder coordinates in a bookkeeping matrix
		integer:: ii ! index in a bookkeeping matrix
        ! potential for cylinder: energy = potdata+potint : probability density = exp(-energy)
        real(rk):: potdata,potint
        real(rk),dimension(1:3):: potarr ! includes potential compontents... useful for testing
		real(rk):: pot_hyp,pot_con ! hypothesis and concentration potentials
        ! some parameters for statistics
        integer:: nr_of_death, nr_of_change, nr_of_change_accept
        integer:: ngal ! nr of galaxies in cylinder
        real(rk):: den ! density in silinder (if all galaxies are equal: den=ngal)
		! %den can be used for photometric redshifts or for galaxy weighting (not yet implemented)
    end type cylinder
	!
	type cylinder_raw
        real(rk),dimension(1:3):: p ! cylinder position, centre
        real(rk):: u,t  ! u:0..1   t:0..2pi -- cylinder orientation
        real(rk):: h,r ! length and radius of a cylinder
		real(rk),dimension(1:3):: eu ! unit vector definition - determined by u,t
		integer:: nr_con ! number of connected cylinders
		real(rk):: potdata
		real(rk):: dh ! this is additional length for cylinder
	end type cylinder_raw
	!
	! this mask is used for fast cylinder finding
	! every cell contains cylinder indexes
	type cylinder_mask
		integer,dimension(:,:,:),allocatable:: ii ! index in a iicyl array
		integer,dimension(:,:),allocatable:: iicyl ! cylinder indexes (1:ncyl_for_cell, 1:c_cylarr_size)
		! iicyl vector includes cyl indexes in every cell, max ncyl_for_cell cylinders in every cell
		real(rk),dimension(1:3):: pmin,pmax,delta
		real(rk),dimension(1:3):: pmin_active,pmax_active ! nendes vahemikes samplitakse silindreid...
		! pmin,pmax define the cylinder centre limits
		integer:: nx,ny,nz
		integer:: nii ! how many grid cells are active
	end type cylinder_mask
	type(cylinder_mask):: cmask
    !
    ! define a global cylinder array
    type(cylinder),dimension(:),allocatable:: cylarr
    integer:: nr_cylarr ! number of used cylinders
    integer:: nr_cylarr_active ! number of active cylinders ! this is different since mcmc deletes cylinders
	integer:: nr_cylarr_fixed ! nr of fixed cylinders in the beginning of cylarr
	integer:: nr_cylarr_n0,nr_cylarr_n1,nr_cylarr_n2
    !
contains
	!
	!
	! initialise fast cylinder finding and initial cylinder arrays
	subroutine init_cylarr(n,pmin,pmax)
		implicit none
		real(rk):: delta
		integer,intent(in):: n ! nr of maximum cylinders
		real(rk),dimension(1:3),intent(in):: pmin,pmax ! minimum and maximum crds, where cylinders are stored
		real(rk):: dum
		integer:: idum
		!
		delta=c_cmask_delta*cdata_cyl_len_min
		! allocate cylarr
		nr_cylarr=0; nr_cylarr_active=0; nr_cylarr_fixed=0
		nr_cylarr_n0=0; nr_cylarr_n1=0; nr_cylarr_n2=0
		allocate(cylarr(1:n))
		cylarr%inuse=.false.
		cylarr%fixed=.false.
		! initialise cylmask for indexes
		! ...Rcon should be added for proper treatment of boundary cylinders
		! ... and the difference between max and min length cylinders
		if (c_cmask_use_boundary_cyls) then
			cmask%pmin_active=pmin-2*cint_cyl_connection_radius-0.5*(cdata_cyl_len_max-cdata_cyl_len_min)
			cmask%pmax_active=pmax+2*cint_cyl_connection_radius+0.5*(cdata_cyl_len_max-cdata_cyl_len_min)
		else
			cmask%pmin_active=pmin
			cmask%pmax_active=pmax
		end if
		cmask%pmin=cmask%pmin_active
		cmask%pmax=cmask%pmax_active
		!
		! add region for boundary cylinders... two cylinders for connection treatment
		if (c_cmask_use_boundary_cyls) then
			cmask%pmin=pmin-2*cdata_cyl_len_max-2*cint_cyl_connection_radius
			cmask%pmax=pmax+2*cdata_cyl_len_max+2*cint_cyl_connection_radius
			print*, "Added region for boundary cylinders: ", real(2*cdata_cyl_len_max+2*cint_cyl_connection_radius)
			dum=2*cint_cyl_connection_radius+0.5*(cdata_cyl_len_max-cdata_cyl_len_min)
			dum=dum+sqrt(0.25*cdata_cyl_len_max**2+cdata_cyl_rad_max**2)+cint_cyl_connection_radius
			print*, "Added region for boundary galaxies:  ", real(dum)
		end if
		!
		idum=ceiling( (cmask%pmax(1)-cmask%pmin(1))/delta )
		cmask%nx=idum+1
		cmask%delta(1)=(cmask%pmax(1)-cmask%pmin(1))/real(cmask%nx-1,kind=rk)
		!
		idum=ceiling( (cmask%pmax(2)-cmask%pmin(2))/delta )
		cmask%ny=idum+1
		cmask%delta(2)=(cmask%pmax(2)-cmask%pmin(2))/real(cmask%ny-1,kind=rk)
		!
		idum=ceiling( (cmask%pmax(3)-cmask%pmin(3))/delta )
		cmask%nz=idum+1
		cmask%delta(3)=(cmask%pmax(3)-cmask%pmin(3))/real(cmask%nz-1,kind=rk)
		!
		allocate(cmask%ii(1:cmask%nx,1:cmask%ny,1:cmask%nz))
		cmask%ii=0
		allocate(cmask%iicyl(1:c_cmask_ncyl_for_cell,1:n))
		cmask%iicyl=0
		!
		cmask%nii=0
	end subroutine init_cylarr
    !
    !---------------------------------------------
    ! below are functions and subroutines which are needed for cylinders
    !---------------------------------------------
    !
    ! calc the cosi between cylinders..
    function get_cosi_between_cyls(cyl1,cyl2) result(res)
        implicit none
        type(cylinder),intent(in):: cyl1,cyl2 ! input cylinders
        real(rk):: res ! angle between these two cylinders (radians)
        res=dot_product(cyl1%eu(1:3),cyl2%eu(1:3))
        res=abs(res)
    end function get_cosi_between_cyls
	! returns the angle with correct sign
    function get_cosi_between_cyls_sign(cyl1,cyl2) result(res)
        implicit none
        type(cylinder),intent(in):: cyl1,cyl2 ! input cylinders
        real(rk):: res ! angle between these two cylinders (radians)
        res=dot_product(cyl1%eu(1:3),cyl2%eu(1:3))
    end function get_cosi_between_cyls_sign
    !
    ! get distance between cylinder centres
    function get_dist_between_cyls(cyl1,cyl2) result(res)
        implicit none
        type(cylinder),intent(in):: cyl1,cyl2 ! input cylinders
        real(rk):: res ! distance between cyls
		real(rk),dimension(1:3):: dum
		dum=cyl1%p-cyl2%p
		res=sqrt( dot_product(dum,dum) )
    end function get_dist_between_cyls
	! return the distance square
    function get_distsq_between_cyls(cyl1,cyl2) result(res)
        implicit none
        type(cylinder),intent(in):: cyl1,cyl2 ! input cylinders
        real(rk):: res ! distance between cyls
		real(rk),dimension(1:3):: dum
		dum=cyl1%p-cyl2%p
		res=dot_product(dum,dum)
    end function get_distsq_between_cyls
    !
    ! return distances between cylinder ending points
    subroutine get_dist_between_cyl_ends(cyl1,cyl2,dist1,dist2,i1,i2)
        implicit none
        type(cylinder),intent(in):: cyl1,cyl2 ! input cylinders
        real(rk),intent(out):: dist1,dist2 ! distance for first,second endpoint for cylinder 1
        integer,intent(out):: i1,i2 ! indicates the endpoint of second cylinder
        real(rk):: dum1,dum2
        ! the first endpoint of cylinder one
		dum1=dot_product(cyl1%endp1-cyl2%endp1,cyl1%endp1-cyl2%endp1)
		dum2=dot_product(cyl1%endp1-cyl2%endp2,cyl1%endp1-cyl2%endp2)
        if (dum1<dum2) then
            dist1=sqrt(dum1); i1=1
        else
            dist1=sqrt(dum2); i1=2
        end if
        ! the second endpoint of cylinder one
		dum1=dot_product(cyl1%endp2-cyl2%endp1,cyl1%endp2-cyl2%endp1)
		dum2=dot_product(cyl1%endp2-cyl2%endp2,cyl1%endp2-cyl2%endp2)
        if (dum1<dum2) then
            dist2=sqrt(dum1); i2=1
        else
            dist2=sqrt(dum2); i2=2
        end if
    end subroutine get_dist_between_cyl_ends
    ! return distance squares between cylinder ending points
    subroutine get_distsq_between_cyl_ends(cyl1,cyl2,dist1,dist2,i1,i2)
        implicit none
        type(cylinder),intent(in):: cyl1,cyl2 ! input cylinders
        real(rk),intent(out):: dist1,dist2 ! distance for first,second endpoint for cylinder 1
        integer,intent(out):: i1,i2 ! indicates the endpoint of second cylinder
        real(rk):: dum1,dum2
        ! the first endpoint of cylinder one
		dum1=dot_product(cyl1%endp1-cyl2%endp1,cyl1%endp1-cyl2%endp1)
		dum2=dot_product(cyl1%endp1-cyl2%endp2,cyl1%endp1-cyl2%endp2)
        if (dum1<dum2) then
            dist1=dum1; i1=1
        else
            dist1=dum2; i1=2
        end if
        ! the second endpoint of cylinder one
		dum1=dot_product(cyl1%endp2-cyl2%endp1,cyl1%endp2-cyl2%endp1)
		dum2=dot_product(cyl1%endp2-cyl2%endp2,cyl1%endp2-cyl2%endp2)
        if (dum1<dum2) then
            dist2=dum1; i2=1
        else
            dist2=dum2; i2=2
        end if
    end subroutine get_distsq_between_cyl_ends
    !
	! distsq_axis:  this is squared! distance square from cylinder axis
	! dist_cnt: distance along cylinder axis from centre...can be negative
    subroutine get_point_dist_from_cyl2(cyl,pts,distsq_axis,dist_cnt)
        implicit none
        type(cylinder),intent(in):: cyl ! input cylinder
        real(rk),dimension(1:3),intent(in):: pts ! input point coordinates
        real(rk),intent(out):: distsq_axis,dist_cnt ! distance from cylinder axis and cylinder centre
        real,dimension(1:3):: dpts ! vector to a point location from cyl centre
		dpts=pts-cyl%p
        ! calculate distance using cosi: distance along cylinder axis
        dist_cnt=dot_product(cyl%eu,dpts)
        distsq_axis=dot_product(dpts,dpts) - dist_cnt*dist_cnt
    end subroutine get_point_dist_from_cyl2
	! dist_axis:  distance from cylinder axis
	! dist_cnt: distance along cylinder axis from centre...can be negative
    subroutine get_point_dist_from_cyl(cyl,pts,dist_axis,dist_cnt)
        implicit none
        type(cylinder),intent(in):: cyl ! input cylinder
        real(rk),dimension(1:3),intent(in):: pts ! input point coordinates
        real(rk),intent(out):: dist_axis,dist_cnt ! distance from cylinder axis and cylinder centre
        real,dimension(1:3):: dpts ! vector to a point location from cyl centre
		dpts=pts-cyl%p
        ! calculate distance using cosi: distance along cylinder axis
        dist_cnt=dot_product(cyl%eu,dpts)
        dist_axis=sqrt( dot_product(dpts,dpts) - dist_cnt*dist_cnt )
    end subroutine get_point_dist_from_cyl
    !
    ! change the value of u,t. Automatically calculates new unit vector
    subroutine set_cyl_tu(cyl,t,u,update)
        implicit none
        type(cylinder),intent(inout):: cyl
        real(rk),intent(in):: u,t ! angle(0..2pi) and cos(theta): (0..1)
        logical,intent(in),optional:: update ! update outomated parameters (default=true)
        ! ... useful for changing many parameters simultaneously
        cyl%t=t; cyl%u=u
		if (present(update) .and. .not.update) return
        call update_cyl_automatic_params(cyl)
    end subroutine set_cyl_tu
    !
    ! change cylinder position and orientation
    subroutine set_cyl_p_tu(cyl,p,t,u,update)
        implicit none
        type(cylinder),intent(inout):: cyl
        real(rk),intent(in):: u,t ! angle(0..2pi) and cos(theta)
        real(rk),dimension(1:3),intent(in):: p
        logical,intent(in),optional:: update ! update outomated parameters (default=true)
        ! ... useful for changing many parameters simultaneously
        cyl%t=t; cyl%u=u
        cyl%p=p
		if (present(update) .and. .not.update) return
		call update_cyl_automatic_params(cyl)
    end subroutine set_cyl_p_tu
    !
    ! update automatic parameters/data for cylinders
    subroutine update_cyl_automatic_params(cyl)
        implicit none
        type(cylinder),intent(inout):: cyl
        real(rk):: dum
        ! update cylinder unit vector
		dum=sqrt(1.0_rk-cyl%u**2)
        cyl%eu(1)=cos(cyl%t)*dum
        cyl%eu(2)=sin(cyl%t)*dum
        cyl%eu(3)=cyl%u
        !
        ! update cylinder endpoint coordinates
        dum=0.5*cyl%h ! half length of a cylinder
		cyl%endp1=cyl%p+dum*cyl%eu
		cyl%endp2=cyl%p-dum*cyl%eu
    end subroutine update_cyl_automatic_params
    !
    !=====================================================
	! below are subroutines that are also needed
    !
    ! find angle (cosi) between cyl and arbitrary vector
    function get_cosi_between_cyl_and_vector(cyl,t,u) result(res)
        implicit none
        type(cylinder),intent(in):: cyl
        real(rk),intent(in):: t,u
        real(rk):: res,dum
        real(rk),dimension(1:3):: eu
        ! generate direction vector
		dum=sqrt(1.0_rk-u**2)
        eu(1)=cos(t)*dum
        eu(2)=sin(t)*dum
        eu(3)=u
        ! calc cosi
        res=dot_product(cyl%eu,eu)
        res=abs(res)
    end function get_cosi_between_cyl_and_vector
    !
    ! find new cyl coordinates based on t,u,h and endpoint iend
	! .. connected cylinder
    subroutine get_xyz_for_new_cylinder(cyl,iend,t,u,h,p)
        implicit none
        type(cylinder),intent(in):: cyl
        real(rk),intent(in):: t,u,h ! h is the length of new cylinder
        integer,intent(in):: iend ! iend is the endpoint of this cylinder
        real(rk),dimension(1:3),intent(out):: p
		real(rk),dimension(1:3):: p0,d
        real(rk):: dum, dist1,dist2
		! calc distance from arbitrary point to new cylinder centre
        dum=0.5*h*sqrt(1.0_rk-u**2)
        p(1)=cos(t)*dum
        p(2)=sin(t)*dum
        p(3)=0.5*h*u
        !
		! set the endpoint based on iend
        if (iend==1) then
			p0=cyl%endp1
        else if (iend==2) then
            p0=cyl%endp2
        else
            stop " ERROR: get coordinates for new cylinder: illegal endpoint"
        end if
        ! check the sign: assume acute angle
		d=p0+p-cyl%p
		dist1=dot_product(d,d)
		d=p0-p-cyl%p
		dist2=dot_product(d,d)
        if (dist1>dist2) then
			p=p0+p
        else
            p=p0-p
        end if
    end subroutine get_xyz_for_new_cylinder
    !
	!
	!=============================================================
	! below are subroutines for global cylinder array
	!=============================================================
    !
	! get the number of connections (for statistics)
    subroutine get_nr_of_connections(n0,n1,n2,pdat,psig,mrad,mlen)
        implicit none
        integer,intent(out):: n0,n1,n2
        real(rk),dimension(1:3),intent(out):: pdat ! data potential 0,1,2 for connected cyls
		real(rk),intent(out),optional:: psig ! sigma of data potential
		real(rk),intent(out),optional:: mrad,mlen ! median radius and length of cylinders
        integer:: i,k,nn
		real(rk),dimension(:),allocatable:: dat
		real(rk):: avg
		if (present(psig)) then
			allocate(dat(1:nr_cylarr))
		end if
		if (present(mrad)) then
			mrad=0.0; mlen=0.0
		end if
        n0=0; n1=0; n2=0; nn=0
        pdat=0.0
        do i=1,nr_cylarr
            if (cylarr(i)%inuse .and. .not.cylarr(i)%fixed) then
				nn=nn+1
                select case (cylarr(i)%nr_con)
                    case (0)
                        n0=n0+1
                        pdat(1)=pdat(1)+cylarr(i)%potdata
                    case (1)
                        n1=n1+1
                        pdat(2)=pdat(2)+cylarr(i)%potdata
                    case (2:)
                        n2=n2+1
                        pdat(3)=pdat(3)+cylarr(i)%potdata
                    case default; stop "Error: wrong place to be... get nr connections"
                end select
				if (present(psig)) then
					dat(nn)=cylarr(i)%potdata
				end if
				if (present(mrad)) then
					mrad=mrad+cylarr(i)%r
					mlen=mlen+cylarr(i)%h
				end if
            end if
        end do
		if (present(psig)) then
			avg=sum(dat(1:nn))/nn
			psig=sum( (dat(1:nn)-avg)**2 )/nn
			deallocate(dat)
		end if
		if (present(mrad)) then
			mrad=mrad/nn
			mlen=mlen/nn
		end if
    end subroutine get_nr_of_connections
end module cylinders







