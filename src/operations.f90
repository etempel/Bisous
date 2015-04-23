!============
! Author: Elmo Tempel
! Date (last changed): 06.02.2015
!============
!
module operations
	use constants
	use galaxies
	use cylinders
	use parameters
	use cubesort
	use data_term
	use interaction_term
	implicit none
contains
	!
	! test point location in a grid...
	! ... sample volume mask should be implemented here...
    function test_sample_mask_location(p) result(res)
        implicit none
        logical:: res
        real(rk),dimension(1:3),intent(in):: p
        res=.true.
		! pxxx_active mask is by default...
		!if ( any(p<cmask%pmin_active) .or. any(p>cmask%pmax_active) ) res=.false.
		! ... apply here other masks if necessary...
    end function test_sample_mask_location
	!
	subroutine update_potential_for_cylinder(cyl,only_datapot)
		implicit none
		type(cylinder),intent(inout):: cyl
		logical,intent(in),optional:: only_datapot
        real(rk),dimension(:,:),allocatable:: ipts
		integer,dimension(:),allocatable:: icyl
		integer:: ngal
		real(rk):: dum,rad,rsearch
		rsearch=sqrt( (cyl%rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cyl%h)**2 ) ! cylinder diagonal
		call extract_galaxy_coordinates(cyl%p,ipts,ngal,rsearch)
		rad=0.5*(cyl%h+cdata_cyl_len_max) + cint_cyl_connection_radius
		call extract_cylinder_indexes(cyl%p,icyl,rad)
		!
		if (present(only_datapot) .and. only_datapot) then
			dum = calc_data_potential_for_cylinder(cyl,ipts,ngal)
			cyl%potint=cint_repulsive
		else
			dum=calc_cylinder_potential(cyl,ipts,ngal,icyl)
		end if
	end subroutine update_potential_for_cylinder
	!
	! update cylinder interaction potentials -- after log_gammas change
	subroutine update_cylinder_interaction_potentials()
		implicit none
		integer:: ii
		real(rk):: dum
		do ii=1,nr_cylarr
			if (cylarr(ii)%inuse) then
				! update interaction potential
				dum=get_cylinder_pot_from_connections(cylarr(ii))
			end if
		end do
	end subroutine update_cylinder_interaction_potentials
	!
    ! calc cylinder potential, it includes potential change due to the connection changes
    function calc_cylinder_potential(cyl,ipts,ngal,icyl) result(res)
        implicit none
        real(rk):: res,dum
        type(cylinder),intent(inout):: cyl
        real(rk),dimension(:,:),intent(in):: ipts
		integer,intent(in):: ngal
        integer,dimension(:),intent(in):: icyl
        !
        ! calculate data term
		res = calc_data_potential_for_cylinder(cyl,ipts,ngal)
		!
        if (res<cdata_pot_undefined) then
			dum=calc_interaction_potential_for_cylinder(cyl,icyl)
			res=res+dum
        else
            cyl%potint=cint_repulsive
        end if
		!print*, "end calc cyl potential"
    end function calc_cylinder_potential
    !
    ! add cylinder to cylinder array
	! ... and add all necessary connections
    subroutine add_cylinder_to_array(cyl,ii,old,fixed)
        implicit none
        type(cylinder),intent(in):: cyl ! added cylinder
        integer,intent(in):: ii ! cylinder location in cylarr
		! ... location ii is found outside of this subroutine...
		logical,intent(in),optional:: old ! is old cylinder?
		logical,intent(in),optional:: fixed ! add fixed cylinder... boundary condition
        integer:: i,j,k,kn,ic
        real(rk):: dum
		real(rk),dimension(1:3):: d
		logical,dimension(1:2):: fend
		logical:: done
		!print*, "start add cylinder to array ", ii,nr_cylarr,nr_cylarr_active+1
        ! increase the number of active cylinders
        nr_cylarr_active=nr_cylarr_active+1
        !
        ! test the length of allocated cylarr
        if (ii>nr_cylarr .or. ii>size(cylarr)) then
            stop "ERROR: cylinder array is too short!!!"
        end if
        !
        cylarr(ii)=cyl
		! if old cylinder, then do not nullify counters
		! ... effective only for cylinder change move in MCMC
		if (.not.(present(old) .and. old)) then
	        cylarr(ii)%nr_of_death=0
	        cylarr(ii)%nr_of_change=0
	        cylarr(ii)%nr_of_change_accept=0
		end if
        !
        ! cylinder is in array
        cylarr(ii)%inuse=.true.
		if (present(fixed)) cylarr(ii)%fixed=fixed
		!
		! find and store the cylinder in mask array (for fast finding)
		! find the location in a grid
		d=(cyl%p-cmask%pmin)/cmask%delta
		i=nint(d(1))+1
		j=nint(d(2))+1
		k=nint(d(3))+1
		! sanity check -- continue program execution
		if (min(i,j,k)<1 .or. i>cmask%nx .or. j>cmask%ny .or. k>cmask%nz) then
			cylarr(ii)%inuse=.false.
			nr_cylarr_active=nr_cylarr_active-1
			return
		end if
		!
		if (present(fixed) .and. fixed) then
			nr_cylarr_active=nr_cylarr_active-1
			nr_cylarr_fixed=nr_cylarr_fixed+1
			if (ii/=nr_cylarr_fixed) stop "Error: fixed index is wrong!!!"
		end if
		!
		! bookkeeping of cylinders for fast cylinder finding
		if (cmask%ii(i,j,k)==0) then
			do ic=1,cmask%nii+1
				if (cmask%iicyl(1,ic)==0) then
					cmask%ii(i,j,k)=ic
					exit
				end if
			end do
			cmask%nii=max(cmask%nii,ic)
			cmask%iicyl(1,ic)=ii
			cylarr(ii)%ix=i; cylarr(ii)%iy=j; cylarr(ii)%iz=k
			cylarr(ii)%ii=1
		else
			cylarr(ii)%ix=i; cylarr(ii)%iy=j; cylarr(ii)%iz=k
			done=.false.
			do ic=1,c_cmask_ncyl_for_cell
				if ( cmask%iicyl(ic,cmask%ii(i,j,k))==0 ) then
					done=.true.
					cmask%iicyl(ic,cmask%ii(i,j,k))=ii
					cylarr(ii)%ii=ic
					exit
				end if
			end do
			if (.not.done) then
				print*, "WARNING!!! WARNING!!! WARNING!!!"
				print*, "... too small array for cmaks_ncyl_for_cell"
				stop "Program terminated!!!!!"
			end if
		end if
        !
        ! add the connection to connected cylinders
        do i=1,cyl%nr_con
            k=cyl%id_con(i) ! interaction with that cylinder
            kn=cylarr(k)%nr_con+1
			cylarr(k)%nr_con=kn
			if (kn>nr_max_con) stop "Error: too many connections, increase array!"
            !
            cylarr(k)%id_con(kn)=ii
            cylarr(k)%id_con_end(kn)=cyl%id_con_here(i)
            cylarr(k)%id_con_here(kn)=cyl%id_con_end(i)
            ! calc/store the new interaction potential
            dum=get_cylinder_pot_from_connections(cylarr(k))
			!
			if (.not.cylarr(k)%fixed) then
				! set n0,n1,n2 global values
				select case (kn)
					case(1)
					nr_cylarr_n1=nr_cylarr_n1+1
					nr_cylarr_n0=nr_cylarr_n0-1
					case(2)
					nr_cylarr_n1=nr_cylarr_n1-1
					nr_cylarr_n2=nr_cylarr_n2+1
				end select
			end if
        end do
		if (.not.cylarr(ii)%fixed) then
			select case (cyl%nr_con)
				case(0); nr_cylarr_n0=nr_cylarr_n0+1
				case(1); nr_cylarr_n1=nr_cylarr_n1+1
				case(2:); nr_cylarr_n2=nr_cylarr_n2+1
			end select
		end if
		!
		!print*, "end add cylinder to array"
    end subroutine add_cylinder_to_array
    !
    ! remove the cylinder and all connections
	! ... clean the cmask
    subroutine remove_cylinder_from_array(id)
        implicit none
        integer,intent(in):: id !cylinder id to be removed
        integer:: i,ii,k,kk,k2,ic
        real(rk):: dum
        logical:: test
        if (.not.cylarr(id)%inuse) stop "Error: remove cylinder... cylinder does not exist!"
        nr_cylarr_active=nr_cylarr_active-1
        cylarr(id)%inuse=.false.
		if (cylarr(id)%fixed) stop "Error: cannot remove fixed cylinder!"
		!
		select case(cylarr(id)%nr_con)
			case(0); nr_cylarr_n0=nr_cylarr_n0-1
			case(1); nr_cylarr_n1=nr_cylarr_n1-1
			case(2:); nr_cylarr_n2=nr_cylarr_n2-1
		end select
        ! remove connections
        do i=1,cylarr(id)%nr_con
            ii=cylarr(id)%id_con(i)
            test=.false.
            indo1: do k=1,cylarr(ii)%nr_con
                if (cylarr(ii)%id_con(k)==id) then
                    kk=k
                    test=.true.
                    exit indo1
                end if
            end do indo1
            if (.not.test) stop "Error: remove cyl, no assumed connectin!"
            ! ii on cyl, millel con tuleb eemaldada ning kk on see connection, mis tuleb eemaldada
            if (kk==cylarr(ii)%nr_con) then
                ! last connection, just decrease the counter
                cylarr(ii)%nr_con=cylarr(ii)%nr_con-1
            else
                ! rearrange connections
                cylarr(ii)%id_con(kk)=cylarr(ii)%id_con(cylarr(ii)%nr_con)
                cylarr(ii)%id_con_end(kk)=cylarr(ii)%id_con_end(cylarr(ii)%nr_con)
                cylarr(ii)%id_con_here(kk)=cylarr(ii)%id_con_here(cylarr(ii)%nr_con)
                ! delete the last connection
                cylarr(ii)%nr_con=cylarr(ii)%nr_con-1
            end if
			if (cylarr(ii)%nr_con<0) stop "Error: remove cyl..connection error!"
            ! calculate/store new interaction term
            dum=get_cylinder_pot_from_connections(cylarr(ii))
			!
			! set n0,n1,n2 global values
			select case (cylarr(ii)%nr_con)
				case(0)
				nr_cylarr_n1=nr_cylarr_n1-1
				nr_cylarr_n0=nr_cylarr_n0+1
				case(1)
				nr_cylarr_n2=nr_cylarr_n2-1
				nr_cylarr_n1=nr_cylarr_n1+1
			end select
        end do
		!
		! remove it from grid mask
		ii=cmask%ii(cylarr(id)%ix,cylarr(id)%iy,cylarr(id)%iz)
		kk=cylarr(id)%ii; cylarr(id)%ii=0
		if (ii==0 .or. kk==0) stop "Error: cylinder not in a cmask array!"
		k2=0
		do ic=c_cmask_ncyl_for_cell,kk+1,-1
			if ( cmask%iicyl(ic,ii)>0 ) then
				k2=ic
				exit
			end if
		end do
		!
		if (k2==0) then
			cmask%iicyl(kk,ii)=0
			if (kk==1) then
				cmask%ii(cylarr(id)%ix,cylarr(id)%iy,cylarr(id)%iz)=0
			end if
		else
			! move the cylinder location
			cylarr(cmask%iicyl(k2,ii))%ii=kk
			cmask%iicyl(kk,ii)=cmask%iicyl(k2,ii)
			cmask%iicyl(k2,ii)=0
		end if
    end subroutine remove_cylinder_from_array
	!
	! extract galaxy coordinates close to the location p
	subroutine extract_galaxy_coordinates(p,ipts,ngal,rad)
		implicit none
		real(rk),dimension(1:3),intent(in):: p
		real(rk),dimension(:,:),allocatable,intent(out):: ipts
		integer,intent(out):: ngal
		real(rk),intent(in):: rad
		integer,dimension(:),allocatable:: ivec
		integer:: i
		call extract_points_raw_index(gal%ii_cubesort,p,rad,ivec,sphere=.true.)
		ngal=size(ivec)
		if (allocated(ipts)) deallocate(ipts)
		allocate(ipts(1:3,1:ngal))
		forall(i=1:ngal) ipts(1:3,i)=gal%p(1:3,ivec(i))
		deallocate(ivec)
	end subroutine extract_galaxy_coordinates
	!
	! extract cylinders close to the location p
	subroutine extract_cylinder_indexes(p,icyl,rad)
		implicit none
		real(rk),dimension(1:3),intent(in):: p
		integer,dimension(:),allocatable,intent(out):: icyl
		real(rk),intent(in):: rad
		integer:: i,j,k,i0,j0,k0,n,ii,c_ext_cyl_dim,ic
		integer,dimension(:),allocatable:: idum
		real(rk),dimension(1:3):: dd
		real(rk):: rad2,dist
		integer:: irad
		! radius in grid
		irad=ceiling(rad/maxval(cmask%delta))
		!
		! array for maxnr cylinders
		c_ext_cyl_dim=(c_cmask_ncyl_for_cell*(1+2*irad)**3)/2
		allocate(idum(1:c_ext_cyl_dim))
		!
		rad2=rad*rad
		! find the location in a grid
		dd=(p-cmask%pmin)/cmask%delta
		i0=nint(dd(1))+1
		j0=nint(dd(2))+1
		k0=nint(dd(3))+1
		!
		! extract cylinders
		n=0
		!
		do k=max(1,k0-irad), min(k0+irad,cmask%nz)
			do j=max(1,j0-irad), min(j0+irad,cmask%ny)
				do i=max(1,i0-irad), min(i0+irad,cmask%nx)
					!
					if (cmask%ii(i,j,k)>0) then
						do ii=1,c_cmask_ncyl_for_cell
							ic=cmask%iicyl(ii,cmask%ii(i,j,k))
							if (ic==0) exit
							if (cylarr( ic )%inuse) then
								dd=abs(p-cylarr(ic)%p)
								if (any(dd>rad)) cycle
								dist=sum( dd*dd )
								if (dist<=rad2) then
									n=n+1
									idum(n)=ic
								end if
							end if
						end do
					end if
					!
				end do
			end do
		end do
		! output the indexes
		if (allocated(icyl)) deallocate(icyl)
		allocate(icyl(1:n))
		icyl=idum(1:n)
		deallocate(idum)
	end subroutine extract_cylinder_indexes
	!
end module operations