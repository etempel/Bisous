!============
! Author: Elmo Tempel
! Date (last changed): 28.04.2016
!============
!
module interaction_term
    use constants
    use cylinders
    use galaxies
    implicit none
    !
contains
    !
    ! connection energy determines the filamentary network... it cannot be too weak
    !
    ! calculate the potential (interaction term) for cylinder interactions
    ! cyl%potint stores the potential of this cylinder
    ! ... res is the global change...takes into account the difference of connected cylinders
    function calc_interaction_potential_for_cylinder(cyl,icyl) result(res)
        implicit none
        real(rk):: res
        type(cylinder),intent(inout):: cyl ! input cylinder
        type(cylinder):: dumcyl ! dummy cylinder
        integer,dimension(:),intent(in):: icyl ! cylinder array of nearby cylinders
        integer:: i,k,nr_icyl,idc
        integer,dimension(1:2):: iend
        integer:: fint
        real(rk):: pot_change,dum
        !
        ! initialise interactions
        cyl%id_con=0; cyl%id_con_end=0; cyl%id_con_here=0; cyl%nr_con=0
		!
		if (ctest_no_interactions) then
			cyl%potint=0.0
			res=0.0
			return
		end if
        !
        ! calc interaction term
        k=0
        nr_icyl=size(icyl) ! close-by cylinders
        !
        do i=1,nr_icyl
            ! negative(-1) -- penalty/repulsive; 0 -- no connection; positive (+1) -- connection 
            ! iend(1:2) annab otsa, kus on interaktsioon: endpoint of 1. cylinder: endpoint of second cylinder
            fint=calc_interaction(cyl,cylarr(icyl(i)),iend)
            if (fint==-1) then ! repulsive cylinder
				! repulsive is forbidden...if allowed in some degree, change here
                res=cint_repulsive
                cyl%potint=res
				cyl%nr_con=0
                return ! no need to look forward
            else ! no connection or connected cylinder
                ! check the connection
                if (iend(1)/=0) then ! connection exists
                    if (fint/=1) stop "Error: calc interaction energy: illegal operation!!!"
                    k=k+1 ! increase the number of connections
                    cyl%id_con(k)=icyl(i)
                    cyl%id_con_end(k)=iend(2)
                    cyl%id_con_here(k)=iend(1)
                end if ! end connected cylinder
                !
            end if
        end do
        ! store the number of connections
        cyl%nr_con=k
        !
		! calculate the connection based on connected cylinders
        res=get_cylinder_pot_from_connections(cyl)
        !
        ! calc the change in potential, when taking into account newly formed connected cylinders
        if (cyl%nr_con>0) then
            pot_change=0.0
            do i=1,cyl%nr_con
                idc=cyl%id_con(i) ! connection with that cylinder
                dumcyl=cylarr(idc)
                !
                ! add the potential connection
                dumcyl%nr_con=dumcyl%nr_con+1
                k=dumcyl%nr_con
                dumcyl%id_con_end(k)=cyl%id_con_here(i)
                dumcyl%id_con_here(k)=cyl%id_con_end(i)
                ! manipulate the potential
                pot_change=pot_change - dumcyl%potint
                dum=get_cylinder_pot_from_connections(dumcyl)
                pot_change=pot_change + dum
            end do
            res=res+pot_change
        end if
    end function calc_interaction_potential_for_cylinder
    !
    ! Calculate potential change for removed cylinder... only hypothetical change (does not change anything)
    function calc_interaction_potential_for_removed_cylinder(cylid) result(res)
        implicit none
        real(rk):: res
        integer,intent(in):: cylid ! id of removed cylinder
        type(cylinder):: dumcyl ! dummy cylinder
        integer:: i,id,k,con
        logical:: test
        if (.not.cylarr(cylid)%inuse) stop "Error: pot for remov cyls, not in inuse!"
        res=0.0
		if (ctest_no_interactions) then
			return
		end if
        ! remove the potential from connections
        do i=1,cylarr(cylid)%nr_con
            id=cylarr(cylid)%id_con(i) ! id of connected cylinder
            dumcyl=cylarr(id)
			!
            res=res-dumcyl%potint ! later will be added if connections are corrected
            !
            test=.false.
            do1: do k=1,dumcyl%nr_con ! find the correct connection
                if (dumcyl%id_con(k)==cylid) then
                    con=k
                    test=.true.
                    exit do1
                end if
            end do do1
            if (.not.test) then
				print*, "ids ", cylid,id
				print*, cylarr(cylid)%nr_con
				print*, cylarr(cylid)%id_con
				print*, "dum ", dumcyl%nr_con
				print*, dumcyl%id_con
				stop "Error: removed cylinder potential... some problems occur!"
			end if
            ! con is the connection, that should be removed from dumcyl array
            ! .. remove this connection
            if (dumcyl%nr_con>1 .and. con/=dumcyl%nr_con) then
                dumcyl%id_con(con)=dumcyl%id_con(dumcyl%nr_con)
                dumcyl%id_con_end(con)=dumcyl%id_con_end(dumcyl%nr_con)
                dumcyl%id_con_here(con)=dumcyl%id_con_here(dumcyl%nr_con)
            end if
            dumcyl%nr_con=dumcyl%nr_con-1
            if (dumcyl%nr_con<0) stop "ERROR: remove cyl potential...no connection (neg)!"
            !
            ! add the potential
            res=res+get_cylinder_pot_from_connections(dumcyl)
        end do
        !
        ! remove the main cylinder potential
        res=res-cylarr(cylid)%potint
    end function calc_interaction_potential_for_removed_cylinder
    !
    ! get silinder potential from connections
    function get_cylinder_pot_from_connections(cyl) result(res)
        implicit none
        real(rk):: res
        type(cylinder),intent(inout):: cyl ! input cylinder
        integer:: i,k
        logical,dimension(1:2):: endpoint ! where is connection
        endpoint=.false. ! noconnections
        res=0.0
        do i=1,cyl%nr_con
            endpoint(cyl%id_con_here(i))=.true.
        end do
        k=count(endpoint) ! nr of edges where are connection
        !
        ! calculate the interaction term
        select case (k)
            case (0)
            res=-cint_conn_0
            case (1)
            res=-cint_conn_1
            case (2)
            res=-cint_conn_2
        end select
        cyl%potint=res
    end function get_cylinder_pot_from_connections
    !
	!
	! negative(-1) -- penalty/repulsive; 0 -- no connection; positive (+1) -- connection 
    ! iend(1:2) :: endpoint of first and second cylinder
    function calc_interaction(cyl1,cyl2,iend) result(res)
        implicit none
        integer:: res ! -1 repulsive, 0 no connection, 1 connection
        integer,dimension(1:2),intent(out):: iend ! 1-cyl1 end; 2- cyl2 end
        type(cylinder),intent(in):: cyl1,cyl2 !
        real(rk):: cosi,dist,dist1,dist2 !! distances for first,second endpoint for cylinder 1
        integer:: i1,i2 !! indicates the endpoint of second cylinder
        res=0; iend=0
        cosi=get_cosi_between_cyls(cyl1,cyl2)
		dist=get_dist_between_cyls(cyl1,cyl2)
        !
        ! check the repulsive ... radius difference cannot be too large
        if (cosi>cint_orthogonal_cosi) then ! not perpendicular...check the distances
			if (max(cyl1%r/cyl2%r,cyl2%r/cyl1%r)>cint_repulsive_radius_dif) then
				! radius difference is too large for repulsion and for connection
				res=0.0
				return
			end if
			! there are slight inaccuracies...but they are negligible (hardly occur)
            if (cyl1%h>cyl2%h) then
                dist2=(cyl1%h+cint_parallel_cosi*cyl2%h)
            else
                dist2=(cint_parallel_cosi*cyl1%h+cyl2%h)
            end if
            if ( dist<( 0.5*dist2 - cint_cyl_connection_radius) ) then
                res=-1
                return
            end if
        end if
        !
        ! calc the real connection
		if (cosi>cint_parallel_cosi .and. dist > 0.5*max(cyl1%h,cyl2%h) .and. max(cyl1%r/cyl2%r,cyl2%r/cyl1%r)<=cint_concyl_radius_dif) then
	        call get_dist_between_cyl_ends(cyl1,cyl2,dist1,dist2,i1,i2)
	        if (i1==i2 .and. dist1<cint_cyl_connection_radius) then
	            res=1
	            iend(1)=1; iend(2)=i1
	        else if (i1==i2 .and. dist2<cint_cyl_connection_radius) then
	            res=1
	            iend(1)=2; iend(2)=i2
	        end if
		end if
    end function calc_interaction
    !
end module interaction_term