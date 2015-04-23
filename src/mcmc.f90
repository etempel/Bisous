!============
! Author: Elmo Tempel
! Date (last changed): 02.04.2015
!============
!
module mcmc_module
    use constants
    use cylinders
    use galaxies
    use operations
    use parameters
    use data_term
    use interaction_term
    implicit none
    !
contains
    !
    ! run one mcmc cycle
    ! ... spacing time pattern extraction number, simulating temperature - these are fixed here
    subroutine run_one_mcmc_cycle(ncycle,temperature)
        implicit none
        integer,intent(in):: ncycle ! spacing time pattern extraction number
        real(rk),intent(in):: temperature ! simulating temperature
		!
        real(rk),dimension(:,:),allocatable:: ipts
		integer:: id,ngal,ii,k,operation
		type(cylinder):: cyl
		integer,dimension(:),allocatable:: icyl
		real(rk):: pot,prob,alpha,potr,rad
		integer:: id_free ! id of free cylinder
		integer:: i
		real(rk):: birth_rate,br0,dum,dumdatapot
		!
        id_free=0 ! initially, there are no free cylinders
		br0=(1.0-cmcmc_prop_connected)/(cmcmc_norm_volume*cmcmc_volume_multiplier) ! precalculate birth rate one component
        !
        ! tsykkel yle pattern extraction numbri
        do k=1,ncycle
            !
            call generate_operation_for_cylinder(operation,cyl,id,ipts,ngal)
			! operation:: 11,12 (birth); 20 (death); 30 (change)
			! cyl:: cylinder to be added
			! id:: id of removed or changed cylinder
			! ipts,ngal:: galaxy coordinates and nr of galaxies around that cylinder
            !
            !apply the operation
            select case (operation)
                !
                !============================================================
                case (11,12) ! cylinder birth_random, birth_connected
                !
				rad=0.5*(cyl%h+cdata_cyl_len_max) + cint_cyl_connection_radius
				call extract_cylinder_indexes(cyl%p,icyl,rad)
                !
                pot=calc_cylinder_potential(cyl,ipts,ngal,icyl)
                !
                ! check the potential
                if (pot<cdata_pot_undefined) then
					! calculate birth_rate: StoiGregMate05.pdf eq. 7.
					birth_rate=br0
					if (nr_cylarr_n0+nr_cylarr_n1>0 .and. .not.ctest_no_interactions) then
						dum=0.0_rk
						do i=1,cyl%nr_con
							if (cylarr(cyl%id_con(i))%nr_con<2) then
								dum=dum+1.0/(2.0-cylarr(cyl%id_con(i))%nr_con)
							end if
						end do
						! boundary correction is ignored...small differences
						! ... actually, boundaries should be corrected during data initialisation
						! ... sampling mask should be inside observed region
						! ... in this case, no boundary correctin is required
						birth_rate=birth_rate + cmcmc_prop_connected*dum/(c_connection_volume*real(nr_cylarr_n0+nr_cylarr_n1))
					end if
					!
                    prob=(cmcmc_adapted_prop_death/cmcmc_adapted_prop_birth)*exp(-pot/temperature)
					prob=prob/( (1+nr_cylarr_active)*birth_rate )
					!
                    ! accept the cylinder?
                    call random_number(alpha)
                    if (alpha<prob) then
						! find the free place for cylinder
                        if (id_free>0) then
                            ii=id_free
                            id_free=0
                        else
                            do ii=1,nr_cylarr+1
                                if (.not.cylarr(ii)%inuse) exit
                            end do
                        end if
                        if (ii>nr_cylarr) nr_cylarr=ii
						! add cylinder to new array
                        call add_cylinder_to_array(cyl,ii)
                        select case (operation)
                            case (11); cstat_nr_cyl_birth_random_accept=cstat_nr_cyl_birth_random_accept+1
                            case (12); cstat_nr_cyl_birth_connected_accept=cstat_nr_cyl_birth_connected_accept+1
                        end select
                    end if
                end if
                !
                !==================================================================
                case (20)  ! cylinder death
                !
                pot=calc_interaction_potential_for_removed_cylinder(id)-cylarr(id)%potdata
                !
				! calculate birth_rate: StoiGregMate05.pdf eq. 7.
				birth_rate=br0
				if (nr_cylarr_n0+nr_cylarr_n1>0) then
					dum=0.0_rk
					do i=1,cylarr(id)%nr_con
						if (cylarr(cylarr(id)%id_con(i))%nr_con-1<2) then
							dum=dum+1.0/(2.0-cylarr(cylarr(id)%id_con(i))%nr_con+1)
							if (cylarr(cylarr(id)%id_con(i))%nr_con==0) stop "Error: some problem with birth rate!!"
						end if
					end do
					birth_rate=birth_rate + cmcmc_prop_connected*dum/(c_connection_volume*real(nr_cylarr_n0+nr_cylarr_n1))
				end if
				!
                prob=(cmcmc_adapted_prop_birth/cmcmc_adapted_prop_death)*exp(-pot/temperature)
				prob=prob*(nr_cylarr_active)*birth_rate
                !
                ! accept the cylinder?
                call random_number(alpha)
                if (alpha<prob) then
                    call remove_cylinder_from_array(id)
                    id_free=id
                    cstat_nr_cyl_death_accept=cstat_nr_cyl_death_accept+1
                    select case (cylarr(id)%nr_con)
                        case (0); cstat_nr_cyl_death_0_accept=cstat_nr_cyl_death_0_accept+1
                        case (1); cstat_nr_cyl_death_1_accept=cstat_nr_cyl_death_1_accept+1
                        case (2:); cstat_nr_cyl_death_2_accept=cstat_nr_cyl_death_2_accept+1
                    end select
                else
                    cylarr(id)%nr_of_death=cylarr(id)%nr_of_death+1
                end if
                !
                !==================================================================
                case (30)   ! cylinder change
				dumdatapot=cdata_coeff_variance*cylarr(id)%pot_con + cdata_coeff_hypothesis*cylarr(id)%pot_hyp
				!
                potr=calc_interaction_potential_for_removed_cylinder(id)-dumdatapot
                call remove_cylinder_from_array(id)
                !
				rad=0.5*(cyl%h+cdata_cyl_len_max) + cint_cyl_connection_radius
				call extract_cylinder_indexes(cyl%p,icyl,rad)
				!
                pot=calc_cylinder_potential(cyl,ipts,ngal,icyl)
                !
                pot=pot+potr
                prob=exp(-pot/temperature)
                !
                ! accept the cylinder?
                call random_number(alpha)
                if (alpha<prob .and. pot<cdata_pot_undefined) then
					! accept the change
                    call add_cylinder_to_array(cyl,id,old=.true.)
                    cylarr(id)%nr_of_change=cyl%nr_of_change+1
                    cylarr(id)%nr_of_change_accept=cyl%nr_of_change_accept+1
					cstat_nr_cyl_change_accept=cstat_nr_cyl_change_accept+1
                else
					! reject the change
                    cyl=cylarr(id)
                    call add_cylinder_to_array(cyl,id,old=.true.)
                    cylarr(id)%nr_of_change=cyl%nr_of_change+1
                end if
                !
                !==================================================================
                case default
                stop "Error: no valid move!!!"
            end select
            !
        end do
        !
    end subroutine run_one_mcmc_cycle
    !
    !
    subroutine generate_operation_for_cylinder(operation,cyl,id,ipts,ngal)
        use constants, only:pi
        implicit none
        integer,parameter:: n_test=100000 ! it must be sufficiently large (does not affect program execution)
        type(cylinder),intent(out):: cyl  ! output cylinder together with parameters
        integer,intent(out):: operation,id  ! operation case and cylinder index
		! operation: 	11	random cylinder (added)
		!				12	random-connected cylinder (added)
		!				20	removed cylinder
		!				30	cylinder change
        real(rk),dimension(:,:),allocatable,intent(out):: ipts
		integer,intent(out):: ngal
		!
		logical:: olemas,test,free_end(1:2)
		real(rk):: rnd,cosi,dum,rr1,rr2
		real(rk),dimension(1:3):: p,vec
		integer:: i,k,irnd,ii,iend,i0
		real(rk),dimension(1:3,1:3):: rx,rz ! rotation matrixes
		real(rk):: rsearch
		!
		rx=0.0_rk; rx(2,2)=1.0_rk; rz=0.0_rk; rz(3,3)=1.0_rk
        olemas=.false.
		id=0
        ! birth, death or change?
        call random_number(rnd)
        if (rnd<cmcmc_adapted_prop_birth) then
            operation=10 ! birth
        else if (rnd>=cmcmc_adapted_prop_birth .and. rnd<cmcmc_adapted_prop_birth+cmcmc_adapted_prop_death) then
            operation=20 ! death
        else
            operation=30 ! change
        end if
        ! if birth... test connected birth
        if (operation==10) then
            call random_number(rnd)
            if (rnd<cmcmc_adapted_prop_connected) then
                operation=12 ! connected cylinder
            else
                operation=11 ! random cylinder
            end if
        end if
        ! if no cylinders...only random birth (for burn in phace)
        if (nr_cylarr_active==0) operation=11
        !------
        !
        ! apply the operation
        select case (operation)
            !
            case (11) ! cylinder birth
            !=== random cylinder
			do i0=1,10 ! try several times and then give up (usually only 1 time is required)
				call get_cylinder_random_rh_adaptive(cyl%h,cyl%rmin,cyl%rmax)
				rsearch=sqrt( (cyl%rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cyl%h)**2 ) ! cylinder diagonal
	            indo11a: do i=1,n_test
	                call random_number(p)
					p=cmask%pmin_active+p*(cmask%pmax_active-cmask%pmin_active)
	                ! check if the position is inside a mask
	                test=test_sample_mask_location(p)
	                if (test) then
						call extract_galaxy_coordinates(p,ipts,ngal,rsearch)
	                    if (ngal<cdata_minpts) cycle indo11a
						! find the orientation
	                    cyl%p=p
						! cycle over orientations helps to speed up calculations a little bit (tiny cheat)
	                    ori11a: do k=1,c_move_n_orientation
							call random_number(cyl%t); cyl%t=cyl%t*twopi
	                        call random_number(cyl%u)
	                        call update_cyl_automatic_params(cyl)
	                        dum=calc_data_potential_for_cylinder(cyl,ipts,ngal,loctest=.true.)
	                        if (dum<cdata_pot_undefined) then
	                            olemas=.true.
	                            cstat_nr_cyl_birth_random=cstat_nr_cyl_birth_random+1
	                            exit indo11a
	                        end if
	                    end do ori11a
	                end if
	            end do indo11a
				if (olemas) exit
			end do
			!
            !
            case (12) ! connected cylinder birth
            !=== connected random cylinder
            indo12: do i=1,min(n_test,10*nr_cylarr,100*nr_cylarr_active) ! to avoid unnecessary searching
                call random_number(rnd); irnd=floor(rnd*nr_cylarr)+1
                if (.not.cylarr(irnd)%inuse) cycle indo12
                free_end(1:2)=.true.
                do k=1,cylarr(irnd)%nr_con
                    free_end(cylarr(irnd)%id_con_here(k))=.false.
                end do
				! fix the free end if possible
                ii=count(free_end)
                select case (ii)
                    case (0)
                    cycle indo12
                    case (1)
                    if (free_end(1)) iend=1
                    if (free_end(2)) iend=2
                    case (2)
                    call random_number(rnd)
                    iend=1
                    if (rnd>0.5) iend=2
                end select
                ! .. irnd is cylinder
                ! .. iend is free end of that cylinder
				! get the new position for that cylinder: find the galaxies in that location
				call get_cylinder_random_rh_adaptive(cyl%h,cyl%rmin,cyl%rmax,cylarr(irnd)%r)
                call get_xyz_for_new_cylinder(cylarr(irnd),iend,cylarr(irnd)%t,cylarr(irnd)%u,cyl%h, p)
				if (any(p-cint_cyl_connection_radius<cmask%pmin_active) .or. any(p+cint_cyl_connection_radius>cmask%pmax_active)) cycle indo12
				! test the grid location
				test = test_sample_mask_location(p)
				if (test) then
	                ! extract galaxies
					rsearch=sqrt( (cyl%rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cyl%h)**2 ) + cint_cyl_connection_radius ! cylinder diag + a_con
					call extract_galaxy_coordinates(p,ipts,ngal,rsearch)
	                if (ngal<cdata_minpts) cycle indo12
	                !
					! initialise rotation matrixes
					rx(1,1)=cylarr(irnd)%u; rx(3,3)=rx(1,1)
					rx(1,3)=sqrt(1.0_rk-cylarr(irnd)%u**2); rx(3,1)=-rx(1,3)
					!
					rz(1,1)=cos(cylarr(irnd)%t); rz(2,2)=rz(1,1)
					rz(1,2)=-sin(cylarr(irnd)%t); rz(2,1)=-rz(1,2)
					! get data_term accepted orientation
	                ori12: do k=1,c_move_n_orientation
	                    ! find the orientation
						call random_number(cyl%u); cyl%u=cint_parallel_cosi+(1.0-cint_parallel_cosi)*cyl%u
						call random_number(cyl%t); cyl%t=cyl%t*twopi
						dum=sqrt(1.0-cyl%u**2)
				        vec(1)=cos(cyl%t)*dum
				        vec(2)=sin(cyl%t)*dum
				        vec(3)=cyl%u
						vec=matmul(rx,vec)
						vec=matmul(rz,vec)
						! vec:xyz to cyl:tu -- take account signs (atan2 definition)
						if (vec(1)>0.0) then
							cyl%t=atan(vec(2)/vec(1))
						else if (vec(2)>=0.0 .and. vec(1)<0.0) then
							cyl%t=atan(vec(2)/vec(1))+pi
						else if (vec(2)<0.0 .and. vec(1)<0.0) then
							cyl%t=atan(vec(2)/vec(1))-pi
						else if (vec(2)>0.0 .and. abs(vec(1))<100*epsilon(1.0) ) then
							cyl%t=pi/2
						else if (vec(2)<0.0 .and. abs(vec(1))<100*epsilon(1.0) ) then
							cyl%t=-pi/2
						else if (abs(vec(2))<100*epsilon(1.0) .and. abs(vec(1))<100*epsilon(1.0) ) then
							cyl%t=0.0
						end if
						if (cyl%t<0.0) cyl%t=twopi+cyl%t
						cyl%u=vec(3)
						if (cyl%u<0.0) then
							cyl%u=-cyl%u
							if (cyl%t<pi) then
								cyl%t=cyl%t+pi
							else
								cyl%t=cyl%t-pi
							end if
						end if
						! done atan2
	                    !
	                    call get_xyz_for_new_cylinder(cylarr(irnd),iend,cyl%t,cyl%u,cyl%h, p)
						! add random location in a sphere
						do
			                call random_number(vec); vec=vec*2-1.0
							if ( sum(vec**2)>1.0 ) cycle
							vec=vec*cint_cyl_connection_radius
							exit
						end do
	                    cyl%p=p+vec
	                    call update_cyl_automatic_params(cyl)
	                    !
	                    dum=calc_data_potential_for_cylinder(cyl,ipts,ngal,loctest=.true.)
	                    if (dum<cdata_pot_undefined) then
	                        olemas=.true.
	                        cstat_nr_cyl_birth_connected=cstat_nr_cyl_birth_connected+1
	                        exit indo12
	                    end if
	                end do ori12
				end if
            end do indo12
            !
            !
            case (20) ! cylinder death
            !=== cylinder death
            indo20: do i=1,min(n_test,5*(nr_cylarr-nr_cylarr_fixed),100*nr_cylarr_active)
	            call random_number(rnd); irnd=floor(rnd*(nr_cylarr-nr_cylarr_fixed))+1+nr_cylarr_fixed
	            if (.not.cylarr(irnd)%inuse .or. cylarr(irnd)%fixed) cycle indo20
                olemas=.true.
				cstat_nr_cyl_death=cstat_nr_cyl_death+1
				select case (cylarr(irnd)%nr_con)
				    case (0); cstat_nr_cyl_death_0=cstat_nr_cyl_death_0+1
				    case (1); cstat_nr_cyl_death_1=cstat_nr_cyl_death_1+1
				    case (2:); cstat_nr_cyl_death_2=cstat_nr_cyl_death_2+1
			    end select
                exit indo20
            end do indo20
			if (olemas) then
	            id=irnd
		        cyl=cylarr(id)
			end if
            !
            !
            case (30) ! cylinder change
            !=== cylinder change
			do i0=1,10 ! try several times and then give up (usually only one is required)
	            indo30a: do i=1,min(n_test,5*(nr_cylarr-nr_cylarr_fixed),10*nr_cylarr_active)
		            call random_number(rnd); irnd=floor(rnd*(nr_cylarr-nr_cylarr_fixed))+1+nr_cylarr_fixed
		            if (.not.cylarr(irnd)%inuse .or. cylarr(irnd)%fixed) cycle indo30a
	                olemas=.true.
	                exit indo30a
	            end do indo30a
	            if (.not.olemas) cycle ! precaution
	            id=irnd
	            cyl=cylarr(id)
	            !
				! get new length and radius
				select case (cyl%nr_con)
					case (0)
					call get_cylinder_random_rh_adaptive_change(cyl%h,cyl%rmin,cyl%rmax,cyl%h,cyl%r)
					case (1)
					rr1=cylarr( cyl%id_con(1) )%r
					call get_cylinder_random_rh_adaptive_change(cyl%h,cyl%rmin,cyl%rmax,cyl%h,rr1)
					case (2:)
					rr1=cylarr( cyl%id_con(1) )%r
					rr2=cylarr( cyl%id_con(2) )%r
					call get_cylinder_random_rh_adaptive_change(cyl%h,cyl%rmin,cyl%rmax,cyl%h,rr1,rr2)
				end select
				!
				rsearch=sqrt( (cyl%rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cyl%h)**2 ) + cmcmc_change_delta_r*cint_cyl_connection_radius
	            call extract_galaxy_coordinates(cyl%p,ipts,ngal,rsearch)
				! initialise rotation matrixes
				rx(1,1)=cylarr(irnd)%u; rx(3,3)=rx(1,1)
				rx(1,3)=sqrt(1.0_rk-cylarr(irnd)%u**2); rx(3,1)=-rx(1,3)
				!
				rz(1,1)=cos(cylarr(irnd)%t); rz(2,2)=rz(1,1)
				rz(1,2)=-sin(cylarr(irnd)%t); rz(2,1)=-rz(1,2)
	            !
				!
				olemas=.false.
	            indo30b: do k=1,c_move_n_change
					! find new point in a sphere mcmc_change_delta_r
					do
		                call random_number(p); p=2*p-1.0
						if (sum(p**2)>1.0) cycle
						p=p*cmcmc_change_delta_r*cint_cyl_connection_radius
						exit
					end do
					cyl%p=cylarr(id)%p + p
	                !
					!---- this is new one...
					call random_number(cyl%u); cyl%u=cmcmc_parallel_cosi_change+(1.0-cmcmc_parallel_cosi_change)*cyl%u
					call random_number(cyl%t); cyl%t=cyl%t*twopi
					dum=sqrt(1.0-cyl%u**2)
			        vec(1)=cos(cyl%t)*dum
			        vec(2)=sin(cyl%t)*dum
			        vec(3)=cyl%u
					vec=matmul(rx,vec)
					vec=matmul(rz,vec)
					! vec:xyz to cyl:tu -- take account signs (atan2 definition)
					if (vec(1)>0.0) then
						cyl%t=atan(vec(2)/vec(1))
					else if (vec(2)>=0.0 .and. vec(1)<0.0) then
						cyl%t=atan(vec(2)/vec(1))+pi
					else if (vec(2)<0.0 .and. vec(1)<0.0) then
						cyl%t=atan(vec(2)/vec(1))-pi
					else if (vec(2)>0.0 .and. abs(vec(1))<100*epsilon(1.0) ) then
						cyl%t=pi/2
					else if (vec(2)<0.0 .and. abs(vec(1))<100*epsilon(1.0) ) then
						cyl%t=-pi/2
					else if (abs(vec(2))<100*epsilon(1.0) .and. abs(vec(1))<100*epsilon(1.0) ) then
						cyl%t=0.0
					end if
					if (cyl%t<0.0) cyl%t=twopi+cyl%t
					cyl%u=vec(3)
					if (cyl%u<0.0) then
						cyl%u=-cyl%u
						if (cyl%t<pi) then
							cyl%t=cyl%t+pi
						else
							cyl%t=cyl%t-pi
						end if
					end if
					! done atan2
	                call update_cyl_automatic_params(cyl)
	                dum=calc_data_potential_for_cylinder(cyl,ipts,ngal,loctest=.true.)
	                if (dum<cdata_pot_undefined) then
	                    olemas=.true.
						cstat_nr_cyl_change=cstat_nr_cyl_change+1
	                    exit indo30b
	                end if
	            end do indo30b
				if (olemas) exit
			end do
			!
        end select ! operation
        !====================================
        !
        ! final test...just in case, something bad happen
        if (.not.olemas) then
			id=0; operation=11
            !=== random cylinder
			do i0=1,100
				call get_cylinder_random_rh_adaptive(cyl%h,cyl%rmin,cyl%rmax)
				rsearch=sqrt( (cyl%rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cyl%h)**2 ) ! cylinder diagonal
	            indo11b: do i=1,n_test*10
		            call random_number(p)
					p=cmask%pmin_active+p*(cmask%pmax_active-cmask%pmin_active)
		            ! check if the position exists
		            test=test_sample_mask_location(p)
		            if (test) then
						call extract_galaxy_coordinates(p,ipts,ngal,rsearch)
		                if (ngal<cdata_minpts) cycle indo11b
						! find the orientation
		                cyl%p=p
						! cycle over orientations helps to speed up calculations a little bit
		                ori11b: do k=1,c_move_n_orientation*5
							call random_number(cyl%t); cyl%t=cyl%t*twopi
		                    call random_number(cyl%u)
		                    call update_cyl_automatic_params(cyl)
		                    dum=calc_data_potential_for_cylinder(cyl,ipts,ngal,loctest=.true.)
		                    if (dum<cdata_pot_undefined) then
		                        olemas=.true.
		                        cstat_nr_cyl_birth_random=cstat_nr_cyl_birth_random+1
		                        exit indo11b
		                    end if
		                end do ori11b
		            end if
	            end do indo11b
				if (olemas) exit
			end do
        end if
        if (.not.olemas) then
            stop "Error: impossible to find new move!"
        end if
    end subroutine generate_operation_for_cylinder
	!
end module mcmc_module