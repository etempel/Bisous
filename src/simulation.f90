!============
! Author: Elmo Tempel
! Date (last changed): 28.04.2016
!============
!
module simulation
    use constants
    use cylinders
    use galaxies
    use operations
    use parameters
    use data_term
    use interaction_term
	use mcmc_module
	use utilities
	use preparation
	implicit none
	integer,parameter,private:: c_max_nr_moves_multiplier=10 ! maksimaalne movede arv yhes iteratsioonis (times ncylsize)
	real(rk),private:: t_delta_max
	type, private :: simtemp_type
		real(rk),dimension(:),allocatable:: stemp ! ST temperature ning data energy
		real(rk),dimension(:,:),allocatable:: en ! ST data energy and sigma
		real(rk),dimension(:,:),allocatable:: nx_cyl ! nr of 0,1,2 connected cylinders (sum of all updates)
		integer,dimension(:),allocatable:: n_en ! nr of energy updates
		real(rk),dimension(:,:),allocatable:: nx_avg ! average nr of cylinders in every energy -- average over all updates
		integer:: ii,ii_min,ii_max ! temperature ladder
		integer:: ncylsize ! average nr of cyls... to determine minimum moves per cycle
		real(rk):: avg_pot,avg_sig
		integer:: avg_cnt
		! parameters for combined (tempering+annealing) scheme
		integer:: tempering_start ! cycle nr where tempering starts
		integer:: annealing_start ! cycle nr where annealing starts
		real(rk):: log_base_const ! logarithm base constant for simulated annealing
	end type simtemp_type
contains
	!
    subroutine start_program_bisous_simulation(resume)
        implicit none
		type(simtemp_type):: st
		logical,intent(in):: resume
		integer:: iunit_stat,iunit2
		real(rk):: tmin1,tmin2,aeg ! time stamps...to track program speed
		real(rk):: tmin1a,tmin2a
		real(rk):: temp
		integer:: ifile,nline ! file index, lines in stat file
		real(rk),dimension(1:3):: pdat ! data potential for statistics
		integer:: i,n0,n1,n2,interrupt, nn_moves
		real(rk):: dum
		character(len=1000):: cdumline
		logical:: suc
		! simulated tempering parameters
		real(rk):: alpha,teg,delta_g,delta_beta,delta_h,dum1,dum2,prob,xdum
		real(rk):: psig,pavg
		integer:: i2,iold,j ! active temperature number
		real(rk):: i3max,i2max,mrad,mlen,psig0
		logical:: st_test
		!
		! initialise simulated tempering
		allocate(st%stemp(1:cmcmc_st_temp_steps),st%en(1:cmcmc_st_temp_steps,1:2),st%n_en(1:cmcmc_st_temp_steps))
		allocate(st%nx_cyl(1:cmcmc_st_temp_steps,1:3),st%nx_avg(1:cmcmc_st_temp_steps,1:3))
		! these must be set because of generalisation
		st%log_base_const=1.0
		st%ii_min=1; st%ii_max=1; st%ii=1; iold=1
		st%annealing_start=cmcmc_cycles +1 ! never starts
		st%tempering_start=cmcmc_cycles
		st%avg_cnt=1; st%avg_pot=1.5; st%avg_sig=0.1 ! some initial values...
		!
		t_delta_max=0.0; interrupt=0; st%ncylsize=0.0
		! set the move probabilities (in principle, can be changed during simulation)
        cmcmc_adapted_prop_birth=cmcmc_prop_birth
		cmcmc_adapted_prop_death=cmcmc_prop_death
		cmcmc_adapted_prop_change=cmcmc_prop_change
		cmcmc_adapted_prop_connected=cmcmc_prop_connected
        !
		if (resume) then
			! resume option should be handled here
			call read_resume_file_single(nline,st,suc)
			if (.false.) then ! here are some test subroutines (for debugging and development)
				!call analyse_adaptive_data_energy()
			end if
		end if
		! prepare statistics file
		if (.not.resume .or. .not.suc) then
			! initialise statistics file
			open(newunit=iunit_stat,file=c_root//'bisous_stats.txt', status='REPLACE')
			write(iunit_stat,fmt='(A)',advance='no') '# cycle n0 n1 n2 temp pot pdat0 pdat1 pdat2'
			write(iunit_stat,fmt='(A)') ' fb fd fc fbcon fbran fd0 fd1 fd2 aeg cint_conn_0 cint_conn_1 cint_conn_2 mrad mlen'
			nline=0
			!
			! set the simulated tempering parameters
			if (cmcmc_cooling_schedule==3 .or. cmcmc_cooling_schedule==4) then ! simulated tempering, set temperature ning initial energies
				st%ii_min=1; st%ii_max=cmcmc_st_temp_steps ! min and max simulated tempering temperature value
				alpha=(cmcmc_initial_temp/cmcmc_final_temp)**(1.0/(cmcmc_st_temp_steps-1))
				st%stemp(1)=cmcmc_final_temp
				do i=1,cmcmc_st_temp_steps-1
					st%stemp(i+1)=st%stemp(i)*alpha
				end do
				call calc_tempering_initial_energies(ncyc=3,nmov=cmcmc_cycle_moves,st=st)
				st%ii=1 ! first temperature value... start with minimum temperature
				st%avg_pot=sum(st%en(:,1))/ sum(st%nx_cyl(:,:))
				st%avg_sig=sum(st%en(:,2))/(1+st%ii_max-st%ii_min)
				!print*, st%avg_pot, sqrt(st%avg_sig)
				!read*
			end if
		else
			! resume options
			print*
			print*, "Resuming previously interrupted simulation: "
			open(newunit=iunit_stat,file=c_root//'bisous_stats.txt')
			read(iunit_stat,fmt=*)
			do i=1,nline
				read(iunit_stat,fmt=*) interrupt
				if (c_resume_from_previous_snapshot .and. interrupt>st%annealing_start) then
					! start simulation from previous snapshot
					j=i
					exit
				end if
			end do
			! set some parameters that are needed for new simulation
			if (c_resume_from_previous_snapshot) then
				nline=j
				print*, "Start new simulation from previous snapshot..."
				st%avg_pot=st%avg_pot/st%avg_cnt
				st%avg_sig=st%avg_sig/st%avg_cnt
				st%avg_cnt=1
			end if
			! -- do something nasty to be compatible with KBFI computing centre
			if (.true.) then
				close(iunit_stat)
				open(newunit=iunit_stat,file=c_root//'bisous_stats.txt')
				open(newunit=iunit2,file=c_root//'stats_tmp.txt')
				read(iunit_stat,fmt=*)
				do i=1,nline
					read(iunit_stat,fmt='(A)') cdumline
					write(iunit2,fmt='(A)') trim(cdumline)
				end do
				close(iunit_stat,status='delete')
				close(iunit2)
				open(newunit=iunit_stat,file=c_root//'bisous_stats.txt')
				open(newunit=iunit2,file=c_root//'stats_tmp.txt')
				write(iunit_stat,fmt='(A)',advance='no') '# cycle n0 n1 n2 temp pot pdat0 pdat1 pdat2'
				write(iunit_stat,fmt='(A)') ' fb fd fc fbcon fbran fd0 fd1 fd2 aeg cint_conn_0 cint_conn_1 cint_conn_2 mrad mlen'
				do i=1,nline
					read(iunit2,fmt='(A)') cdumline
					write(iunit_stat,fmt='(A)') trim(cdumline)
				end do
				close(iunit2,status='delete')
			end if
		end if
        !
        ! cycle over
		call get_time_min(tmin1); tmin1a=tmin1
		print*, "start sampling (res,suc,inter) ", resume,suc, interrupt
		!
		call reset_counting_parameters(initial=.true.)
        do i=1,cmcmc_cycles
			if (resume .and. suc .and. i<=interrupt) cycle
			if (i==1 .or. mod(i,cmcmc_every_cycles_to_file)==0) then
				call reset_counting_parameters(initial=.false.)
				call random_seed()
			end if
			!
            ! set the temperature
			select case (cmcmc_cooling_schedule)
				case(1) ! constant temperature
				temp=cmcmc_initial_temp
				case(2) ! simulated annealing
				dum=(cmcmc_initial_temp/cmcmc_final_temp-1.0)/log(real(cmcmc_cycles))
				temp=cmcmc_initial_temp/(dum*log(real(i))+1.0)
				case(3) ! simulated tempering
				temp=st%stemp(st%ii)
				case(4) ! tempering + annealing
				if (i<=st%annealing_start) then
					temp=st%stemp(st%ii)
				else
					temp=st%stemp(st%ii_max)/(st%log_base_const * log(real(i-st%annealing_start))+1.0)
				end if
				case default
				stop "Error: cooling schedule not defined!"
			end select
			!
			! use prior for data term coefficients
			call random_number(dum)
			cdata_coeff_variance=cdata_coeff_variance_min + dum*(cdata_coeff_variance_max-cdata_coeff_variance_min)
			call random_number(dum)
			cdata_coeff_hypothesis=cdata_coeff_hypothesis_min + dum*(cdata_coeff_hypothesis_max-cdata_coeff_hypothesis_min)
			!
			! set new values for interaction gammas (based on data energy)
			if (cint_conn_use_self_regulation) then
				pavg=st%avg_pot/st%avg_cnt
				psig=sqrt(st%avg_sig/st%avg_cnt)
				psig0=psig
				!
				cint_conn_0_min=cint_self_conn_0_mean*pavg-psig*cint_self_conn_0_disp
				cint_conn_1_min=cint_self_conn_1_mean*pavg-psig*cint_self_conn_1_disp
				cint_conn_2_min=cint_self_conn_2_mean*pavg-psig*cint_self_conn_2_disp
				cint_conn_0_max=cint_self_conn_0_mean*pavg+psig*cint_self_conn_0_disp
				cint_conn_1_max=cint_self_conn_1_mean*pavg+psig*cint_self_conn_1_disp
				cint_conn_2_max=cint_self_conn_2_mean*pavg+psig*cint_self_conn_2_disp
				! for constant values
				cint_conn_0=cint_self_conn_0_mean*pavg
				cint_conn_1=cint_self_conn_1_mean*pavg
				cint_conn_2=cint_self_conn_2_mean*pavg
			end if
			!
			if (cint_conn_use_interval) then
				call random_number(dum); cint_conn_0=cint_conn_0_min + dum*(cint_conn_0_max-cint_conn_0_min)
				call random_number(dum); cint_conn_1=cint_conn_1_min + dum*(cint_conn_1_max-cint_conn_1_min)
				call random_number(dum); cint_conn_2=cint_conn_2_min + dum*(cint_conn_2_max-cint_conn_2_min)
			end if
			!
			call update_cylinder_interaction_potentials()
			!
			! run one mcmc cyle
			nn_moves=max(cmcmc_cycle_moves,int( st%ncylsize*cmcmc_nr_moves_min_frac_con ))
			nn_moves=min(nn_moves,st%ncylsize*c_max_nr_moves_multiplier)
			nn_moves=max(nn_moves,cmcmc_cycle_moves/10) ! vahemalt 1/10 esialgsetest movedest
            call run_one_mcmc_cycle(nn_moves,temp)
			call get_nr_of_connections(n0,n1,n2,pdat,psig=psig) ! needed for set temperatures and automatic variables
			!
			! set scheme from tempering to annealing
			if (cmcmc_cooling_schedule==4 .and. i>st%tempering_start+cmcmc_nr_cycle_for_tempering .and. st%ii==st%ii_max) then
				st%annealing_start=i
				! set the log base for annealing
				dum=(cmcmc_cycles-i) !*cmcmc_stsa_frac_for_min_temp
				st%log_base_const=(st%stemp(st%ii_max)/st%stemp(st%ii_min)-1.0)/log(dum)
				st%ii=1; iold=1
			end if
			if (cint_conn_use_self_regulation) then
				if (n0+n1+n2>0) then
					st%avg_pot=st%avg_pot+sum(pdat)/(n0+n1+n2)
					st%avg_sig=st%avg_sig+psig
					st%avg_cnt=st%avg_cnt+1
					st%ncylsize=n1+n2
				end if
			end if
			!
			! update simulated tempering temperature (st%ii value)
			if ((cmcmc_cooling_schedule==3 .or. cmcmc_cooling_schedule==4) .and. i<st%annealing_start) then
				! update energies
				if (n0+n1+n2>0) then
					st%en(st%ii,1)=st%en(st%ii,1)+sum(pdat(1:3))
					st%en(st%ii,2)=st%en(st%ii,2)+psig
					st%nx_cyl(st%ii,1)=st%nx_cyl(st%ii,1)+n0
					st%nx_cyl(st%ii,2)=st%nx_cyl(st%ii,2)+n1
					st%nx_cyl(st%ii,3)=st%nx_cyl(st%ii,3)+n2
					st%n_en(st%ii)=st%n_en(st%ii)+1
				end if
				!print*, "nrcyl  ", n0,n1,n2, st%ii
				!print*, real(st%nx_cyl(st%ii,1:3)),st%n_en(st%ii)
				!
				! calculate average nr of cylinders (for automatic adjusting temperature ladder)
				st%ncylsize=0.0; i3max=0.0; i2max=0.0
				do j=st%ii_min,st%ii_max
					st%ncylsize=st%ncylsize + real(st%nx_cyl(j,2)+st%nx_cyl(j,3))/real(st%n_en(j))
					st%nx_avg(j,1)=real(st%nx_cyl(j,1))/real(st%n_en(j))
					st%nx_avg(j,2)=real(st%nx_cyl(j,2))/real(st%n_en(j))
					st%nx_avg(j,3)=real(st%nx_cyl(j,3))/real(st%n_en(j))
					if (st%nx_avg(j,3)>i3max) i3max=st%nx_avg(j,3)
					if (st%nx_avg(j,3)+st%nx_avg(j,2)>i2max) i2max=st%nx_avg(j,3)+st%nx_avg(j,2)
				end do
				st%ncylsize=st%ncylsize/real(st%ii_max-st%ii_min+1) ! needed for temperature ladder changes
				!print*, real(st%nx_avg(st%ii,1:3))
				!
				! adjust temperature ladder --  remove bad temperatures
				! ... only for combined simulating schedule
				if (cmcmc_cooling_schedule==4 .and. st%ii>st%ii_min .and. st%ii<st%ii_max .and. i>cmcmc_nr_cycle_for_burnin &
				.and. cmcmc_nr_cycle_for_burnin>0) then
					st_test=.true.
					if (st%nx_avg(st%ii_max,3)<cmcmc_st_adjust_temp_ladder_coef*i3max .or. &
					st%nx_avg(st%ii_max,3)+st%nx_avg(st%ii_max,2)<cmcmc_st_adjust_temp_ladder_coef*i2max .or. &
					st%nx_avg(st%ii_max,1)>cmcmc_st_adjust_temp_ladder_coef_free_max*i3max .or.&
					st%nx_avg(st%ii_max,1)<cmcmc_st_adjust_temp_ladder_coef_free_min*i3max) then
						st_test=.false.
						st%ii_max=st%ii_max-1
						cmcmc_nr_cycle_for_burnin=cmcmc_nr_cycle_for_burnin*1.1
					end if
					if (st_test .and. (st%nx_avg(st%ii_min,3)<cmcmc_st_adjust_temp_ladder_coef*i3max .or. &
					st%nx_avg(st%ii_min,3)+st%nx_avg(st%ii_min,2)<cmcmc_st_adjust_temp_ladder_coef*i2max .or.&
					st%nx_avg(st%ii_min,1)>cmcmc_st_adjust_temp_ladder_coef_free_max*i3max .or.&
					st%nx_avg(st%ii_min,1)<cmcmc_st_adjust_temp_ladder_coef_free_min*i3max)) then
						st_test=.false.
						st%ii_min=st%ii_min+1
						cmcmc_nr_cycle_for_burnin=cmcmc_nr_cycle_for_burnin*1.1
					end if
					st%tempering_start=cmcmc_nr_cycle_for_burnin+1 ! last change for temperature ladder
				end if
				if (cmcmc_cooling_schedule==4 .and. cmcmc_nr_cycle_for_burnin<=0) then
					st%tempering_start=1 ! no burn-in phase
				end if
				!
				! try temperature update, i2 is potential new temperature ladder
				if (st%ii==st%ii_min) then
					i2=st%ii_min+1; teg=0.5
				else if (st%ii==st%ii_max) then
					i2=st%ii_max-1; teg=0.5
				else
					call random_number(dum)
					i2=st%ii+1; teg=1.0
					if (dum<0.5) i2=st%ii-1
				end if
				iold=st%ii
				! probability
				delta_beta=1.0/st%stemp(i2)-1.0/st%stemp(st%ii)
				dum1=st%en(st%ii,1)/st%n_en(st%ii) + st%nx_cyl(st%ii,1)*(-cint_conn_0)/st%n_en(st%ii) + &
				st%nx_cyl(st%ii,2)*(-cint_conn_1)/st%n_en(st%ii) + st%nx_cyl(st%ii,3)*(-cint_conn_2)/st%n_en(st%ii)
				dum2=st%en(i2,1)/st%n_en(i2) + st%nx_cyl(i2,1)*(-cint_conn_0)/st%n_en(i2) + &
				st%nx_cyl(i2,2)*(-cint_conn_1)/st%n_en(i2) + st%nx_cyl(i2,3)*(-cint_conn_2)/st%n_en(i2)
				!
				dum=sum(pdat(1:3)) + n0*(-cint_conn_0) + n1*(-cint_conn_1) + n2*(-cint_conn_2)
				dum=dum/(st%ncylsize*cmcmc_st_expected_nrcyl_size)
				dum1=dum1/(st%ncylsize*cmcmc_st_expected_nrcyl_size)
				dum2=dum2/(st%ncylsize*cmcmc_st_expected_nrcyl_size)
				!
				! add second order correction... the Variance (I have to store the individual potentials...)
				!
				delta_g=delta_beta*0.5*(dum1+dum2)
				delta_h=delta_beta*dum-delta_g
				prob=teg*exp(-delta_h)
				prob=max(0.05,prob) ! there is always 5% probability for temperature change (usually not efective)
				if (i<cmcmc_nr_cycle_for_burnin) then
					! increase the probability for better burning in
					dum=maxval(st%n_en(st%ii_min:st%ii_max))
					if (st%n_en(i2)<0.5*dum) prob=max(0.5,prob)
				end if
				!
				! set new temperature
				call random_number(alpha)
				if (alpha<prob) then
					st%ii=i2
				end if
				!
			end if ! --- end of simulated tempering temperature value update
            !
			!==============================================
			! write statistics file
			!==============================================
            if (mod(i,cmcmc_every_cycles_to_stats)==0 .or. i==cmcmc_cycles) then
				!
				call get_time_min(tmin2)
				aeg=tmin2-tmin1; tmin1=tmin2
				!
				nline=nline+1 ! lines in output file
				if (c_verbose) then
					print*
		            print*, "MCMC cycle, aeg: ", i, nr_cylarr, count(cylarr(:)%inuse),real(aeg)
					if (.true.) then ! print extended output
						call get_nr_of_connections(n0,n1,n2,pdat,psig)
						print*, "nrcyl(0,1,2): ",n0,n1,n2,real(temp)
						print*, "moves, ex. cyls:: ", nn_moves, st%ncylsize
						!
						xdum=sum(st%nx_cyl(iold,1:3))/real(st%n_en(iold))
						print*, "ST, pot,sig:: ", real(temp),real((st%en(iold,1)/st%n_en(iold))/xdum),real(sqrt(st%en(iold,2)/st%n_en(iold)))
						print*, "energies (avg,sig):: ", pavg,psig0
						do j=st%ii_min,st%ii_max
							print*, real(st%stemp(j)),st%n_en(j)," : ",nint(st%nx_avg(j,1:3))
						end do
						print*, "ii min and max: ", st%ii_min,st%ii_max
						print*
					end if
				end if
				!
                !=== write statistics to file:
                call get_nr_of_connections(n0,n1,n2,pdat,mrad=mrad,mlen=mlen) ! nr of 0,1,2 connected cylinders
                if (n0+n1+n2/=nr_cylarr_active) then
                    stop "Error: number of connected cylinders is wrong...write stats..."
                end if
				dum=sum(pdat(1:3)) + n0*(-cint_conn_0) + n1*(-cint_conn_1) + n2*(-cint_conn_2)
				!
                write(iunit_stat,fmt="(I10,3I7,F9.4,4F10.5, 8F9.5,F12.6,3F10.4,2F7.3)") i,n0,n1,n2, temp, &
				dum/nr_cylarr_active, pdat(1)/n0,pdat(2)/n1,pdat(3)/n2, &
				real(cstat_nr_cyl_birth_connected_accept+cstat_nr_cyl_birth_random_accept)/real(cstat_nr_cyl_birth_random+cstat_nr_cyl_birth_connected),&
				real(cstat_nr_cyl_death_accept)/real(cstat_nr_cyl_death),real(cstat_nr_cyl_change_accept)/real(cstat_nr_cyl_change),&
				real(cstat_nr_cyl_birth_connected_accept)/real(cstat_nr_cyl_birth_connected),&
				real(cstat_nr_cyl_birth_random_accept)/real(cstat_nr_cyl_birth_random),&
				real(cstat_nr_cyl_death_0_accept)/real(cstat_nr_cyl_death_0),real(cstat_nr_cyl_death_1_accept)/real(cstat_nr_cyl_death_1),&
				real(cstat_nr_cyl_death_2_accept)/real(cstat_nr_cyl_death_2),aeg, cint_conn_0,cint_conn_1,cint_conn_2, mrad,mlen
                !
            end if
			!
            !===========================================
            ! write the output file and resume file
            !===========================================
			if (mod(i,cmcmc_every_cycles_to_file)==0 .or. i==cmcmc_cycles) then
				!
                ifile=i/cmcmc_every_cycles_to_file
                call write_cylinders_to_file_binary(c_root//'zyls_'//get_file_counter(ifile, cfile_number_spaces)//'.cyl',nline,st)
				! write resume file - for convenience (file name do not change)
				call write_cylinders_to_file_binary(c_root//'resume.bin',nline,st)
				!
				call get_time_min(tmin2a)
				aeg=tmin2a-tmin1a; tmin1a=tmin2a
				if (aeg/60.0>t_delta_max) t_delta_max=aeg/60.0 ! aeg tundides kahe kirjutamise vahel
				aeg=tmin2a-prog_start_time; aeg=aeg/60.0 ! aeg tundides programmi algusest
				! lopetame programmi too, kui ei ole piisavalt aega yle
				if (ctime_max_time-aeg<t_delta_max*1.5) then
					print*, "Stoping the program because of the time limit!"
					print*, "... program can be resumed from resume file..."
					stop "WARNING: Program terminated because of the time limit!!!"
				end if
				!
			end if
            !
        end do
    end subroutine start_program_bisous_simulation
	!
	subroutine calc_tempering_initial_energies(ncyc,nmov,st)
		implicit none
		type(simtemp_type),intent(inout):: st
		integer,intent(in):: ncyc,nmov ! ncyc x nmov steps for energy calculation
		integer:: ii,i,k,n0,n1,n2
		real(rk):: dum,psig,pavg
		real(rk),dimension(1:3):: pdat
		!
		print*, "Calc tempering initial energies..."
		st%en=0.0
		st%nx_cyl=0.0
		pavg=1.5
		psig=0.2
		!
		do ii=size(st%stemp),1,-1 ! cycle over temperatures, in decreasing order
			do k=1,ncyc ! cycle over cycles...change potential every time
				! use prior for data term coefficients
				call random_number(dum)
				cdata_coeff_variance=cdata_coeff_variance_min + dum*(cdata_coeff_variance_max-cdata_coeff_variance_min)
				call random_number(dum)
				cdata_coeff_hypothesis=cdata_coeff_hypothesis_min + dum*(cdata_coeff_hypothesis_max-cdata_coeff_hypothesis_min)
				!
				cint_conn_0_min=cint_self_conn_0_mean*pavg-psig*cint_self_conn_0_disp
				cint_conn_1_min=cint_self_conn_1_mean*pavg-psig*cint_self_conn_1_disp
				cint_conn_2_min=cint_self_conn_2_mean*pavg-psig*cint_self_conn_2_disp
				cint_conn_0_max=cint_self_conn_0_mean*pavg+psig*cint_self_conn_0_disp
				cint_conn_1_max=cint_self_conn_1_mean*pavg+psig*cint_self_conn_1_disp
				cint_conn_2_max=cint_self_conn_2_mean*pavg+psig*cint_self_conn_2_disp
				call random_number(dum); cint_conn_0=cint_conn_0_min + dum*(cint_conn_0_max-cint_conn_0_min)
				call random_number(dum); cint_conn_1=cint_conn_1_min + dum*(cint_conn_1_max-cint_conn_1_min)
				call random_number(dum); cint_conn_2=cint_conn_2_min + dum*(cint_conn_2_max-cint_conn_2_min)
				!
				call update_cylinder_interaction_potentials()
				! run one mcmc cyle
	            call run_one_mcmc_cycle(nmov,st%stemp(ii))
				!
				call get_nr_of_connections(n0,n1,n2,pdat,psig)
				pavg=sum(pdat)/(n0+n1+n2)
				!
				st%en(ii,1)=st%en(ii,1)+sum(pdat(1:3))
				st%en(ii,2)=st%en(ii,2)+psig
				st%nx_cyl(ii,1)=st%nx_cyl(ii,1)+n0
				st%nx_cyl(ii,2)=st%nx_cyl(ii,2)+n1
				st%nx_cyl(ii,3)=st%nx_cyl(ii,3)+n2
				if (k==ncyc) then
					print*, "T,en,sig: ", real(st%stemp(ii)),real(st%en(ii,1)/(sum(st%nx_cyl(ii,:)))),real(st%en(ii,2)/k)
					print*, "  Ncyl: 0,1,2:: ", n0,n1,n2
				end if
			end do
			st%en(ii,1)=st%en(ii,1)/ncyc
			st%en(ii,2)=st%en(ii,2)/ncyc
			st%nx_cyl(ii,1)=st%nx_cyl(ii,1)/ncyc
			st%nx_cyl(ii,2)=st%nx_cyl(ii,2)/ncyc
			st%nx_cyl(ii,3)=st%nx_cyl(ii,3)/ncyc
			st%n_en(ii)=1
		end do
		print*, "Done -- tempering initial energies"
	end subroutine calc_tempering_initial_energies
	!
	subroutine reset_counting_parameters(initial)
		implicit none
		logical,intent(in):: initial
		if (initial) then
			cstat_nr_cyl_birth_random=0
			cstat_nr_cyl_birth_connected=0
			cstat_nr_cyl_death=0
			cstat_nr_cyl_death_0=0
			cstat_nr_cyl_death_1=0
			cstat_nr_cyl_death_2=0
			cstat_nr_cyl_change=0
			cstat_nr_cyl_birth_random_accept=0
			cstat_nr_cyl_birth_connected_accept=0
			cstat_nr_cyl_death_accept=0
			cstat_nr_cyl_change_accept=0
			cstat_nr_cyl_death_0_accept=0
			cstat_nr_cyl_death_1_accept=0
			cstat_nr_cyl_death_2_accept=0
		else
			cstat_nr_cyl_birth_random=cstat_nr_cyl_birth_random/cmcmc_every_cycles_to_file
			cstat_nr_cyl_birth_connected=cstat_nr_cyl_birth_connected/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death=cstat_nr_cyl_death/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_0=cstat_nr_cyl_death_0/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_1=cstat_nr_cyl_death_1/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_2=cstat_nr_cyl_death_2/cmcmc_every_cycles_to_file
			cstat_nr_cyl_change=cstat_nr_cyl_change/cmcmc_every_cycles_to_file
			cstat_nr_cyl_birth_random_accept=cstat_nr_cyl_birth_random_accept/cmcmc_every_cycles_to_file
			cstat_nr_cyl_birth_connected_accept=cstat_nr_cyl_birth_connected_accept/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_accept=cstat_nr_cyl_death_accept/cmcmc_every_cycles_to_file
			cstat_nr_cyl_change_accept=cstat_nr_cyl_change_accept/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_0_accept=cstat_nr_cyl_death_0_accept/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_1_accept=cstat_nr_cyl_death_1_accept/cmcmc_every_cycles_to_file
			cstat_nr_cyl_death_2_accept=cstat_nr_cyl_death_2_accept/cmcmc_every_cycles_to_file
		end if
	end subroutine reset_counting_parameters
	!
	!------------------------------------------------------------
    subroutine write_cylinders_to_file_binary(filename,nline,st)
        implicit none
        character(len=*),intent(in):: filename
		integer,intent(in):: nline
		type(simtemp_type),intent(in):: st
        integer:: iunit,nn,i
        ! write cylinders to file
		open(newunit=iunit,file=filename,form='unformatted')
		close(iunit,status='delete')
        open(newunit=iunit,file=filename,form='unformatted',status='REPLACE')
		write(iunit) nr_cylarr_active,nr_cylarr_fixed,nline
		nn=count(cylarr(:)%inuse)
		if (nn/=nr_cylarr_active+nr_cylarr_fixed) stop "Err: wrong nr of cyls, write binary!"
		nn=0
        do i=1,nr_cylarr
            if (cylarr(i)%inuse) then
				nn=nn+1
                write(iunit) cylarr(i)%fixed,cylarr(i)%p,cylarr(i)%t,cylarr(i)%u,cylarr(i)%h,cylarr(i)%r,&
				cylarr(i)%rmin,cylarr(i)%rmax,&
                cylarr(i)%nr_con,cylarr(i)%nr_of_death, cylarr(i)%nr_of_change, cylarr(i)%nr_of_change_accept
            end if
        end do
		!
		! write simulated tempering parameters
		write(iunit) cmcmc_st_temp_steps,st%ii,st%ii_min,st%ii_max
		write(iunit) st%stemp
		write(iunit) st%en
		write(iunit) st%nx_cyl
		write(iunit) st%n_en
		write(iunit) st%nx_avg
		write(iunit) st%ncylsize,st%tempering_start,st%annealing_start,st%log_base_const
		write(iunit) st%avg_pot,st%avg_sig,st%avg_cnt
		!
        close(iunit)
		if (nn/=nr_cylarr_active+nr_cylarr_fixed) stop "Err: write cyls to binary... wrong #ncyls!"
    end subroutine write_cylinders_to_file_binary
	!
	! read resume file
	subroutine read_resume_file_single(nline,st,suc)
		implicit none
		character(len=:),allocatable:: filename
		integer,intent(out):: nline
		integer:: iunit,nn_active,nn_fixed,i,n,ii
		logical,intent(out):: suc
		integer:: ierr,stepold
		type(cylinder):: cyl
		type(simtemp_type),intent(inout):: st
		!
		print*, "read resume file..."
		!
		stepold=cmcmc_st_temp_steps
		filename=c_root//'resume.bin'
		open(newunit=iunit,file=filename,form='unformatted',action='READ',iostat=ierr)
		if (ierr/=0) then
			suc=.false.
			return
		end if
		suc=.true.
		read(iunit) nn_active,nn_fixed,nline
		n=0; i=0
		!
		do ii=1,nn_active+nn_fixed
            read(iunit) cyl%fixed, cyl%p,cyl%t,cyl%u,cyl%h,cyl%r,cyl%rmin,cyl%rmax,&
            cyl%nr_con,cyl%nr_of_death, cyl%nr_of_change, cyl%nr_of_change_accept
			!
			call update_cyl_automatic_params(cyl)
			call update_potential_for_cylinder(cyl)
			!print*, i,cyl%nr_con,cyl%potdata,cyl%potint,cyl%ngal
			if (cyl%potdata<cdata_pot_undefined) then
				i=i+1
				nr_cylarr=i
				call add_cylinder_to_array(cyl,i,old=.true.,fixed=cyl%fixed)
				if (.not.cylarr(i)%inuse) stop "Err: resume failed from single file!"
			end if
			!
			n=n+1
		end do
		if (n/=nn_active+nn_fixed) stop "Error: resume unsuccessful (0)!"
		if (.not.c_resume_from_previous_snapshot) then ! see test ei ole alati vajalik
			!if (nn_active/=nr_cylarr_active .or. nn_fixed/=nr_cylarr_fixed) stop "Error: resume unsuccessful (1)!"
		end if
		!
		! read simulated tempering parameters
		read(iunit) cmcmc_st_temp_steps,st%ii,st%ii_min,st%ii_max
		if (cmcmc_st_temp_steps/=stepold) stop "ERROR: resume: wrong or modified parameter file!!!"
		read(iunit) st%stemp
		read(iunit) st%en
		read(iunit) st%nx_cyl
		read(iunit) st%n_en
		read(iunit) st%nx_avg
		read(iunit) st%ncylsize,st%tempering_start,st%annealing_start,st%log_base_const
		read(iunit) st%avg_pot,st%avg_sig,st%avg_cnt
		!
		close(iunit)
		!
		call update_cylinder_interaction_potentials()
	end subroutine read_resume_file_single
end module simulation




