!============
! Author: Elmo Tempel
! Date (last changed): 10.11.2015
!============
!
module parameters
	use constants
	use config
	implicit none
	!
	! real parameters - in principle, can be set in config file (as hidden parameters) - not yet implemented
	integer,parameter:: cfile_number_spaces=5 ! up to 99999 files to output
	! parameters to control fast cylinder finding -- these numbers cannot be exceeded
	integer,parameter:: c_cmask_ncyl_for_cell=8 ! maximum number of cylinders for each cell in cmask
	real(rk),parameter:: c_cmask_delta=0.75 ! cmask cell size in units of cylinder length (minimum length)
	!
	! general input and output options
	character(len=:),allocatable:: c_root ! root of output
	!--- scale factor... can be used for multiscale filaments or for rescaling of datapoints
	real(rk):: c_scale_factor
	integer:: c_cylarr_size ! max size of cylarr. Memory is allocated for these cylinders
	!
	integer:: c_input_type, c_main_option ! input filetype and action to be carried out
	integer:: c_analyse_option ! input option for analyse (c_main_option=2)
	character(len=:),allocatable:: fname_data
	real(rk):: c_dummy_1
	logical:: c_resume ! resume previously interrupted simulation (from resume file)
	logical:: c_resume_from_previous_snapshot ! start simulation using previous snapshot resume file
	logical:: c_verbose ! output
	!
	! mcmc parameters
	integer:: cmcmc_cycles ! nr of mcmc cycles (temperature is updated after every cycle)
	integer:: cmcmc_cycle_moves ! moves in every cycle (mcmc run with constant temperature)
	real(rk):: cmcmc_nr_moves_min_frac_con ! min nr of moves wrt nr of connected cylinders
	integer:: cmcmc_every_cycles_to_stats ! every nth cycle goes to stats and output: mod(stats,file)==0
	integer:: cmcmc_every_cycles_to_file ! every nth cycle goes to output file (and resume file)
	integer:: cmcmc_cooling_schedule ! cooling schedule type
	!
	real(rk):: cmcmc_initial_temp, cmcmc_final_temp ! temperatures for simulation
	integer:: cmcmc_st_temp_steps ! simulated tempering nr of temperatures
 	integer:: cmcmc_nr_cycle_for_burnin,cmcmc_nr_cycle_for_tempering ! combined temperature scheme
	real(rk):: cmcmc_st_expected_nrcyl_size ! in units of n1+n2 cyls (for simulated tempering)
	real(rk):: cmcmc_st_adjust_temp_ladder_coef
	real(rk):: cmcmc_st_adjust_temp_ladder_coef_free_min,cmcmc_st_adjust_temp_ladder_coef_free_max
	!
	real(rk):: cmcmc_prop_birth,cmcmc_prop_death,cmcmc_prop_change
	real(rk):: cmcmc_prop_connected ! fraction of connected births in birth proposal
	real(rk):: cmcmc_change_delta_r ! delta r for change..cylinder center (physical units)
	real(rk):: cmcmc_parallel_cosi_change
	real(rk):: cmcmc_volume_fraction
	real(rk):: cmcmc_volume_multiplier,cmcmc_auto_volume_multiplier
	logical:: cmcmc_auto_sampling_volume ! use automatic sampling volume
	real(rk):: cmcmc_auto_samvol_acc,cmcmc_auto_samvol_maxtime
	! internal parameters
	real(rk):: cmcmc_adapted_prop_birth,cmcmc_adapted_prop_death,cmcmc_adapted_prop_change
	real(rk):: cmcmc_adapted_prop_connected ! fraction of connected births in birth proposal
	real(rk):: cmcmc_norm_volume ! normalised volume for MH
	!
	!data potential parameters
	integer:: cdata_minpts  ! min nr of points in a cylinder
	logical:: cdata_fixed_radius
	real(rk):: cdata_cyl_rad_min, cdata_cyl_rad_max ! cylinder radius (min and max values)
	real(rk):: cdata_cyl_len_min,cdata_cyl_len_max ! cyl min and max length (in units of cylinder diameter)
	real(rk):: cdata_cyl_rel_lvsd_min,cdata_cyl_rel_lvsd_max ! length vs diameter min and max values
	real(rk):: cdata_cyl_shadow ! shadow radius compared with cylinder radius
	real(rk):: cdata_p_dens_cut
	real(rk):: cdata_p_uniform_cut
	real(rk):: cdata_uniform_dens ! density along a cylinder
	real(rk):: cdata_local_dens ! local contrast
	real(rk):: cdata_coeff_variance, cdata_coeff_hypothesis
	real(rk):: cdata_coeff_variance_min, cdata_coeff_hypothesis_min
	real(rk):: cdata_coeff_variance_max, cdata_coeff_hypothesis_max
	integer:: cdata_nr_of_rad_samples
	! internal parameters
	real(rk):: cdata_pot_undefined
	real(rk):: cdata_shadow_dvol ! internal parameter
	real(rk):: cdata_locdens_dens ! local contrast.. takes account cyl and shadow volumes
	!
	! interaction potential parameters
	real(rk):: cint_conn_0,cint_conn_0_min,cint_conn_0_max
	real(rk):: cint_conn_1,cint_conn_1_min,cint_conn_1_max
	real(rk):: cint_conn_2,cint_conn_2_min,cint_conn_2_max
	logical:: cint_conn_use_interval, cint_conn_use_self_regulation
	real(rk):: cint_orthogonal_cosi
	real(rk):: cint_parallel_cosi
	real(rk):: cint_cyl_connection_radius
	real(rk):: cint_concyl_radius_dif,cint_repulsive_radius_dif
	! self-regulating interaction potential values
	real(rk):: cint_self_conn_0_mean,cint_self_conn_1_mean,cint_self_conn_2_mean ! mean value in units of mean data potential
	real(rk):: cint_self_conn_0_disp,cint_self_conn_1_disp,cint_self_conn_2_disp ! interaction dispersion in units of data pot sigma
	! internal
	real(rk):: cint_repulsive
	!
	logical:: c_cmask_plim_use ! sampling limits are defined in config file. if false, galaxy data will determine them
	logical:: c_cmask_use_boundary_cyls ! use boundary cylinders
	real(rk),dimension(1:3):: c_cmask_pmin,c_cmask_pmax ! min and max values for new cylinder locations
	! --- c_cmask_pmin/pmax does not define mask, only allocatable region
	!... actually, it defines the mask if data is in a cube (simulation)
	!... e.g. useful for periodic box... however, does not use real periodiciy
	!
	! parameters for post processing
	character(len=:),allocatable:: craw_root_dir, craw_file_prefix, craw_root_dir_run
	logical:: craw_multiple_runs, craw_old_format
	integer:: craw_run_first,craw_run_last,craw_rs_first,craw_rs_last,craw_rs_every
	real(rk):: craw_weight_0,craw_weight_1,craw_weight_2 ! weights for 0,1,2 connected cylinders
	real(rk):: craw_rsmooth
	! parameters for spine extraction
	real(rk):: craw_spine_initial_dx ! initial grid size in input units
	real(rk):: craw_spine_min_visitmap_value ! min visitmap value 0..1
	real(rk):: craw_spine_location_acc ! spine location accuracy in input units
	integer:: craw_spine_location_acc_tune ! final accuracy is better
	real(rk):: craw_spine_min_ori_strength ! min orientation strength 0..1
	integer:: craw_max_nr_of_fil_spines, craw_max_nr_of_fil_points
	real(rk):: craw_spine_spacing_rad_frac, craw_spine_spacing_max
	character(len=:),allocatable:: c_output_spines, c_output_points, c_output_general, c_input_points ! root and file beginning of output spines
	!
	! testing parameters
	logical:: ctest_constant_dataterm ! use constant dataterm
	real(rk):: ctest_constant_dataterm_value
	logical:: ctest_no_interactions
	!
	! internal parameters
	real(rk):: c_cubesort_rsearch
	real(rk):: c_connection_volume
	integer:: c_move_n_orientation
	integer:: c_move_n_change
	!
	! statistics about mcmc dynamics
	integer*8:: cstat_nr_cyl_birth_random ! nr of random birth cylinders
	integer*8:: cstat_nr_cyl_birth_connected ! nr of connected random birth cylinders
	integer*8:: cstat_nr_cyl_death ! nr of deaths
	integer*8:: cstat_nr_cyl_death_0,cstat_nr_cyl_death_1,cstat_nr_cyl_death_2
	integer*8:: cstat_nr_cyl_change ! nr of change
	integer*8:: cstat_nr_cyl_birth_random_accept,cstat_nr_cyl_birth_connected_accept
	integer*8:: cstat_nr_cyl_death_accept, cstat_nr_cyl_change_accept
	integer*8:: cstat_nr_cyl_death_0_accept,cstat_nr_cyl_death_1_accept,cstat_nr_cyl_death_2_accept
	!
	real(rk):: prog_start_time ! programmi algusaeg
	real(rk):: ctime_max_time ! maksimaalne programmi aeg tundides
	!
contains
	!
	subroutine read_input_parameters(filename)
		implicit none
		character(len=*),intent(in):: filename
		call init_parameters()
		call set_conf_file(filename,warn=.false.)
		!
		call get_conf('general','main_option',c_main_option,1)
		call get_conf('general','verbose',c_verbose,.true.)
		call get_conf('general','resume',c_resume,.false.)
		call get_conf('general','resume_from_previous_snapshot',c_resume_from_previous_snapshot,.false.)
		call get_conf('general','max_time_hour',ctime_max_time,24.0_rk)
		call get_conf('general','input_type',c_input_type,1)
		call get_conf('general','file_data',fname_data)
		call get_conf('general','dummy_1',c_dummy_1,0.0_rk)
		call get_conf('general','root',c_root,'bm_')
		call get_conf('general','scale_factor',c_scale_factor,1.0_rk)
		call get_conf('general','cylarr_size',c_cylarr_size,100000)
		!
		call get_conf('mcmc','nr_cycles',cmcmc_cycles,10000)
		call get_conf('mcmc','nr_moves',cmcmc_cycle_moves,10000)
		call get_conf('mcmc','nr_moves_min_frac_con',cmcmc_nr_moves_min_frac_con,1.0_rk)
		call get_conf('mcmc','every_cycle_to_stat',cmcmc_every_cycles_to_stats,10)
		call get_conf('mcmc','every_cycle_to_output',cmcmc_every_cycles_to_file,100)
		call get_conf('mcmc','cooling_schedule',cmcmc_cooling_schedule,4)
		call get_conf('mcmc','temp_initial',cmcmc_initial_temp,2.0_rk)
		call get_conf('mcmc','temp_final',cmcmc_final_temp,0.5_rk)
		call get_conf('mcmc','st_temp_steps', cmcmc_st_temp_steps,10)
		call get_conf('mcmc','st_expected_nrcyl_size', cmcmc_st_expected_nrcyl_size,0.5_rk)
		call get_conf('mcmc','nr_cycle_for_tempering',cmcmc_nr_cycle_for_tempering,3000)
		call get_conf('mcmc','nr_cycle_for_burnin',cmcmc_nr_cycle_for_burnin,1000)
		call get_conf('mcmc','st_adjust_temp_ladder_coef', cmcmc_st_adjust_temp_ladder_coef,0.4_rk)
		call get_conf('mcmc','st_adjust_temp_ladder_coef_free_min', cmcmc_st_adjust_temp_ladder_coef_free_min,0.05_rk)
		call get_conf('mcmc','st_adjust_temp_ladder_coef_free_max', cmcmc_st_adjust_temp_ladder_coef_free_max,1.0_rk)
		!
		call get_conf('mcmc','prop_birth',cmcmc_prop_birth,0.5_rk)
		call get_conf('mcmc','prop_death',cmcmc_prop_death,0.3_rk)
		call get_conf('mcmc','prop_change',cmcmc_prop_change,0.2_rk)
		call get_conf('mcmc','prop_birth_connected',cmcmc_prop_connected,0.8_rk)
		call get_conf('mcmc','change_delta_r',cmcmc_change_delta_r,0.2_rk)
		call get_conf('mcmc','change_delta_cosi',cmcmc_parallel_cosi_change,0.95_rk)
		call get_conf('mcmc','auto_sampling_volume',cmcmc_auto_sampling_volume,.true.)
		call get_conf('mcmc','auto_volume_multiplier',cmcmc_auto_volume_multiplier,1.0_rk)
		call get_conf('mcmc','auto_samvol_acc',cmcmc_auto_samvol_acc,0.05_rk)
		call get_conf('mcmc','auto_samvol_maxtime',cmcmc_auto_samvol_maxtime,1.0_rk)
		call get_conf('mcmc','volume_fraction',cmcmc_volume_fraction,0.1_rk)
		call get_conf('mcmc','volume_multiplier',cmcmc_volume_multiplier,10.0_rk)
		!
		call get_conf('data_term','min_pts',cdata_minpts,3)
		call get_conf('data_term','use_fixed_radius',cdata_fixed_radius,.false.)
		call get_conf('data_term','cyl_rad_min',cdata_cyl_rad_min,0.4_rk)
		call get_conf('data_term','cyl_rad_max',cdata_cyl_rad_max,1.2_rk)
		call get_conf('data_term','cyl_len_min',cdata_cyl_len_min,3.0_rk)
		call get_conf('data_term','cyl_len_max',cdata_cyl_len_max,10.0_rk)
		call get_conf('data_term','cyl_rel_lvsd_min',cdata_cyl_rel_lvsd_min,3.0_rk)
		call get_conf('data_term','cyl_rel_lvsd_max',cdata_cyl_rel_lvsd_max,10.0_rk)
		call get_conf('data_term','cyl_shadow_rad',cdata_cyl_shadow,1.0_rk)
		call get_conf('data_term','hyptest_uniform_den',cdata_uniform_dens,4.0_rk)
		call get_conf('data_term','hyptest_local_den',cdata_local_dens,3.0_rk)
		call get_conf('data_term','hyptest_uniform_p_cut',cdata_p_uniform_cut,0.001_rk)
		call get_conf('data_term','hyptest_local_p_cut',cdata_p_dens_cut,0.001_rk)
		call get_conf('data_term','variance_coeff_min',cdata_coeff_variance_min,0.5_rk)
		call get_conf('data_term','variance_coeff_max',cdata_coeff_variance_max,0.5_rk)
		call get_conf('data_term','hypothesis_coeff_min',cdata_coeff_hypothesis_min,1.0_rk)
		call get_conf('data_term','hypothesis_coeff_max',cdata_coeff_hypothesis_max,1.0_rk)
		call get_conf('data_term','nr_of_rad_samples',cdata_nr_of_rad_samples,10)
		!
		call get_conf('interaction_term','rad_connection',cint_cyl_connection_radius,0.5_rk)
		call get_conf('interaction_term','cos_orthogonal',cint_orthogonal_cosi,0.3_rk)
		call get_conf('interaction_term','cos_parallel',cint_parallel_cosi,0.85_rk)
		call get_conf('interaction_term','concyl_radius_difference',cint_concyl_radius_dif,1.5_rk)
		call get_conf('interaction_term','repulsive_radius_difference',cint_repulsive_radius_dif,3.0_rk)
		call get_conf('interaction_term','lg_gamma_use_self_regulation',cint_conn_use_self_regulation,.true.)
		call get_conf('interaction_term','lg_gamma_0_mean',cint_self_conn_0_mean,-1.5_rk)
		call get_conf('interaction_term','lg_gamma_1_mean',cint_self_conn_1_mean,-0.2_rk)
		call get_conf('interaction_term','lg_gamma_2_mean',cint_self_conn_2_mean,1.0_rk)
		call get_conf('interaction_term','lg_gamma_0_disp',cint_self_conn_0_disp,0.8_rk)
		call get_conf('interaction_term','lg_gamma_1_disp',cint_self_conn_1_disp,0.8_rk)
		call get_conf('interaction_term','lg_gamma_2_disp',cint_self_conn_2_disp,0.8_rk)
		!
		call get_conf('interaction_term','lg_gamma_use_interval',cint_conn_use_interval,.true.)
		call get_conf('interaction_term','lg_gamma_0_min',cint_conn_0_min,-1.0_rk)
		call get_conf('interaction_term','lg_gamma_1_min',cint_conn_1_min,-0.3_rk)
		call get_conf('interaction_term','lg_gamma_2_min',cint_conn_2_min,1.5_rk)
		call get_conf('interaction_term','lg_gamma_0_max',cint_conn_0_max,-1.0_rk)
		call get_conf('interaction_term','lg_gamma_1_max',cint_conn_1_max,-0.3_rk)
		call get_conf('interaction_term','lg_gamma_2_max',cint_conn_2_max,1.5_rk)
		call get_conf('interaction_term','lg_gamma_0',cint_conn_0,-1.0_rk)
		call get_conf('interaction_term','lg_gamma_1',cint_conn_1,-0.3_rk)
		call get_conf('interaction_term','lg_gamma_2',cint_conn_2,1.5_rk)
		!
		call get_conf('cmask','c_cmask_use_boundary_cyls',c_cmask_use_boundary_cyls,.false.)
		call get_conf('cmask','c_cmask_plim_use',c_cmask_plim_use,.false.)
		call get_conf('cmask','c_cmask_pmin_x',c_cmask_pmin(1),0.0_rk)
		call get_conf('cmask','c_cmask_pmin_y',c_cmask_pmin(2),0.0_rk)
		call get_conf('cmask','c_cmask_pmin_z',c_cmask_pmin(3),0.0_rk)
		call get_conf('cmask','c_cmask_pmax_x',c_cmask_pmax(1),0.0_rk)
		call get_conf('cmask','c_cmask_pmax_y',c_cmask_pmax(2),0.0_rk)
		call get_conf('cmask','c_cmask_pmax_z',c_cmask_pmax(3),0.0_rk)
		!
		call get_conf('raw_data','root_dir',craw_root_dir,'results')
		call get_conf('raw_data','file_prefix',craw_file_prefix,'')
		call get_conf('raw_data','multiple_runs',craw_multiple_runs,.false.)
		call get_conf('raw_data','root_dir_run',craw_root_dir_run,'run_000')
		call get_conf('raw_data','run_first',craw_run_first,1)
		call get_conf('raw_data','run_last',craw_run_last,1)
		call get_conf('raw_data','rs_first',craw_rs_first,1)
		call get_conf('raw_data','rs_last',craw_rs_last,1)
		call get_conf('raw_data','rs_every',craw_rs_every,1)
		call get_conf('raw_data','weight_0',craw_weight_0,0.0_rk)
		call get_conf('raw_data','weight_1',craw_weight_1,0.5_rk)
		call get_conf('raw_data','weight_2',craw_weight_2,1.0_rk)
		call get_conf('raw_data','r_smooth',craw_rsmooth,1.0_rk)
		call get_conf('raw_data','old_format',craw_old_format,.false.)
		!
		call get_conf('analyse','analyse_option',c_analyse_option,0)
		call get_conf('analyse','output_general',c_output_general,'temp_output_')
		call get_conf('analyse','output_spines',c_output_spines,'bisousfil_')
		call get_conf('analyse','input_points',c_input_points,'nofile.txt')
		call get_conf('analyse','output_points',c_output_points,'bisous_stat_for_pts.txt')
		call get_conf('analyse','initial_grid_size',craw_spine_initial_dx,1.0_rk)
		call get_conf('analyse','minimim_visitmap_value',craw_spine_min_visitmap_value,0.05_rk)
		call get_conf('analyse','spine_location_accuracy',craw_spine_location_acc,0.1_rk)
		call get_conf('analyse','spine_location_accuracy_tuning',craw_spine_location_acc_tune,5)
		call get_conf('analyse','minimum_orientation_strength',craw_spine_min_ori_strength,0.5_rk)
		call get_conf('analyse','spine_spacing_rad_frac',craw_spine_spacing_rad_frac,1.0_rk)
		call get_conf('analyse','spine_spacing_max',craw_spine_spacing_max,0.5_rk)
		call get_conf('analyse','max_nr_of_fil_spines',craw_max_nr_of_fil_spines,100000)
		call get_conf('analyse','max_nr_of_fil_points',craw_max_nr_of_fil_points,100000)
		!
		call get_conf('testing','constant_dataterm',ctest_constant_dataterm,.false.)
		call get_conf('testing','constant_dataterm_value',ctest_constant_dataterm_value,0.0_rk)
		call get_conf('testing','no_interactions',ctest_no_interactions,.false.)
		!call get_conf('','',)
		!
		call fix_input_parameters()
	end subroutine read_input_parameters
	!
	! fix and calculate input values for some parameters
	subroutine fix_input_parameters()
		implicit none
		real(rk):: svol ! shadow volume / cylinder volume
		!
		! [general parameters]
		if (ctime_max_time<0.0) ctime_max_time=0.0
		if (c_scale_factor<=0.0) stop "Err: scale factor unphysical !!!"
		if (c_cylarr_size<=1) stop "Err: cylinder array size unrealistic!"
		!
		! [mcmc]
		! test mcmc cycles... important for resume
		if (mod(cmcmc_cycles,cmcmc_every_cycles_to_stats)/=0) stop "Error1: nr of moves incorrect (mod/=0)"
		if (mod(cmcmc_cycles,cmcmc_every_cycles_to_file)/=0) stop "Error2: nr of moves incorrect (mod/=0)"
		if (mod(cmcmc_every_cycles_to_file,cmcmc_every_cycles_to_stats)/=0) stop "Error3: nr of moves incorrect (mod/=0)"
		!
		if (cmcmc_initial_temp<cmcmc_final_temp) stop "ERR: simulation temperatures unphysical!!!"
		cmcmc_st_temp_steps=max(3,cmcmc_st_temp_steps) ! three is minimum, because of implementation
		if (cmcmc_prop_birth+cmcmc_prop_death+cmcmc_prop_change>1.0+tiny(1.0)) stop "Err: birth, death, change proposals are invalid!"
		!
		! [data term]
		! minpts>=3
		if (cdata_minpts<3) stop "Error: minimum nr of galaxies must be >=3 !!!"
		! parameters that depend on other parameters
		svol=(1.0_rk+cdata_cyl_shadow)**2-1.0_rk ! shadow volume compared to cylinder volume
		cdata_shadow_dvol=(1.0+cdata_cyl_shadow)**2
		cdata_locdens_dens=cdata_local_dens/svol
		! big value if bad potential
		cdata_pot_undefined=max(100.0_rk,-5*log(cdata_p_uniform_cut*cdata_p_dens_cut)) ! precaution
		cint_repulsive=cdata_pot_undefined
		! initialise coefficients
		cdata_coeff_hypothesis=0.5*(cdata_coeff_hypothesis_min+cdata_coeff_hypothesis_max)
		cdata_coeff_variance=0.5*(cdata_coeff_variance_min+cdata_coeff_variance_max)
		!
		! [interaction term]
		if (cint_cyl_connection_radius<=0.0) stop "Err: connection radius too small!!"
		! repulsive radius difference
		cint_concyl_radius_dif=max(1.0_rk,cint_concyl_radius_dif)
		cint_repulsive_radius_dif=max(2.0_rk,cint_repulsive_radius_dif, 1.05*cint_concyl_radius_dif**2)
		! gammas check
		if (cint_conn_use_interval) then
			if (cint_conn_0_min>cint_conn_0_max) stop "gamma 0 bounds are wrong!"
			if (cint_conn_1_min>cint_conn_1_max) stop "gamma 1 bounds are wrong!"
			if (cint_conn_2_min>cint_conn_2_max) stop "gamma 2 bounds are wrong!"
		end if
		!
		! [other parameters]
		! galaxy search radius for cubesort (initial sorting for galaxies)
		c_cubesort_rsearch=sqrt( (cdata_cyl_rad_max*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cdata_cyl_len_min)**2 ) ! cylinder diagonal
		!
		c_connection_volume=cmcmc_volume_multiplier*(1.0-cint_parallel_cosi)*(4.0*pi/3.0)*cint_cyl_connection_radius**3
		!
		! post processing parameters
		! [raw_data]
		if (.not.craw_multiple_runs) then
			craw_run_first=1; craw_run_last=1
		end if
	end subroutine fix_input_parameters
	!
	subroutine init_parameters()
		implicit none
		!
		! set values for input parameters before read from config file
		c_move_n_orientation=10
		c_move_n_change=20
		!
		cstat_nr_cyl_birth_random=0
		cstat_nr_cyl_birth_connected=0
		cstat_nr_cyl_death=0
		cstat_nr_cyl_death_0=0
		cstat_nr_cyl_death_1=0
		cstat_nr_cyl_death_2=0
		cstat_nr_cyl_change=0
		cstat_nr_cyl_birth_random_accept=0
		cstat_nr_cyl_birth_connected_accept=0
		cstat_nr_cyl_death_accept=0; cstat_nr_cyl_change_accept=0
		cstat_nr_cyl_death_0_accept=0
		cstat_nr_cyl_death_1_accept=0
		cstat_nr_cyl_death_2_accept=0
	end subroutine init_parameters
	!
	subroutine write_input_parameters(filename)
		implicit none
		character(len=*),intent(in):: filename
		integer:: iunit
		open(newunit=iunit,file=filename)
		close(iunit,status='delete')
		open(newunit=iunit,file=filename)
		!
		write(iunit,fmt='(A)') '[general]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",g0)') 'main_option',c_main_option
		write(iunit,fmt='(A," = ",L)') 'verbose',c_verbose
		write(iunit,fmt='(A," = ",L)') 'resume',c_resume
		write(iunit,fmt='(A," = ",L)') 'resume_from_previous_snapshot',c_resume_from_previous_snapshot
		write(iunit,fmt='(A," = ",g0)') 'max_time_hour',real(ctime_max_time)
		write(iunit,fmt='(A," = ",g0)') 'input_type',c_input_type
		write(iunit,fmt='(A," = ",A)') 'file_data',fname_data
		write(iunit,fmt='(A," = ",g0)') 'dummy_1',real(c_dummy_1)
		write(iunit,fmt='(A," = ",A)') 'root',c_root
		write(iunit,fmt='(A," = ",g0)') 'scale_factor',real(c_scale_factor)
		write(iunit,fmt='(A," = ",g0)') 'cylarr_size',c_cylarr_size
		write(iunit,fmt=*); write(iunit,fmt=*)
		!
		write(iunit,fmt='(A)') '[mcmc]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",g0)') 'nr_cycles',cmcmc_cycles
		write(iunit,fmt='(A," = ",g0)') 'nr_moves',cmcmc_cycle_moves
		write(iunit,fmt='(A," = ",g0)') 'nr_moves_min_frac_con',real(cmcmc_nr_moves_min_frac_con)
		write(iunit,fmt='(A," = ",g0)') 'every_cycle_to_stat',cmcmc_every_cycles_to_stats
		write(iunit,fmt='(A," = ",g0)') 'every_cycle_to_output',cmcmc_every_cycles_to_file
		write(iunit,fmt='(A," = ",g0)') 'cooling_schedule',cmcmc_cooling_schedule
		write(iunit,fmt='(A," = ",g0)') 'temp_initial',real(cmcmc_initial_temp)
		write(iunit,fmt='(A," = ",g0)') 'temp_final',real(cmcmc_final_temp)
		write(iunit,fmt='(A," = ",g0)') 'st_temp_steps',cmcmc_st_temp_steps
		write(iunit,fmt='(A," = ",g0)') 'st_expected_nrcyl_size', real(cmcmc_st_expected_nrcyl_size)
		write(iunit,fmt='(A," = ",g0)') 'nr_cycle_for_tempering',cmcmc_nr_cycle_for_tempering
		write(iunit,fmt='(A," = ",g0)') 'nr_cycle_for_burnin',cmcmc_nr_cycle_for_burnin
		write(iunit,fmt='(A," = ",g0)') 'st_adjust_temp_ladder_coef', real(cmcmc_st_adjust_temp_ladder_coef)
		write(iunit,fmt='(A," = ",g0)') 'st_adjust_temp_ladder_coef_free_min', real(cmcmc_st_adjust_temp_ladder_coef_free_min)
		write(iunit,fmt='(A," = ",g0)') 'st_adjust_temp_ladder_coef_free_max', real(cmcmc_st_adjust_temp_ladder_coef_free_max)
		!
		write(iunit,fmt='(A," = ",g0)') 'prop_birth',real(cmcmc_prop_birth)
		write(iunit,fmt='(A," = ",g0)') 'prop_death',real(cmcmc_prop_death)
		write(iunit,fmt='(A," = ",g0)') 'prop_change',real(cmcmc_prop_change)
		write(iunit,fmt='(A," = ",g0)') 'prop_birth_connected',real(cmcmc_prop_connected)
		write(iunit,fmt='(A," = ",g0)') 'change_delta_r',real(cmcmc_change_delta_r)
		write(iunit,fmt='(A," = ",g0)') 'change_delta_cosi',real(cmcmc_parallel_cosi_change)
		write(iunit,fmt='(A," = ",L)') 'auto_sampling_volume',cmcmc_auto_sampling_volume
		write(iunit,fmt='(A," = ",g0)') 'auto_volume_multiplier',real(cmcmc_auto_volume_multiplier)
		write(iunit,fmt='(A," = ",g0)') 'auto_samvol_acc',real(cmcmc_auto_samvol_acc)
		write(iunit,fmt='(A," = ",g0)') 'auto_samvol_maxtime',real(cmcmc_auto_samvol_maxtime)
		write(iunit,fmt='(A," = ",g0)') 'volume_fraction',real(cmcmc_volume_fraction)
		write(iunit,fmt='(A," = ",g0)') 'volume_multiplier',real(cmcmc_volume_multiplier)
		write(iunit,fmt=*); write(iunit,fmt=*)
		!
		write(iunit,fmt='(A)') '[data_term]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",g0)') 'min_pts',cdata_minpts
		write(iunit,fmt='(A," = ",L)') 'use_fixed_radius',cdata_fixed_radius
		write(iunit,fmt='(A," = ",g0)') 'cyl_rad_min',real(cdata_cyl_rad_min)
		write(iunit,fmt='(A," = ",g0)') 'cyl_rad_max',real(cdata_cyl_rad_max)
		write(iunit,fmt='(A," = ",g0)') 'cyl_len_min',real(cdata_cyl_len_min)
		write(iunit,fmt='(A," = ",g0)') 'cyl_len_max',real(cdata_cyl_len_max)
		write(iunit,fmt='(A," = ",g0)') 'cyl_rel_lvsd_min',real(cdata_cyl_rel_lvsd_min)
		write(iunit,fmt='(A," = ",g0)') 'cyl_rel_lvsd_max',real(cdata_cyl_rel_lvsd_max)
		write(iunit,fmt='(A," = ",g0)') 'cyl_shadow_rad',real(cdata_cyl_shadow)
		write(iunit,fmt='(A," = ",g0)') 'hyptest_uniform_den',real(cdata_uniform_dens)
		write(iunit,fmt='(A," = ",g0)') 'hyptest_local_den',real(cdata_local_dens)
		write(iunit,fmt='(A," = ",g0)') 'hyptest_uniform_p_cut',real(cdata_p_uniform_cut)
		write(iunit,fmt='(A," = ",g0)') 'hyptest_local_p_cut',real(cdata_p_dens_cut)
		write(iunit,fmt='(A," = ",g0)') 'variance_coeff_min',real(cdata_coeff_variance_min)
		write(iunit,fmt='(A," = ",g0)') 'variance_coeff_max',real(cdata_coeff_variance_max)
		write(iunit,fmt='(A," = ",g0)') 'hypothesis_coeff_min',real(cdata_coeff_hypothesis_min)
		write(iunit,fmt='(A," = ",g0)') 'hypothesis_coeff_max',real(cdata_coeff_hypothesis_max)
		write(iunit,fmt='(A," = ",g0)') 'nr_of_rad_samples',cdata_nr_of_rad_samples
		write(iunit,fmt=*); write(iunit,fmt=*)
		!
		write(iunit,fmt='(A)') '[interaction_term]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",g0)') 'rad_connection',real(cint_cyl_connection_radius)
		write(iunit,fmt='(A," = ",g0)') 'cos_orthogonal',real(cint_orthogonal_cosi)
		write(iunit,fmt='(A," = ",g0)') 'cos_parallel',real(cint_parallel_cosi)
		write(iunit,fmt='(A," = ",g0)') 'concyl_radius_difference',real(cint_concyl_radius_dif)
		write(iunit,fmt='(A," = ",g0)') 'repulsive_radius_difference',real(cint_repulsive_radius_dif)
		!
		write(iunit,fmt='(A," = ",L)') 'lg_gamma_use_self_regulation',cint_conn_use_self_regulation
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_0_mean',real(cint_self_conn_0_mean)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_1_mean',real(cint_self_conn_1_mean)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_2_mean',real(cint_self_conn_2_mean)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_0_disp',real(cint_self_conn_0_disp)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_1_disp',real(cint_self_conn_1_disp)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_2_disp',real(cint_self_conn_2_disp)
		!
		write(iunit,fmt='(A," = ",L)') 'lg_gamma_use_interval',cint_conn_use_interval
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_0_min',real(cint_conn_0_min)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_1_min',real(cint_conn_1_min)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_2_min',real(cint_conn_2_min)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_0_max',real(cint_conn_0_max)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_1_max',real(cint_conn_1_max)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_2_max',real(cint_conn_2_max)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_0',real(cint_conn_0)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_1',real(cint_conn_1)
		write(iunit,fmt='(A," = ",g0)') 'lg_gamma_2',real(cint_conn_2)
		write(iunit,fmt=*); write(iunit,fmt=*)
		!
		write(iunit,fmt='(A)') '[cmask]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",L)') 'c_cmask_use_boundary_cyls',c_cmask_use_boundary_cyls
		write(iunit,fmt='(A," = ",L)') 'c_cmask_plim_use',c_cmask_plim_use
		write(iunit,fmt='(A," = ",g0)') 'c_cmask_pmin_x',real(c_cmask_pmin(1))
		write(iunit,fmt='(A," = ",g0)') 'c_cmask_pmin_y',real(c_cmask_pmin(2))
		write(iunit,fmt='(A," = ",g0)') 'c_cmask_pmin_z',real(c_cmask_pmin(3))
		write(iunit,fmt='(A," = ",g0)') 'c_cmask_pmax_x',real(c_cmask_pmax(1))
		write(iunit,fmt='(A," = ",g0)') 'c_cmask_pmax_y',real(c_cmask_pmax(2))
		write(iunit,fmt='(A," = ",g0)') 'c_cmask_pmax_z',real(c_cmask_pmax(3))
		!
		write(iunit,fmt='(A)') '[raw_data]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",A)') 'root_dir',craw_root_dir
		write(iunit,fmt='(A," = ",A)') 'file_prefix',craw_file_prefix
		write(iunit,fmt='(A," = ",L)') 'multiple_runs',craw_multiple_runs
		write(iunit,fmt='(A," = ",A)') 'root_dir_run',craw_root_dir_run
		write(iunit,fmt='(A," = ",g0)') 'run_first',craw_run_first
		write(iunit,fmt='(A," = ",g0)') 'run_last',craw_run_last
		write(iunit,fmt='(A," = ",g0)') 'rs_first',craw_rs_first
		write(iunit,fmt='(A," = ",g0)') 'rs_last',craw_rs_last
		write(iunit,fmt='(A," = ",g0)') 'rs_every',craw_rs_every
		write(iunit,fmt='(A," = ",g0)') 'weight_0',real(craw_weight_0)
		write(iunit,fmt='(A," = ",g0)') 'weight_1',real(craw_weight_1)
		write(iunit,fmt='(A," = ",g0)') 'weight_2',real(craw_weight_2)
		write(iunit,fmt='(A," = ",g0)') 'r_smooth',craw_rsmooth
		write(iunit,fmt='(A," = ",L)') 'old_format',craw_old_format
		!
		write(iunit,fmt='(A)') '[analyse]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",A)') 'analyse_option',c_analyse_option
		write(iunit,fmt='(A," = ",A)') 'output_general',c_output_general
		write(iunit,fmt='(A," = ",A)') 'output_spines',c_output_spines
		write(iunit,fmt='(A," = ",A)') 'input_points',c_input_points
		write(iunit,fmt='(A," = ",A)') 'output_points',c_output_points
		write(iunit,fmt='(A," = ",g0)') 'initial_grid_size',real(craw_spine_initial_dx)
		write(iunit,fmt='(A," = ",g0)') 'minimim_visitmap_value',real(craw_spine_min_visitmap_value)
		write(iunit,fmt='(A," = ",g0)') 'spine_location_accuracy',real(craw_spine_location_acc)
		write(iunit,fmt='(A," = ",g0)') 'spine_location_accuracy_tuning',craw_spine_location_acc_tune
		write(iunit,fmt='(A," = ",g0)') 'minimum_orientation_strength',real(craw_spine_min_ori_strength)
		write(iunit,fmt='(A," = ",g0)') 'spine_spacing_rad_frac',real(craw_spine_spacing_rad_frac)
		write(iunit,fmt='(A," = ",g0)') 'spine_spacing_max',real(craw_spine_spacing_max)
		write(iunit,fmt='(A," = ",g0)') 'max_nr_of_fil_spines',craw_max_nr_of_fil_spines
		write(iunit,fmt='(A," = ",g0)') 'max_nr_of_fil_points',craw_max_nr_of_fil_points
		!
		write(iunit,fmt='(A)') '[testing]'
		write(iunit,fmt=*)
		write(iunit,fmt='(A," = ",L)') 'constant_dataterm',ctest_constant_dataterm
		write(iunit,fmt='(A," = ",g0)') 'constant_dataterm_value',real(ctest_constant_dataterm_value)
		write(iunit,fmt='(A," = ",L)') 'no_interactions',ctest_no_interactions
		write(iunit,fmt=*); write(iunit,fmt=*)
		!
		close(iunit)
	end subroutine write_input_parameters
	!
    subroutine get_cylinder_random_rh_adaptive(h,rmin,rmax,rcon)
        implicit none
        real(rk),intent(out):: h,rmin,rmax ! length and min and max radius for adaptive data term
		real(rk),intent(in),optional:: rcon ! connected cylinder radius
        real(rk):: rnd,hmin,hmax
		!
		if (present(rcon)) then
			! connected cylinder
			rmin=max(cdata_cyl_rad_min, rcon/cint_concyl_radius_dif)
			rmax=min(cdata_cyl_rad_max, rcon*cint_concyl_radius_dif)
			hmin=min(cdata_cyl_len_max, max(cdata_cyl_len_min, rmin*2*cdata_cyl_rel_lvsd_min) )
			hmax=max(cdata_cyl_len_min, min(cdata_cyl_len_max, rmax*2*cdata_cyl_rel_lvsd_max) )
			call random_number(rnd)
			h=hmin + rnd*(hmax-hmin)
			rmin=max(rmin, h/(2*cdata_cyl_rel_lvsd_max))
			rmax=min(rmax, h/(2*cdata_cyl_rel_lvsd_min))
		else
			! free cylinder
			call random_number(rnd)
			h=cdata_cyl_len_min + rnd*(cdata_cyl_len_max-cdata_cyl_len_min)
			rmin=max(cdata_cyl_rad_min, h/(2*cdata_cyl_rel_lvsd_max))
			rmax=min(cdata_cyl_rad_max, h/(2*cdata_cyl_rel_lvsd_min))
		end if
    end subroutine get_cylinder_random_rh_adaptive
	!
    subroutine get_cylinder_random_rh_adaptive_change(h,rmin,rmax,hori,rcon1,rcon2)
        implicit none
        real(rk),intent(out):: h,rmin,rmax ! length and min and max radius for adaptive data term
		real(rk),intent(in):: hori,rcon1 ! original length connected cylinder radius or original radius
		real(rk),intent(in),optional:: rcon2 ! second connected cylinder
        real(rk):: rnd,hmin,hmax
		!
		if (present(rcon2)) then
			rmin=max(cdata_cyl_rad_min, max(rcon1,rcon2)/cint_concyl_radius_dif)
			rmax=min(cdata_cyl_rad_max, min(rcon1,rcon2)*cint_concyl_radius_dif)
		else
			rmin=max(cdata_cyl_rad_min, rcon1/cint_concyl_radius_dif)
			rmax=min(cdata_cyl_rad_max, rcon1*cint_concyl_radius_dif)
		end if
		hmin=max(cdata_cyl_len_min,hori-cint_cyl_connection_radius)
		hmax=min(cdata_cyl_len_max,hori+cint_cyl_connection_radius)		
		hmin=min(cdata_cyl_len_max, max(hmin, rmin*2*cdata_cyl_rel_lvsd_min) )
		hmax=max(cdata_cyl_len_min, min(hmax, rmax*2*cdata_cyl_rel_lvsd_max) )
		call random_number(rnd)
		h=hmin + rnd*(hmax-hmin)
		rmin=max(rmin, h/(2*cdata_cyl_rel_lvsd_max))
		rmax=min(rmax, h/(2*cdata_cyl_rel_lvsd_min))
    end subroutine get_cylinder_random_rh_adaptive_change
	!
end module parameters