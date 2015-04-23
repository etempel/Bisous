!============
! Author: Elmo Tempel
! Date (last changed): 22.04.2015
!============
!
program bisous_model
	use constants
	use parameters
	use utilities
	use cylinders
	use galaxies
	use data_term
	use operations
	use mcmc_module
	use simulation
	use preparation
	use read_cylinders
	use post_processing
	use fil_cat
    implicit none
	character(len=:),allocatable:: infile ! local variable inside this file
	!
	call initialise_program()
	!
	select case (c_main_option)
		!==============
		case (1) ! normal MCMC sampling
		call initialise_for_mcmc()
		call start_program_bisous_simulation(c_resume)
		!==============
		case (2) ! post processing
		!call convert_cyl_files_to_ascii_v2()
		call read_cylinders_from_realisations()
		!call test_post_processing()
		call extract_filament_spines()
		call write_basic_catalogue_to_file('test_filament_spines.txt')
		!==============
		case default; stop "ERROR: program action is not defined!"
	end select
    !
contains
    !
	! read input file and initialise program
	subroutine initialise_program()
		implicit none
		character(len=300):: c_input
		integer:: clen,status
		call random_seed()
		call get_time_min(prog_start_time)
		!
        ! read the input file from command line argument (first argument)
        call get_command_argument(1, c_input, clen, status)
        if (status .ne. 0) then
            print*, "No input file to read!"
            stop
        else
            infile=c_input(1:clen)
			print*, "Reading input file ", infile
        end if
		call read_input_parameters(infile)
		!
		! read the root of output (from command line if present)
        call get_command_argument (2, c_input, clen, status)
        if (status .ne. 0) then
            print*, "Output root from input file: ", c_root
        else
			c_root=c_input(1:clen)
			print*, "Output root from command line: ", c_root
        end if
		!
		select case (c_input_type)
			case (1) ! regular x,y,z data
			call read_points_xyzlum(fname_data)
			case (2) ! regular x,y,z data + distance limit
			call read_points_xyzlum(fname_data,c_dummy_1)
			case default
			stop "Error: input type not defined!"
		end select
		!
		! set sampling volume by hand... or from galaxy distribution
		if (c_cmask_plim_use) then
			call init_cylarr(c_cylarr_size,c_cmask_pmin,c_cmask_pmax)
		else
			call init_cylarr(c_cylarr_size,minval(gal%p,dim=2),maxval(gal%p,dim=2))
		end if
	end subroutine initialise_program
	subroutine initialise_for_mcmc()
		implicit none
		if (cmcmc_auto_sampling_volume) then
			call set_sampling_volume_fraction_for_mcmc(cmcmc_auto_samvol_acc,cmcmc_auto_samvol_maxtime)
			call fix_input_parameters()
		end if
		!
		! initialise mcmc normalising volume
		cmcmc_norm_volume=cmcmc_norm_volume*cmcmc_volume_fraction
		print*, "MCMC normalising volume: ", real(cmcmc_norm_volume), real(cmcmc_norm_volume*cmcmc_volume_multiplier)
		!
		call write_input_parameters(c_root//'par_'//infile)
	end subroutine initialise_for_mcmc
	!
end program bisous_model