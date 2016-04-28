!============
! Author: Elmo Tempel
! Date: 28.04.2016
!============
!
module read_cylinders
	use constants
	use cylinders
	use utilities
	use operations
	implicit none
    type realisation ! one realisation - all cylinders from one realisation
        integer:: nrcyl
        type(cylinder_raw),dimension(:),allocatable:: cyl
    end type realisation
	type(realisation),dimension(:),allocatable:: run
	integer:: nrun ! number of active runs
	!
	type iifast_type
		integer:: n
		integer,dimension(:),allocatable:: irun,icyl
	end type iifast_type
	type(iifast_type),dimension(:),allocatable:: iifast
	!
	integer,private:: total_cyls
	real(rk):: csearch_min_rad
	real(rk),dimension(1:3),private:: potarr,potarr2 ! data potential array for statistics
	integer,private:: npotarr2
	!
contains
	!
	! subroutine to write binary data for plotting in processing
	subroutine write_cylinders_information_for_3d_plotting(dirout)
		implicit none
		character(len=*),intent(in):: dirout
		integer:: i,iunit,k
		print*, "Write data for plotting: ", dirout
		do i=1,nrun
			open(newunit=iunit,file=dirout//'/cyls_bin_'//get_file_counter(i, 5)//'.cyl', access='stream', form='unformatted')
			!open(55,file='testcyl.txt')
			write(iunit) run(i)%nrcyl
			!write(55,fmt=*) run(i)%nrcyl
			do k=1,run(i)%nrcyl
				write(iunit) real(run(i)%cyl(k)%p(1:3)),real(run(i)%cyl(k)%u),real(run(i)%cyl(k)%t),&
				real(run(i)%cyl(k)%h),real(run(i)%cyl(k)%r),run(i)%cyl(k)%nr_con
				!write(55,fmt='(7F8.3,I3)') run(i)%cyl(k)%p(1:3),run(i)%cyl(k)%u,run(i)%cyl(k)%t,run(i)%cyl(k)%h,run(i)%cyl(k)%r,run(i)%cyl(k)%nr_con
			end do
			!close(55)
			close(iunit)
			!print*, "done yks"; read*
		end do
		print*, "done writing data for output!", nrun
	end subroutine write_cylinders_information_for_3d_plotting
	!
	subroutine write_cylinders_information_for_validation(fnameout)
		implicit none
		character(len=*),intent(in):: fnameout
		integer:: i,iunit,k,cnt
		type(cylinder):: cyl0
		print*, "Write data for validation: ", fnameout
		open(newunit=iunit,file=fnameout)
		write(iunit,fmt='(A)') '# x y z u t h r nrcon pot_hyp pot_con pot_hyp_local pot_hyp_uniform ngal'
		cnt=0
		do i=1,nrun
			do k=1,run(i)%nrcyl
				cyl0%p=run(i)%cyl(k)%p
				cyl0%t=run(i)%cyl(k)%t; cyl0%u=run(i)%cyl(k)%u
				cyl0%h=run(i)%cyl(k)%h; cyl0%r=run(i)%cyl(k)%r
				cyl0%rmin=run(i)%cyl(k)%r; cyl0%rmax=run(i)%cyl(k)%r
				call update_cyl_automatic_params(cyl0)
				call update_potential_for_cylinder(cyl0,only_datapot=.true.)
				write(iunit,fmt='(3F13.5,4F9.5,I3,4F10.5,I6)') run(i)%cyl(k)%p(1:3),run(i)%cyl(k)%u,run(i)%cyl(k)%t,run(i)%cyl(k)%h,run(i)%cyl(k)%r,&
				run(i)%cyl(k)%nr_con,cyl0%pot_hyp*cdata_coeff_hypothesis,cyl0%pot_con*cdata_coeff_variance,cyl0%potarr(1:2)*cdata_coeff_hypothesis,cyl0%ngal
				cnt=cnt+1
			end do
		end do
		close(iunit)
		print*, "Total nr of cylinders: ", cnt
		print*, "done writing data for validation output!", nrun
	end subroutine write_cylinders_information_for_validation
	!
	!
    subroutine read_cylinders_from_realisations()
        implicit none
		integer:: i,k,n,nc
		character(len=:),allocatable:: dir,filename
		logical:: suc
		!
		print*, "=== Read cylinder data from realisations ==="
		potarr=0.0; potarr2=0.0; npotarr2=0
		!
		nrun=(craw_run_last-craw_run_first+1)*( (craw_rs_last-craw_rs_first)/craw_rs_every +1)
		!
		if (allocated(run)) stop "ERR: run is already allocated, cannot read cylinders!"
        allocate(run(1:nrun))
		n=1
		nc=0
		! cycle over different runs
        do k=craw_run_first,craw_run_last
			if (craw_multiple_runs) then
				print*, "tot: ", craw_run_last, " run: ", k
				select case(k)
					case (1:9); dir=craw_root_dir//'/'//craw_root_dir_run(1:len(craw_root_dir_run)-1)//get_file_counter(k, 1)//'/'//craw_file_prefix//'zyls_'
					case (10:99); dir=craw_root_dir//'/'//craw_root_dir_run(1:len(craw_root_dir_run)-2)//get_file_counter(k, 2)//'/'//craw_file_prefix//'zyls_'
					case (100:999); dir=craw_root_dir//'/'//craw_root_dir_run(1:len(craw_root_dir_run)-3)//get_file_counter(k, 3)//'/'//craw_file_prefix//'zyls_'
					case (1000:9999); dir=craw_root_dir//'/'//craw_root_dir_run(1:len(craw_root_dir_run)-4)//get_file_counter(k, 4)//'/'//craw_file_prefix//'zyls_'
					case default; stop "ERR: run number too large!"
				end select
			else
				dir=craw_root_dir//'/'//craw_file_prefix//'zyls_'
			end if
			!
			! cycle over realisations in this run
			do i=craw_rs_first,craw_rs_last,craw_rs_every
				if (craw_old_format) then ! old format for cylinders
					dir(len_trim(dir)-4:len_trim(dir)-4)='c'
					filename=dir//get_file_counter(i, cfile_number_spaces)//'.cyl'
					call read_cylinders_from_old_binary_file(run(n)%cyl,filename,suc)
				else ! new format for cylinders
					filename=dir//get_file_counter(i, cfile_number_spaces)//'.cyl'
					call read_cylinders_from_binary_file(run(n)%cyl,filename,suc)
				end if
				if (suc) then
					run(n)%nrcyl=size(run(n)%cyl)
					nc=nc+run(n)%nrcyl
					n=n+1
				end if
			end do
        end do
		nrun=n-1 ! actual size of run array
		print*, "Runs: ", craw_run_last-craw_run_first+1, " RS: ", (craw_rs_last-craw_rs_first)/craw_rs_every+1
		print*, "Assumed size: ", (craw_run_last-craw_run_first+1)*( (craw_rs_last-craw_rs_first)/craw_rs_every +1), " Actual: ", nrun
		print*, "Total nr of cylinders: ", nc, " M: ", real(nc)/1E6
		print*, "Expected memory usage (GB): ", (real(nc)/1E6)*0.147
		total_cyls=nc
		! cylinder datapot statistics
		print*
		print*, "--- All cylinders ---"
		print*, "DP-hyptest: local,uniform: ", real(potarr(1:2)/total_cyls)
		print*, "DP: hyptest - consentrat:  ", real(sum(potarr(1:2))/total_cyls),real(potarr(3)/total_cyls)
		print*, "DP: hypt - cons (w coef):  ", real(sum(potarr(1:2))/total_cyls)*real(cdata_coeff_hypothesis),&
		real(potarr(3)/total_cyls)*real(cdata_coeff_variance)
		print*
		print*, "--- 2*ngal_min cylinders ---"
		print*, "DP2-hyptest: local,uniform: ", real(potarr2(1:2)/npotarr2)
		print*, "DP2: hyptest - consentrat:  ", real(sum(potarr2(1:2))/npotarr2),real(potarr2(3)/npotarr2)
		print*, "DP2: hypt - cons (w coef):  ", real(sum(potarr2(1:2))/npotarr2)*real(cdata_coeff_hypothesis),&
		real(potarr2(3)/npotarr2)*real(cdata_coeff_variance)
		!
		call init_fast_cyl_finder()
    end subroutine read_cylinders_from_realisations
	!
	subroutine init_fast_cyl_finder()
		implicit none
		integer,parameter:: nmax_sort = 4500000 ! how many cylinders per one sort
		! 5.6 miljon silindrit on liiga palju ... 4.5 on juba ok!
		real(rk),dimension(:),allocatable:: x,y,z
		real(rk):: rsearch,dum
		integer:: i,n,k,ii
		integer:: idelta=80 ! kui mitu realisatsiooni on yhes sorteerimises
		integer:: nnfast
		!
		! estimate nr of realisations per one sort
		idelta=max(1,floor(real(nmax_sort)/(1.5*(real(total_cyls)/real(nrun))))) ! 1.5 is 50% fluctuations (safety)
		if (idelta>=nrun) then
			nnfast=1
		else
			nnfast=ceiling(real(nrun)/real(idelta))
		end if
		allocate(iifast(1:nnfast))
		!
		print*, "-------- initialise fast cylinder finder --------"
		! calculate cube size for fast search
		dum=0.5*cint_cyl_connection_radius/cos(0.5*acos(cint_parallel_cosi)) + cdata_cyl_rad_max*tan(0.5*acos(cint_parallel_cosi))
		rsearch=sqrt( (cdata_cyl_rad_max)**2 + (0.5*cdata_cyl_len_max+dum)**2 )*1.01 ! cylinder diagonal + small correction (safety)
		csearch_min_rad=rsearch ! minimum search radius
		!
		ii=0
		do
			ii=ii+1
			if (1+(ii-1)*idelta>nrun) exit
			if (ii>nnfast) stop "ERR: something is wrong (init fast cyl finder)!!!"
			n=0
			do i=1+(ii-1)*idelta, ii*idelta
				if (i>nrun) exit
				n=n+run(i)%nrcyl
			end do
			if (n>nmax_sort) stop "ERR: init fast cyl finder... 50% fluctuation is not enough!" 
			allocate(x(1:n),y(1:n),z(1:n))
			allocate(iifast(ii)%irun(1:n),iifast(ii)%icyl(1:n))
			print*, "Init cyl: ", nnfast,'/', ii, real(n)/1E6,'M cyl'
			n=0
			do i=1+(ii-1)*idelta, ii*idelta
				if (i>nrun) exit
				do k=1,run(i)%nrcyl
					n=n+1
					x(n)=run(i)%cyl(k)%p(1)
					y(n)=run(i)%cyl(k)%p(2)
					z(n)=run(i)%cyl(k)%p(3)
					iifast(ii)%irun(n)=i
					iifast(ii)%icyl(n)=k
				end do
			end do
			!
			call init_cubesort(x,y,z,rsearch,iifast(ii)%n)
			deallocate(x,y,z)
		end do
	end subroutine init_fast_cyl_finder
	!
	subroutine get_nearest_cyls(pts,ivec,rsearch)
		use quicksort
		implicit none
		real(rk),dimension(1:3),intent(in):: pts ! input point
		type(iifast_type),intent(out):: ivec
		real(rk),intent(in):: rsearch
		real(rk):: dist1,dist2
		integer:: i,n,k,nn,ii
		logical,dimension(:),allocatable:: mask,good
		type testtype
			integer,dimension(:),allocatable:: ivec
		end type testtype
		type(testtype),dimension(1:100):: tt
		!
		n=0
		do ii=1,size(iifast)
			call extract_points_index(iifast(ii)%n,pts,rsearch,tt(ii)%ivec)
			n=n+size(tt(ii)%ivec)
		end do
		ivec%n=n
		if (n==0) return
		allocate(ivec%irun(1:n),ivec%icyl(1:n))		
		n=0
		do ii=1,size(iifast)
			do i=1,size(tt(ii)%ivec)
				n=n+1
				ivec%irun(n)=iifast(ii)%irun(tt(ii)%ivec(i))
				ivec%icyl(n)=iifast(ii)%icyl(tt(ii)%ivec(i))
			end do
		end do
		!
		call qsort(ivec%irun,ivec%icyl)
	end subroutine get_nearest_cyls
	!
	!------------------------------------------------------------
    subroutine read_cylinders_from_binary_file(cyl,filename,suc)
        implicit none
        character(len=*),intent(in):: filename
		type(cylinder):: cyl0
		type(cylinder_raw),dimension(:),allocatable,intent(inout):: cyl
		logical,intent(out):: suc
        integer:: iunit,i,ierr
		integer:: nn_active,nn_fixed,nline,nn
        ! read cylinders from file
		!print*, filename
		open(newunit=iunit,file=filename,form='unformatted',action='READ',iostat=ierr)
		if (ierr/=0) then
			suc=.false.
			return
		end if
		suc=.true.
		read(iunit,iostat=ierr) nn_active,nn_fixed,nline
		if (ierr/=0) then
			suc=.false.
			return
		end if
		nn=nn_active+nn_fixed
		allocate(cyl(1:nn))
        do i=1,nn
            read(iunit,iostat=ierr) cyl0%fixed, cyl0%p,cyl0%t,cyl0%u,cyl0%h,cyl0%r,cyl0%rmin,cyl0%rmax,&
            cyl0%nr_con,cyl0%nr_of_death, cyl0%nr_of_change, cyl0%nr_of_change_accept
			if (ierr/=0) then
				suc=.false.
				deallocate(cyl)
				return
			end if
			!
			call update_cyl_automatic_params(cyl0)
			call update_potential_for_cylinder(cyl0,only_datapot=.true.)
			! update statistics array
			potarr(1)=potarr(1)+cyl0%potarr(1)
			potarr(2)=potarr(2)+cyl0%potarr(2)
			potarr(3)=potarr(3)+cyl0%pot_con
			if (cyl0%ngal>=2*cdata_minpts) then
				potarr2(1)=potarr2(1)+cyl0%potarr(1)
				potarr2(2)=potarr2(2)+cyl0%potarr(2)
				potarr2(3)=potarr2(3)+cyl0%pot_con
				npotarr2=npotarr2+1
			end if
			!
			cyl(i)%p=cyl0%p
			cyl(i)%u=cyl0%u; cyl(i)%t=cyl0%t
			cyl(i)%h=cyl0%h; cyl(i)%r=cyl0%r
			cyl(i)%eu=cyl0%eu
			cyl(i)%nr_con=cyl0%nr_con
			cyl(i)%potdata=cyl0%potdata
			!
			cyl(i)%r=cyl(i)%r*craw_rsmooth
			! increase cylinder length to take into account connection region
			cyl(i)%dh=0.5*cint_cyl_connection_radius/cos(0.5*acos(cint_parallel_cosi)) + cyl(i)%r*tan(0.5*acos(cint_parallel_cosi))
        end do
        close(iunit)
    end subroutine read_cylinders_from_binary_file
	!
	! subroutine to read cylinders from old format...
	! ... needed to reanalyse the old datasets
	subroutine read_cylinders_from_old_binary_file(cyl,filename,suc)
		implicit none
        character(len=*),intent(in):: filename
		type(cylinder):: cyl0
		type(cylinder_raw),dimension(:),allocatable,intent(inout):: cyl
		logical,intent(out):: suc
        integer:: iunit,i,ierr
		integer:: nn
        ! read cylinders from file
		open(newunit=iunit,file=filename,form='unformatted',action='READ',iostat=ierr)
		if (ierr/=0) then
			suc=.false.
			return
		end if
		suc=.true.
		read(iunit,iostat=ierr) nn
		if (ierr/=0) then
			suc=.false.
			return
		end if
		allocate(cyl(1:nn))
        do i=1,nn
            read(iunit,iostat=ierr) cyl0%p,cyl0%t,cyl0%u,cyl0%h,cyl0%r,&
            cyl0%nr_con,cyl0%ngal,cyl0%potdata,cyl0%potint, cyl0%potarr, &
			cyl0%nr_of_death, cyl0%nr_of_change, cyl0%nr_of_change_accept
			if (ierr/=0) then
				suc=.false.
				deallocate(cyl)
				return
			end if
			call update_cyl_automatic_params(cyl0)
			cyl0%rmin=cyl0%r; cyl0%rmax=cyl0%r
			call update_potential_for_cylinder(cyl0,only_datapot=.true.)
			! update statistics array
			potarr(1)=potarr(1)+cyl0%potarr(1)
			potarr(2)=potarr(2)+cyl0%potarr(2)
			potarr(3)=potarr(3)+cyl0%pot_con
			if (cyl0%ngal>=2*cdata_minpts) then
				potarr2(1)=potarr2(1)+cyl0%potarr(1)
				potarr2(2)=potarr2(2)+cyl0%potarr(2)
				potarr2(3)=potarr2(3)+cyl0%pot_con
				npotarr2=npotarr2+1
			end if
			!
			cyl(i)%p=cyl0%p
			cyl(i)%u=cyl0%u; cyl(i)%t=cyl0%t
			cyl(i)%h=cyl0%h; cyl(i)%r=cyl0%r
			cyl(i)%eu=cyl0%eu
			cyl(i)%nr_con=cyl0%nr_con
			cyl(i)%potdata=cyl0%potdata
			!
			cyl(i)%r=cyl(i)%r*craw_rsmooth
			cyl(i)%dh=0.5*cint_cyl_connection_radius/cos(0.5*acos(cint_parallel_cosi)) + cyl(i)%r*tan(0.5*acos(cint_parallel_cosi))
        end do
        close(iunit)
	end subroutine read_cylinders_from_old_binary_file
	!
end module read_cylinders