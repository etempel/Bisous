!============
! Author: Elmo Tempel
! Date (last changed): 11.02.2015
!============
!
module preparation
    use constants
    use cylinders
    use galaxies
    use operations
    use parameters
    use data_term
    use interaction_term
	use utilities
	implicit none
contains
	!
	! read and write cylinder files
	subroutine convert_cyl_files_to_ascii()
		implicit none
        character(len=:),allocatable:: fileout,filein
		!character(len=*),parameter:: root='/Volumes/wrk/aai/elmo/Saienss/mcmcfils/sdss_r05b/run050/'
		!character(len=*),parameter:: root='/Volumes/wrk/aai/elmo/Saienss/twomass/bisous/obs1/'
		character(len=*),parameter:: root='/Users/elmo/Saienss/boris_clusters/results/'
		!character(len=*),parameter:: root='/Users/elmo/Temp/'
        integer:: iunit,nn,i,k,iunit2
		type(cylinder):: cyl
		PRINT*, "Convert cyl files to ascii files..."
        ! write cylinders to file
		do k=200,500
			!
			filein=root//'bm_cyls_'//get_file_counter(k, cfile_number_spaces)//'.cyl'
			fileout=root//'acyls_'//get_file_counter(k, cfile_number_spaces)//'.txt'
			!
	        open(newunit=iunit,file=filein,form='unformatted')
			open(newunit=iunit2,file=fileout)
			write(iunit2,fmt='(A)') "# x y z t u h pot pot_con pot_hyp ncon ngal pu ph"
			read(iunit) nn,nn,nn
	        do i=1,nn
	            read(iunit) cyl%fixed, cyl%p,cyl%t,cyl%u,cyl%h,cyl%r,cyl%rmin,cyl%rmax,&
	            cyl%nr_con,cyl%nr_of_death, cyl%nr_of_change, cyl%nr_of_change_accept
				!
				write(iunit2, fmt='(9F15.5,I3,I6,2F15.5)') cyl%p,cyl%t,cyl%u,cyl%h,cyl%potdata,cyl%potarr(2:3),&
				cyl%nr_con,cyl%ngal,cyl%potarr(1), exp(-cyl%pot_hyp)/cyl%potarr(1)
	        end do
	        close(iunit)
			close(iunit2)
		end do
	end subroutine convert_cyl_files_to_ascii
	!
	subroutine convert_cyl_files_to_ascii_v2()
		implicit none
		logical,parameter:: makedirs=.false.
        character(len=:),allocatable:: fileout,filein
		! jargnevad kaks rida on kaspari jaoks
		character(len=*),parameter:: root='/Volumes/wrk/bolshoi/test/'
		character(len=*),parameter:: output='/Volumes/wrk/bolshoi/visual/vmap/'
		!character(len=*),parameter:: root='/Users/elmo/Saienss/bisous/test_near2/'
		!character(len=*),parameter:: output='/Users/elmo/Saienss/bisous/test_near2/'
		!
		character(len=:),allocatable:: dir,dirout
        integer:: iunit,nn,i,k,iunit2,irun
		integer:: nn_active,nn_fixed,nline,nrcon
		type(cylinder):: cyl
		real(rk),dimension(1:3):: p1,p2
		real(rk):: rdum
		PRINT*, "Convert cyl files to ascii files for processing..."
        ! write cylinders to file
		!if (makedirs) call system('mkdir '//output)
		do irun=1,1!50
			!dir=root//'run'//get_file_counter(irun, 3)//'/bm_cyls_'
			dir=root//'/bm_zyls_'
			!dirout=output//'run'//get_file_counter(irun, 3)
			dirout=output
			!if (makedirs) call system('mkdir '//dirout)
			print*, "generating run: ", irun
			!
			do k=1,3,1
				!
				filein=dir//get_file_counter(k, cfile_number_spaces)//'.cyl'
				fileout=dirout//'/cyls_'//get_file_counter(k, cfile_number_spaces)//'.txt'
				!
		        open(newunit=iunit,file=filein,form='unformatted')
				open(newunit=iunit2,file=fileout)
				write(iunit2,fmt='(A)') "# x1 y1 z1 x2 y2 z2 x0 y0 z0 t u h pot ncon ngal rad pot_hyp pot_con"
				read(iunit) nn_active,nn_fixed,nline
				nn=nn_active+nn_fixed
		        do i=1,nn
		            read(iunit) cyl%fixed, cyl%p,cyl%t,cyl%u,cyl%h,cyl%r,cyl%rmin,cyl%rmax,&
		            cyl%nr_con,cyl%nr_of_death, cyl%nr_of_change, cyl%nr_of_change_accept
					!
					p1(1)=cyl%p(1)-0.5*cyl%h*cos(cyl%t)*sin(acos(cyl%u))
					p1(2)=cyl%p(2)-0.5*cyl%h*sin(cyl%t)*sin(acos(cyl%u))
					p1(3)=cyl%p(3)-0.5*cyl%h*cyl%u
					p2(1)=cyl%p(1)+0.5*cyl%h*cos(cyl%t)*sin(acos(cyl%u))
					p2(2)=cyl%p(2)+0.5*cyl%h*sin(cyl%t)*sin(acos(cyl%u))
					p2(3)=cyl%p(3)+0.5*cyl%h*cyl%u
					cyl%ngal=1
					!
					!if (cyl%nr_con<2) cycle
					nrcon=cyl%nr_con
					cyl%rmin=cyl%r*0.999; cyl%rmax=cyl%r*1.001
					rdum=cyl%r
					call update_cyl_automatic_params(cyl)
					call update_potential_for_cylinder(cyl)
					!print*, "end potential calculation.... ", rdum
					!if (nrcon>1) read*
					!
					cyl%nr_con=nrcon
					!
					write(iunit2, fmt='(13F15.5,I3,I6,3F15.5)') p1,p2, cyl%p,cyl%t,cyl%u,cyl%h,cyl%potdata,&
					cyl%nr_con,cyl%ngal,cyl%r,cyl%pot_hyp,cyl%pot_con
		        end do
		        close(iunit)
				close(iunit2)
			end do
		end do
		stop "end write data for processing..."
	end subroutine convert_cyl_files_to_ascii_v2
	!
	! set sampling volume fraction for mcmc sampling
	subroutine set_sampling_volume_fraction_for_mcmc(acc,maxtime)
		implicit none
		real(rk),parameter:: vff=0.05 ! assumed volume filling fraction, just for output...
		real(rk),intent(in):: acc ! required accuracy, e.g. 0.05
		real(rk),intent(in):: maxtime ! maximum time in minutes
		!
		integer,parameter:: n_test=20
		integer,parameter:: n_moves0=5000
		integer,parameter:: n_sample0=500
		!
		integer:: i,k,j,nsample,nmoves
		real(rk),dimension(:,:),allocatable:: ipts
		real(rk),dimension(1:3):: p
		integer:: ngal,ncounter,ncounter1,nmovestot
		integer:: ncounter_tot,ncounter1_tot,iout,nnn
		type(cylinder):: cyl
		logical:: test,test2
		real(rk):: dum,apot,apot1,avgpot,avgpot1,xtegold,rsearch,voltot
		real(rk):: frac,f2,fncyl,vcyl,xteg1,t1,t2,vol,vol2,rmin,rmax,vcyl2
		integer*8:: ncnt_all,ncnt_pass
		!
		print*, "set sampling volume fraction: "
		!
		ncnt_all=0; ncnt_pass=0
		ncounter_tot=0; ncounter1_tot=0
		apot=0.0; apot1=0.0
		nmovestot=0
		iout=0; xtegold=0.0
		vol=0.0; vol2=0.0
		call get_time_min(t1)
		do k=1,n_test
			nsample=n_sample0*k
			nmoves=n_moves0*k
			ncounter=0; ncounter1=0
			nmovestot=nmovestot+nmoves
			do i=1,nmoves
				call get_cylinder_random_rh_adaptive(cyl%h,rmin,rmax)
				rsearch=sqrt( (rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cyl%h)**2 ) ! cylinder diagonal
				infdo: do
		            call random_number(p)
					p=cmask%pmin_active+p*(cmask%pmax_active-cmask%pmin_active)
					test=test_sample_mask_location(p)
					ncnt_all=ncnt_all+1
					if (test) then
						call extract_galaxy_coordinates(p,ipts,ngal,rsearch)
						if (ngal>=cdata_minpts) then
							ncnt_pass=ncnt_pass+1
							exit infdo
						end if
					end if
					!
				end do infdo
				!
				cyl%p=p
				do j=1,nsample
					call random_number(cyl%t); cyl%t=cyl%t*twopi
                    call random_number(cyl%u)
                    call update_cyl_automatic_params(cyl)
					cyl%rmin=rmin; cyl%rmax=rmax
					dum=calc_data_potential_for_cylinder(cyl,ipts,ngal)
					if (j==1 .and. dum<cdata_pot_undefined) then
						ncounter1=ncounter1+1
						ncounter1_tot=ncounter1_tot+1
						apot1=apot1+dum
						!
						vol=vol + pi*cyl%h*cyl%r**2
						!vol=vol + 4* pi*cdata_cyl_len_min*cdata_cyl_rad_min**2
						!vol=vol+10.2835/3
					end if
                    if (dum<cdata_pot_undefined) then
						ncounter=ncounter+1
						ncounter_tot=ncounter_tot+1
						apot=apot+dum
						vol2=vol2 + pi*cyl%h*cyl%r**2
                        exit
                    end if
				end do
			end do
			!
			f2=real(ncounter_tot,kind=rk)/real(nmovestot,kind=rk)
			frac=real(ncounter1_tot,kind=rk)/real(nmovestot,kind=rk)
			avgpot=apot/real(ncounter_tot,kind=rk)
			avgpot1=apot1/real(ncounter1_tot,kind=rk)
			vcyl=vol/real(ncounter1_tot,kind=rk)
			vcyl2=vol2/real(ncounter_tot,kind=rk)
			!
			voltot=product(cmask%pmax_active-cmask%pmin_active)*real(ncnt_pass,kind=rk)/real(ncnt_all,kind=rk)
			fncyl=voltot*f2*vff!/vcyl
			!xteg1=(f2/(frac*vcyl))*exp(avgpot1)
			xteg1=(f2/(frac))*exp(avgpot1)
			!
			!print*, "mov:sam:frac ", nmoves,nsample,real(ncounter)/real(nmoves),real(ncounter1)/real(nmoves)
			!print*, "sum    :frac ", nmoves,nsample,real(f2),real(frac)
			!print*, "potential: ", k,real(avgpot),real(avgpot1)
			!print*, "ncyl, xteg: ", nint(fncyl),real(xteg1)
			!print*
			!
			call get_time_min(t2)
			if (t2-t1>maxtime) exit
			if (abs(xteg1-xtegold)/xteg1<acc) then
				iout=iout+1
			else
				iout=0
			end if
			if (iout>1) exit
			xtegold=xteg1
		end do
		!
		print*, "frac, mult:: ", real(frac), real(xteg1*cmcmc_auto_volume_multiplier)
		print*, "Expected nr of cyls:", nint(fncyl)
		print*, "Average potentials: ", real(avgpot1),real(avgpot),real(exp(avgpot1))
		print*, "f2,vcyl,vcyl2 ", real(f2),real(vcyl),real(vcyl2)
		print*, "Volume and frac: ", real(product(cmask%pmax_active-cmask%pmin_active)),real(ncnt_pass)/real(ncnt_all)
		!
		! xteg1,frac,nfcyl on parameetrid, mida vaja valjundis
		cmcmc_volume_fraction=frac
		cmcmc_volume_multiplier=xteg1*cmcmc_auto_volume_multiplier
		cmcmc_norm_volume=voltot
	end subroutine set_sampling_volume_fraction_for_mcmc
	!
	!
	! this is testing subroutine to speed up internal calculations
	subroutine test_move_n_orientation_time()
		implicit none
		integer,parameter:: n_test=5000
		integer,parameter:: n_max_moves=40
		real(rk):: t1,t2
		integer:: i,n,j
		real(rk),dimension(:,:),allocatable:: ipts
		real(rk),dimension(1:3):: p
		integer:: ngal,ntot
		logical:: test,olemas
		type(cylinder):: cyl
		real(rk):: dum
		print*, "Start testing move n orientation time..."
		!
		do i=1,n_max_moves
			n=0; ntot=0
			call get_time_min(t1)
			do
				do
		            call random_number(p)
					p(1)=gal%xmin+p(1)*(gal%xmax-gal%xmin)
		            p(2)=gal%ymin+p(2)*(gal%ymax-gal%ymin)
		            p(3)=gal%zmin+p(3)*(gal%zmax-gal%zmin)
					test=test_sample_mask_location(p)
					if (test) exit
				end do
				!!call extract_galaxy_coordinates(p,ipts,ngal,c_ext_gal_fix_center)
				if (ngal<cdata_minpts) cycle
				ntot=ntot+1
				cyl%p=p
				olemas=.false.
				do j=1,i
                    !call get_cylinder_random_rh(cyl%r,cyl%h)
					call random_number(cyl%t); cyl%t=cyl%t*twopi
                    call random_number(cyl%u)
                    call update_cyl_automatic_params(cyl)
					cyl%rmin=0.5; cyl%rmax=0.5
                    dum=calc_data_potential_for_cylinder(cyl,ipts,ngal,loctest=.true.)
                    if (dum<cdata_pot_undefined) then
                        olemas=.true.
						n=n+1
                        exit
                    end if
				end do
				if (n==n_test) exit
			end do
			call get_time_min(t2)
			print*, "nmoves,ntot,time: ", i,ntot,t2-t1
			t1=t2
		end do
	end subroutine test_move_n_orientation_time
	!
	subroutine analyse_adaptive_data_energy()
		implicit none
		real(rk),parameter:: rmin=0.1
		real(rk),parameter:: rmax=2.0
		integer,parameter:: nrp=600
		real(rk),parameter:: teg=1.00
		integer:: i,k,n,ngal
		real(rk),dimension(:,:),allocatable:: ipts
		real(rk):: diag,dum,dum1,dum2,pmin
		integer:: iu,i2,n2,n1,ki(1:1),n0,j,kmin
		real(rk),dimension(1:nrp):: arr,rr
		logical:: test,test2
		!
		print*, "Start adaptive data energy analyse..."
		!
		cdata_cyl_len_max=15.0
		!
		diag=sqrt( (rmax*(1.0_rk+cdata_cyl_shadow))**2 + (0.5*cdata_cyl_len_max)**2 ) ! cylinder diagonal
		!
		open(newunit=iu,file='test.txt')
		open(newunit=i2,file='test2.txt')
		n=0; n2=0; n1=0; n0=0
		do i=1,nr_cylarr
			if (cylarr(i)%h<3.0 .or. cylarr(i)%nr_con<2 .or. cylarr(i)%ngal<3 .or.&
			cylarr(i)%nr_of_death<0 .or. cylarr(i)%nr_of_change<0) cycle
			n=n+1
			!if (n>100) exit
			!
			cylarr(i)%h=cdata_cyl_len_max
			!call random_number(cylarr(i)%t); cylarr(i)%t=cylarr(i)%t*twopi
            !call random_number(cylarr(i)%u)
			!
			call extract_galaxy_coordinates(cylarr(i)%p,ipts,ngal,diag)
			do k=1,nrp
				cylarr(i)%r=rmin+real(k-1)*(rmax-rmin)/real(nrp-1)
				call update_cyl_automatic_params(cylarr(i))
				cylarr(i)%rmin=0.5; cylarr(i)%rmax=0.5
				dum=calc_data_potential_for_cylinder(cylarr(i),ipts,ngal)
				if (dum<100) write(iu,fmt=*) cylarr(i)%r,dum
				rr(k)=cylarr(i)%r
				arr(k)=dum
			end do
			write(iu,fmt=*); write(iu,fmt=*)
			!
			!test=.false.; test2=.false.
			!do k=1,nrp-1
			!	if (.not.test .and. arr(k)>arr(k+1)*teg .and. arr(k)<100) then
			!		!write(i2,fmt=*) rr(k),arr(k)
			!		test=.true.
			!		n1=n1+1
			!		!exit
			!	end if
			!	if (test .and. arr(k)*teg<arr(k+1) .and. arr(k+1)<100) then
			!		write(i2,fmt=*) rr(k),arr(k),cylarr(i)%h
			!		n2=n2+1
			!		test2=.true.
			!		exit
			!	end if
			!end do
			!if (.not.test2 .or. .false.) then
			!	ki=minloc(arr); k=ki(1)
			!	write(i2,fmt=*) rr(k),arr(k),cylarr(i)%h
			!	n0=n0+1
			!end if
			!
			ki=minloc(arr); k=ki(1)
			do j=k,nrp
				if (arr(j)>99) exit
			end do
			arr(j:)=101.0
			!
			!ki=maxloc(arr,mask=arr<100); j=ki(1)
			test=.false.
			do j=1,nrp-1
				if (arr(j)<100 .and. arr(j)>arr(j+1)) then
					test=.true.
					ki(1)=j
				end if
			end do
			if (.not.test) then
				do j=k,nrp
					if (arr(j)<100) then
						ki(1)=j
					end if
					if (arr(j)>99) exit
				end do
			end if
			j=ki(1)
			!
			pmin=arr(j)
			test=.false.
			do k=j,nrp-1
				if (arr(k)<pmin .and. arr(k)<100) then
					pmin=arr(k)
					kmin=k
					test=.true.
				end if
				if (arr(k)<arr(k+1)) exit
			end do
			!
			if (.not.test) then
				do k=j-1,2,-1
					if (arr(k)<pmin .and. arr(k)<100) then
						pmin=arr(k)
						kmin=k
						test=.true.
					end if
					if (arr(k)<arr(k-1)) exit
				end do
			end if
			!
			if (test) then
				n1=n1+1
				ki=minloc(arr); kmin=ki(1)
				write(i2,fmt=*) rr(kmin),arr(kmin),i!cylarr(i)%ngal
			end if
			!
			!do k=nrp-1,1,-1
			!	if (arr(k)>arr(k+1)) then
			!		write(i2,fmt=*) rr(k+1),arr(k+1)
			!		exit
			!	end if
			!end do
		end do
		close(iu); close(i2)
		!
		print*, "Kokku kasulikke silindreid: ", n,n1,n2
		print*, "                            ", n2+n0,n0
		stop "... end test adaptive data energy ..."
	end subroutine analyse_adaptive_data_energy
	!
	! generate points for periodic boxes
	subroutine generate_points_for_periodic_boxes()
		implicit none
		character(len=*),parameter:: fname='example_halos.txt'
		character(len=*),parameter:: fout='example_halos_periodic.txt'
		integer,parameter:: nn=292289
		real(rk),parameter:: rmax=64.0
		real(rk),parameter:: radd=5.1
		integer:: ii,i,j,k
		real(rk),dimension(1:3):: p,p0
		print*, "generate points for periodic box"
		open(1,file=fname)
		open(2,file=fout)
		read(1,fmt=*)
		do ii=1,nn
			read(1,fmt=*) p0
			!
			do i=-1,1,1
				p(1)=p0(1)+i*rmax
				do j=-1,1,1
					p(2)=p0(2)+j*rmax
					do k=-1,1,1
						p(3)=p0(3)+k*rmax
						if (all(p>-radd) .and. all(p<rmax+radd)) then
							write(2,fmt=*) p
						end if
					end do
				end do
			end do
		end do
		close(1)
		close(2)
	end subroutine generate_points_for_periodic_boxes
	!
	! generate Poisson sample for testing
	subroutine generate_poisson_sample(xsize,ngal,filename)
		implicit none
		real(rk),intent(in):: xsize
		integer,intent(in):: ngal
		character(len=*),intent(in):: filename
		integer:: i,iunit
		real(rk),dimension(1:3):: p
		call random_seed()
		open(newunit=iunit,file=filename)
		do i=1,ngal
			call random_number(p)
			p=p*xsize
			write(iunit,fmt=*) real(p)
		end do
		close(iunit)
	end subroutine generate_poisson_sample
end module preparation


