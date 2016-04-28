!============
! Author: Elmo Tempel
! Date: 28.04.2016
!============
!
module post_processing
	use constants
	use parameters
	use cylinders
	use galaxies
	use cubesort
	use read_cylinders
	!use fits_op
	use utilities
	use quicksort
	use fil_cat
	!
	! parameters for spine extraction
	type(fil_spine),private:: fdum
	integer,private:: n1_fdum,n2_fdum ! fdum real size
	real(rk),dimension(:,:),allocatable,private:: filpts ! temporary filament points for cubesort
	integer,dimension(:,:),allocatable,private:: ifilpts ! index for filament points
	integer,private:: nfilpts ! nr of temporary filament points
	integer,private:: ii_filpts ! cubesort index
	!
contains
	!----- several test, specific and temporary subroutines
	!
	subroutine test_bolshoi_visitmap_movie()
		implicit none
		real(rk),dimension(1:3):: p,dori
		real(rk):: vmap,wmap,dstrength
		integer:: i,j,k
		real(rk),dimension(1:5):: dum
		real(rk),dimension(:,:,:),allocatable:: fmap,fstr
		!
		print*, "calc boslhoi visitmap movie..."
		allocate(fmap(1:1250,1:1250,1:11),fstr(1:1250,1:1250,1:11))
		!
		!print*, "bolshoi visitmap test...."
		open(11,file=c_output_spines//'movie.bin',form='unformatted')
		do i=1,size(fmap,dim=1)
			!print*, "progress: ", i
			p(1)=0.1 + (i-1)*0.2_rk
			do j=1,size(fmap,dim=2)
				p(2)=0.1 + (j-1)*0.2_rk
				do k=1,size(fmap,dim=3)
					p(3)=124.0+(k-1)*0.2_rk
					!
					call get_vmap_and_ori_for_point(p,vmap,wmap,dori,dstrength,only_vmap=.false.)
					!dum(k)=vmap
					fmap(i,j,k)=vmap/nrun
					fstr(i,j,k)=dstrength
				end do
				!dum=dum/nrun
				!write(11,fmt='(2F8.2,3F8.5)') p(1:2), dum(3), sum(dum(2:4))/3.0, sum(dum(1:5))/5.0
			end do
			!write(11,fmt=*)
		end do
		!
		write(11) size(fmap,dim=1), size(fmap,dim=2), size(fmap,dim=3)
		write(11) fmap
		write(11) fstr
		close(11)
	end subroutine test_bolshoi_visitmap_movie
	!
	subroutine test_post_processing()
		implicit none
		real(rk),dimension(1:3):: p,dori
		real(rk),dimension(:,:,:),allocatable:: vmap,wmap
		type(iifast_type),dimension(:,:,:),allocatable:: ic
		integer:: i,n
		real(rk):: t1,t2,vm,wm,dst
		real(rk):: vm0,wm0,dst0,wmax
		print*, "start testing post processing..."
		call test_post_processing_sdss()
		call cpu_time(t1)
		n=0
		vm0=0.0; wm0=0.0; dst0=0.0
		open(1,file='test.txt')
		open(2,file='test2.txt')
		wmax=0.0
		do i=1,1000
			call random_number(p)
			p=cmask%pmin_active+p*(cmask%pmax_active-cmask%pmin_active)
			!
			call get_vmap_and_ori_for_point(p,vm,wm,dori,dst)
			if (vm>0.0 .and. dst>0.0) then
				n=n+1
				vm0=vm0+vm; wm0=wm0+wm; dst0=dst0+dst
				write(1,fmt=*) vm,wm,dst
				write(2,fmt=*) dori(1:3)
				if (vm>wmax) then
					wmax=vm
				end if
			end if
			!
			!call calc_visit_and_orientation_cubes(p,25,5.0_rk,vmap,wmap,ic)
			!print*, "point ", real(p)
			!print*, size(ic(1,1,1)%irun), size(ic(1,1,1)%icyl), ic(1,1,1)%n
			!read*
			!if (size(ic(1,1,1)%irun)>0) n=n+1
			!if (size(ic(1,1,1)%irun)>500) then
			!	print*, "write fits file..."
			!	where(vmap<5) wmap=0.0
			!	call write_cube_to_fits(vmap,'test_vmap.fit')
			!	call write_cube_to_fits(wmap,'test_wmap.fit')
			!	read*
			!end if
		end do
		close(1); close(2)
		call cpu_time(t2)
		print*, "time ", t2-t1,n
		print*, vm0/n,wm0/n,dst0/n
		print*, "vmax ", wmax
	end subroutine test_post_processing
	!
	subroutine test_post_processing_sdss()
		implicit none
		real(rk),dimension(1:3):: p,dori
		real(rk),dimension(:,:,:),allocatable:: vmap,wmap
		type(iifast_type),dimension(:,:,:),allocatable:: ic
		integer:: i,n
		real(rk):: t1,t2,vm,wm,dst
		real(rk):: vm0,wm0,dst0,wmax
		print*, "start testing post processing..."
		call cpu_time(t1)
		n=0
		vm0=0.0; wm0=0.0; dst0=0.0
		open(1,file='sdss_xyz.txt')
		open(2,file='sdss_bis3.txt')
		write(2,fmt='(A)') '# vis wden dst'
		read(1,fmt=*)
		wmax=0.0
		do i=1,588193
			read(1,fmt=*) p
			!
			call get_vmap_and_ori_for_point(p,vm,wm,dori,dst)
			write(2,fmt=*) vm/404.0,wm,dst
		end do
		close(1); close(2)
		call cpu_time(t2)
		print*, "time ", t2-t1,n
		stop "end sdss testing"
	end subroutine test_post_processing_sdss
	subroutine test_post_processing_pks()
		implicit none
		real(rk),dimension(1:3):: p,dori
		real(rk),dimension(:,:,:),allocatable:: vmap,wmap
		type(iifast_type),dimension(:,:,:),allocatable:: ic
		integer:: i,n
		real(rk):: t1,t2,vm,wm,dst
		real(rk):: vm0,wm0,dst0,wmax
		print*, "start testing post processing (pks)..."
		call cpu_time(t1)
		n=0
		vm0=0.0; wm0=0.0; dst0=0.0
		open(1,file='210215_2dF_PKS2155-304_xyz.txt')
		open(2,file='pks_bisous_galaxies.txt')
		write(2,fmt='(A)') '# vis wden dst'
		!read(1,fmt=*)
		wmax=0.0
		do i=1,21471
			read(1,fmt=*) p
			!
			call get_vmap_and_ori_for_point(p,vm,wm,dori,dst)
			write(2,fmt=*) vm/1000.0,wm/1000.0,dst
		end do
		close(1); close(2)
		call cpu_time(t2)
		print*, "time ", t2-t1,n
		stop "end sdss testing"
	end subroutine test_post_processing_pks
	!
	subroutine test_post_processing_bolshoi()
		implicit none
		real(rk),dimension(1:3):: p,dori
		real(rk),dimension(:,:,:),allocatable:: vmap,wmap
		type(iifast_type),dimension(:,:,:),allocatable:: ic
		integer:: i,n
		real(rk):: t1,t2,vm,wm,dst
		real(rk):: vm0,wm0,dst0,wmax
		print*, "start testing post processing (bolshoi)..."
		call cpu_time(t1)
		n=0
		vm0=0.0; wm0=0.0; dst0=0.0
		open(1,file='halos_xyz_sub/halos_sn_1.txt')
		read(1,fmt=*)
		open(2,file='bolshoi_bisous_galaxies.txt')
		write(2,fmt='(A)') '# vis wden dst'
		wmax=0.0
		do i=1,1033718
			read(1,fmt=*) p
			!
			call get_vmap_and_ori_for_point(p,vm,wm,dori,dst)
			write(2,fmt=*) vm/nrun,wm/nrun,dst
		end do
		close(1); close(2)
		call cpu_time(t2)
		print*, "time ", t2-t1,n
		stop "end bolshoi testing"
	end subroutine test_post_processing_bolshoi
	!
	!============================================================
	! ===== siit edasi on subroutined, mis on olulised ning ei ole yhekordsed!
	!============================================================
	!
	subroutine calc_filament_statistics_for_input_points(file_in,file_out)
		implicit none
		character(len=*),intent(in):: file_in,file_out
		integer:: uin,uout ! unit in and out
		integer:: i,nin,ncnt,ierr
		real(rk),dimension(1:3):: crd,dori
		real(rk):: vm,wm,dst
		print*
		print*, "Calculate filament statistics for points..."
		print*, "   Filename: ", file_in
		open(newunit=uin,file=file_in, status='old', iostat=ierr, action='read')
		if (ierr/=0) then
			print*, "ERR: calc_filament_statistics_for_input_points"
			print*, "Input file does not exist or is corrupted: ", file_in
			stop "Stopping the program..."
		end if
		open(newunit=uout,file=file_out)
		write(uout,fmt='(A)') '# x y z  vmap wmap d_str  dx dy dz'
		nin=0; ncnt=1
		do
			301 continue
			if (nin>0) ncnt=ncnt+1 ! jatame valja faili paise
			read(uin,fmt=*,err=301,end=302) crd
			nin=nin+1 ! valid xyz coordinates
			!
			call get_vmap_and_ori_for_point(crd,vm,wm,dori,dst)
			vm=vm/nrun; wm=wm/nrun
			!
			write(uout,fmt='(3F13.5, 3F8.5, 3F12.8)') crd(1:3), vm,wm,dst, dori(1:3)
		end do
		302 continue; ncnt=ncnt-1
		!
		if (nin/=ncnt) then
			print*, "WARNING: number of input and output points do not match!"
			print*, "in, out: ", nin, ncnt
		end if
		close(uin)
		close(uout)
	end subroutine calc_filament_statistics_for_input_points
	!
	subroutine extract_filament_spines(pmin,pmax)
		implicit none
		real(rk),parameter:: rad_acc=0.01 ! relative radius accuracy
		real(rk),parameter:: vmap_safety=0.8 ! safety tegur for vmap
		integer,parameter:: nsm=3000000 ! sort maximum array size
		real(rk),dimension(1:3),intent(in),optional:: pmin,pmax ! minimum and maximum crds for spine cube
		real(rk),dimension(:),allocatable:: xx,yy,zz ! coordinates (grid steps) in xyz axis
		real(rk),dimension(1:3):: cmin,cmax,gp,p,dori,p2
		real,dimension(:,:,:),allocatable:: vmap0 ! initial visitmap for spine extraction
		real(rk),dimension(:,:,:),allocatable:: vmapdum ! dummy cube
		integer:: n,i,k,ii,nrad,cnt,kk,nn,idum,j,kk2
		integer:: irun,icyl,k1,k2,k3
		real(rk):: radmax,vdum,wdum,dstrength,frad,frad0,dist
		integer,dimension(1:3):: imin,imax,nmax
		logical:: inside, testfpt
		integer:: cntold
		type vmap_index
			integer:: ix,iy,iz
			real(rk):: vis
		end type vmap_index
		type(vmap_index),dimension(:),allocatable:: vmi
		type(vmap_index),dimension(:,:),allocatable:: vmidum
		real(rk),dimension(:),allocatable:: vval
		integer,dimension(:),allocatable:: vind
		integer,dimension(1:1):: id
		!
		type(iifast_type),dimension(:,:,:),allocatable:: ic
		type(iifast_type),dimension(:,:),allocatable:: ic2
		real(rk),dimension(:,:,:),allocatable:: vm,wm,crd2
		real(rk),dimension(:,:),allocatable:: vm2,wm2
		real(rk),dimension(:,:,:,:),allocatable:: crd
		real(rk):: step,cosi,tac,aa,bb,cc,rr
		integer,dimension(-1:1):: ecode
		!
		! set cube limits and gridding
		if (present(pmin).and.present(pmax)) then
			cmin=pmin; cmax=pmax
		else
			cmin=cmask%pmin; cmax=cmask%pmax
		end if
		print*
		print*, "===== Extract filament spines ====="
		print*, "Coordinate limits for spine cube:"
		print*, "   MIN:: ", real(cmin)
		print*, "   MAX:: ", real(cmax)
		print*
		n=ceiling((cmax(1)-cmin(1))/craw_spine_initial_dx)+1; allocate(xx(1:n))
		n=ceiling((cmax(2)-cmin(2))/craw_spine_initial_dx)+1; allocate(yy(1:n))
		n=ceiling((cmax(3)-cmin(3))/craw_spine_initial_dx)+1; allocate(zz(1:n))
		nmax(1)=size(xx); nmax(2)=size(yy); nmax(3)=size(zz)
		do i=1,nmax(1)
			xx(i)=cmin(1) + (i-1)*craw_spine_initial_dx
		end do
		do i=1,nmax(2)
			yy(i)=cmin(2) + (i-1)*craw_spine_initial_dx
		end do
		do i=1,nmax(3)
			zz(i)=cmin(3) + (i-1)*craw_spine_initial_dx
		end do
		print*, "Expected memory size (GB): ", real(size(xx))*real(size(yy))*real(size(zz))*4.0/1024.0**3
		allocate(vmap0(1:size(xx),1:size(yy),1:size(zz))); vmap0=0.0
		!
		tac=max(0.5, tan(acos(cint_parallel_cosi)) )
		!
		! generate initial vmap0 for spine extraction
		print*; print*, "Generate initial visitmap for spine extraction..."
		do irun=1,nrun ! cycle over runs
			if (mod(irun,nrun/10)==0) print*, "tot: ", nrun, " run: ", irun
			do icyl=1,run(irun)%nrcyl
				radmax=sqrt( ( run(irun)%cyl(icyl)%r )**2 + (0.5*run(irun)%cyl(icyl)%h+run(irun)%cyl(icyl)%dh)**2 )
				imin=max(1, ceiling((run(irun)%cyl(icyl)%p-radmax - cmin)/craw_spine_initial_dx)+1 )
				imax=min(nmax, ceiling((run(irun)%cyl(icyl)%p+radmax - cmin)/craw_spine_initial_dx) )
				if (any(imin>nmax) .or. any(imax<1)) cycle ! cylinder is outside of box
				!
				!--- cycle over gridpoints
				do k3=imin(3),imax(3)
					gp(3)=zz(k3)
					do k2=imin(2),imax(2)
						gp(2)=yy(k2)
						do k1=imin(1),imax(1)
							gp(1)=xx(k1)
							! gp(1:3) are gridpoint coordinates
							call get_visitmap_value_from_cyl(run(irun)%cyl(icyl),gp,vdum,wdum,inside,rad_min=craw_spine_initial_dx)
							if (inside) then
								vmap0(k1,k2,k3)=vmap0(k1,k2,k3) + vdum
							end if
						end do
					end do
				end do
				!---
			end do
		end do
		! normalise
		vmap0=vmap0/nrun
		!
		print*, "Total grid cells: ", size(vmap0), " visits: ", count(vmap0>vmap_safety*craw_spine_min_visitmap_value), &
		count(vmap0>vmap_safety*craw_spine_min_visitmap_value)/1E6
		! next lines for testing
		!allocate(vmapdum(1:size(xx),1:size(yy),1:size(zz))); vmapdum=vmap0
		!call write_cube_to_fits(vmapdum,'test_vmap.fit')
		!
		! sort gridcells for spine detection:: nsm==3000000 is maximum size for sort
		n=count(vmap0>vmap_safety*craw_spine_min_visitmap_value) ! ===== here is minimum visitmap value !!!!!! =====
		k=floor(real(n,kind=rk)/real(nsm,kind=rk))+1
		allocate(vmidum(1:nsm,1:k)); vmidum(:,:)%vis=0.0
		allocate(vmi(1:n)); i=0; k=1; ii=0
		do k3=1,nmax(3)
			do k2=1,nmax(2)
				do k1=1,nmax(1)
					if (vmap0(k1,k2,k3)>vmap_safety*craw_spine_min_visitmap_value) then
						i=i+1
						if (i>n) stop "ERR: visitmap initialising index is wrong!"
						ii=ii+1
						if (ii>nsm) then
							ii=1
							k=k+1
						end if
						vmidum(ii,k)%ix=k1
						vmidum(ii,k)%iy=k2
						vmidum(ii,k)%iz=k3
						vmidum(ii,k)%vis=vmap0(k1,k2,k3)
					end if
				end do
			end do
		end do
		if (i/=n) stop "ERR: visitmap initialising index is wrong (2)!"
		do ii=1,k
			call qsort(vmidum(:,ii)%vis,vmidum(:,ii)%ix,vmidum(:,ii)%iy,vmidum(:,ii)%iz)
		end do
		allocate(vval(1:k),vind(1:k)); vind=nsm
		vval(1:k)=vmidum(nsm,1:k)%vis
		i=0
		do
			id=maxloc(vval); k=id(1)
			if (vval(k)==0.0) exit
			i=i+1
			if (i>n) then
				print*, i,id
				print*, vind
				print*, real(vval)
				print*, vmidum(vind(k),k)%vis
				stop "ERR: visitmap initialising index is wrong (3)!"
			end if
			vmi(i)%ix=vmidum(vind(k),k)%ix
			vmi(i)%iy=vmidum(vind(k),k)%iy
			vmi(i)%iz=vmidum(vind(k),k)%iz
			vmi(i)%vis=vmidum(vind(k),k)%vis
			vind(k)=vind(k)-1
			if (vind(k)>0) then
				vval(k)=vmidum(vind(k),k)%vis
			else
				vind(k)=vind(k)+1
				vval(k)=0.0
			end if
		end do
		if (i/=n) stop "ERR: visitmap initialising index is wrong (4)!"
		deallocate(vval,vind,vmidum)
		!
		!=================================================================
		!== preparation is finished, start spine extraction
		!==  vmi(1:n) sisaldab endas sorteeritud visitmapi
		!=================================================================
		print*; print*, "----- Start spine extraction -----"
		radmax=0.7*craw_spine_initial_dx ! 0.5 is minimum, 0.7 is to ensure safety
		nrad=ceiling(radmax/craw_spine_location_acc)
		! preparation
		call prepare_fil_spine_arrays()
		! spine extraction starts here
		cnt=0; cntold=0 ! count good starting points
		do i=1,n ! cycle over gridcells...sorted based on visitmaps
			if (mod(i,1000)==0) then
				!print*, "progress: ", n, " / ", i, cnt
				if (mod(i,10000)==0 .and. cnt/=cntold) then
					!call write_basic_catalogue_to_file('test_filament_spines.txt')
					call write_basic_catalogue_to_file(c_output_spines//'points.txt')
					cntold=cnt
				end if
			end if
			!
			k1=vmi(i)%ix; k2=vmi(i)%iy; k3=vmi(i)%iz
			if ( vmap0(k1,k2,k3)<=0.0 ) cycle ! gridpoint is already visited or spine is nearby
			vmap0(k1,k2,k3)=0.0 ! nullify
			gp(1)=xx(k1); gp(2)=yy(k2); gp(3)=zz(k3)
			! gp is gridpoint, where filament starts
			!
			! find the best location around that point
			nrad=ceiling(radmax/craw_spine_location_acc)
			call calc_visit_and_orientation_cubes(gp,nrad,radmax,vm,wm,ic)
			imax=maxloc(wm)
			if ( vm(imax(1),imax(2),imax(3))/nrun<craw_spine_min_visitmap_value ) cycle ! too low visitmap value
			call get_vmap_cubes_coordinates(gp,nrad,radmax,crd)
			p=crd(:, imax(1),imax(2),imax(3))
			if (any( abs(gp-p)>0.5*craw_spine_initial_dx )) cycle ! gridpoint gp is not the nearest gridpoint
			! fine tune the location
			if (craw_spine_location_acc_tune>1) then
				call calc_visit_and_orientation_cubes(p,craw_spine_location_acc_tune,craw_spine_location_acc,vm,wm,ic)
				imax=maxloc(wm)
				if ( vm(imax(1),imax(2),imax(3))/nrun<craw_spine_min_visitmap_value ) cycle ! too low visitmap value
				call get_vmap_cubes_coordinates(p,craw_spine_location_acc_tune,craw_spine_location_acc,crd)
				p=crd(:, imax(1),imax(2),imax(3))
			end if
			vdum=vm(imax(1),imax(2),imax(3)); wdum=wm(imax(1),imax(2),imax(3))
			if (vdum/nrun<craw_spine_min_visitmap_value) cycle
			! best location is p - filament starting point
			!
			! get the orientation info for point p
			call get_best_orientation(ic(imax(1),imax(2),imax(3)),dori,dstrength,frad); frad0=frad
			kk=0
			do
				call get_best_orientation(ic(imax(1),imax(2),imax(3)),dori,dstrength,frad,rad_prior=.true.)
				if (frad<=0.0 .or. abs(frad-frad0)/frad<rad_acc) exit
				frad0=frad
				kk=kk+1
				if (kk>10) stop "ERR: extract spines... orientation from radius is invalid!!!"
			end do
			if (dstrength<craw_spine_min_ori_strength) cycle ! orientation is not defined
			!
			! test, whether other filament spine points are nearby
			testfpt=test_nearby_other_fil_points(p,dori,frad)
			if (testfpt) cycle ! other filament point is nearby (dublicate filament)
			!
			! good starting point for spine have been found
			! p - location :: dori - orientation
			n1_fdum=0; n2_fdum=0
			fdum%fpt(0)%p=p
			fdum%fpt(0)%dvis=dori
			fdum%fpt(0)%vis=vdum/nrun
			fdum%fpt(0)%wvis=wdum
			fdum%fpt(0)%dstr=dstrength
			fdum%fpt(0)%rad=frad
			!---------------------------------------
			!-- first filament point exists... look to both sides
			!---------------------------------------
			do kk2=1,2
				select case (kk2)
					case(1); kk=-1
					case(2); kk=1
				end select
				nn=0
				!print*, "algus ", kk2,kk,nn
				infdo: do !-- infinite loop
					step=min(fdum%fpt(nn)%rad*craw_spine_spacing_rad_frac, craw_spine_spacing_max)/1.4 ! 1.4=sqrt(2)...for safety
					gp=fdum%fpt(nn)%p + kk*step*fdum%fpt(nn)%dvis
					! gp is new potential point
					! dori is new orientation
					nrad=ceiling(2*step*tac*1.2_rk/craw_spine_location_acc)
					p2=gp; idum=0; dori=fdum%fpt(nn)%dvis
					frad=fdum%fpt(nn)%rad
					do
						idum=idum+1
						call calc_visit_and_orientation_in_a_plane(p2,dori,nrad,step*tac*1.2_rk,vm2,wm2,ic2)
						imax(1:2)=maxloc(wm2)
						call get_perpendicular_plane_coordinates(p2,dori,nrad,step*tac*1.2_rk,crd2)
						p=crd2(:, imax(1),imax(2))
						dist=sqrt( sum((p2-p)**2) )
						call get_best_orientation(ic2(imax(1),imax(2)),dori,dstrength,frad,rad_prior=.true.)
						if (dist<0.5*craw_spine_location_acc .or. dstrength==0.0 .or. nrad>100) exit
						p2=p
						if (idum>10) then 
							dist=sqrt( sum((gp-p)**2) )
							if (dist>step*tac) exit
							nrad=nrad*2
							idum=0
							!stop "ERR: coordinates for new point impossible to find!!!"
						end if
					end do
					! fine tuning the location
					if (craw_spine_location_acc_tune>1 .and. dstrength>0.0) then
						call calc_visit_and_orientation_in_a_plane(p2,dori,craw_spine_location_acc_tune,craw_spine_location_acc,vm2,wm2,ic2)
						imax(1:2)=maxloc(wm2)
						call get_perpendicular_plane_coordinates(p2,dori,craw_spine_location_acc_tune,craw_spine_location_acc,crd2)
						p=crd2(:, imax(1),imax(2))
						call get_best_orientation(ic2(imax(1),imax(2)),dori,dstrength,frad,rad_prior=.true.)
					end if
					vdum=vm2(imax(1),imax(2))/nrun; wdum=wm2(imax(1),imax(2))
					!
					if ( vdum<craw_spine_min_visitmap_value ) then
						ecode(kk)=1 !!! too low visitmap value
						exit infdo
 					end if
					if (dstrength<craw_spine_min_ori_strength) then
						ecode(kk)=2 !!! orientation strength is too low
						exit infdo
					end if
					dist=sqrt( sum((gp-p)**2) )
					if (dist>step*tac) then ! must be smaller than step size
						ecode(kk)=3 !!! next point is too much off from previous point
						exit infdo
					end if
					if (frad<fdum%fpt(nn)%rad/cint_concyl_radius_dif .or. frad>fdum%fpt(nn)%rad*cint_concyl_radius_dif) then
						ecode(kk)=4 !!! radius is too much different
						exit infdo
					end if
					cosi=dot_product(dori,fdum%fpt(nn)%dvis)
					if (abs(cosi)<cint_parallel_cosi) then
						ecode(kk)=5 !!! orientation is too different between two points
						exit infdo
					end if
					if (cosi<0.0) dori=-dori ! normalise the vector
					if (abs(nn)>0 .or. n1_fdum<0) then
						aa=sqrt( (p(1)-fdum%fpt(nn)%p(1))**2 + (p(2)-fdum%fpt(nn)%p(2))**2 + (p(3)-fdum%fpt(nn)%p(3))**2 )
						bb=sqrt( (fdum%fpt(nn-1*kk)%p(1)-fdum%fpt(nn)%p(1))**2 + (fdum%fpt(nn-1*kk)%p(2)-fdum%fpt(nn)%p(2))**2 +&
						 (fdum%fpt(nn-1*kk)%p(3)-fdum%fpt(nn)%p(3))**2 )
						cc=sqrt( (fdum%fpt(nn-1*kk)%p(1)-p(1))**2 + (fdum%fpt(nn-1*kk)%p(2)-p(2))**2 + (fdum%fpt(nn-1*kk)%p(3)-p(3))**2 )
						cosi=(aa+bb+cc)/2.0_rk
						rr=4*sqrt(cosi*(cosi-aa)*(cosi-bb)*(cosi-cc))/(aa*bb*cc)
						if (rr>1.0/(frad*1.0) .or. cc<max(aa,bb)) then
							ecode(kk)=6 !!! curvature too large
							!open(1,file='test2.txt')
							!print*,  "#rr ", rr,kk
							!print*, n1_fdum,n2_fdum,nn,nn-1*kk
							!write(1,fmt=*) p
							!write(1,fmt=*) fdum%fpt(nn)%p
							!write(1,fmt=*) fdum%fpt(nn-1*kk)%p
							!write(1,fmt=*); write(1,fmt=*)
							!do j=n1_fdum,n2_fdum
							!	write(1,fmt=*) fdum%fpt(j)%p
							!end do
							!close(1)
							!print*, "print results"
							!read*
							exit infdo
						end if
					end if
					testfpt=test_nearby_other_fil_points(p,dori,frad)
					if (testfpt) then
						 ecode(kk)=7 !!! other filament point is nearby (dublicate filament)
						 exit infdo
					end if
					! check, wheter infinite loop inside current filament
					if (abs(nn)>0 .or. n1_fdum<0) then
						if (kk<0) then
							imin(1)=n1_fdum+1; imin(2)=n2_fdum
						else
							imin(1)=n1_fdum; imin(2)=n2_fdum-1
						end if
						rr=sum( (p-fdum%fpt(nn)%p)**2 )
						testfpt=.false.
						do j=imin(1),imin(2)
							dist=sum( (p-fdum%fpt(j)%p)**2 )
							if (dist<rr) then
								testfpt=.true.
								exit infdo
							end if
						end do
						if (testfpt) then
							 ecode(kk)=8 !!! other current filament point is nearby
							 exit infdo
						end if
					end if
					!
					!----- good filament point exists, add to array -----
					nn=nn+1*kk
					!print*, "pts ", nn,kk
					if (abs(nn)>craw_max_nr_of_fil_points) stop "ERR: too few points in temporary array!!!"
					if (kk<0) then
						n1_fdum=nn
					else
						n2_fdum=nn
					end if
					fdum%fpt(nn)%p=p
					fdum%fpt(nn)%dvis=dori
					fdum%fpt(nn)%vis=vdum
					fdum%fpt(nn)%wvis=wdum
					fdum%fpt(nn)%dstr=dstrength
					fdum%fpt(nn)%rad=frad
				end do infdo !-- end infinite loop
				!print*, "loppi ", kk,nn
			end do ! look to both sides
			!---------------------------------------------------------
			!-- filament exists, tune values and update main array
			!---------------------------------------------------------
			!
			!---------------------------------------------------------
			! remove gridpoints from array
			do j=n1_fdum,n2_fdum
				imin = max(1, nint((fdum%fpt(j)%p - cmin)/craw_spine_initial_dx)+1 )
				imin = min(imin,nmax)
				vmap0(imin(1),imin(2),imin(3))=0.0 ! nullify
			end do
			! add fdum to main array
			if (n2_fdum-n1_fdum+1>4) then
				cnt=cnt+1 ! increase the counter for good filaments
				fdum%ecode(1)=ecode(-1)
				fdum%ecode(2)=ecode(1)
				call add_temporary_spine_to_main_fil()
			end if
		end do ! cycle over gridcells in visitmap
		!
		!print*, "FINAL CATALOGUE TO FILE..."
		!call write_basic_catalogue_to_file('test_filament_spines.txt')
	end subroutine extract_filament_spines
	!
	subroutine prepare_fil_spine_arrays()
		implicit none
		allocate(fil(1:craw_max_nr_of_fil_spines))
		nfil=0
		allocate(fdum%fpt(-craw_max_nr_of_fil_points:craw_max_nr_of_fil_points))
		n1_fdum=0; n2_fdum=0
		allocate(filpts(1:3000000,1:3),ifilpts(1:3000000,1:2)) ! temporary array for filament points
		nfilpts=0
		ii_filpts=0
	end subroutine prepare_fil_spine_arrays
	subroutine add_temporary_spine_to_main_fil()
		implicit none
		integer:: n,i
		n=n2_fdum-n1_fdum+1
		nfil=nfil+1
		allocate(fil(nfil)%fpt(1:n))
		fil(nfil)%npts=n
		n=0
		do i=n1_fdum,n2_fdum
			n=n+1
			fil(nfil)%fpt(n)=fdum%fpt(i)
			fil(nfil)%ecode=fdum%ecode
			nfilpts=nfilpts+1
			filpts(nfilpts,1:3)=fdum%fpt(i)%p
			ifilpts(nfilpts,1)=nfil
			ifilpts(nfilpts,2)=n
		end do
		if (ii_filpts>0) then
			call free_cubesort(ii_filpts); ii_filpts=0
		end if
		if (nfilpts==0) return
		call init_cubesort(filpts(1:nfilpts,1),filpts(1:nfilpts,2),filpts(1:nfilpts,3),craw_spine_initial_dx,ii_filpts)
	end subroutine add_temporary_spine_to_main_fil
	function test_nearby_other_fil_points(p,d,rad) result(res)
		implicit none
		logical:: res
		real(rk),dimension(1:3),intent(in):: p,d ! point and orientation
		real(rk),intent(in):: rad ! radius
		integer,dimension(:),allocatable:: ivec
		real(rk):: rsearch
		integer:: i
		real(rk),dimension(1:3):: dpts,d2
		real(rk):: dist_cnt,dist_axis,teg,rad2,cosi
		!
		res=.false.
		if (ii_filpts<1) return
		teg=0.5*(cint_repulsive_radius_dif+cint_concyl_radius_dif)
		rsearch=1.415*rad*teg
		call extract_points_index(ii_filpts,p,rsearch,ivec,sphere=.true.)
		if (size(ivec)<1) return
		do i=1,size(ivec)
			dpts=p-filpts(ivec(i),:)
	        dist_cnt=dot_product(d,dpts)
	        dist_axis=sqrt( dot_product(dpts,dpts) - dist_cnt*dist_cnt )
			if (dist_axis<=rad*teg .and. dist_cnt<rad*teg) then
				! punkt on filamendi raadiuse sees...edasised testid
				rad2=fil( ifilpts(ivec(i),1) )%fpt( ifilpts(ivec(i),2) )%rad
				d2=fil( ifilpts(ivec(i),1) )%fpt( ifilpts(ivec(i),2) )%dvis
				if (rad2<rad/teg .or. rad2>rad*teg) cycle ! radius too different
				! same scale filament exists -- final test
				rad2=max(rad2,rad)
				cosi=abs( dot_product(d,d2) )
				if (cosi>cint_orthogonal_cosi .and. max(dist_axis,dist_cnt)<rad2) then
					res=.true.
					return ! similar filament exists
				end if
			end if
		end do
	end function test_nearby_other_fil_points
	!
	! this subroutine determines visitmap calculation details
	subroutine get_visitmap_value_from_cyl(cyl,p,vmap,wmap,inside,rad_min)
		implicit none
		type(cylinder_raw),intent(in):: cyl
		real(rk),dimension(1:3),intent(in):: p
		logical,intent(out):: inside
		real(rk),intent(out):: vmap,wmap
		real(rk),intent(in),optional:: rad_min ! minimum radius (needed for spine initialisation)
		real(rk),dimension(1:3):: dpts
		real(rk):: dist_cnt,dist_axis
		real(rk):: w,rad
		!
		if (cyl%nr_con==0 .and. craw_weight_0<=0.0) return ! ignore 0-connected cylinders
		if (cyl%nr_con==1 .and. craw_weight_1<=0.0) return ! ignore 1-connected cylinders
		!
		dpts=p-cyl%p
        dist_cnt=dot_product(cyl%eu,dpts)
        dist_axis=sqrt(dot_product(dpts,dpts) - dist_cnt*dist_cnt)
		!
		rad=cyl%r
		if (present(rad_min)) rad=max(rad_min,cyl%r)
		vmap=0.0; wmap=0.0; inside=.false.
		if (abs(dist_cnt)<=0.5*cyl%h+cyl%dh .and. dist_axis<=rad) then
			inside=.true.
			select case (cyl%nr_con)
				case (0); w=craw_weight_0
				case (1); w=craw_weight_1
				case (2:); w=craw_weight_2
			end select
			vmap=1.0*w
			if (present(rad_min)) return ! no need to calculate further (spine initialisation)
			wmap=exp(-cyl%potdata)*(1.0-dist_axis/cyl%r)*w
		end if
	end subroutine get_visitmap_value_from_cyl
	!
	! get the best orientation... it is one of the input cylinder orientation
	subroutine get_best_orientation(ic,dori,dstrength,frad,rad_prior)
		implicit none
		type(iifast_type),intent(in):: ic ! input list of cylinders
		real(rk),dimension(1:3),intent(out):: dori ! output best orinetation
		real(rk),intent(out):: dstrength ! orientation strength
		real(rk),intent(inout),optional:: frad ! filament average radius
		logical,intent(in),optional:: rad_prior
		integer:: i,k,idum(1:1)
		real(rk):: wsum,w0,cosi,rad,teg
		logical,dimension(:),allocatable:: mask
		real(rk),dimension(:),allocatable:: w,den
		!
		if (ic%n<1) then
			dori=0.0
			dstrength=0.0
			if (present(frad)) frad=0.0
			return
		end if
		allocate(mask(1:ic%n),w(1:ic%n),den(1:ic%n)); mask=.false.
		wsum=0.0; w=0.0; rad=0.0
		teg=0.5*(cint_repulsive_radius_dif+cint_concyl_radius_dif)
		do i=1,ic%n
			if (run(ic%irun(i))%cyl(ic%icyl(i))%nr_con==0 .and. craw_weight_0==0.0) cycle
			if (run(ic%irun(i))%cyl(ic%icyl(i))%nr_con==1 .and. craw_weight_1==0.0) cycle
			if (present(rad_prior) .and. rad_prior) then
				! radius is limited with input frad
				if (run(ic%irun(i))%cyl(ic%icyl(i))%r<frad/teg .or. run(ic%irun(i))%cyl(ic%icyl(i))%r>frad*teg) cycle
			end if
			mask(i)=.true.
			select case (run(ic%irun(i))%cyl(ic%icyl(i))%nr_con)
				case (0); w0=craw_weight_0
				case (1); w0=craw_weight_1
				case (2:); w0=craw_weight_2
			end select
			w(i)=exp(-run(ic%irun(i))%cyl(ic%icyl(i))%potdata)*w0
			wsum=wsum+w(i)
			rad=rad + run(ic%irun(i))%cyl(ic%icyl(i))%r * w(i)
		end do
		if (present(frad)) frad=rad/wsum
		if (wsum==0.0) then
			dori=0.0
			dstrength=0.0
			if (present(frad)) frad=0.0
			return
		end if
		!
		den=w ! cosi with self
		do i=1,ic%n-1
			if (.not.mask(i)) cycle
			do k=i+1,ic%n
				if (.not.mask(k)) cycle
		        cosi=dot_product(run(ic%irun(i))%cyl(ic%icyl(i))%eu(1:3),run(ic%irun(k))%cyl(ic%icyl(k))%eu(1:3))
				den(i)=den(i) + w(k)*(cosi**2)
				den(k)=den(k) + w(i)*(cosi**2)
			end do
		end do
		den=den/wsum
		idum=maxloc(den,mask=mask)
		dstrength=den(idum(1))
		dori=run(ic%irun(idum(1)))%cyl(ic%icyl(idum(1)))%eu
		!
		deallocate(mask,den,w)
	end subroutine get_best_orientation
	!
	subroutine get_vmap_and_ori_for_point(p,vmap,wmap,dori,dstrength,only_vmap)
		implicit none
		real(rk),dimension(1:3),intent(in):: p
		real(rk),dimension(1:3),intent(out):: dori
		real(rk),intent(out):: vmap,wmap,dstrength
		type(iifast_type),dimension(:,:,:),allocatable:: ic
		real(rk),dimension(:,:,:),allocatable:: vm,wm
		logical,optional:: only_vmap
		call calc_visit_and_orientation_cubes(p,0,0.0_rk,vm,wm,ic)
		vmap=vm(1,1,1); wmap=wm(1,1,1)
		if (present(only_vmap) .and. only_vmap) then
			dori=0.0; dstrength=0.0
			return
		end if
		call get_best_orientation(ic(1,1,1),dori,dstrength)
	end subroutine get_vmap_and_ori_for_point
	!
	! calculate vistmaps in a plane, perpendicular to eu vector
	subroutine calc_visit_and_orientation_in_a_plane(p,dd,nrad,rmax,vmap,wmap,ic)
		implicit none
		real(rk),dimension(1:3),intent(in):: p,dd ! input point and orientation
		integer,intent(in):: nrad ! number of point from centre: 0-only in central point
		real(rk),intent(in):: rmax ! maximum radius
		real(rk),dimension(:,:),allocatable,intent(out):: vmap ! pure visitmap
		real(rk),dimension(:,:),allocatable,intent(out):: wmap ! weighted visitmap for analysis
		type(iifast_type),dimension(:,:),allocatable,intent(out):: ic ! array, storing cylinder indexes for each cell
		type(iifast_type):: ivec ! vector of cylinder indexes
		integer:: n,i,k1,k2,irun,icyl,run0
		real(rk):: rsearch,step,radmax,vdum,wdum
		real(rk),dimension(:,:,:),allocatable:: crd
		real(rk),dimension(:,:),allocatable:: vmap0,wmap0
		integer,dimension(:,:),allocatable:: cnt
		logical:: inside
		!
		if (allocated(vmap)) deallocate(vmap)
		if (allocated(wmap)) deallocate(wmap)
		if (allocated(ic)) deallocate(ic)
		rsearch=csearch_min_rad+rmax
		call get_nearest_cyls(p,ivec,rsearch)
		! set the grid and allocate memory
		n=1+2*nrad
		allocate(vmap(1:n,1:n),wmap(1:n,1:n),ic(1:n,1:n))
		allocate(vmap0(1:n,1:n),wmap0(1:n,1:n),cnt(1:n,1:n))
		vmap=0.0_rk; wmap=0.0_rk
		ic%n=0
		if (ivec%n==0) return
		do k2=1,n
			do k1=1,n
				allocate(ic(k1,k2)%irun(1:ivec%n),ic(k1,k2)%icyl(1:ivec%n))
			end do
		end do
		!
		if (nrad==0) then
			step=0.0_rk
			if (rmax>0.0) stop "Err: unrealistic parameters....!"
		else
			step=rmax/real(nrad,kind=rk)
		end if
		!
		! prepare crd coordinates... crd(1:3,1:n,1:n)
		call get_perpendicular_plane_coordinates(p,dd,nrad,rmax,crd)
		!
		run0=0
		do i=1,ivec%n ! cycle over cylinders
			irun=ivec%irun(i)
			icyl=ivec%icyl(i)
			!
			if (irun/=run0) then
				cnt=0; vmap0=0.0; wmap0=0.0
				run0=irun
			end if
			! cycle over gridpoints
			do k2=1,n
				do k1=1,n
					call get_visitmap_value_from_cyl(run(irun)%cyl(icyl),crd(1:3,k1,k2),vdum,wdum,inside)
					if (inside) then
						ic(k1,k2)%n =  ic(k1,k2)%n+1
						ic(k1,k2)%irun(ic(k1,k2)%n) = irun
						ic(k1,k2)%icyl(ic(k1,k2)%n) = icyl
						!
						vmap0(k1,k2)=vmap0(k1,k2) + vdum
						wmap0(k1,k2)=wmap0(k1,k2) + wdum
						cnt(k1,k2)=cnt(k1,k2)+1
					end if
				end do
			end do
			! every realisation have same weight
			! ... average overlapping cylinders
			if (i==ivec%n .or. ivec%irun(i+1)>run0) then
				where(cnt>0)
					vmap=vmap + vmap0/cnt
					wmap=wmap + wmap0/cnt
				end where
			end if
		end do
		deallocate(crd,vmap0,wmap0,cnt)
	end subroutine calc_visit_and_orientation_in_a_plane
	subroutine get_perpendicular_plane_coordinates(p,dd,nrad,rmax,crd)
		implicit none
		real(rk),dimension(1:3),intent(in):: p,dd
		integer,intent(in):: nrad
		real(rk),intent(in):: rmax
		real(rk),dimension(:,:,:),allocatable,intent(out):: crd
		real(rk),dimension(1:3):: ex,ey
		integer:: n,k1,k2
		real(rk):: step,rx,ry
		!
		if (allocated(crd)) deallocate(crd)
		n=1+2*nrad
		allocate(crd(1:3,1:n,1:n))
		if (nrad==0) then
			step=0.0_rk
			if (rmax>0.0) stop "Err: unrealistic parameters....!"
		else
			step=rmax/real(nrad,kind=rk)
		end if
		call get_perpendicular_vectors(dd,ex,ey)
		do k2=1,n
			ry=(k2-1)*step-rmax
			do k1=1,n
				rx=(k1-1)*step-rmax
				crd(1:3,k1,k2)=rx*ex+ry*ey + p
			end do
		end do
	end subroutine get_perpendicular_plane_coordinates
	!
	! calculate visitmap in a given location
	! store the cylinder indexes
	! corners are inaccurate... only spherical volume rmax is correct
	subroutine calc_visit_and_orientation_cubes(p,nrad,rmax,vmap,wmap,ic)
		implicit none
		real(rk),dimension(1:3),intent(in):: p ! input point
		integer,intent(in):: nrad ! number of point from centre: 0-only in central point
		real(rk),intent(in):: rmax ! maximum radius
		real(rk),dimension(:,:,:),allocatable,intent(out):: vmap ! pure visitmap
		real(rk),dimension(:,:,:),allocatable,intent(out):: wmap ! weighted visitmap for analysis
		type(iifast_type),dimension(:,:,:),allocatable,intent(out):: ic ! array, storing cylinder indexes for each cell
		type(iifast_type):: ivec ! vector of cylinder indexes
		integer:: n,i,k1,k2,k3,irun,icyl,run0
		real(rk):: rsearch,step,radmax,vdum,wdum
		real(rk),dimension(1:3):: gp
		integer,dimension(1:3):: imin,imax
		real(rk),dimension(:,:),allocatable:: crd
		real(rk),dimension(:,:,:),allocatable:: vmap0,wmap0
		integer,dimension(:,:,:),allocatable:: cnt
		logical:: inside
		!
		if (allocated(vmap)) deallocate(vmap)
		if (allocated(wmap)) deallocate(wmap)
		if (allocated(ic)) deallocate(ic)
		rsearch=csearch_min_rad+rmax
		call get_nearest_cyls(p,ivec,rsearch)
		! set the grid and allocate memory
		n=1+2*nrad
		allocate(vmap(1:n,1:n,1:n),wmap(1:n,1:n,1:n),ic(1:n,1:n,1:n))
		allocate(vmap0(1:n,1:n,1:n),wmap0(1:n,1:n,1:n),cnt(1:n,1:n,1:n))
		vmap=0.0_rk; wmap=0.0_rk
		ic%n=0
		if (ivec%n==0) return
		do k3=1,n
			do k2=1,n
				do k1=1,n
					allocate(ic(k1,k2,k3)%irun(1:ivec%n),ic(k1,k2,k3)%icyl(1:ivec%n))
				end do
			end do
		end do
		!
		!
		if (nrad==0) then
			step=0.0_rk
			imin=1; imax=1
			if (rmax>0.0) stop "Err: unrealistic parameters....!"
		else
			step=rmax/real(nrad,kind=rk)
		end if
		!
		allocate(crd(1:3,1:n))
		do i=1,n
			crd(1,i)=p(1)+step*(i-nrad-1)
			crd(2,i)=p(2)+step*(i-nrad-1)
			crd(3,i)=p(3)+step*(i-nrad-1)
		end do
		!
		run0=0
		do i=1,ivec%n ! cycle over cylinders
			irun=ivec%irun(i)
			icyl=ivec%icyl(i)
			radmax=sqrt( ( run(irun)%cyl(icyl)%r )**2 + (0.5*run(irun)%cyl(icyl)%h+run(irun)%cyl(icyl)%dh)**2 )
			if (nrad>0) then
				imin=max(1, ceiling((run(irun)%cyl(icyl)%p-radmax - p+rmax)/step)+1 )
				imax=min(n, ceiling((run(irun)%cyl(icyl)%p+radmax - p+rmax)/step) )
			end if
			if (any(imin>n) .or. any(imax<1)) cycle ! cylinder is outside of box
			!
			if (irun/=run0) then
				cnt=0; vmap0=0.0; wmap0=0.0
				run0=irun
			end if
			!
			! cycle over gridpoints
			do k3=imin(3),imax(3)
				gp(3)=crd(3,k3)
				do k2=imin(2),imax(2)
					gp(2)=crd(2,k2)
					do k1=imin(1),imax(1)
						gp(1)=crd(1,k1)
						! gp(1:3) are gridpoint coordinates
						call get_visitmap_value_from_cyl(run(irun)%cyl(icyl),gp,vdum,wdum,inside)
						if (inside) then
							ic(k1,k2,k3)%n = ic(k1,k2,k3)%n+1
							ic(k1,k2,k3)%irun(ic(k1,k2,k3)%n) = irun
							ic(k1,k2,k3)%icyl(ic(k1,k2,k3)%n) = icyl
							!
							vmap0(k1,k2,k3)=vmap0(k1,k2,k3) + vdum
							wmap0(k1,k2,k3)=wmap0(k1,k2,k3) + wdum
							cnt(k1,k2,k3)=cnt(k1,k2,k3)+1
						end if
					end do
				end do
			end do
			! every realisation have same weight
			! ... average overlapping cylinders
			if (i==ivec%n .or. ivec%irun(i+1)>run0) then
				where(cnt>0)
					vmap=vmap + vmap0/cnt
					wmap=wmap + wmap0/cnt
				end where
			end if
		end do
		deallocate(crd,vmap0,wmap0,cnt)
	end subroutine calc_visit_and_orientation_cubes
	subroutine get_vmap_cubes_coordinates(p,nrad,rmax,crd)
		implicit none
		real(rk),dimension(1:3),intent(in):: p ! input point
		integer,intent(in):: nrad ! number of point from centre: 0-only in central point
		real(rk),intent(in):: rmax ! maximum radius
		real(rk),dimension(:,:,:,:),allocatable,intent(out):: crd ! (1:3, :,:,:)
		integer:: n,k1,k2,k3
		real(rk):: step
		!
		if (allocated(crd)) deallocate(crd)
		n=1+2*nrad
		allocate(crd(1:3,1:n,1:n,1:n))
		if (nrad==0) then
			step=0.0_rk
			if (rmax>0.0) stop "Err: unrealistic parameters....!"
		else
			step=rmax/real(nrad,kind=rk)
		end if
		do k3=1,n
			do k2=1,n
				do k1=1,n
					crd(1,k1,k2,k3)=p(1)+step*(k1-nrad-1)
					crd(2,k1,k2,k3)=p(2)+step*(k2-nrad-1)
					crd(3,k1,k2,k3)=p(3)+step*(k3-nrad-1)
				end do
			end do
		end do
	end subroutine get_vmap_cubes_coordinates
end module post_processing