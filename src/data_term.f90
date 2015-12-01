!============
! Author: Elmo Tempel
! Date (last changed): 08.04.2015
!============
!
module data_term
    use constants
    use cylinders
    use galaxies
	use parameters
    implicit none
	! internal parameters to speed up data term hypothesis testing
	integer,private,parameter:: cmem_uniform_ntot=500 ! it cannot be larger than 1000 (NaN will appear)
	real(rk),dimension(:,:),allocatable,private:: cmem_uniform_pval
	logical,private:: cmem_uniform_init=.true.
	integer,private,parameter:: cmem_locdens_ntot=300 ! it cannot be larger than 1000 (NaN will appear)
	real(rk),dimension(:,:),allocatable,private:: cmem_locdens_pval
	logical,private:: cmem_locdens_init=.true.
	integer,private:: nndum1=0,nndum2=0
    !
    ! Probability density is exp(-U) = exp(-[U_data+U_int])
contains
	!
    function calc_data_potential_for_cylinder(cyl,ipts,ngal,loctest) result(res)
        implicit none
        type(cylinder),intent(inout):: cyl ! input cylinder
        real(rk),dimension(:,:),intent(in):: ipts ! input points (it was faster to pass along points than indexes)
        logical,intent(in),optional:: loctest ! test only location (can neglect some calculations)
		integer,intent(in):: ngal ! number of point
        real(rk):: res,rad
        integer:: i,ii,iii(1:1)
        integer:: n1,n2,n3,ns1,ns2,ns3 ! nr of galaxies in cylinder and shadow sections
        real(rk),dimension(:),allocatable:: dist_axis,dist_cnt,dpot1,dpot2,dens,dgrad
        integer:: nn,ns
        real(rk):: dum,dum2,r2,r2s,dummin,step
        ! initialise for bad cylinder
        res=cdata_pot_undefined
        cyl%potdata=res; cyl%potarr=cdata_pot_undefined
		! set constant dataterm for testing
		if (ctest_constant_dataterm) then
			res=ctest_constant_dataterm_value
			cyl%potdata=res
			return
		end if
        !
        allocate(dist_axis(1:ngal),dist_cnt(1:ngal), dens(1:cdata_nr_of_rad_samples+2),dgrad(1:cdata_nr_of_rad_samples+2))
		allocate(dpot1(1:cdata_nr_of_rad_samples+2),dpot2(1:cdata_nr_of_rad_samples+2))
		dpot1=cdata_pot_undefined; dpot2=cdata_pot_undefined
		dens=0.0; dgrad=0.0
        do i=1,ngal
			! dist_axis is squared: dist_cnt is normal
            call get_point_dist_from_cyl2(cyl,ipts(1:3,i),dist_axis(i),dist_cnt(i))
        end do
		!
		step=(cyl%rmax-cyl%rmin)/real(cdata_nr_of_rad_samples-1,kind=rk)
		do ii=1,cdata_nr_of_rad_samples+2
			if (cdata_fixed_radius) then
				rad=cdata_cyl_rad_min
			else
				rad=cyl%rmin-step+(ii-1)*(cyl%rmax-cyl%rmin+2*step)/real(cdata_nr_of_rad_samples+1)
			end if
	        !---- apply the rules for p-value (initial selection)
	        ! minimum number of galaxies inside a cylinder
	        r2=rad*rad
	        nn=count(dist_axis<r2 .and. abs(dist_cnt)<0.5*cyl%h )
	        if (nn<cdata_minpts) then
	            cycle
	        end if
	        !
	        n2=count(dist_axis<r2 .and. abs(dist_cnt)<=cyl%h/6.0)
	        n1=count(dist_axis<r2 .and. abs(dist_cnt)<0.5*cyl%h .and. dist_cnt>cyl%h/6.0)
	        n3=count(dist_axis<r2 .and. abs(dist_cnt)<0.5*cyl%h .and. dist_cnt<-cyl%h/6.0)
			!
	        dum2=get_uniform_pvalue(n1,n2,n3)
	        if (dum2<cdata_p_uniform_cut) then
	            cycle
	        end if
			!
	        ! find shadow galaxies
	        r2s=r2*cdata_shadow_dvol
	        ns1=count(dist_axis<r2s .and. abs(dist_cnt)<0.5*cyl%h .and. dist_axis>=r2 .and. dist_cnt>cyl%h/6.0)
	        ns2=count(dist_axis<r2s .and. abs(dist_cnt)<=cyl%h/6.0 .and. dist_axis>=r2)
	        ns3=count(dist_axis<r2s .and. abs(dist_cnt)<0.5*cyl%h .and. dist_axis>=r2 .and. dist_cnt<-cyl%h/6.0)
	        ns=ns1+ns2+ns3
	        !
	        ! calc locally high density and uniform hypotheis testing and apply the probability cuts
	        dum=get_locdens_pvalue(nn,ns)
	        if (dum<cdata_p_dens_cut) then
	            cycle
	        end if
			!
	        ! all tests are successful... calculate the data term
			!
			! loctest asjad...
	        if (present(loctest).and.loctest) then
	            res=0.0
	            deallocate(dist_axis,dist_cnt,dpot1,dpot2,dens,dgrad)
	            return
	        end if
	        !
	        !---- all tests are passed ---
			dens(ii)=real(nn)/( r2 )
			dpot1(ii)=-log( dum ) ! hypothesis testing potential: locally high density
			dpot2(ii)=-log( dum2 ) ! hypothesis testing potential: uniform along cylinder
			!dgrad(ii)=dum ! locally high density (it does not work very well...)
			if (cdata_fixed_radius) exit
		end do ! cycle over rad samples
		!
		if (cdata_fixed_radius) then
			ii=1
		else
			if (.true.) then ! tiheduse gradient... if false, then locally high density
				! find the density gradient
				dgrad=0.0
				do ii=2,cdata_nr_of_rad_samples+1
					if (dens(ii)==0.0) cycle
					if (dens(ii+1)==0.0) then
						dgrad(ii)=(dens(ii-1)-dens(ii))*2
					else
						dgrad(ii)=dens(ii-1)-dens(ii+1)
					end if
				end do
			end if
			! find the best radius and set that potential
			iii=maxloc(dgrad); ii=iii(1)
		end if
		!
		if (dpot1(ii)+dpot2(ii)<cdata_pot_undefined) then
			if (cdata_fixed_radius) then
				rad=cdata_cyl_rad_min
			else
				rad=cyl%rmin-step+(ii-1)*(cyl%rmax-cyl%rmin+2*step)/real(cdata_nr_of_rad_samples+1)
			end if
			cyl%r=rad
			r2=rad*rad
	        ! calculate the cylinder variance
	        res=0.0; nn=0
	        do i=1,ngal
	            if (dist_axis(i)<r2 .and. abs(dist_cnt(i))<0.5*cyl%h) then
	                res=res+dist_axis(i)
					nn=nn+1
	            end if
	        end do
			!
			cyl%ngal=nn
			cyl%pot_hyp=dpot1(ii)+dpot2(ii)
			cyl%pot_con=res/( (nn-2)*r2 )
			! store potential values to testing array
			cyl%potarr(1)=dpot1(ii)
			cyl%potarr(2)=dpot2(ii)
			!
			res=cdata_coeff_variance*cyl%pot_con + cdata_coeff_hypothesis*cyl%pot_hyp
			cyl%potdata=res
		else
	        res=cdata_pot_undefined
	        cyl%potdata=res
			cyl%r=0.0
			cyl%ngal=0
		end if
        !
        deallocate(dist_axis,dist_cnt,dpot1,dpot2,dens,dgrad)
    end function calc_data_potential_for_cylinder
    !
	!
	!--- testime hypothesis varki
	subroutine test_hypothesis()
		implicit none
		integer:: n1,n2,nmax,ntot,n3
		cdata_uniform_dens=3.0
		cdata_locdens_dens=1.0
		do
			print*, "sisesta n1,n2: "
			read*, n1,n2
			ntot=n1+n2
			!
			if (.true.) then
				nmax=max(n1,n2)
				n3=(n1+n2)/2
				cdata_uniform_dens=3.0
				print*, get_uniform_pvalue_internal(ntot,nmax)
				cdata_uniform_dens=1.0/cdata_uniform_dens
				print*, get_uniform_pvalue_internal2(ntot,nmax)
				cdata_uniform_dens=1.0/cdata_uniform_dens
				print*, get_uniform_pvalue(n1,n2,n3)/6
				!
				print*
				nmax=min(n1,n2)
				n3=(n1+n2)/2
				cdata_uniform_dens=3.0
				print*, get_uniform_pvalue_internal(ntot,nmax)
				cdata_uniform_dens=1.0/cdata_uniform_dens
				print*, get_uniform_pvalue_internal2(ntot,nmax)
				cdata_uniform_dens=1.0/cdata_uniform_dens
				print*, get_uniform_pvalue(n1,n2,n3)/6
			end if
			!
			if (.false.) then
				print*
				print*, get_locdens_pvalue_internal(ntot,n1)
				print*, get_locdens_pvalue_internal2(ntot,n1)
				print*, get_locdens_pvalue(n1,n2)
			end if
		end do
	end subroutine test_hypothesis
    !
	!
	function get_uniform_pvalue(n1,n2,n3) result(pvalue)
		implicit none
		integer,intent(in):: n1,n2,n3
		real(rk):: pvalue,dum
		integer:: nmin,nmax,ntot
		! find the min and max values and correct if necessary
        nmin=min(n1,n2,n3)
        nmax=max(n1,n2,n3)
		ntot=nmin+nmax
		if (ntot>cmem_uniform_ntot) then
			dum=real(ntot)/real(cmem_uniform_ntot)
			nmax=floor(nmax/dum)
			nmin=floor(nmin/dum)
			ntot=nmin+nmax
		end if
		!
		! calculate using incomplete beta function
		!dum=1.0-betainc(real(nmin+1,kind=rk),real(nmax,kind=rk),1.0_rk-cdata_uniform_dens/(1.0_rk+cdata_uniform_dens))
		pvalue=get_uniform_pvalue_memory(ntot,nmax)
	end function get_uniform_pvalue
	function get_uniform_pvalue_real(n1,n2,n3) result(pvalue)
		implicit none
		real(rk),intent(in):: n1,n2,n3
		real(rk):: pvalue
		real(rk):: nmin,nmax
		! find the min and max values and correct if necessary
        nmin=min(n1,n2,n3)
        nmax=max(n1,n2,n3)
		! calculate using incomplete beta function
		pvalue=1.0-betainc(nmin+1.0_rk,nmax,1.0_rk-cdata_uniform_dens/(1.0_rk+cdata_uniform_dens))
	end function get_uniform_pvalue_real
	function get_uniform_pvalue_memory(ntot,nmax) result(pval)
		implicit none
		integer,intent(in):: ntot,nmax
		real(rk):: pval
		integer:: i,k
		if (cmem_uniform_init) then
			allocate(cmem_uniform_pval(1:cmem_uniform_ntot,1:cmem_uniform_ntot))
			cmem_uniform_pval=0.0
			do i=1,cmem_uniform_ntot
				do k=max(1,floor(i*0.5)),i
					cmem_uniform_pval(i,k)=get_uniform_pvalue_internal(i,k)
				end do
			end do
			cmem_uniform_init=.false.
		end if
		if (2*nmax<ntot) stop "Error: uniform pvalue from memory: incorrect numbers!"
		if (ntot>cmem_uniform_ntot) stop "Error: uniform pvalue mem not big enough!"
		pval=cmem_uniform_pval(ntot,nmax)
	end function get_uniform_pvalue_memory
	function get_uniform_pvalue_internal(ntot,nmax) result(pval)
		implicit none
		integer,intent(in):: ntot,nmax
		real(rk):: pval,dum,p_uniform
		integer:: i,j
		p_uniform=cdata_uniform_dens/(1.0_rk+cdata_uniform_dens)
		pval=0.0
        do i=nmax,ntot
			dum=1.0_rk
			! calculate binomial factor
			do j=1,i
				dum=dum*(1.0_rk+real(ntot-i,kind=rk)/real(j,kind=rk))
			end do
            dum=dum*(p_uniform**i)*((1.0_rk-p_uniform)**(ntot-i))
            pval=pval+dum
        end do
		if (isnan(pval)) stop "Error: get uniform pvalue NaN!!!"
	end function get_uniform_pvalue_internal
	function get_uniform_pvalue_internal2(ntot,nmax) result(pval)
		implicit none
		integer,intent(in):: ntot,nmax
		real(rk):: pval,dum,p_uniform
		integer:: i,j
		p_uniform=cdata_uniform_dens/(1.0_rk+cdata_uniform_dens)
		pval=0.0
        do i=0,nmax!nmax,ntot
			dum=1.0_rk
			! calculate binomial factor
			do j=1,i
				dum=dum*(1.0_rk+real(ntot-i,kind=rk)/real(j,kind=rk))
			end do
            dum=dum*(p_uniform**i)*((1.0_rk-p_uniform)**(ntot-i))
            pval=pval+dum
        end do
		if (isnan(pval)) stop "Error: get uniform pvalue NaN!!!"
	end function get_uniform_pvalue_internal2
	!
	! get local density pvalue
	function get_locdens_pvalue(ncyl,nshadow) result(pvalue)
		implicit none
		integer,intent(in):: ncyl,nshadow
		real(rk):: pvalue,dum
		integer:: ncyl2,ntot
		! correct numbers if necessary
		ntot=ncyl+nshadow
		ncyl2=ncyl
		if (ntot>cmem_locdens_ntot) then
			dum=real(ntot)/real(cmem_locdens_ntot)
			ncyl2=floor(ncyl/dum)
			ntot=floor(nshadow/dum)
			ntot=ncyl2+ntot
		end if
		!
		! calculate using incomplete beta function - for real values
		!dum=betainc(real(nshadow,kind=rk),real(ncyl+1,kind=rk),1.0_rk-cdata_locdens_dens/(1.0_rk+cdata_locdens_dens))
		pvalue=get_locdens_pvalue_memory(ntot,ncyl2)
	end function get_locdens_pvalue
	function get_locdens_pvalue_real(ncyl,nshadow) result(pvalue)
		implicit none
		real(rk),intent(in):: ncyl,nshadow
		real(rk):: pvalue
		! calculate using incomplete beta function - for real values
		pvalue=betainc(nshadow,ncyl+1.0_rk,1.0_rk-cdata_locdens_dens/(1.0_rk+cdata_locdens_dens))
	end function get_locdens_pvalue_real
	function get_locdens_pvalue_memory(ntot,ncyl) result(pval)
		implicit none
		integer,intent(in):: ntot,ncyl
		real(rk):: pval
		integer:: i,k
		if (cmem_locdens_init) then
			allocate(cmem_locdens_pval(1:cmem_locdens_ntot,1:cmem_locdens_ntot))
			cmem_locdens_pval=0.0
			do i=1,cmem_locdens_ntot
				do k=1,cmem_locdens_ntot
					cmem_locdens_pval(i,k)=get_locdens_pvalue_internal(i,k)
				end do
			end do
			cmem_locdens_init=.false.
		end if
		if (ntot>cmem_locdens_ntot) stop "Error: locdens pvalue mem not big enough!"
		pval=cmem_locdens_pval(ntot,ncyl)
	end function get_locdens_pvalue_memory
	function get_locdens_pvalue_internal(ntot,ncyl) result(pval)
		implicit none
		integer,intent(in):: ntot,ncyl
		real(rk):: pval,dum,p_uniform
		integer:: i,j
		p_uniform=cdata_locdens_dens/(1.0_rk+cdata_locdens_dens)
		pval=0.0
		! next do is effectively 1,ncyl
        do i=ncyl+1,ntot ! this +1 is necessary... to be compatible with reversed probability
			dum=1.0_rk
			! calculate binomial factor
			do j=1,i
				dum=dum*(1.0_rk+real(ntot-i,kind=rk)/real(j,kind=rk))
			end do
            dum=dum*(p_uniform**i)*((1.0_rk-p_uniform)**(ntot-i))
            pval=pval+dum
        end do
		if (isnan(pval)) stop "Error: get locdens pvalue NaN!!!"
		pval=1.0_rk-pval
	end function get_locdens_pvalue_internal
	function get_locdens_pvalue_internal2(ntot,ncyl) result(pval)
		implicit none
		integer,intent(in):: ntot,ncyl
		real(rk):: pval,dum,p_uniform
		integer:: i,j
		p_uniform=cdata_locdens_dens/(1.0_rk+cdata_locdens_dens)
		pval=0.0
		! next do is effectively 1,ncyl
        do i=0,ncyl!ncyl+1,ntot ! this +1 is necessary... to be compatible with reversed probability
			dum=1.0_rk
			! calculate binomial factor
			do j=1,i
				dum=dum*(1.0_rk+real(ntot-i,kind=rk)/real(j,kind=rk))
			end do
            dum=dum*(p_uniform**i)*((1.0_rk-p_uniform)**(ntot-i))
            pval=pval+dum
        end do
		if (isnan(pval)) stop "Error: get locdens pvalue NaN!!!"
		pval=pval
	end function get_locdens_pvalue_internal2
    !
    ! calc factorial - not needed anymore
    recursive function factorial(p) result(res)
        implicit none
        integer,intent(in):: p
        real(rk):: res
        if (p==1) then
            res=1.0_rk
        else if (p==0) then
            res=1.0_rk
        else
            res=p*factorial(p-1)
        end if
    end function factorial
end module data_term