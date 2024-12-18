!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!this fortran codes perform the 3d structural order analysis for a binary LJ mixture (bljm) and 
!silica-based network system (nsx) as documented in the paper of Zhang and Kob (see ref at the end). 
!The configuration (in xyz format) to be analyzed is by default 3d-periodic, but non-periodic 
!configs can easily be adapted to this analysis. 
!
!the analysis are basically two folds: 
!
!0) local coordinate system: 
! this lies at the core of the proposed four-point correlation. the general rule 
! is that one takes a relatively stable local structure motif, formed by a triplet of atoms,
! as reference. as such, one first needs some basical knowledge of the structure of 
! the system under investigation, such as g(r)
! specifically, for hard-shere-like systems, one can use the triplet of particles that are 
! nearest neighbor to each other to construct the reference coordinate system; 
! for open-network systems with a well-defined local structure, e.g. SiO4 in sio2-based glasses,
! one can use the O-Si-O linkage to construct the local reference frame
! in any case, one should play a bit around with the choices of local reference frames
! to see their performance
!
!1) spherical density distribution (aka sff2 in the code): 
! project the spatial distribution of the particles on a sphere with respect to
! the chosen local coordinate system
!
!2) quantitative characterization of the 3d structure (aka sff1 in the code):
! decompose the signal as seen from the density distribution by spherical harmonics (SH), 
! for which different modes (l) catch different types of symmetries and order
! it is more accurate to calculate the SH for each of the point projected on the sphere 
! rather than using the density distribution obtained from 1)
! note that the larger the mode (l), the more time consuming.
!
!Reference: Zhang, Z., & Kob, W. (2020). Revealing the three-dimensional structure of liquids
!using four-point correlation functions. PNAS, 117(25), 14032-14037.
!
!usage of the code:
!gfortran -o nsx.out 3d-structure-analysis-nsx.f90 && ./nsx.out sff1 22 2.5 5.5 0.05 0.5 2000 3 0.1
!see below for the definition of each of the loaded parameters
!
!by Zhen Zhang, zhen.zhang1991@hotmail.com, Jun 6, 2022
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program orderparameter 
!~    USE MPI
   implicit none
      
   !predefiend parameters
   real(kind=8), parameter :: min1_gr=2.0  !the first min of gr: 2.0 for SiO2, 1.4 for BLJM
   real(kind=8), parameter :: pi=4*atan(1.0) 
   integer,parameter :: n_min=100000, nmax=30  !max nb of atoms allowed in each grid_cell
   integer,parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
   
   !variables
   character (len=200) :: filename4,filename0
   character (len=10) :: seqstring,ana 
   character (len=128) :: label, correlation
   
   real(kind=8) :: lbox(3), cellsize, distance,P_nn  
   real(kind=8) :: S_ff,r_sp,phi_sp,theta_sp,n_real,n_imag!,n_real_test,n_imag_test!,whatever
   real(kind=8) :: dr, rmax, rmin, d_area,tmp2(3),delta0(3),vect_z(3),&
                   & vect_y(3),vect_x(3),vect_tmp(3),vect_N(n_min,3), vect_tar(n_min,3), minr, maxr
   real(kind=8) :: rad,delta,theta_del, frac !ang_min_1st,
   real(kind=8) :: goldenratio,theta_fib,phi_fib,coord_fib(3),s_cell
           
   integer :: l_value !for y_lm, the modes of SH to be analyzed 
   integer :: nline,n1,n2,i,ierror,n,j,num1,num2  
   integer :: n_atom,num_total,num_max,n_area,k,ntotal_atom !num_1,    
   integer :: m1,n_3d, num_totalPC, num_maxpc,num_minpc,ntotal_1cfg!,num_min !,num_totalsp
   integer :: real_totalbond,l,m_l,n_shell, nshell, j0, nn
   integer :: ntype(3),ncells(3), d(3), d0(3), d1(3),n_neighbor,i0,k0, nn1,nn2,numN,&
              & corr, nshells, n_fib !=4610 !nb of points be to evently spaced on the sphere  
   integer :: num_args, ix, i_fib, nfrac !number of a fraction of atoms to be analyzed   
   
   complex(wp) :: ylm_numeric
   complex(wp), external :: ylm 
   
   !allocable arrays
   integer, dimension(:), allocatable :: ntotal
   integer, dimension(:), allocatable ::typat,npoints
   integer, dimension(:,:,:), allocatable :: nb
   real(kind=8),dimension(:,:),allocatable :: wrapped,vect_neighbor
   real(kind=8),dimension(:,:,:,:,:),allocatable :: codcell
   real(kind=8), dimension (:,:), allocatable :: integral_imag,integral_real,area_coordinate 
   real(kind=8), dimension (:), allocatable :: s_f
   character(len=200), dimension(:), allocatable :: args

   type :: atom    !define type variable for identifying and saving atom information
      character (len=5) :: label
      real(kind=8) :: x,y,z
   end type atom  
   
   !-----command_argument to be loaded while executing the code-----
   num_args = command_argument_count()
   allocate(args(num_args))
   do ix = 1, num_args
      call get_command_argument(ix,args(ix))
      print*, trim(args(ix))
   end do

   print*,"COMMAND_ARGUMENT_COUNT",command_argument_count()

   ana=args(1) !sff1 or sff2
   correlation=args(2) !atyp_ref-atyp_target
   read(args(3),'(f5.2)')minr !min value of r to be analyzed
   read(args(4),'(f5.2)')maxr !max vale of r 
   read(args(5),'(f5.2)')dr  !0.05, define the resoultion of the S_rho curve (only active for sff1)
   read(args(6),'(f5.2)')delta  !0.5, define the layer [r-delta, r+delta] to be looked in practice
   read(args(7),'(i6)') n_fib  !nb of cells on the sphere for bining
   read(args(8),'(i3)') l_value  !the mode of SH to be calculated: 3 for SiO2, 6 for BLJM 
   read(args(9),'(f6.4)')frac !fraction of the config to be analyzed (in case the config is large)
   
   !-----------------------------------------------------------------
   if(correlation=="21") then  !atyp2_atyp1
      corr=1 
   elseif(correlation=="22") then  !atyp2_atyp2
      corr=2
   elseif(correlation=="23") then  !atyp2_atyp2
      corr=3
   endif
   
   !---------read the config file, in xyz format------------------------
   open(unit=12, file='config-nsx.xyz',status="old",action="read",iostat=ierror)
   print*,"NOTE: by default the origin of the box in at [0,0,0], otherwise translate the box!"

   read(12,*,end=80) nline !nb. of atoms
   read(12,*) lbox(1:3) ! box lengths in x,y,z
   print *,"n_atoms, Lbox:",nline, lbox(1:3) !tmp0,tmp2(:)  
   
   allocate (typat(nline),wrapped(nline,3))
   ntype=0
   do j=1,nline
    read(12,*) typat(j), wrapped(j,1:3)
    ntype(typat(j))=ntype(typat(j))+1  !nsx(O=1, Si=2, Na=3), lj(A=1,B=2)
   enddo
   
   80 continue 
   write(*,*) "n_atom of different typs= ", ntype      
   close(12)   
   print *,"-- readin coordinates done! --"
!-----------------------------------------------------------------

!------------bin atoms into uniformly sized cells-----------------
!--a procedure that makes efficient the analysis of large configs-
  
  cellsize=min1_gr
  ncells(:)=int(lbox(:)/cellsize)+1
  print *,"bin size:",cellsize,"n_cells in x,y,z: ",ncells
  
  allocate (codcell(ncells(1),ncells(2),ncells(3),nmax,5))
  allocate (nb(ncells(1),ncells(2),ncells(3)))
  nb=0  !array saving the index of the atoms in each cell
  
  do j=1,nline
    d(:)=int(wrapped(j,1:3)/cellsize)+1  !cell identity
    nb(d(1),d(2),d(3))=nb(d(1),d(2),d(3))+1  ! index of the atoms in cell[d(1),d(2),d(3)]
    codcell(d(1),d(2),d(3),nb(d(1),d(2),d(3)),1:3)=wrapped(j,1:3)
    codcell(d(1),d(2),d(3),nb(d(1),d(2),d(3)),4)=0!atomid(i,j)
    codcell(d(1),d(2),d(3),nb(d(1),d(2),d(3)),5)=typat(j)!atomtyp(i,j)
  enddo

  if (maxval(nb)>nmax) then
    print *,"Too many atoms in a cell! adjust nmax or cellsize!"; call exit
  endif  
  !call sleep(30)  
  print *,"-- bin atoms into cells done! --"         
!-------------------------------------------------------------  

!-------------------------------------------------------------  
  if (ana=="sff2") then   
  !----get the coordinates of the cells on the sphere-----------
  !this is needed for latter bining of the 3d distributions of the particles 
  !see more: http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
  
  !calculate the theta_del to be used for latter bining of particles in space
  !s_cell=2*pi*(1-cos(theta_del/180*pi)) 
  s_cell=3*4*pi/n_fib  !coarse-grain to make sure 3 times the spherical areas are covered in total
  theta_del=180/pi*acos(1-s_cell/(2*pi))
  print*,"theta_del:",theta_del
  
  n_area=n_fib
  print*,"nb of Fibonacci lattices evenly distributed on sphere (for sff2):", n_fib
  goldenratio=(1 + 5**0.5)/2
  allocate(area_coordinate(n_fib,3),npoints(n_fib))
  npoints=0
  
  !--canonical Fibonacci lattice----
  do i_fib=0, n_fib-1
   theta_fib = 2 *pi * i_fib / goldenratio
   phi_fib = acos(1 - 2*(i_fib+0.5)/n_fib)
   coord_fib(1) = cos(theta_fib) * sin(phi_fib)
   coord_fib(2) = sin(theta_fib) * sin(phi_fib)
   coord_fib(3) = cos(phi_fib)
   area_coordinate(i_fib+1,1:3)=coord_fib(1:3)
  enddo
  endif !ana=sff2
!-----------------------------------------------------------

!-----------------------------------------------------------  
  if (ana=="sff1") then  
   print*,"the mode of spherical harmonics to be calculated: l=", l_value 
!---------create files to write for sff1 analysis-----------
   write (seqstring,'(i0.2)') l_value 
   call system ('mkdir -p structure')  !create a folder to save the outputs
   filename4='Srho-l'//trim(seqstring)//'-nsx-corr'//trim(correlation)//'.dat' !S_rho
   !append to an existing file, otherwise "create and write" will overwrite the existing file.
   open(unit=1000+l_value,file='structure/'//filename4,action='write') !,position='append')  
  endif !ana=sff1
!-----------------------------------------------------------

!-----------------------------------------------------------  
  print*,"-loop over all distances (shells) with dr as increment-"
!the analysis of the strcutue is done in a stepwise manner, i.e., layer by layer
  
  !allocate arrays for the latter SH analysis
  allocate(integral_imag(l_value+1,2*l_value+1),integral_real(l_value+1,2*l_value+1),s_f(l_value+1))
  allocate(vect_neighbor(nmax,3)) !at most nmax nearest neighbors are allow
   
  n_shell=int((maxr-minr)/dr)+1
  shells: do nshell=1,n_shell
  
  rad=minr+(nshell-1)*dr  !define the shell: [r-delta,r+delta]
  rmax=rad+delta
  rmin=rad-delta
  print '(a,2f6.2,a,i5)',"r_range:",rmin,rmax, " n_shell:",nshell

 !------INITIALIZATION--------
  ntotal_atom=0; num_max=0;n_3d=0;num_total=0;num_totalPC=0;real_totalbond=0
  integral_real=0.0; integral_imag=0.0      
      
!------start the algorithm for 3d structural order----
  print*,"-- starts loop over all ref atoms --" !i.e., atom #1 in the local coordinate system
  !in case the config is too large, one may one need a fraction of the atoms to be analyzed
  nfrac=int(frac*nline) !frac -> a fraction of the atoms to be considered
  print*,"nfrac",nfrac
  
  do i=1,nfrac !nline
   if (typat(i)==2) then !type_2 as ref, i.e. particle #1 
   d(:)=int(wrapped(i,1:3)/cellsize)+1 !find the hosting cell for atom_ref
   
   n_neighbor=0 ; vect_neighbor=0  !array to save the list of nearest neighbors
   do i0=d(1)-2,d(1)+2 !search 2 cells out in each direction
   do j0=d(2)-2,d(2)+2
   do k0=d(3)-2,d(3)+2
	d0=(/ i0, j0, k0 /)
	d1(:)=d0(:)-ncells(:)*floor(real(d0(:)-1)/ncells(:))  !PBC   3D
	
	do k=1,nb(d1(1),d1(2),d1(3))
	 if(int(codcell(d1(1),d1(2),d1(3),k,5))==1) then !atype_1
	 delta0(:)= codcell(d1(1),d1(2),d1(3),k,1:3)-wrapped(i,:)
	 delta0(:)= delta0(:)-lbox(:)*anint(delta0(:)/lbox(:))     !PBC in 3D				  
	 if (norm2(delta0)>0.01 .and. norm2(delta0)<=min1_gr) then  !typ2_typ1
		 n_neighbor=n_neighbor+1
		 vect_neighbor(n_neighbor,:)=delta0  !saving the nearest neighbors
	 endif              
	 endif
	enddo
   
   enddo
   enddo
   enddo !O
    
    numN=0 ; vect_N=0  !saving the vectors from the target particle to the ref particle  
    nn1=int(rmax/cellsize)+2  !outer_bound

    do i0=d(1)-nn1,d(1)+nn1 !search cells only within the outer bounds
    do j0=d(2)-nn1,d(2)+nn1
    do k0=d(3)-nn1,d(3)+nn1
	 
	 d0=(/ i0, j0, k0 /)
	 d1(:)=d0(:)-ncells(:)*floor(real(d0(:)-1)/ncells(:))  !PBC 3D
	 
	 do k=1,nb(d1(1),d1(2),d1(3))
	  if(int(codcell(d1(1),d1(2),d1(3),k,5))==corr) then   !atyp_corr as target particles 
	  delta0(:)= codcell(d1(1),d1(2),d1(3),k,1:3)-wrapped(i,:)
	  delta0(:)= delta0(:)-lbox(:)*anint(delta0(:)/lbox(:))  !PBC in 3D
	  distance=norm2(delta0)  
	  if (distance<rmax .and. distance>rmin) then
		 numN=numN+1
		 vect_N(numN,:)=delta0  !saving the target particles
	  endif              
	  endif
	 enddo
 
    enddo
    enddo
    enddo

    if (numN>n_min) then
      print *,"Too many atoms in the shell, change size of the vect_N array!"; call exit
    endif          

    !loop over the nearest neighbors for constrcuting local coordinate system    
    do n1=1,n_neighbor-1   
     do n2=n1+1,n_neighbor
      
      !the 1st bond define unit vect in z
      vect_z(:)=vect_neighbor(n1,:)/(norm2(vect_neighbor(n1,:)))  
      !define unit vector in y
      vect_tmp(1)=vect_z(2)*vect_neighbor(n2,3)-vect_z(3)*vect_neighbor(n2,2)
      vect_tmp(2)=vect_z(3)*vect_neighbor(n2,1)-vect_z(1)*vect_neighbor(n2,3)
      vect_tmp(3)=vect_z(1)*vect_neighbor(n2,2)-vect_z(2)*vect_neighbor(n2,1)
      vect_y(:)=vect_tmp(:)/(norm2(vect_tmp(:)))
      !define unit vector in x
      vect_x(1)=vect_y(2)*vect_z(3)-vect_y(3)*vect_z(2)
      vect_x(2)=vect_y(3)*vect_z(1)-vect_y(1)*vect_z(3)
      vect_x(3)=vect_y(1)*vect_z(2)-vect_y(2)*vect_z(1)
    
    !now rotate all the vectors formed by the ref and target particles 
    !with respect to the local coordinate system
    do j=1,numN
     vect_tar(j,:)=vect_N(j,:)
     call newcoordinate(vect_x,vect_y,vect_z,vect_tar(j,:))  
     ntotal_atom=ntotal_atom+1
    enddo
    
    !-------------------------------------------------------------
    !---calculate spherical harmonics---(sff1)
    if (ana=="sff1") then 
    do j=1,numN
     call cart2sphere(vect_tar(j,:),r_sp,phi_sp,theta_sp)   !rad rather than degree
      l=l_value
      do m_l=-l,l
     
      ylm_numeric = ylm (l, m_l, phi_sp, theta_sp)
	  !the realpart of ylm with precision defined for the complex !realpart(ylm_numeric)
	  n_real= real(ylm_numeric)  !realpart
	  !the imagpart of ylm with precision defined for the compleximagpart(ylm_numeric)
	  n_imag= aimag(ylm_numeric) !imagpart
	  integral_real(l+1,m_l+l_value+1)=integral_real(l+1,m_l+l_value+1)+n_real ! l+1,m+l+1 means l=l & m=m
	  integral_imag(l+1,m_l+l_value+1)=integral_imag(l+1,m_l+l_value+1)-n_imag !
    
      enddo
     enddo 
    endif !ana=sff1
    !---------------------------------------------------------
    
    !---------------------------------------------------------
    !---for para: bin atoms into cells--(sff2)
    if (ana=="sff2") then
    do j=1,numN
     do k=1,n_area
       call angleof2vector(area_coordinate(k,:),vect_tar(j,:),P_nn)
       
       if (P_nn<theta_del) then    !search atom in the circle
         npoints(k)=npoints(k)+1     
         real_totalbond=real_totalbond+1  !take into account the recounted bonds
       end if 
       
     enddo
    enddo
    endif !ana=sff2
    !---------------------------------------------------------
    
     enddo !n2
    enddo  !n1
   
   endif !ref_atom
  enddo
  
  print *,"nb_atoms on the sphere; actual(sff1), corase-grained(sff2): ",ntotal_atom,real_totalbond
  
  
  print *,"------- output results to files ------------"
  
  !--------sff2----------------------------------------
  if (ana=="sff2") then
  call system ('mkdir -p structure')
  write (seqstring,'(i0.2,f3.2)') int(rad),rad-int(rad)
  filename0='spherical-density-r'//trim(seqstring)//'-nsx-corr'//trim(correlation)//'.xyz'
  open(unit=600,file='structure/'//filename0,status='replace',action='write')
  
  write (600,*) n_area
  write (600,*) "spatial density distribution"
  d_area=4*pi/n_area  !area of each cell
  do i=1, n_area
   !type,x,y,z,rho_normalized
   write (600,'(i2, 3f7.3,f9.5)') 1, area_coordinate(i,:),real(npoints(i))/sum(npoints(:))/d_area 
  enddo
  close(600)
  
  endif !ana=sff2
  !----------------------------------------------------
  
  !-----------sff1------------------------------------- 
  if (ana=="sff1") then
    l=l_value; S_ff=0.0
			
	over_mm: do m_l=-l,l
	S_ff=S_ff+(integral_real(l+1,m_l+l_value+1)/dble(ntotal_atom))**2        
	end do over_mm
		
	S_ff=S_ff/dble(2*l+1)  !angular power spectrum (for Schmidt semi-normalized harmonics)
	s_f(l+1)=S_ff  !array to save S_ff at each l
	
	if (s_f(l+1)>=0) then
	write (1000+l,*) rad, sqrt(s_f(l+1))  !r, s_rho
	endif 
  endif !ana=sff1
  !---------------------------------------------------      
      
  enddo shells  
  
  if (ana=="sff1") then  
   close(1000+l_value)      
  endif !ana=sff1
  
end program orderparameter


!============subroutines and functions====================
subroutine newcoordinate(u,v,w,s)

  implicit none 
  real(kind=8) :: a,b,c
  real(kind=8), dimension (3) :: u,v,w,s
  
  a=s(1)*u(1)+s(2)*u(2)+s(3)*u(3)
  b=s(1)*v(1)+s(2)*v(2)+s(3)*v(3)
  c=s(1)*w(1)+s(2)*w(2)+s(3)*w(3)
  
  s(1)=a
  s(2)=b
  s(3)=c

end subroutine newcoordinate

subroutine angleof2vector (v,u,ang)
   
   implicit none
   real(kind=8) :: u(3),v(3),ang
   real(kind=8), parameter :: pi=4*atan(1.0)
  
  ang=(v(1)*u(1)+v(2)*u(2)+v(3)*u(3))/(norm2(v)*norm2(u)) 
  
  if (ang==1.00000012) then
     ang=1.0
  end if

  ang=acos(ang)/pi*180
  return
end subroutine angleof2vector

  
  subroutine cart2sphere (v,r_sp,phi_sp,theta_sp)             
  !originally from Vicente M. Reyes,!return radians
    implicit none
    
    character*31 tag
    real(kind=8) :: x, y, z, pi, r_sp, phi_sp, theta_sp, S, v(3)
    
    x=v(1); y=v(2); z=v(3)
    pi=4*atan(1.0)
    tag = "SPHER(RHO,PHI,THETA) in radians" 
    
    r_sp = sqrt(x**2 + y**2 + z**2) 
    phi_sp =acos(z/r_sp)  !acos(z/r_sp)*180.0/pi
    
    S = sqrt(x**2 + y**2)
    if ((x.gt.0.00).and.(y.ge.0.00)) then
      theta_sp =(ASIN(y/S)) !((ASIN(y/S))*(180.0/pi))
    elseif ((x.le.0.00).and.(y.gt.0.00)) then
      theta_sp =pi - ASIN(y/S) !((pi - ASIN(y/S))*(180.0/pi))
    elseif ((x.lt.0.00).and.(y.le.0.00)) then
      theta_sp =pi - ASIN(y/S) !((pi - ASIN(y/S))*(180.0/pi))
    elseif ((x.ge.0.00).and.(y.lt.0.00)) then
      theta_sp =2*pi + ASIN(y/S) !((2*pi + ASIN(y/S))*(180.0/pi))
    endif 
    
    !----“spher2cart------------
    !x=rsin(φ)cos(θ)  0 ≤ θ ≤ 2π.
    !y=rsin(φ)sin(θ)  0 ≤ φ ≤ π.
    !z=rcos(φ)
    !---------------------------
    
  end subroutine cart2sphere


  function ylm (l, m, thrad, phirad)
  implicit none
!--------------------------------------------------------------------
! Computes the spherical harmonic Y_lm (theta,phi) using the
! reduced rotation matrix d^l_{m 0} (theta) and using the
! external function fac10(n) = factorial(n)/10**n
!--------------------------------------------------------------------
! input: angular momentum quantum numbers l, m (integers)
!        angles theta and phi (radian)
! -------------------------------------------------------------------
! Reference: D.M. Brink and G.R. Satchler, Angular Momentum,
!            second edition, Oxford University Press, p.22 and p. 145
! -------------------------------------------------------------------
  integer, parameter  :: wp = kind(1.0d0)  ! working precision = double (portable)
!--------------------------------------------------------------------
!   local constants
!--------------------------------------------------------------------
  real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
  complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
!--------------------------------------------------------------------
!   formal arguments
!--------------------------------------------------------------------
  integer, intent(in)  :: l, m
  real(wp), intent(in) :: thrad, phirad
!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
  integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, it, iphase, &
             ia, ib, ic
  real(wp) :: sqrt_fac, sumt, denom, term, dlm0, const, cosb2, sinb2
  complex(wp) :: ylm, exphi
!--------------------------------------------------------------------
!   external function
!--------------------------------------------------------------------
  real(wp), external :: fac10
!--------------------------------------------------------------------
!  program starts here
!  first calculate d^l_{m 0} (theta)
!--------------------------------------------------------------------
  cosb2 = cos(thrad/2.0_wp)
  sinb2 = sin(thrad/2.0_wp)
!--------------------------------------------------------------------
! determine lower and upper limits for summation index it; these
! are derived from the requirement that all factorials n! in the
! denominator are restricted to values with n >=0.
!--------------------------------------------------------------------
  itmin1 = 0
  itmin2 = m
  itmin = max(itmin1,itmin2)
  itmax1 = l+m
  itmax2 = l
  itmax = min(itmax1,itmax2)
!  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
  sqrt_fac = sqrt( fac10(l+m) * fac10(l-m) * fac10(l) * fac10(l) )
!
  sumt = 0.0_wp
  do it = itmin, itmax
     iphase = (-1)**it
     ia = l + m - it
     ib = l - it
     ic = it - m
!     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
     denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
     term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
     sumt = sumt + term
  end do
  dlm0 = sqrt_fac * sumt
!--------------------------------------------------------------------
!  now compute Y_{l m} (theta,phi) from d^l_{m 0} (theta)
!--------------------------------------------------------------------
  const = sqrt( (2.0_wp *l + 1.0_wp) / (4.0_wp * pi) )
  exphi = exp( eye * m * phirad )
  ylm = const * exphi * dlm0
!
  return
  end function ylm

  function fac10 (n)
  implicit none
! -----------------------------------------------
! function fac10(n) calculates factorial(n)/10**n
! -----------------------------------------------
! input: integer n >= 0 (you may want to check this
!        in the program calling this function)
! -----------------------------------------------
  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!------------------------------------------------
!      formal arguments
!------------------------------------------------
  integer, intent(in) :: n
!------------------------------------------------
!      local variables
!------------------------------------------------
  integer :: i
  real(wp) :: fac10, q
! -----------------------------------------------
  if (n == 0) then
     fac10 = 1.0_wp
  else
     fac10 = 1.0_wp
     q = 1.0_wp
     do i = 1, n
        fac10 = fac10 * q / 10.0_wp
        q = q + 1.0_wp
     end do
  endif
!
  return
  end function fac10
