!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!this fortran code decompose the signal on the sphere into spherical harmonics of various modes
!the orignal spherical density rho gives S_rho, and the normalized density eta gives S_eta 
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
   character (len=10) :: seqstring,ana,dummy 
   character (len=128) :: label, correlation
   
   real(kind=8) :: lbox(3), cellsize, distance,P_nn  
   real(kind=8) :: S_ff,r_sp,phi_sp,theta_sp,n_real,n_imag!,n_real_test,n_imag_test!,whatever
   real(kind=8) :: dr, rmax, rmin, d_area,tmp2(3),delta0(3),vect_z(3),&
                   & vect_y(3),vect_x(3),vect_tmp(3),vect_N(n_min,3), vect_tar(n_min,3), minr, maxr
   real(kind=8) :: rad,delta,theta_del, frac !ang_min_1st,
   real(kind=8) :: goldenratio,theta_fib,phi_fib,coord_fib(3),s_cell,totweight
           
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
   real(kind=8),dimension(:,:),allocatable :: wrapped,vect_neighbor,coord_cell
   real(kind=8),dimension(:,:,:,:,:),allocatable :: codcell
   real(kind=8), dimension (:,:), allocatable :: integral_imag,integral_real,area_coordinate 
   real(kind=8), dimension (:), allocatable :: s_f, rho_cell
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

   read(args(1),'(i3)') l_value  !the mode of SH to be calculated: 3 for SiO2, 6 for BLJM 
   read(args(2),'(i6)') nline  !nb of cells on the sphere
   
   !---------read the config file, in xyz format------------------------
   open(unit=12, file='normalized-density.dat',status="old",action="read",iostat=ierror)
     
   allocate (coord_cell(nline,3),rho_cell(nline))
   
   totweight=0
   do j=1,nline
    read(12,*) coord_cell(j,1:3),rho_cell(j)
    totweight=totweight+rho_cell(j)
   enddo
   
  !-------S_rho----------------------------------------
  !calculate spherical harmonics from the density distribution
  allocate(integral_imag(l_value+1,2*l_value+1),integral_real(l_value+1,2*l_value+1),s_f(l_value+1))
  integral_real=0.0; integral_imag=0.0 
  
   print*,"the mode of spherical harmonics to be calculated: l=", l_value 
   open(unit=1000,file='S_eta.dat',action='write') !,position='append')  

    do i=1,nline
     call cart2sphere(coord_cell(i,1:3),r_sp,phi_sp,theta_sp)   !rad rather than degree
      
      l=l_value
      do m_l=-l,l
     
      ylm_numeric = ylm (l, m_l, phi_sp, theta_sp)
	  !the realpart of ylm with precision defined for the complex !realpart(ylm_numeric)
	  n_real= real(ylm_numeric)  !realpart
	  !the imagpart of ylm with precision defined for the compleximagpart(ylm_numeric)
	  n_imag= aimag(ylm_numeric) !imagpart
	  integral_real(l+1,m_l+l_value+1)=integral_real(l+1,m_l+l_value+1)+n_real*rho_cell(i) ! l+1,m+l+1 means l=l & m=m
	  integral_imag(l+1,m_l+l_value+1)=integral_imag(l+1,m_l+l_value+1)-n_imag*rho_cell(i) !
    
      enddo
     enddo 
 
	 l=l_value; S_ff=0.0
				
	 over_mm: do m_l=-l,l
	  S_ff=S_ff+(integral_real(l+1,m_l+l_value+1)/totweight)**2        
	 end do over_mm
			
	 S_ff=S_ff/dble(2*l+1)  !angular power spectrum (for Schmidt semi-normalized harmonics)
	 s_f(l+1)=S_ff  !array to save S_ff at each l
		
	 if (s_f(l+1)>=0) then
	  write (1000,*) sqrt(s_f(l+1))  !r, s_rho
	 endif 
	 close(1000)  
  !--------------------------------------------------- 
  
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
