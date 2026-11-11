! $Id$

!  Bohm quantum-potential force in 3D Cartesian geometry.
!  Adds (hbar^2/2 m^2) grad(lap(sqrt(n))/sqrt(n)) to the momentum equation.
!
!  With f = ln(n) and c = hbar^2/(4 m^2), the force is:
!    F_i = c (d_i d_j d_j f + (d_j f) (d_i d_j f)), summed over j.
!  The mixed third derivatives f_ijj (i /= j) are implemented here.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module

! CPARAM logical, parameter :: lspecial = .true.

! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1

!***************************************************************
module Special

  use Cdata
  use Quiet
  use Messages

  implicit none

  include '../special.h'

  logical :: lbohm_term = .true.
  real(KIND=rkind8) :: m_particle_cgs = 0.0d0

  real(KIND=rkind8) :: m_particle = 0.0d0
  real :: qhd_bohm_coeff = 0.0
  integer :: ibohm_pressure = 0

  integer :: bohm_stencil_order = 0

  logical, dimension(3) :: lactive = .false.

! Central-difference weights of the mixed stencil on an equidistant grid, set from the order.
! d1_weight(a) multiplies f(+a) - f(-a); d2_weight(|a|) multiplies f(+a) + f(-a), a=0 once.
! stencil_reach is the number of points on each side, half the order.
  integer :: stencil_reach = 0
  real, dimension(:), allocatable :: d1_weight, d2_weight

  namelist /special_init_pars/ lbohm_term, m_particle_cgs, bohm_stencil_order

  namelist /special_run_pars/ lbohm_term, m_particle_cgs, bohm_stencil_order

  contains
!***********************************************************************
    subroutine register_special()

      use FArrayManager

! The Bohm term only adds to existing hydro/density equations.
! The auxiliary bohm_pressure P_Q = -c rho lap(ln rho) is stored for output.
      call svn_id("$Id$")

      call farray_register_auxiliary('bohm_pressure',ibohm_pressure)

      if (naux+naux_com < maux+maux_com) aux_var(aux_count)=',bohm_pressure $'
      if (naux+naux_com == maux+maux_com) aux_var(aux_count)=',bohm_pressure'
      aux_count=aux_count+1

      if (lroot) write(15,*) 'bohm_pressure = fltarr(mx,my,mz)*one'

    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)

    real(KIND=rkind8) :: m_particle = 0.0d0

      real, dimension(mx,my,mz,mfarray) :: f

      m_particle = m_particle_cgs/unit_mass
      lactive = (/ nxgrid/=1, nygrid/=1, nzgrid/=1 /)

      if (lbohm_term) then
        if (.not.ldensity .or. .not.lhydro) &
            call fatal_error('initialize_special','Bohm term requires DENSITY and HYDRO')
        if (ldensity_nolog) &
            call fatal_error('initialize_special','Bohm term requires logarithmic density')
        if (lcylindrical_coords .or. lspherical_coords) &
            call fatal_error('initialize_special','Bohm term requires Cartesian coordinates')
        if (any(lactive .and. .not.lequidist)) &
            call fatal_error('initialize_special','Bohm term requires equidistant grids')
        if (nghost<2) &
            call fatal_error('initialize_special','Bohm term requires a DERIV module with nghost>=2')
        if (trim(der2_type)/='standard') &
            call fatal_error('initialize_special',"Bohm term requires der2_type='standard'")
        if (m_particle<=0.0d0) &
            call fatal_error('initialize_special', 'set m_particle_cgs>0 (special_init_pars/run_pars)')

        call initialize_stencils()

        qhd_bohm_coeff = real(0.25d0*hbar**2/m_particle**2)
      else
        qhd_bohm_coeff = 0.0
      endif

      call keep_compiler_quiet(f)

    endsubroutine initialize_special
!***********************************************************************
    subroutine initialize_stencils()

! Set the order and the weights of the mixed stencils.

      if (bohm_stencil_order == 0) then
        bohm_stencil_order = 2*nghost
      endif

      if (bohm_stencil_order < 2 .or. bohm_stencil_order > 2*nghost .or. &
          mod(bohm_stencil_order, 2) /= 0) then
        call fatal_error('initialize_stencils', &
            'bohm_stencil_order must be even and between 2 and 2*nghost (0 selects 2*nghost)')
      endif

      stencil_reach = bohm_stencil_order/2
      if (allocated(d1_weight)) then
        deallocate(d1_weight, d2_weight)
      endif
      allocate(d1_weight(stencil_reach), d2_weight(0:stencil_reach))
      d1_weight = first_derivative_weights(stencil_reach)
      d2_weight = second_derivative_weights(stencil_reach)

      if (lroot) then
        print*, 'initialize_stencils: mixed stencil of order', bohm_stencil_order, &
            ' (nghost =', nghost, ')'
      endif

    endsubroutine initialize_stencils
!***********************************************************************
    subroutine init_special(f)

      real, dimension(mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()

!  The Bohm term needs grad(ln(n)) and its Hessian.
!  n is taken proportional to the fluid mass density rho.
      if (lbohm_term) then
        lpenc_requested(i_glnrho) = .true.
        lpenc_requested(i_hlnrho) = .true.
        lpenc_requested(i_del2lnrho) = .true.
        lpenc_requested(i_rho) = .true.
      endif

    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)

      logical, dimension(npencils) :: lpencil_in

      call keep_compiler_quiet(lpencil_in)

    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)

      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)

! Add the Bohm quantum-potential force to the momentum equation.

      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

      integer :: i

      if (lbohm_term) then
        f(l1:l2,m,n,ibohm_pressure) = -qhd_bohm_coeff*p%rho*p%del2lnrho

        do i = 1, 3
          if (.not.lactive(i)) cycle
          df(l1:l2,m,n,iux+i-1) = df(l1:l2,m,n,iux+i-1) + qhd_bohm_coeff*( &
              grad_laplacian_lnrho(f,i) + glnrho_dot_hessian(p,i))
        enddo

! The linear Bohm wave has omega_B^2 = coeff*k^4.
! Bound the time step at the largest grid wavevector, k^2 = pi^2 sum_i dx_i^-2.
        if (lupdate_courant_dt) &
            advec_cs2 = max(advec_cs2, qhd_bohm_coeff*nyquist_k2(dx_1(l1:l2),dy_1(m),dz_1(n))**2)
      endif

    endsubroutine special_calc_hydro
!***********************************************************************
    function glnrho_dot_hessian(p,i) result(g)

! sum_j (d_j ln(n)) (d_i d_j ln(n))

      type (pencil_case), intent(in) :: p
      integer, intent(in) :: i
      real, dimension(nx) :: g

      integer :: j

      g = 0.0
      do j = 1, 3
        g = g + p%glnrho(:,j)*p%hlnrho(:,i,j)
      enddo

    endfunction glnrho_dot_hessian
!***********************************************************************
    function grad_laplacian_lnrho(f,i) result(g)

! d_i lap(ln(n)) = sum_j d_i d_j d_j ln(n)

      use Deriv, only: der3

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: i
      real, dimension(nx) :: g

      real, dimension(nx) :: d3lnrho
      integer :: j

      g = 0.0
      do j = 1, 3
        if (.not.lactive(j)) cycle
        if (j==i) then
          call der3(f,ilnrho,d3lnrho,i)
          g = g + d3lnrho
        else
          g = g + mixed_third_derivative_lnrho(f,i,j)
        endif
      enddo

    endfunction grad_laplacian_lnrho
!***********************************************************************
    function mixed_third_derivative_lnrho(f,i,j) result(d3)

! d_i d_j d_j ln(n) at the current pencil (m,n), for i /= j.

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: i, j
      real, dimension(nx) :: d3

      real, dimension(3) :: inv_dx
      integer, dimension(3) :: offset
      real :: w1
      integer :: a, b

      inv_dx = (/ dx_1(l1), dy_1(m), dz_1(n) /)

      d3 = 0.0
      do a = -stencil_reach, stencil_reach
        if (a==0) cycle
        w1 = sign(1,a)*d1_weight(abs(a))
        do b = -stencil_reach, stencil_reach
          offset = 0
          offset(i) = a
          offset(j) = b
          d3 = d3 + w1*d2_weight(abs(b))*f(l1+offset(1):l2+offset(1), &
                                          m+offset(2),n+offset(3),ilnrho)
        enddo
      enddo
      d3 = d3*inv_dx(i)*inv_dx(j)**2

    endfunction mixed_third_derivative_lnrho
!***********************************************************************
    function first_derivative_weights(p) result(w)

! Weights of the central first derivative of order 2p, in units of 1/h:
!   f'(0) = sum_{k=1}^p w(k) (f(+k) - f(-k)).
! Closed-form solution of the Taylor moment conditions.

      integer, intent(in) :: p
      real, dimension(p) :: w

      integer :: k

      do k = 1, p
        w(k) = (-1.0)**(k+1)*factorial(p)**2/(k*factorial(p-k)*factorial(p+k))
      enddo

    endfunction first_derivative_weights
!***********************************************************************
    function second_derivative_weights(p) result(w)

! Weights of the central second derivative of order 2p, in units of 1/h^2:
!   f''(0) = w(0) f(0) + sum_{k=1}^p w(k) (f(+k) + f(-k)).

      integer, intent(in) :: p
      real, dimension(0:p) :: w

      integer :: k

      do k = 1, p
        w(k) = 2.0*(-1.0)**(k+1)*factorial(p)**2/(k**2*factorial(p-k)*factorial(p+k))
      enddo
      w(0) = -2.0*sum(w(1:p))

    endfunction second_derivative_weights
!***********************************************************************
    function factorial(n) result(fact)

      integer, intent(in) :: n
      real :: fact

      integer :: i

      fact = 1.0
      do i = 2, n
        fact = fact*i
      enddo

    endfunction factorial
!***********************************************************************
    function nyquist_k2(dx_1_x,dy_1_m,dz_1_n) result(k2)

! (largest wavevector)^2 on the grid, over the active directions.

      real, dimension(nx), intent(in) :: dx_1_x
      real, intent(in) :: dy_1_m, dz_1_n
      real, dimension(nx) :: k2

      k2 = 0.0
      if (lactive(1)) k2 = k2 + dx_1_x**2
      if (lactive(2)) k2 = k2 + dy_1_m**2
      if (lactive(3)) k2 = k2 + dz_1_n**2
      k2 = pi**2*k2

    endfunction nyquist_k2
!***********************************************************************
    subroutine read_special_init_pars(iomsg)

      use File_io, only: parallel_unit

      character(LEN=iomsglen), intent(out) :: iomsg
      integer :: iostat

      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""

    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)

      integer, intent(in) :: unit

      write(unit, NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iomsg)

      use File_io, only: parallel_unit

      character(LEN=iomsglen), intent(out) :: iomsg
      integer :: iostat

      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""

    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)

      integer, intent(in) :: unit

      write(unit, NML=special_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************

!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
