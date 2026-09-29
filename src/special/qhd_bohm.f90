! $Id$

!  Bohm quantum-potential force.
!  Adds (hbar^2/2 m^2) grad(lap(sqrt(n))/sqrt(n)) to the momentum equation.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module

! CPARAM logical, parameter :: lspecial = .true.

! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0

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
  real(KIND=rkind8) :: qhd_bohm_coeff = 0.0d0

  namelist /special_init_pars/ lbohm_term, m_particle_cgs

  namelist /special_run_pars/ lbohm_term, m_particle_cgs

  contains
!***********************************************************************
    subroutine register_special()

! The Bohm term only adds to existing hydro/density equations.
      call svn_id("$Id$")

    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)

      real, dimension(mx,my,mz,mfarray) :: f

      m_particle = m_particle_cgs/unit_mass

      if (lbohm_term) then
        if (.not.ldensity .or. .not.lhydro) &
            call fatal_error('initialize_special','Bohm term requires density and hydro')
        if (ldensity_nolog) &
            call fatal_error('initialize_special','Bohm term requires logarithmic density')
        if (nygrid/=1.or.nzgrid/=1 .or. lcylindrical_coords .or. lspherical_coords) &
            call fatal_error('initialize_special','Bohm term currently supports 1D Cartesian x only')
        if (.not.lequidist(1)) &
            call fatal_error('initialize_special','Bohm term requires an equidistant x grid')
        if (m_particle<=0.0d0 .or. hbar<=0.) &
            call fatal_error('initialize_special','set m_particle_cgs>0 (special_init_pars/run_pars) and hbar must be positive in code units')
        qhd_bohm_coeff = 0.25d0*hbar**2/m_particle**2
      else
        qhd_bohm_coeff = 0.0d0
      endif

      call keep_compiler_quiet(f)

    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)

      real, dimension(mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()

!  The Bohm term needs grad(ln n) and lap(ln n).
!  n is taken proportional to the fluid mass density rho.
      if (lbohm_term) then
        lpenc_requested(i_glnrho) = .true.
        lpenc_requested(i_del2lnrho) = .true.
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

      use Deriv, only: der3

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

      real, dimension(nx) :: d3lnrho

      if (lbohm_term) then
        call der3(f,ilnrho,d3lnrho,1)

        df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + &
            qhd_bohm_coeff*(d3lnrho + p%glnrho(:,1)*p%del2lnrho)
      endif

      call keep_compiler_quiet(f)

    endsubroutine special_calc_hydro
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
