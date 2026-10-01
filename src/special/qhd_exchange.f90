! $Id$

!  Dirac exchange force of a zero-temperature electron gas.
!  Adds -(1/rho) grad(P_X) to the momentum equation with
!  P_X = -(3/pi)^(1/3) e^2 n^(4/3) / 4, i.e. the acceleration
!  (3/pi)^(1/3) e^2 n^(1/3) / (3 m) grad(ln rho).
!  The auxiliary variable exchange_pressure holds P_X.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of qhd_exchange_dummies.inc) the number of f array
! variables and auxiliary variables added by this module

! CPARAM logical, parameter :: lspecial = .true.

! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1

!***************************************************************
module qhd_exchange

  use Cdata
  use Quiet
  use Messages

  implicit none

  include '../special.h'

  logical :: lexchange_term = .true.
  real(KIND=rkind8) :: m_particle_cgs = 0.0d0

  real(KIND=rkind8) :: m_particle = 0.0d0
  real :: qhd_exchange_coeff = 0.0
  real :: inv_m_particle = 0.0
  integer :: iexchange_pressure = 0

  namelist /qhd_exchange_init_pars/ lexchange_term, m_particle_cgs

  namelist /qhd_exchange_run_pars/ lexchange_term, m_particle_cgs

  contains
!***********************************************************************
    subroutine register_special()

      use FArrayManager

! The exchange term only adds to existing hydro/density equations.
      call svn_id("$Id$")

      call farray_register_auxiliary('exchange_pressure',iexchange_pressure)

      if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=',exchange_pressure $'
      if (naux+naux_com  == maux+maux_com) aux_var(aux_count)=',exchange_pressure'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'exchange_pressure = fltarr(mx,my,mz)*one'

    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)

      real, dimension(mx,my,mz,mfarray) :: f

      real(KIND=rkind8) :: e2_code

      m_particle = m_particle_cgs/unit_mass

      if (lexchange_term) then
        if (.not.ldensity .or. .not.lhydro) &
            call fatal_error('initialize_special','Exchange term requires DENSITY and HYDRO')
        if (unit_system/='cgs') &
            call fatal_error('initialize_special','Exchange term requires unit_system=cgs')
        if (m_particle<=0.0d0) &
            call fatal_error('initialize_special', 'set m_particle_cgs>0 (special_init_pars/run_pars)')

        e2_code = e_cgs**2/(dble(unit_energy)*dble(unit_length))
        qhd_exchange_coeff = real((3.0d0/(4.0d0*atan(1.0d0)))**(1.0d0/3.0d0)*e2_code/(3.0d0*m_particle))
        inv_m_particle = real(1.0d0/m_particle)
      else
        qhd_exchange_coeff = 0.0
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

!  The exchange force needs n = rho/m, grad(ln rho) and, for the
!  stability check, the sound speed.
      if (lexchange_term) then
        lpenc_requested(i_rho) = .true.
        lpenc_requested(i_glnrho) = .true.
        lpenc_requested(i_cs2) = .true.
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

! Add the exchange force to the momentum equation.
! The exchange pressure lowers the effective sound speed squared by
! coeff*n^(1/3); the total must stay positive.

      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

      real, dimension(nx) :: exchange_cs2
      integer :: j

      if (lexchange_term) then
        exchange_cs2 = qhd_exchange_coeff*(p%rho*inv_m_particle)**(1.0/3.0)
        f(l1:l2,m,n,iexchange_pressure) = -0.75*real(m_particle)*exchange_cs2*p%rho*inv_m_particle

        if (lupdate_courant_dt) then
          if (any(p%cs2 <= exchange_cs2)) &
              call fatal_error('special_calc_hydro', &
              'exchange pressure exceeds the degenerate pressure: sound speed squared <= 0')
        endif

        do j = 1, 3
          df(l1:l2,m,n,iuu+j-1) = df(l1:l2,m,n,iuu+j-1) + exchange_cs2*p%glnrho(:,j)
        enddo
      endif

      call keep_compiler_quiet(f)

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine read_special_init_pars(iomsg)

      use File_io, only: parallel_unit

      character(LEN=iomsglen), intent(out) :: iomsg
      integer :: iostat

      read(parallel_unit, NML=qhd_exchange_init_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""

    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)

      integer, intent(in) :: unit

      write(unit, NML=qhd_exchange_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iomsg)

      use File_io, only: parallel_unit

      character(LEN=iomsglen), intent(out) :: iomsg
      integer :: iostat

      read(parallel_unit, NML=qhd_exchange_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""

    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)

      integer, intent(in) :: unit

      write(unit, NML=qhd_exchange_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************

!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../qhd_exchange_dummies.inc'
!********************************************************************
endmodule qhd_exchange
