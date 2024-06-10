!! Taken from Octopus (2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch)
!! Adapted to oncvpsp by A. Castaneda M. (2019)
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! generic interfaces
module parser
  implicit none
  interface parse_variable
     module procedure real_parser ! interface for a module
     module procedure logical_parser ! procedure is implicit
  end interface parse_variable
contains
  subroutine real_parser(varval, var)
    real(8), intent(in) :: varval
    real(8), intent(out) :: var
    var = varval
  end subroutine real_parser
  subroutine logical_parser(varval, var)
    logical, intent(in) :: varval
    logical, intent(out) :: var
    var = varval
  end subroutine logical_parser
end module parser

! For messaging
module xc_messages
  implicit none
  ! Messages
  character(1200) :: message(2)
contains
  subroutine messages_info(lmess, iunit)
    integer, intent(in) :: lmess, iunit
    integer :: i
    do i=1,lmess
       write(iunit,'("LXC ",A)') trim(message(i))
    end do
  end subroutine messages_info

  subroutine messages_fatal(lmess)
    integer, intent(in) :: lmess
    integer :: i
    do i=1,lmess
       write(*,'("Fatal ",A)') trim(message(i))
    end do
    stop
  end subroutine messages_fatal

  subroutine messages_warning(lmess)
    integer, intent(in) :: lmess
    integer :: i
    do i=1,lmess
       write(*,'("Warning ",A)') trim(message(i))
    end do
  end subroutine messages_warning

  subroutine messages_input_error(str1,str2)
    character(*), intent(in) :: str1,str2
    write(*,'("Input error ",A,A)') trim(str1),trim(str2)
  end subroutine messages_input_error

  subroutine messages_experimental(str)
    character(*), intent(in) :: str
    write(*,'("Experimental xc ",A)') trim(str)
  end subroutine messages_experimental

  subroutine messages_not_implemented(str)
    character(*), intent(in) :: str
    write(*,'("Not implemented ",A)') trim(str)
  end subroutine messages_not_implemented
end module xc_messages

! The main module
module functionals_m
  use parser
  use xc_messages
  use xc_f90_lib_m !libxc
  use iso_c_binding
#include "xc_version.h"

  implicit none

  ! Parameters
  integer, parameter :: iunit = 6
  real(8), parameter :: pi = 3.141592653589793238462643383279502884197d0
  integer(c_size_t), parameter, private :: xc_one = 1

  private
  public ::                     &
       xc_functl_t,                &
       xc_functl_init_functl,      &
       xc_functl_end,              &
       xc_functl_write_info, &
       xc_functl_get_vxc

  integer, public, parameter :: XC_DERIV_NUMERICAL = 1,XC_DERIV_ANALYTICAL = 2

  type xc_functl_t
     integer         :: family            !< LDA, GGA, etc.
     integer         :: type              !< exchange, correlation, or exchange-correlation
     integer         :: id                !< identifier

     integer         :: nspin             !< XC_UNPOLARIZED | XC_POLARIZED
     integer         :: flags             !< XC_FLAGS_HAVE_EXC + XC_FLAGS_HAVE_VXC + ...

     type(xc_f90_func_t) :: conf          !< the pointer used to call the library
     type(xc_f90_func_info_t) :: info     !< information about the functional

     integer         :: LB94_modified     !< should I use a special version of LB94 that
     real(8)         :: LB94_threshold    !< needs to be handled specially

     integer  :: deriv_method
  end type xc_functl_t

contains

  ! ---------------------------------------------------------
  subroutine xc_functl_init(functl, nspin, deriv_method)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: nspin, deriv_method


    functl%family = 0
    functl%type   = 0
    functl%id     = 0
    functl%flags  = 0
    functl%nspin = nspin
    functl%deriv_method = deriv_method

  end subroutine xc_functl_init

  ! ---------------------------------------------------------

  subroutine xc_functl_init_functl(functl, id, ndim,nel, nspin, deriv_method)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: id
    integer,           intent(in)  :: ndim
    real(8),             intent(in)  :: nel
    integer,           intent(in)  :: nspin
    integer,           intent(in)  :: deriv_method

    real(8)   :: alpha
    real(8)   :: parameters(2)
    logical :: ok, lb94_modified

#if XC_MAJOR_VERSION<5
    call messages_input_error('LibXC version', 'at least v5 is now required')
#endif

    ! initialize structure
    call xc_functl_init(functl, nspin, deriv_method)

    functl%id = id

    if(functl%id == 0) then
       functl%family = XC_FAMILY_NONE
    else
       ! get the family of the functional
       functl%family = xc_f90_family_from_id(functl%id)
       ! this also ensures it is actually a functional defined by the linked version of libxc

       if(functl%family == XC_FAMILY_UNKNOWN) then
          call messages_input_error('XCFunctional', 'Unknown functional')
       end if
    end if

    if(functl%family == XC_FAMILY_OEP) then
       functl%type = XC_EXCHANGE

    else if(functl%family  ==  XC_FAMILY_NONE) then
       functl%type = -1
       functl%flags = 0
    else ! handled by libxc
       ! initialize
       call xc_f90_func_init(functl%conf, functl%id, nspin)
       functl%info     = xc_f90_func_get_info(functl%conf)
       functl%type     = xc_f90_func_info_get_kind(functl%info)
       functl%flags    = xc_f90_func_info_get_flags(functl%info)

       ! FIXME: no need to say this for kernel
       if(iand(functl%flags, XC_FLAGS_HAVE_EXC) == 0) then
          message(1) = 'Specified functional does not have total energy available.'
          message(2) = 'Corresponding component of energy will just be left as zero.'
          call messages_warning(2)
       end if

       if(iand(functl%flags, XC_FLAGS_HAVE_VXC) == 0) then
          message(1) = 'Specified functional does not have XC potential available.'
          message(2) = 'Cannot run calculations. Choose another XCFunctional.'
          call messages_fatal(2)
       end if
    end if

    !    XC_NON_RELATIVISTIC     =   0
    !    XC_RELATIVISTIC         =   1

    ! FIXME: aren`t there other parameters that can or should be set?
    ! special parameters that have to be configured
    select case(functl%id)
    case(XC_LDA_C_XALPHA)
       ! FIXME: aren`t there other Xalpha functionals?
       ! Variable Xalpha
       ! The parameter of the Slater X<math>\alpha</math> functional
       call parse_variable(1.0d0, alpha)
       parameters(1) = alpha
       call xc_f90_func_set_ext_params(functl%conf, parameters(1))
    case(XC_GGA_X_LB)
       ! FIXME: libxc has XC_GGA_X_LBM, isn`t that the modified one?
       ! Whether to use a modified form of the LB94 functional
       call parse_variable(.false., lb94_modified)
       if(lb94_modified) then
          functl%LB94_modified = 1
       else
          functl%LB94_modified = 0
       end if
       ! FIXME: libxc seems to have 1e-32 as a threshold, should we not use that?
       ! A threshold for the LB94 functional
       call parse_variable(1.0d-6, functl%LB94_threshold)
    end select

  end subroutine xc_functl_init_functl


  ! ---------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_t), intent(inout) :: functl


    if(functl%family /= XC_FAMILY_NONE .and. functl%family /= XC_FAMILY_OEP)  then
       call xc_f90_func_end(functl%conf)
    end if

  end subroutine xc_functl_end


  ! ---------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit)
    type(xc_functl_t), intent(in) :: functl
    integer,           intent(in) :: iunit

    character(len=1000) :: s1, s2
    integer :: ii
    type(xc_f90_func_reference_t) :: xc_ref


    if(functl%family /= XC_FAMILY_NONE) then ! all the other families
       select case(functl%type)
       case(XC_EXCHANGE)
          write(message(1), '(2x,a)') 'Exchange'
       case(XC_CORRELATION)
          write(message(1), '(2x,a)') 'Correlation'
       case(XC_EXCHANGE_CORRELATION)
          write(message(1), '(2x,a)') 'Exchange-correlation'
       case(XC_KINETIC)
          call messages_not_implemented("kinetic-energy functionals")
       case default
          write(message(1), '(a,i6,a,i6)') "Unknown functional type ", functl%type, ' for functional ', functl%id
          call messages_fatal(1)
       end select

       s1 = xc_f90_func_info_get_name(functl%info)
       select case(functl%family)
       case (XC_FAMILY_LDA);       write(s2,'(a)') "LDA"
       case (XC_FAMILY_GGA);       write(s2,'(a)') "GGA"
       case (XC_FAMILY_HYB_GGA);   write(s2,'(a)') "Hybrid GGA"
       case (XC_FAMILY_HYB_MGGA);  write(s2,'(a)') "Hybrid MGGA"
       case (XC_FAMILY_MGGA);      write(s2,'(a)') "MGGA"
       end select
       write(message(2), '(4x,4a)') trim(s1), ' (', trim(s2), ')'
       call messages_info(2, iunit)

       ii = 0
       xc_ref = xc_f90_func_info_get_references(functl%info, ii)
       s1 = xc_f90_func_reference_get_ref(xc_ref)
       do while(ii >= 0)
          write(message(1), '(4x,a,i1,2a)') '[', ii, '] ', trim(s1)
          call messages_info(1, iunit)
          xc_ref = xc_f90_func_info_get_references(functl%info, ii)
          s1 = xc_f90_func_reference_get_ref(xc_ref)
       end do
    end if
  end subroutine xc_functl_write_info

  subroutine xc_functl_get_vxc(functl, np, al, rr, rho, rho_grad, rho_lapl, tau, ip, v, e, vtau)
    !-----------------------------------------------------------------------!
    ! Given a density, computes the corresponding exchange/correlation      !
    ! potentials and energies.                                              !
    !                                                                       !
    !  functl   - functional                                                !
    !  np       - number of mesh points                                     !
    !  al       - grid parameter                                            !
    !  rr       - radial points                                             !
    !  rho      - electronic radial density                                 !
    !  rho_grad - gradient of the electronic radial density                 !
    !  rho_lapl - laplacian of the electronic radial density                !
    !  tau      - radial kinetic energy density                             !
    !  ip       - ionization potential                                      !
    !  v        - potential                                                 !
    !  e        - energy per-volume                                         !
    !  vtau     - extra term arising from MGGA potential                    !
    !-----------------------------------------------------------------------!
    type(xc_functl_t), intent(inout) :: functl
    integer     ,       intent(in)    :: np
    real(8),           intent(in)    :: al
    real(8),           intent(in)    :: rr(np)
    real(8),           intent(in)    :: rho(np, functl%nspin)
    real(8),           intent(in)    :: rho_grad(np, functl%nspin)
    real(8),           intent(in)    :: rho_lapl(np, functl%nspin)
    real(8),           intent(in)    :: tau(np, functl%nspin)
    real(8),           intent(in)    :: ip(functl%nspin)
    real(8),           intent(out)   :: v(np, functl%nspin), e(np)
    real(8),           intent(out)   :: vtau(np, functl%nspin)

    integer  :: i, is, nspin
    real(8) :: a, b, c
    real(8), parameter   :: alpha = -0.012d0, beta = 1.023d0

    ! Global variables
    real(8), allocatable :: dedrho(:,:), dedgrad(:,:), dedlapl(:,:), dedtau(:,:)
    real(8), allocatable :: d2edrhodgrad(:,:), d2edgrad2(:,:)

    ! Local variables
    real(8), allocatable :: n(:), s(:), l(:), t(:)
    real(8), allocatable :: dedn(:), deds(:), dedl(:), dedt(:)
    real(8), allocatable :: d2edn2(:), d2eds2(:), d2ednds(:)

    real(8), allocatable :: dpr(:), dppr(:), dlap(:)
    real(8)   :: parameters(2)


    if (.not. (size(v, dim=2) == functl%nspin)) stop 'functionals: ERROR bad nspin definition'

    ! Initialize all output quantities to zero
    v =0.0d0 ; e =0.0d0 ; vtau =0.0d0

    ! If the functional is not set, there is nothing to be done
    if (functl%family == 0) then
       return
    end if

    ! Shortcut
    nspin = functl%nspin


    ! Compute c parameter of the TB09 functional
    if (functl%id == XC_MGGA_X_TB09) then
       if (maxval(ip) ==0.0d0) then
          c = 1.0d0
       else
          a =0.0d0
          do
             c = alpha + beta*sqrt(2.0d0*sqrt(2.0d0*(maxval(ip) + a)))
             b = (3.0d0*c - 2.0d0)/pi*sqrt(5.0d0/6.0d0*(maxval(ip) + a))
             if (abs(a - b) < 1.0e-8) exit
             a = b
          end do
       end if
       parameters(1) = c
       call xc_f90_func_set_ext_params(functl%conf, parameters(1))
    end if


    !---Allocate work arrays---!

    ! LDA
    allocate(n(nspin), dedn(nspin))
    allocate(dedrho(np, nspin))
    n =0.0d0; dedn =0.0d0
    dedrho =0.0d0

    if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
       if (nspin == 1) then
          allocate(d2edn2(1))
       else
          allocate(d2edn2(3))
       end if
       d2edn2 =0.0d0
    end if

    ! GGA
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
       if (nspin == 1) then
          allocate(s(1), deds(1))
       else
          allocate(s(3), deds(3))
       end if
       allocate(dedgrad(np, nspin))
       s =0.0d0; deds =0.0d0
       dedgrad =0.0d0

       if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          if (nspin == 1) then
             allocate(d2eds2(1), d2ednds(1))
             allocate(d2edgrad2(np, 1))
             allocate(d2edrhodgrad(np, 1))
          else
             allocate(d2eds2(6), d2ednds(6))
             allocate(d2edgrad2(np, 3))
             allocate(d2edrhodgrad(np, 4))
          end if
          d2eds2 =0.0d0; d2ednds =0.0d0
          d2edgrad2 =0.0d0; d2edrhodgrad =0.0d0
       end if

    end if

    ! MGGA
    if (functl%family == XC_FAMILY_MGGA) then
       allocate(l(nspin), dedl(nspin))
       allocate(t(nspin), dedt(nspin))
       allocate(dedlapl(np, nspin))
       allocate(dedtau(np, nspin))
       l =0.0d0; dedl =0.0d0
       t =0.0d0; dedt =0.0d0
       dedlapl =0.0d0
       dedtau =0.0d0

       if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          !Not yet implemented
       end if

    end if


    !---Space loop---!

    do i = 1, np
       ! make a local copy with the correct memory order
       n(1:nspin) = rho(i, 1:nspin)
       if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
          s(1) = rho_grad(i, 1)**2
          if(nspin == 2) then
             s(2) = rho_grad(i, 1)*rho_grad(i, 2)
             s(3) = rho_grad(i, 2)**2
          end if
       end if
       if (functl%family == XC_FAMILY_MGGA) then
          t(1:nspin) = tau(i, 1:nspin)/2.0d0
          l(1:nspin) = rho_lapl(i, 1:nspin)
       end if

       if (iand(xc_f90_func_info_get_flags(functl%info), XC_FLAGS_HAVE_EXC) .ne. 0) then

          select case(functl%family)
          case(XC_FAMILY_LDA)
             call xc_f90_lda_exc_vxc(functl%conf, xc_one, n(1), e(i), dedn(1))
          case(XC_FAMILY_GGA)
             call xc_f90_gga_exc_vxc(functl%conf, xc_one, n(1), s(1), e(i), dedn(1), deds(1))
          case(XC_FAMILY_MGGA)
             call xc_f90_mgga_exc_vxc(functl%conf, xc_one, n(1), s(1), l(1), t(1), e(i), &
                  dedn(1), deds(1), dedl(1), dedt(1))
          end select

       else !Just get the potential

          select case(functl%family)
          case(XC_FAMILY_LDA)
             call xc_f90_lda_vxc(functl%conf, xc_one, n(1), dedn(1))
          case(XC_FAMILY_GGA)
             call xc_f90_gga_vxc(functl%conf, xc_one, n(1), s(1), dedn(1), deds(1))
          case(XC_FAMILY_MGGA)
             call xc_f90_mgga_vxc(functl%conf, xc_one, n(1), s(1), l(1), t(1), &
                  dedn(1), deds(1), dedl(1), dedt(1))
          end select
          e(i) =0.0d0

       end if

       e(i) = e(i)*sum(n)
       dedrho(i, :) = dedn(:)
       if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
          if (nspin == 1) then
             dedgrad(i, 1) = 2.0d0*deds(1)*rho_grad(i, 1)
          else
             dedgrad(i, 1) = 2.0d0*deds(1)*rho_grad(i, 1) + deds(2)*rho_grad(i, 2)
             dedgrad(i, 2) = 2.0d0*deds(3)*rho_grad(i, 2) + deds(2)*rho_grad(i, 1)
          end if
       end if

       if(functl%family == XC_FAMILY_MGGA) then
          dedlapl(i, 1:nspin) = dedl(1:nspin)
          dedtau(i, 1:nspin) = dedt(1:nspin)/2.0d0
       end if

       if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          !Evaluate second-order derivatives
          if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
             call xc_f90_gga_fxc(functl%conf, xc_one, n(1), s(1), d2edn2(1), d2ednds(1), d2eds2(1))

             if (nspin == 1) then
                d2edrhodgrad(i, 1) = 2.0d0*rho_grad(i, 1)*d2ednds(1)
                d2edgrad2(i, 1) = 2.0d0*deds(1) + 4.0d0*s(1)*d2eds2(1)
             else
                d2edrhodgrad(i, 1) = 2.0d0*rho_grad(i, 1)*d2ednds(1) + rho_grad(i, 2)*d2ednds(2)
                d2edrhodgrad(i, 2) = 2.0d0*rho_grad(i, 1)*d2ednds(4) + rho_grad(i, 2)*d2ednds(5)
                d2edrhodgrad(i, 3) = 2.0d0*rho_grad(i, 2)*d2ednds(3) + rho_grad(i, 1)*d2ednds(2)
                d2edrhodgrad(i, 4) = 2.0d0*rho_grad(i, 2)*d2ednds(6) + rho_grad(i, 1)*d2ednds(5)

                d2edgrad2(i, 1) = 2.0d0*deds(1) + 4.0d0*s(1)*d2eds2(1) + 4.0d0*s(2)*d2eds2(2) + s(3)*d2eds2(4)
                d2edgrad2(i, 2) = deds(2) + 4.0d0*s(2)*d2eds2(3) + 2.0d0*s(1)*d2eds2(2) + 2.0d0*s(3)*d2eds2(5) + s(2)*d2eds2(4)
                d2edgrad2(i, 3) = 2.0d0*deds(3) + 4.0d0*s(3)*d2eds2(6) + 4.0d0*s(2)*d2eds2(5) + s(1)*d2eds2(4)
             end if

          else if (functl%family == XC_FAMILY_MGGA) then
             !Not yet implemented
          end if
       end if

    end do ! loop over points i


    !---Compute potentials---!

    ! LDA contribution
    v = dedrho

    ! GGA contribution
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
       if (functl%deriv_method == XC_DERIV_NUMERICAL) then
          allocate (dpr(np))
          allocate (dppr(np))
          allocate (dlap(np))
          do is = 1, nspin
             call derivs(np, dedgrad(:,is), al, rr, dpr, dppr, dlap)
             !          v(:, is) = v(:, is) - mesh_divergence(m, dedgrad(:, is))
             v(:, is) = v(:, is) - (2.0d0/rr(:)*dedgrad(:,is) + dpr(:))
          end do
          deallocate (dpr, dppr, dlap)
       elseif (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          stop 'functionals: ERROR - no analytical derivatives coded'
       end if
    end if

    ! MGGA contribution
    if (functl%family == XC_FAMILY_MGGA) then
       stop 'functrionals: ERROOR : meta GGA not coded yet '
       if (functl%deriv_method == XC_DERIV_NUMERICAL) then
          allocate (dpr(np))
          allocate (dppr(np))
          allocate (dlap(np))
          do is = 1, nspin
             call derivs(np, dedlapl(:,is), al, rr, dpr, dppr, dlap)
             !          v(:, is) = v(:, is) + mesh_laplacian(m, dedlapl(:, is))
             v(:, is) = v(:, is) + dlap(:)
          end do
          deallocate (dpr, dppr, dlap)
       elseif (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          !Not yet implemented
       end if

       vtau = dedtau
    end if

    !Shift potentials that do not go to zero at infinity
    do is = 1, nspin
       select case (functl%id)
       case (XC_MGGA_X_BJ06)
          a = sqrt(5.0d0/6.0d0)/pi
       case (XC_MGGA_X_TB09)
          a = (3.0d0*c - 2.0d0)*sqrt(5.0d0/6.0d0)/pi
       case (XC_GGA_X_AK13)
          a = sqrt(2.0d0)*(1.0d0/(54.0d0*pi) + 2.0d0/15.0d0)
       case default
          a =0.0d0
       end select
       do i = np, 1, -1
          if (v(i, is) /=0.0d0) then
             v(1:i, is) = v(1:i, is) - a*sqrt(ip(is))
             exit
          end if
       end do
    end do

    !---Deallocate arrays---!

    ! LDA
    deallocate(n, dedn, dedrho)
    if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
       deallocate(d2edn2)
    end if

    ! GGA
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
       deallocate(s, deds, dedgrad)
       if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          deallocate(d2eds2, d2ednds, d2edgrad2, d2edrhodgrad)
       end if
    end if

    ! MGGA
    if (functl%family == XC_FAMILY_MGGA) then
       deallocate(l, dedl, dedlapl)
       deallocate(t, dedt, dedtau)
       if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
          message(1) = 'Not implemented analytical derivatives'
          call messages_fatal(1)
       end if
    end if

  end subroutine xc_functl_get_vxc

  ! Not used for the moment
  subroutine functional_get_tau(functl, np, rho, rho_grad, rho_lapl, tau)
    !-----------------------------------------------------------------------!
    ! Computes the approximated kinetic energy density.                     !
    !                                                                       !
    !  functl   - functional                                                !
    !  m        - mesh                                                      !
    !  rho      - electronic radial density                                 !
    !  rho_grad - gradient of the electronic radial density                 !
    !  rho_lapl - laplacian of the electronic radial density                !
    !  tau      - radial kinetic energy density                             !
    !-----------------------------------------------------------------------!
    type(xc_functl_t), intent(in)  :: functl
    integer     ,       intent(in)  :: np
    real(8),           intent(in)  :: rho(np, functl%nspin)
    real(8),           intent(in)  :: rho_grad(np, functl%nspin)
    real(8),           intent(in)  :: rho_lapl(np, functl%nspin)
    real(8),           intent(out) :: tau(np, functl%nspin)

    integer  :: i, is, nspin
    real(8), allocatable :: n(:), s(:), l(:), t(:)

    nspin = functl%nspin

    !Allocate work arrays
    allocate(n(nspin), t(nspin))
    n =0.0d0; t =0.0d0
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
      if (nspin == 1) then
        allocate(s(1))
      else
        allocate(s(3))
      end if
      s =0.0d0
    end if
    if (functl%family == XC_FAMILY_MGGA) then
      allocate(l(nspin))
      l =0.0d0
    end if

    !Spin loop
    do is = 1, nspin
      !Space loop
      do i = 1, np
        ! make a local copy with the correct memory order
        n(is) = rho(i, is)

        if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
          s(1) = rho_grad(i, 1)**2
          if(nspin == 2) then
            s(2) = rho_grad(i, 1)*rho_grad(i, 2)
            s(3) = rho_grad(i, 2)**2
          end if
        end if
        if (functl%family == XC_FAMILY_MGGA) then
          l(is) = rho_lapl(i, is)
        end if

        select case(functl%family)
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(functl%conf, xc_one, n(1), t(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(functl%conf, xc_one, n(1), s(1), t(1))
        end select

        tau(i, is) = 2.0d0*t(1)*n(is)
      end do
    end do

    !Deallocate arrays
    deallocate(n, t)
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) deallocate(s)
    if (functl%family == XC_FAMILY_MGGA) deallocate(l)

  end subroutine functional_get_tau

end module functionals_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
