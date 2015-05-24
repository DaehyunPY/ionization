module gsl_special
    use iso_c_binding
    implicit none
    
!
! Kind and length parameters are default integer
!
  integer, parameter, public :: fgsl_double = c_double
  integer, parameter, public :: fgsl_double_complex = c_double_complex
!  integer, parameter, public :: fgsl_extended = selected_real_kind(18)
  integer, parameter, public :: fgsl_extended = selected_real_kind(13)
! FIXME - c_long_double unsupported, selected_real_kind(30) unsupported in g95
  integer, parameter, public :: fgsl_float = c_float
  integer, parameter, public :: fgsl_int = c_int
  integer, parameter, public :: fgsl_long = c_long
  integer, parameter, public :: fgsl_size_t = c_size_t
  integer, parameter, public :: fgsl_char = c_char
  integer, parameter, public :: fgsl_strmax = 128
  integer, parameter, public :: fgsl_pathmax = 2048

!
! Types: Special Functions
!
  type, public :: fgsl_sf_result
     real(fgsl_double) :: val, err
  end type fgsl_sf_result
! FIXME ifort refuses = overload if not public
  type, public, bind(c) :: gsl_sf_result
     real(c_double) :: val, err
  end type
  type, public :: fgsl_sf_result_e10
     real(fgsl_double) :: val, err
     integer(fgsl_int) :: e10
  end type fgsl_sf_result_e10
! FIXME ifort refuses = overload if not public
  type, public, bind(c) :: gsl_sf_result_e10
     real(c_double) :: val, err
     integer(c_int) :: e10
  end type 
  type, public :: fgsl_mode_t
     private
     integer(c_int) :: gsl_mode = 0
  end type fgsl_mode_t
  type(fgsl_mode_t), parameter, public :: &
       fgsl_prec_double = fgsl_mode_t(0), &
       fgsl_prec_single = fgsl_mode_t(1), &
       fgsl_prec_approx = fgsl_mode_t(2)
    
    interface 
     function gsl_sf_bessel_jsl(n,x) bind(c, name='gsl_sf_bessel_jl')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_jsl
     end function gsl_sf_bessel_jsl
     function gsl_sf_bessel_jsl_array(lmax, x, result) bind(c, name='gsl_sf_bessel_jl_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jsl_array
     end function gsl_sf_bessel_jsl_array
     function gsl_sf_bessel_jsl_steed_array(lmax, x, result) bind(c, name='gsl_sf_bessel_jl_steed_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jsl_steed_array
     end function gsl_sf_bessel_jsl_steed_array


     function gsl_sf_bessel_ysl(n,x) bind(c, name='gsl_sf_bessel_yl')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ysl
     end function gsl_sf_bessel_ysl
     function gsl_sf_bessel_ysl_array(lmax, x, result) bind(c, name='gsl_sf_bessel_yl_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ysl_array
     end function gsl_sf_bessel_ysl_array


     function gsl_sf_bessel_isl_scaled(n,x) bind(c, name='gsl_sf_bessel_il_scaled')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_isl_scaled
     end function gsl_sf_bessel_isl_scaled
     function gsl_sf_bessel_isl_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_il_scaled_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_isl_scaled_e
     end function gsl_sf_bessel_isl_scaled_e
     function gsl_sf_bessel_isl_scaled_array(lmax, x, result) bind(c, name='gsl_sf_bessel_il_scaled_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_isl_scaled_array
     end function gsl_sf_bessel_isl_scaled_array


     function gsl_sf_bessel_ksl_scaled(n,x) bind(c, name='gsl_sf_bessel_kl_scaled')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ksl_scaled
     end function gsl_sf_bessel_ksl_scaled
     function gsl_sf_bessel_ksl_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_kl_scaled_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ksl_scaled_e
     end function gsl_sf_bessel_ksl_scaled_e
     function gsl_sf_bessel_ksl_scaled_array(lmax, x, result) bind(c, name='gsl_sf_bessel_kl_scaled_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ksl_scaled_array
     end function gsl_sf_bessel_ksl_scaled_array


     function gsl_sf_legendre_pl(l, x) bind(c, name='gsl_sf_legendre_Pl')
       import
       integer(c_int), value :: l
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_pl
     end function gsl_sf_legendre_pl
     function gsl_sf_legendre_pl_array(lmax, x, res_arr) bind(c, name='gsl_sf_legendre_Pl_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr
       integer(c_int) :: gsl_sf_legendre_pl_array
     end function gsl_sf_legendre_pl_array
     function gsl_sf_legendre_pl_deriv_array(lmax, x, res_arr, der_arr) &
          bind(c, name='gsl_sf_legendre_Pl_deriv_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr, der_arr
       integer(c_int) :: gsl_sf_legendre_pl_deriv_array
     end function gsl_sf_legendre_pl_deriv_array
   end interface 
end module gsl_special





! program test
!     use gsl_special
!     implicit none
!     integer :: callc
!     real(8) :: result(5)
    
!     write(*, *) gsl_sf_bessel_jsl(1, 1.d0)
!     write(*, *)
!     callc = gsl_sf_bessel_jsl_array(5, 1.d0, result)
!     write(*, *) callc, result 
!     write(*, *)
!     callc = gsl_sf_bessel_jsl_steed_array(5, 1.d0, result)
!     write(*, *) callc, result 
!     write(*, *)
! end program test
