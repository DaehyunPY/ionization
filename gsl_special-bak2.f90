module gsl_special
    use iso_c_binding
    implicit none
    
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
