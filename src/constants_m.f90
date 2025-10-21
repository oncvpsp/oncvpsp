module constants_m
   implicit none
   private
   integer, parameter, public :: MAX_NUM_PROJ = 5
      !! Maximum number of projectors per angular momentum
   integer, parameter, public :: MAX_NUM_STATE = 30
      !! Maximum number of states in a configuration
   integer, parameter, public :: MAX_NUM_TEST = 5
      !! Maximum number of test configurations
   integer, parameter, public :: MAX_NUM_ELL = 6
      !! Maximum number of angular momentum channels
   integer, parameter, public :: LLOC_POLY_EXTRAP = 4
      !! Special value of lloc indicating construction via polynomial
      !! extrapolation of the all-electron potential to r=0.
end module constants_m
