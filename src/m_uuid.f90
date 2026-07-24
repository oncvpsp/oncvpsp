module mtprng

!---------------------------------------------------------------------
! From the Algorithmic Conjurings of Scott Robert Ladd comes...
!---------------------------------------------------------------------
!
!  mtprng.f90 (a Fortran 95 module)
!
!  An implementation of the Mersenne Twister algorithm for generating
!  psuedo-random sequences.
!
!  History
!  -------
!   1.0.0   Initial release
!
!   1.1.0   6 February 2002
!           Updated to support algorithm revisions posted
!           by Matsumoto and Nishimura on 26 January 2002
!
!   1.5.0   12 December 2003
!           Added to hypatia project
!           Minor style changes
!           Tightened code
!           Now state based; no static variables
!           Removed mtprng_rand_real53
!
!   2.0.0   4 January 2004
!           Corrected erroneous unsigned bit manipulations
!           Doubled resolution by using 64-bit math
!           Added mtprng_rand64

!  Version for distribution with FoX <http://uszla.me.uk/FoX>
!  Very small cosmetic changes to fit FoX naming scheme and
!  avoid additional dependencies.
!  Toby White <tow@uszla.me.uk>, 2007

!
!  ORIGINAL ALGORITHM COPYRIGHT
!  ============================
!  Copyright (C) 1997,2002 Makoto Matsumoto and Takuji Nishimura.
!  Any feedback is very welcome. For any question, comments, see
!  http://www.math.keio.ac.jp/matumoto/emt.html or email
!  matumoto@math.keio.ac.jp
!---------------------------------------------------------------------
!
!  COPYRIGHT NOTICE, DISCLAIMER, and LICENSE:
!
!  This notice applies *only* to this specific expression of this
!  algorithm, and does not imply ownership or invention of the
!  implemented algorithm.
!  
!  If you modify this file, you may insert additional notices
!  immediately following this sentence.
!  
!  Copyright 2001, 2002, 2004 Scott Robert Ladd.
!  All rights reserved, except as noted herein.
!
!  This computer program source file is supplied "AS IS". Scott Robert
!  Ladd (hereinafter referred to as "Author") disclaims all warranties,
!  expressed or implied, including, without limitation, the warranties
!  of merchantability and of fitness for any purpose. The Author
!  assumes no liability for direct, indirect, incidental, special,
!  exemplary, or consequential damages, which may result from the use
!  of this software, even if advised of the possibility of such damage.
!  
!  The Author hereby grants anyone permission to use, copy, modify, and
!  distribute this source code, or portions hereof, for any purpose,
!  without fee, subject to the following restrictions:
!  
!      1. The origin of this source code must not be misrepresented.
!  
!      2. Altered versions must be plainly marked as such and must not
!         be misrepresented as being the original source.
!  
!      3. This Copyright notice may not be removed or altered from any
!         source or altered source distribution.
!  
!  The Author specifically permits (without fee) and encourages the use
!  of this source code for entertainment, education, or decoration. If
!  you use this source code in a product, acknowledgment is not required
!  but would be appreciated.
!  
!  Acknowledgement:
!      This license is based on the wonderful simple license that
!      accompanies libpng.
!
!-----------------------------------------------------------------------
!
!  For more information on this software package, please visit
!  Scott's web site, Coyote Gulch Productions, at:
!
!      http://www.coyotegulch.com
!  
!-----------------------------------------------------------------------

    implicit none

    ! Kind types for 64-, 32-, 16-, and 8-bit signed integers
    integer, parameter :: INT64 = selected_int_kind(18)
    integer, parameter :: INT32 = selected_int_kind(9)
    integer, parameter :: INT16 = selected_int_kind(4)
    integer, parameter :: INT08 = selected_int_kind(2)

    ! Kind types for IEEE 754/IEC 60559 single- and double-precision reals
    integer, parameter :: IEEE32 = selected_real_kind(  6,  37 )
    integer, parameter :: IEEE64 = selected_real_kind( 15, 307 )

    !------------------------------------------------------------------------------
    ! Everything is private unless explicitly made public
    private

    public :: mtprng_state, &
              mtprng_init, mtprng_init_by_array, &
              mtprng_rand64, mtprng_rand, mtprng_rand_range, &
              mtprng_rand_real1, mtprng_rand_real2, mtprng_rand_real3

    !------------------------------------------------------------------------------
    ! Constants
    integer(INT32), parameter :: N = 624_INT32
    integer(INT32), parameter :: M = 397_INT32

    !------------------------------------------------------------------------------
    ! types
    type mtprng_state
        integer(INT32)                   :: mti = -1
        integer(INT64), dimension(0:N-1) :: mt
    end type 

contains
    !--------------------------------------------------------------------------
    !  Initializes the generator with "seed"
    subroutine mtprng_init(seed, state)
    
        ! arguments
        integer(INT32),     intent(in)  :: seed
        type(mtprng_state), intent(out) :: state
        
        ! working storage
        integer :: i

        ! save seed        
        state%mt(0) = seed
        
        ! Set the seed using values suggested by Matsumoto & Nishimura, using
        !   a generator by Knuth. See original source for details.
        do i = 1, N - 1
            state%mt(i) = iand(4294967295_INT64,1812433253_INT64 * ieor(state%mt(i-1),ishft(state%mt(i-1),-30_INT64)) + i)
        end do
        
        state%mti = N

    end subroutine mtprng_init
    
    !--------------------------------------------------------------------------
    ! Initialize with an array of seeds
    subroutine mtprng_init_by_array(init_key, state)
    
        ! arguments
        integer(INT32), dimension(:), intent(in) :: init_key
        type(mtprng_state), intent(out) :: state
        
        ! working storage
        integer :: key_length
        integer :: i
        integer :: j
        integer :: k
        
        call mtprng_init(19650218_INT32,state)
        
        i = 1
        j = 0
        key_length = size(init_key)
        
        do k = max(N,key_length), 0, -1
            state%mt(i) = ieor(state%mt(i),(ieor(state%mt(i-1),ishft(state%mt(i-1),-30_INT64) * 1664525_INT64))) + init_key(j) + j
            
            i = i + 1
            j = j + 1
            
            if (i >= N) then
                state%mt(0) = state%mt(N-1)
                i = 1
            end if
            
            if (j >= key_length) j = 0
        end do
        
        do k = N-1, 0, -1
            state%mt(i) = ieor(state%mt(i),(ieor(state%mt(i-1),ishft(state%mt(i-1),-30_INT64) * 1566083941_INT64))) - i
            
            i = i + 1
            
            if (i>=N) then
                state%mt(0) = state%mt(N-1)
                i = 1
            end if
        end do

        state%mt(0) = 1073741824_INT64 ! 0x40000000, assuring non-zero initial array 
        
    end subroutine mtprng_init_by_array
    
    !--------------------------------------------------------------------------
    !   Obtain the next 32-bit integer in the psuedo-random sequence
    function mtprng_rand64(state) result(r)
    
        ! arguments
        type(mtprng_state), intent(inout) :: state
    
        !return type
        integer(INT64) :: r

        ! internal constants
        integer(INT64), dimension(0:1), parameter :: mag01 = (/ 0_INT64, -1727483681_INT64 /)

        ! Period parameters
        integer(INT64), parameter :: UPPER_MASK =  2147483648_INT64
        integer(INT64), parameter :: LOWER_MASK =  2147483647_INT64

        ! Tempering parameters
        integer(INT64), parameter :: TEMPERING_B = -1658038656_INT64
        integer(INT64), parameter :: TEMPERING_C =  -272236544_INT64
        
        ! Note: variable names match those in original example
        integer(INT32) :: kk
        
        ! Generate N words at a time
        if (state%mti >= N) then
            ! The value -1 acts as a flag saying that the seed has not been set.
            if (state%mti == -1) call mtprng_init(4357_INT32,state)
            
            ! Fill the mt array
            do kk = 0, N - M - 1
                r = ior(iand(state%mt(kk),UPPER_MASK),iand(state%mt(kk+1),LOWER_MASK))
                state%mt(kk) = ieor(ieor(state%mt(kk + M),ishft(r,-1_INT64)),mag01(iand(r,1_INT64)))
            end do
            
            do kk = N - M, N - 2
                r = ior(iand(state%mt(kk),UPPER_MASK),iand(state%mt(kk+1),LOWER_MASK))
                state%mt(kk) = ieor(ieor(state%mt(kk + (M - N)),ishft(r,-1_INT64)),mag01(iand(r,1_INT64)))
            end do
            
            r = ior(iand(state%mt(N-1),UPPER_MASK),iand(state%mt(0),LOWER_MASK))
            state%mt(N-1) = ieor(ieor(state%mt(M-1),ishft(r,-1)),mag01(iand(r,1_INT64)))
            
            ! Start using the array from first element
            state%mti = 0
        end if
        
        ! Here is where we actually calculate the number with a series of
        !   transformations 
        r = state%mt(state%mti)
        state%mti = state%mti + 1
        
        r = ieor(r,ishft(r,-11))
        r = iand(4294967295_INT64,ieor(r,iand(ishft(r, 7),TEMPERING_B)))
        r = iand(4294967295_INT64,ieor(r,iand(ishft(r,15),TEMPERING_C)))
        r = ieor(r,ishft(r,-18))
        
    end function mtprng_rand64
    
    !--------------------------------------------------------------------------
    !   Obtain the next 32-bit integer in the psuedo-random sequence
    function mtprng_rand(state) result(r)
    
        ! arguments
        type(mtprng_state), intent(inout) :: state
    
        !return type
        integer(INT32) :: r
        
        ! working storage
        integer(INT64) :: x
        
        ! done
        x = mtprng_rand64(state)
        
        if (x > 2147483647_INT64) then
            r = x - 4294967296_INT64
        else
            r = x
        end if
        
    end function mtprng_rand
    
    !---------------------------------------------------------------------------
    !   Obtain a psuedorandom integer in the range [lo,hi]
    function mtprng_rand_range(state, lo, hi) result(r)

        ! arguments
        type(mtprng_state), intent(inout) :: state
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        
        ! return type
        integer(INT32) :: r
        
        ! Use real value to caluclate range
        r = lo + floor((hi - lo + 1.0_IEEE64) * mtprng_rand_real2(state))
        
    end function mtprng_rand_range

    !--------------------------------------------------------------------------
    !   Obtain a psuedorandom real number in the range [0,1], i.e., a number
    !   greater than or equal to 0 and less than or equal to 1.
    function mtprng_rand_real1(state) result(r)

        ! arguments
        type(mtprng_state), intent(inout) :: state
    
        ! return type
        real(IEEE64) :: r
        
        ! Local constant; precalculated to avoid division below
        real(IEEE64), parameter :: factor = 1.0_IEEE64 / 4294967295.0_IEEE64
        
        ! compute
        r = real(mtprng_rand64(state),IEEE64) * factor
        
    end function mtprng_rand_real1

    !--------------------------------------------------------------------------
    !   Obtain a psuedorandom real number in the range [0,1), i.e., a number
    !   greater than or equal to 0 and less than 1.
    function mtprng_rand_real2(state) result(r)

        ! arguments
        type(mtprng_state), intent(inout) :: state
    
        ! return type
        real(IEEE64) :: r
        
        ! Local constant; precalculated to avoid division below
        real(IEEE64), parameter :: factor = 1.0_IEEE64 / 4294967296.0_IEEE64
        
        ! compute
        r = real(mtprng_rand64(state),IEEE64) * factor
        
    end function mtprng_rand_real2

    !--------------------------------------------------------------------------
    !   Obtain a psuedorandom real number in the range (0,1), i.e., a number
    !   greater than 0 and less than 1.
    function mtprng_rand_real3(state) result(r)

        ! arguments
        type(mtprng_state), intent(inout) :: state
    
        ! return type
        real(IEEE64) :: r
        
        ! Local constant; precalculated to avoid division below
        real(IEEE64), parameter :: factor = 1.0_IEEE64 / 4294967296.0_IEEE64
        
        r = (real(mtprng_rand64(state),IEEE64) + 0.5_IEEE64) * factor
        
    end function mtprng_rand_real3

end module mtprng

module m_uuid
  ! Fortran module to generate UUIDs
  ! Extracted from FoX (work by Toby White and Andrew Walker)


  !This generates UUIDs according to RFC 4122

  ! Only types 1 (time-based) and 4 (pseudo-RNG-based) are implemented.

  use mtprng, only : mtprng_state, mtprng_init, mtprng_rand64

  implicit none
  private
  
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i8b = selected_int_kind(18)
  
  character, parameter :: hexdigits(0:15) = &
    (/'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'/)
  
  type(mtprng_state), save :: rng_state
  logical, save :: initialized = .false.
  integer, save :: values_save ! must be default for date_and_time
  integer(kind=i4b), save :: hires_count = 0

! clock-seq holds a random number constant for the lifetime of the program
! using this module. That's the best we can do per S 4.1.5
  integer, save :: clock_seq = 0

  public :: generate_uuid
  
contains
  
  function generate_uuid(version) result(uuid)
    integer, intent(in), optional :: version
    character(len=36) :: uuid

    integer(kind=i8b) :: timestamp, node
    integer(kind=i4b) :: clock_sequence

    integer(kind=i4b) :: time_low, time_mid, time_hi_and_version
    integer(kind=i4b) :: clk_seq_hi_res, clk_seq_low

    integer :: values(8) ! must be default for date_and_time
    integer(kind=i4b) :: variant, v


    if (.not.initialized) then
      ! Use the current date and time to init mtprng
      ! but this gives limited varaibility, so mix 
      ! the result up.  Can we do better? In any
      ! case, this gets passed through a quick 
      ! generator inside mtprng_init.
      call date_and_time(values=values)
      values(7) = values(7)*1000+values(5)*100+values(3)*10+values(1)
      values(8) = values(2)*1000+values(4)*100+values(6)*10+values(8)
      call mtprng_init(int(values(7)*10000+values(8), i4b), rng_state)
      clock_seq = int(mtprng_rand64(rng_state), i4b)
      initialized = .true.
    endif

    variant = 1

    if (present(version)) then
      v = version
    else
      v = 4
    endif

    select case (v)
    case (0)
      ! Nil UUID  - S 4.1.7
      uuid = repeat('0',8)//'-'//repeat('0',4)//'-'//repeat('0',4)// &
        '-'//repeat('0',4)//'-'//repeat('0',12)
      return
    case(1)
      call date_and_time(values=values)
      ! In case of too-frequent requests, we will replace time_low
      ! with the count below ...
      if (all(values==values_save)) then
        hires_count = hires_count + 1
      else
        hires_count = 0
      endif
    case(2-3)
      !Unimplemented
      uuid = ''
      return
    case(4)
      continue
    case(5)
      !Unimplemented
      uuid = ''
      return
    case default
      !Unspecified
      uuid = ''
      return
    end select

!4.1.4 Timestamp

    select case(v)
    case(1)
      timestamp = get_utc_since_1582(values)
    case(4)
      timestamp = ior(mtprng_rand64(rng_state), ishft(mtprng_rand64(rng_state), 28))
    end select

!4.1.5 Clock Sequence
    ! 14 bits
    select case(v)
    case(1)
      clock_sequence = clock_seq
    case(4)
      clock_sequence = int(mtprng_rand64(rng_state), i4b)
    end select

!4.1.6 Node
    ! 48 bits
    select case(v)
    case(1)
      node = ior(mtprng_rand64(rng_state), ishft(mtprng_rand64(rng_state), 16))
      ! No MAC address accessible - see section 4.5 !FIXME
    case(4)
      node = ior(mtprng_rand64(rng_state), ishft(mtprng_rand64(rng_state), 16))
    end select

    time_low = ibits(timestamp, 0, 32)
    time_mid = ibits(timestamp, 32, 16)
    if (hires_count==0) then
      time_hi_and_version = ior(int(ibits(timestamp, 48, 12), i4b), ishft(v, 12))
    else
      time_hi_and_version = ior(hires_count, ishft(v, 12))
    endif

    clk_seq_low = ibits(clock_sequence, 0, 8)
    clk_seq_hi_res = ior(ibits(clock_sequence, 8, 6), ishft(variant, 6))

    uuid = int32ToHexOctets(time_low, 4)//"-"// &
      int32ToHexOctets(time_mid, 2)//"-"// &
      int32ToHexOctets(time_hi_and_version, 2)//"-"// & 
      int32ToHexOctets(clk_seq_hi_res, 1)// &
      int32ToHexOctets(clk_seq_low, 1)//"-"// &
      int64ToHexOctets(node, 6)

  contains

    function int32ToHexOctets(b, n) result(s)
      integer(i4b), intent(in) :: b
      integer, intent(in) :: n ! number of octets to print
      character(len=2*n) :: s
      
      integer :: i
      
      do i = 0, 2*n-1
        s(2*n-i:2*n-i) = hexdigits(ibits(b, i*4, 4))
      enddo
      
    end function int32ToHexOctets
    function int64ToHexOctets(b, n) result(s)
      integer(i8b), intent(in) :: b
      integer, intent(in) :: n ! number of octets to print
      character(len=2*n) :: s
      
      integer :: i
      
      do i = 0, 2*n-1
        s(2*n-i:2*n-i) = hexdigits(ibits(b, i*4, 4))
      enddo
      
    end function int64ToHexOctets

  end function generate_uuid

  function get_utc_since_1582(values) result(ns)
    ! This subroutine is a little broken. It only works
    ! for times after 1/1/2006 and takes no account
    ! of any future leapseconds. It ought to serve regardless.

    ! It returns the number of 100-ns intervals since 1582-10-15-00-00-00

    integer, dimension(8), intent(in) :: values
    integer(kind=i8b) :: ns

    integer :: days
    integer :: years

    integer, parameter :: days_in_normal_year(12) = &
      (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    ns = 23_i8b * 1000_i8b * 1000_i8b * 10_i8b ! 23 leap seconds until 24:00:00 31/12/2005

    ! A count of the 100-nanosecond intervals since the
    ! beginning of the day.
    ns = ns &
      ! milliseconds
      + int(values(8), i8b)             * 10_i8b * 1000_i8b &
      ! seconds
      + int(values(7), i8b)             * 10_i8b * 1000_i8b * 1000_i8b &
      ! minutes (with timezone adjustment)
      + int(values(6) + values(4), i8b) * 10_i8b * 1000_i8b * 1000_i8b * 60_i8b &
      ! hours
      + int(values(5), i8b)             * 10_i8b * 1000_i8b * 1000_i8b * 60_i8b * 60_i8b

    ! Number of days this year:
    days = sum(days_in_normal_year(:values(2)-1))
    days = days + values(3) - 1 !add days in current month
    if (values(2)>2 .and. isLeapYear(values(1))) then
      days = days + 1
    endif
    !That's all the time since the turn of this year

    days = days + 78 ! From the start of 15th Oct to the end of 31st Dec in 1582
    !That's the additional time before the turn of the year 1583

    days = days + 102  ! 102 leap years from 1584 to 2000 inclusive
    ! That's all the intercalataed days until 2000
    
    years = values(1) - 2000 - 1 ! years since 2000 - not including this year

    days = days + years/4 - years/100 + years/400 !Add extra leap days to this total:
    ! That's all out intercalated days - remaining years are all 365 days long.

    years = years + 418 ! Add the years from 1583-2000 inclusive back on.

    ! Multiply by number of time units in one day & add to today's total.
    ns = ns + 864000000000_i8b * (int(days,i8b) + 365_i8b * int(years,i8b))

  contains
    function isLeapYear(y) result(p)
      integer, intent(in) :: y
      logical :: p
      p = (mod(y,4)==0 .and. .not.mod(y,100)==0 .or. mod(y,400)==0)
    end function isLeapYear

  end function get_utc_since_1582

end module m_uuid

subroutine get_uuid(uuid)
  use m_uuid, only: generate_uuid

    character(len=36), intent(out) :: uuid

    uuid = generate_uuid(version=1)
end subroutine get_uuid
