!
! To the best of my knowledge, no restrictions have been
! placed on the use and re-distribution of this module.
!
module m_getopts
! Assume now a f2003 compliant compiler...
!use f2kcli
public :: getopts

CONTAINS

subroutine getopts(optionstring,name,optarg,optind,exitcode)
!$$$  Subprogram Documentation Block
!
! Subprogram:  Getopts    Process command line arguments for valid options
!   Prgmmr: Iredell       Org: W/NP23        Date: 2000-08-22
!
! Abstract: This subprogram processes command-line arguments for valid options.
!           It is the Fortran equivalent of the built-in shell command getopts.
!           Options on the command line come before the positional arguments.
!           Options are preceded by a - (minus sign) or a + (plus sign).
!           Options are single case-sensitive alphanumeric characters.
!           Options either do or do not have an expected argument.
!           Options without an argument may be immediately succeeded by
!           further options without the accompanying - or + prefix.
!           Options with an argument may be separated from their argument
!           by zero or more blanks.  The argument cannot end with a blank.
!           Options end when not preceded by a - or a + or after -- or ++.
!           This subprogram processes one option per invocation.
!           This subprogram is not thread-safe.
!
! Program History Log:
!   2007-05-05  Alberto Garcia: Use f2k-compliant f2kcli module
!   2000-08-22  Iredell
!
! Usage:    call getopts(optionstring,name,optarg,optind,exitcode)
!
!   Input Argument List:
!     optionstring
!       character string containing a list of all valid options;
!       options succeeded by a : require an argument
!
!   Input and Output Argument List:
!     optind
!       integer index of the next argument to be processed;
!       set to 0 before initial call or to restart processing
!
!   Output Argument List:
!     name
!       character string containing the name of the next option
!       or ? if no option or an unknown option is found
!       or : if an option had a missing required argument;
!       a + is prepended to the value in name if the option begins with a +
!     optarg
!       character string containing the option argument if required
!       or null if not required or not found;
!       optarg contains the option found if name is ? or :.
!     exitcode
!       integer return code (0 if an option was found, 1 if end of options)
!     
! Subprograms Called:
!   command_argument_count
!     Retrieve number of command-line arguments
!   get_command_argument
!     Retrieve a command-line argument
!   index 
!     Retrieve the starting position of a substring within a string
!     
! Remarks:
!   Here is an example of how to use this subprogram.
!     -----------------------------------------------------------
!     implicit none
!     character(len=20) copt,carg,cb,cpos
!     integer ia,ib,iopt,iret,narg,npos,ipos
!
!     ia=0     ! Programmer's flag for option a
!     ib=0     !  "                           b
!
!     iopt=0
!     do
!       call getopts('ab:',copt,carg,iopt,iret)
!       if(iret.ne.0) exit
!       select case(copt)
!       case('a','+a')
!         ia=1
!       case('b','+b')
!         ib=1
!         cb=carg
!       case('?',':')
!         print *,'invalid option ',carg(1:1)
!         stop 1
!       end select
!     enddo
!     if(ia.eq.1) print *,'option a selected'
!     if(ib.eq.1) print *,'option b selected; argument=',cb
!
!     Now process positional parameters
!     narg=command_argument_count()
!     npos=narg-iopt+1
!     do ipos=1,npos
!       call get_command_argument(number=ipos+iopt-1,value=cpos)
!       print *,'positional argument ',ipos,' is ',cpos
!     enddo
!     end
!     -----------------------------------------------------------
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none

!  Passed data
  character(len=*),intent(in)      :: optionstring
  character(len=*),intent(out)     :: name
  character(len=*),intent(out)     :: optarg
  integer,         intent(inout)   :: optind
  integer,         intent(out)     :: exitcode

!  Saved data
  character(len=256), save  :: carg
  character(len=1) ,  save  :: cone
  integer,save              :: narg,larg,lcur

!  Local data
  character(len=1) ::  copt
  integer          ::  lname,lopt

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Initially set saved data.
  if(optind.le.0) then
    optind=0
    narg=command_argument_count()
    carg=''
    cone=''
    larg=0
    lcur=1
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Retrieve next command-line argument if necessary;
!  exit if at end of options
  if(lcur.gt.larg) then
    optind=optind+1
    if(optind.gt.narg) then
      name='?'
      optarg=''
      exitcode=1
      RETURN
    endif
    call get_command_argument(number=optind,value=carg)
    cone=carg(1:1)
    larg=len_trim(carg)
    lcur=2
    if(larg.eq.1.or.(cone.ne.'-'.and.cone.ne.'+')) then
      name='?'
      optarg=''
      exitcode=1
      RETURN
    elseif(larg.eq.2.and.carg(2:2).eq.cone) then
      optind=optind+1
      name='?'
      optarg=''
      exitcode=1
      RETURN
    endif
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Find next option in the list; exit if option is unknown
  exitcode=0
  copt=carg(lcur:lcur)
  lcur=lcur+1
  lopt=index(optionstring,copt)
  if(lopt.eq.0) then
    name='?'
    optarg=copt
    RETURN
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Option found; retrieve its argument if requested
    if(cone.eq.'-') then
      name=""
      lname=1
    else
      name="+"
      lname=2
    endif
  name(lname:lname)=copt
  optarg=''
  if(lopt.lt.len(optionstring)) then
     if (optionstring(lopt+1:lopt+1).eq.':') then
        if(lcur.gt.larg) then
           optind=optind+1
           if(optind.gt.narg) then
              name=':'
              optarg=copt
              RETURN
           endif
           call get_command_argument(number=optind,value=carg)
           larg=len_trim(carg)
           lcur=1
        endif
        optarg=carg(lcur:larg)
        lcur=larg+1
     endif
  endif
end subroutine
end module m_getopts
