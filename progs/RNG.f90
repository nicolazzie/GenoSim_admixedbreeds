MODULE RNG
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! HISTORY                                                                                                                  
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Original program: Hossein Jorjani, 2009-2011                                                                             
! Updates:                                                                                                                 
!                                                                                                                          
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
IMPLICIT NONE

CONTAINS


! CAN YOU CHANGE RAN TO RAN1, RAN1 TO RAN2? ####################################
!--------1---------2---------3---------4---------5---------6---------7---------8
!------------------------------------- Function ran1 ---------------------------
!--------1---------2---------3---------4---------5---------6---------7---------8

! I CHANGED ran to ran1.   Hossein Jorjani   2009 08 04

FUNCTION ran1(idum)
! Uniform random number generator
! Copied from Numerical Recipes in Fortran 90, Second edition, Volume 2 of
! Fortran Numerical Recipes, Chapter B7, p 1142
! Calling: x=ran(idum)
! Call with idum a negative integer to initialize (See below)
IMPLICIT NONE

INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL :: ran1

! "Minimal" random number generator of Park and Miller combined with a
! Marsaglia shift sequence. Returns a uniform random deviate between 0.0 and
! 1.0 (exclusive of the endpoint values). This fully portable, scalar generator
! has the "traditional" (not Fortran 90) calling sequence with a random deviate
! as the returned function value: call with idum a negative integer to
! initialize; thereafter, do not alter idum except to reinitialize. The period
! of this generator is about 3.1 × 10e18.

INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k

IF (idum <= 0 .or. iy < 0) THEN           ! Initialize
  am=nearest(1.0,-1.0)/IM
  iy=ior(ieor(888889999,abs(idum)),1)
  ix=ieor(777755555,abs(idum))
  idum=abs(idum)+1                        ! Set idum positive.
END IF

ix=ieor(ix,ishft(ix,13))        ! Marsaglia shift sequence with period 2e32-1.
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ
iy=IA*(iy-k*IQ)-IR*k            ! Park-Miller sequence by Schrage's method,
IF (iy < 0) iy=iy+IM            ! period 2e31-2.
ran1=am*ior(iand(IM,ieor(ix,iy)),1) ! Combines the two generators with masking
                                   ! to ensure nonzero values.
END FUNCTION ran1
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------------------------------- Function ran2 -------------------------------
!--------1---------2---------3---------4---------5---------6---------7---------8

! I CHANGED ran1 to ran2.   Hossein Jorjani   2009 08 04

FUNCTION ran2(idum)
INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
REAL ran2,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

! "Minimal" random number generator of Park and Miller with Bays-Durham shuffle
! and added safeguards. Returns a uniform random deviate between 0.0 and 1.0
! (exclusive of the endpoint values). Call with idum a negative integer to
! initialize; thereafter, do not alter idum between successive deviates in a
! sequence. RNMX should approximate the largest floating value that is less
! than 1.

INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
if (idum.le.0.or.iy.eq.0) THEN     ! Initialize.
  idum=max(-idum,1)                ! Be sure to prevent idum = 0.
  DO j=NTAB+8,1,-1                 ! Load the shuffle table (after 8 warm-ups).
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF (idum.LT.0) idum=idum+IM
    IF (j.LE.NTAB) iv(j)=idum
  END DO
  iy=iv(1)
END IF
k=idum/IQ                          ! Start here when not initializing.
idum=IA*(idum-k*IQ)-IR*k           ! Compute idum=mod(IA*idum,IM) without
IF (idum.LT.0) idum=idum+IM        ! overflows by Schrage's method.
j=1+iy/NDIV                        ! Will be in the range 1:NTAB.
iy=iv(j)                           ! Output previously stored value and refill
iv(j)=idum                         ! the shuffle table.
ran2=min(AM*iy,RNMX)               ! Because users don't expect endpoint values.
RETURN

END FUNCTION ran2
!--------1---------2---------3---------4---------5---------6---------7---------8
!------------------------------ Function gasdev --------------------------------
!--------1---------2---------3---------4---------5---------6---------7---------8

! I CHANGED ran1 to ran2.   Hossein Jorjani   2009 08 04

FUNCTION gasdev(idum)
! Normal random number generator
! Copied from Numerical Recipes in Fortran 77, The art of Scientific computing,
! Second edition, Volume 1 of Fortran Numerical Recipes, Chapter 7, p 280
! Calling: x=gasdev(idum)
! Call with idum a negative integer to initialize ?

INTEGER idum
REAL gasdev

! Uses ran2
! Returns a normally distributed deviate with zero mean and unit variance,
! using ran2(idum) as the source of uniform deviates.

INTEGER iset
REAL fac,gset,rsq,v1,v2
SAVE iset,gset
DATA iset/0/

IF (idum.LT.0) iset=0                ! Reinitialize.
IF (iset.eq.0) THEN                  ! We don't have an extra deviate handy, so
1  v1=2.*ran2(idum)-1.               ! pick two uniform numbers in the square
   v2=2.*ran2(idum)-1.               ! extending from -1 to +1 in each
   rsq=v1**2+v2**2                   ! direction, see if they are in the unit
   IF(rsq.GT.1..OR.rsq.EQ.0.) GOTO 1 ! circle, and if they are not, try again.
   fac=sqrt(-2.*log(rsq)/rsq)        ! Now make the Box-Muller transformation
   gset=v1*fac                       ! to get two normal deviates. Return one
   gasdev=v2*fac                     ! and save the other for next time.
   iset=1                            ! Set flag.
ELSE                                 ! We have an extra deviate handy,
   gasdev=gset                       ! so return it,
   iset=0                            ! and unset the flag.
END IF
RETURN

END FUNCTION gasdev
!--------1---------2---------3---------4---------5---------6---------7---------8
!-------------------------------------------------------------------------------
!--------1---------2---------3---------4---------5---------6---------7---------8






!--------1---------2---------3---------4---------5---------6---------7---------8     
!------------------------------ Function gamma --------------------------------                                   
!--------1---------2---------3---------4---------5---------6---------7---------8     
FUNCTION random_gamma(s,first) RESULT(fn_val)
!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).
!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).

REAL, INTENT(IN)    :: s
LOGICAL,INTENT(IN)  :: first
REAL                :: fn_val

IF (s <= 0.D0) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s > 1.D0) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < 1.D0) THEN
  fn_val = random_gamma2(s, first)
ELSE
  fn_val = random_exponential()
END IF

RETURN
END FUNCTION random_gamma


FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.
USE Prog_Param, ONLY : idum1
REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val, half

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (first) THEN
  d = s - 1./3.
  c = 1./SQRT(9.0*d)
END IF
half = 0.5D+00

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.


  DO
    x = gasdev(idum1)!random_normal()
    v = (1. + c*x)**3
    IF (v > 0.) EXIT
  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < 1. - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(1. - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1



FUNCTION random_gamma2(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

!     Local variables
REAL       :: r, x, w, vsmall
REAL, SAVE :: a, p, c, uf, vr, d

vsmall = TINY(1.0)

IF (s <= 0. .OR. s >= 1.) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = 1. - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = 1./s
  uf = p*(vsmall/a)**s
  vr = 1. - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((1. - r)/(1. - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = 0.
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (1.-r <= w .AND. r > 0.) THEN
    IF (r*(w + 1.) >= 1.) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2




FUNCTION random_exponential() RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

REAL  :: fn_val

!     Local variable
REAL  :: r

DO
  CALL RANDOM_NUMBER(r)
  IF (r > 0.) EXIT
END DO

fn_val = -LOG(r)
RETURN

END FUNCTION random_exponential

subroutine init_random_seed()
  ! ----- variables for portable seed setting -----
  INTEGER :: i_seed
  INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
  INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----
  REAL :: r

  ! ----- Set up random seed portably -----
  CALL RANDOM_SEED(size=i_seed)
  ALLOCATE(a_seed(1:i_seed))
  CALL RANDOM_SEED(get=a_seed)
  CALL DATE_AND_TIME(values=dt_seed)
  a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
  CALL RANDOM_SEED(put=a_seed)
  DEALLOCATE(a_seed)
  ! ----- Done setting up random seed -----

  CALL RANDOM_NUMBER(r)
!  WRITE(6,*) 'random number is ',r
end subroutine init_random_seed


!subroutine init_random_seed()
!  use iso_fortran_env, only: int64
!  use Numeric_kinds
!  implicit none

!  integer, allocatable :: seed(:)
!  integer :: i, n, un, istat, dt(8), pid
!  integer(i4) :: t

!  call random_seed(size = n)
!  allocate(seed(n))
!  ! First try if the OS provides a random number generator
!  open(4134, file="here", access="stream", &
!       form="unformatted", action="read", status="old", iostat=istat)
!  if (istat == 0) then
 !    read(un) seed
 !    close(un)
 ! else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
 !    call system_clock(t)
 !    if (t == 0) then
 !       call date_and_time(values=dt)
 !       t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
 !            + dt(2) * 31 * 24 * 60 * 60 * 1000 &
 !            + dt(3) * 24 * 60 * 60 * 1000 &
 !            + dt(5) * 60 * 60 * 1000 &
 !            + dt(6) * 60 * 1000 + dt(7) * 1000 &
 !            + dt(8)
 !    end if
 !    pid = getpid()
 !    t = ieor(t, int(pid, kind(t)))
 !    do i = 1, n
 !       seed(i) = lcg(t)
 !    end do
 ! end if
 ! call random_seed(put=seed)
!contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
 ! function lcg(s)
 !   integer :: lcg
 !   integer(i4) :: s
 !   if (s == 0) then
 !      s = 104729
 !   else
 !      s = mod(s, 4297296)
 !   end if
 !   s = mod(s * 279470273, 429496791)
 !   lcg = int(mod(s, int(huge(0), i4)), kind(0))
 ! end function lcg
!end subroutine init_random_seed







END MODULE RNG
