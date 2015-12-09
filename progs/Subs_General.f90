MODULE Subs_General

USE Numeric_Kinds

IMPLICIT NONE

CONTAINS


! The three sort subroutines: sort_rec, sort_ascending, and sort_descending seem
! to be obsolete after inclusion of Ezequiel's sort subroutine 'sort_updo'.
! However, I will keep them active until further testing.


!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine sort_updo ----------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
! History:
!  - Original Program:                                 EL Nicolazzi - Sep   2008
!  - Update for sorting up/down                               - ELN - Nov   2009
!  - Adapted for Hossein's sim                                - ELN - Jun   2012
!  - Added "ord(j)" to allow non-conscutive arrays            - ELN - Jan   2013
!  - Modified CYCLE by EXIT to contemplate the HIGHLY improbable possibility of
!    2 random numbers equal, that would result in an animal 0 - ELN - Feb   2013
!--------1---------2---------3---------4---------5---------6---------7---------8

SUBROUTINE sort_updo (VAL, ORD, N, DIR)
!
! Description of parameters:
!   VAL - array (IN-OUT) of values to be sorted
!   ORD - array with positions (IN-OUT). In should be an array of 1 to N.
!    N  - length of the array (N anims to be sorted)
!   DIR - direction of sorting: 'high_to_low' or 'low_to_high'
!--------1---------2---------3---------4---------5---------6---------7---------8

  IMPLICIT NONE

  REAL(dp)                          :: matchz
  REAL(dp),INTENT(INOUT)            :: val(*)
  REAL(dp),ALLOCATABLE,DIMENSION(:) :: tempval, outval
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: tempord
  INTEGER,INTENT(INOUT)             :: ord(*)
  INTEGER,INTENT(IN)                :: N
  INTEGER                           :: i, j
  CHARACTER*11                      :: dir

!--------1---------2---------3---------4---------5---------6---------7---------8

ALLOCATE(tempval(N),tempord(N),outval(N))
outval=0.D0;matchz=0.D0
tempord=0
tempval=val(1:N)

IF(dir/='high_to_low' .AND. dir/='low_to_high')THEN
   PRINT*,'ERROR: Sort options are high_to_low - low_to_high';PRINT*,'EXITING SORT...'
   STOP
ENDIF

DO i=1,N
  IF (dir == 'high_to_low') matchz = MAXVAL(tempval)     
  IF (dir == 'low_to_high') matchz = MINVAL(tempval)     
  DO j=1,N
        IF (val(j) .EQ. matchz) THEN
           tempord(i) = ord(j)
           IF (dir == 'high_to_low') tempval(j) = -99999999.
           IF (dir == 'low_to_high') tempval(j) =  99999999.
           EXIT
        ENDIF
     ENDDO
  outval(i) = matchz
ENDDO

val(1:N) = outval
ord(1:N) = tempord

DEALLOCATE(tempval,tempord,outval)

!--------1---------2---------3---------4---------5---------6---------7---------8

END SUBROUTINE sort_updo

!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8




!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine sort_int_updo ------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
! History:
!  - Original Program:   EL Nicolazzi - Jan 2013
!  - Modified CYCLE by EXIT to contemplate the HIGHLY improbable possibility of
!    2 random numbers equal, that would result in an animal 0     (ELN Feb 2013)
!--------1---------2---------3---------4---------5---------6---------7---------8

SUBROUTINE sort_int_updo (VAL, ORD, N, DIR)
!
! Description of parameters:
!   VAL - array (IN-OUT) of values to be sorted (INTEGERS)
!   ORD - array with positions (IN-OUT). In should be an array of 1 to N.
!    N  - length of the array (N anims to be sorted)
!   DIR - direction of sorting: 'high_to_low' or 'low_to_high'
!--------1---------2---------3---------4---------5---------6---------7---------8

  IMPLICIT NONE

  INTEGER,INTENT(INOUT)             :: val(*), ord(*)
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: tempval, tempord, outval
  INTEGER,INTENT(IN)                :: N
  INTEGER                           :: i, j,  matchz
  CHARACTER*11                      :: dir

!--------1---------2---------3---------4---------5---------6---------7---------8

ALLOCATE(tempval(N),tempord(N),outval(N))
outval=0;matchz=0
tempord=0
tempval=val(1:N)

IF(dir/='high_to_low' .AND. dir/='low_to_high')THEN
   PRINT*,'ERROR: Sort options are high_to_low - low_to_high';PRINT*,'EXITING SORT...'
   STOP
ENDIF

DO i=1,N
  IF (dir == 'high_to_low') matchz = MAXVAL(tempval)     
  IF (dir == 'low_to_high') matchz = MINVAL(tempval)     
  DO j=1,N
        IF (val(j) .EQ. matchz) THEN
           tempord(i) = j
           IF (dir == 'high_to_low') tempval(j) = -99999999
           IF (dir == 'low_to_high') tempval(j) =  99999999
           EXIT
        ENDIF
     ENDDO
  outval(i) = matchz
ENDDO

val(1:N) = outval
ord(1:N) = tempord

DEALLOCATE(tempval,tempord,outval)

!--------1---------2---------3---------4---------5---------6---------7---------8

END SUBROUTINE sort_int_updo

!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8



!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine samp_WR       ------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
! History:
!  - Original Program:   EL Nicolazzi - Jan 2013
!--------1---------2---------3---------4---------5---------6---------7---------8
SUBROUTINE samp_WR (FIRST, LAST, NSAMP, VEC)
USE Numeric_Kinds
USE RNG
! This program samples without replacement from a field N=last-first, allowing to sample 
! Nsamp samples (lower or equal to N). Produces a vector VEC of dimension NSAMP.
! Description of parameters:
! FIRST - Starting value of sampling
!  LAST - Ending value of sampling
!   VEC - length of the array (N anims to be sorted)
!--------1---------2---------3---------4---------5---------6---------7---------8

  IMPLICIT NONE

  REAL(dp),ALLOCATABLE, DIMENSION(:) :: tempval,rte
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ite,disc
  INTEGER, INTENT(IN)                :: first,last, nsamp
  INTEGER, INTENT(OUT)               :: vec(nsamp)
  INTEGER                            :: i, j, cont, N, time(3)
  REAL                               :: thres, scart

!--------1---------2---------3---------4---------5---------6---------7---------8
N=1+last-first
call itime(time)
ALLOCATE(tempval(N),rte(Nsamp*2),ite(Nsamp*2))
vec=0;cont=0;rte=0.D0;rte=0.D0;ite=1

thres=Nsamp*1.5/(last-first)       ! To reduce computational burden, assuring to get all QTLs I want

DO i=1,time(3)               ! Not the most efficient method, but used to get real pseud-random QTL positions.
   call random_number(scart)
ENDDO

DO i=1,N
!   tempval(i)=ran1(idum1)
   call random_number(tempval(i))
   IF(tempval(i).LT.thres) THEN
      cont=cont+1;
      IF(cont.EQ.Nsamp*2)EXIT      ! To avoid segfaults
      ite(cont)=i
      rte(cont)=tempval(i)
   ENDIF
ENDDO

IF (cont.LT.Nsamp)THEN             ! Just to be sure I have all required QTLs (too much precaution?)
   thres=Nsamp*2/(last-first)      ! To reduce computational burden, assuring to get all QTLs I want
   DO i=1,N
      IF(tempval(i).LT.thres)THEN
         cont=cont+1
         IF(cont.EQ.Nsamp*2)EXIT
         ite(cont)=i;rte(cont)=tempval(i)
      ENDIF
   ENDDO
ENDIF
call sort_updo(rte,ite,cont,'high_to_low')    ! sorting a (reduced) random n to obtain a vec of positions 

ALLOCATE(disc(Nsamp))
FORALL(j=1:Nsamp)disc(j)=j
vec=ite(1:Nsamp)
call sort_int_updo(vec,disc,Nsamp,'low_to_high') ! re-sorting to use them in an easier way
DEALLOCATE(tempval)

!--------1---------2---------3---------4---------5---------6---------7---------8

END SUBROUTINE samp_WR

!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8


SUBROUTINE Cpoisson(mean,output,leng)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! History:
!!  - Original Program:   EL Nicolazzi - July      2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculates cumulative poisson distribution,creating a vector of values as output
!! P(x)=((h**x)*(e**-h))/x!, whereas h= expected recombination events (user defined)
!!                                   x= number of effects (used to create a distribution)
!! In this subroutine: one=(h**x),two=(e**-h),three=x! and output[x]=output[x-1]+one*two/three                      
!! Output should be an allocatable (non allocated) vector of real*8, leng should be an integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER     , PARAMETER    :: i4 = SELECTED_INT_KIND(9)           
INTEGER     , PARAMETER    :: dp  = SELECTED_REAL_KIND( 15,  307 )
INTEGER     , PARAMETER    :: largenum=999999 ! This should do for "real" simulations
REAL    (dp), ALLOCATABLE  :: temp (:), output (:)
REAL    (sp), INTENT(in)   :: mean
REAL    (dp)               :: one, two
INTEGER (i4)               :: i, leng, three

ALLOCATE (temp(100+INT(mean)))
leng=0;temp=0

DO i = 0,largenum
   one = mean ** i      !! h**x
   two = EXP(-mean)     !! e**-h
   three = i
   CALL factor(three)   !! x!
   IF (i == 0) THEN;
      temp(i+1)=(one*two)/1.
   ELSE;
      temp(i+1)=temp(i)+(one*two)/real(three)
   ENDIF
   IF (temp(i+1) .GE. .99999)EXIT
   leng = leng + 1
ENDDO

ALLOCATE (output(leng))
DO i = 1,leng
   output(i) = temp(i)
ENDDO
END SUBROUTINE CPoisson
!--------1---------2---------3---------4---------5---------6---------7---------8

!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine Cpoisson  ----------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8


!SUBROUTINE Cpoisson(mean,leng,type)
!USE Prog_Param, ONLY: recomb_distr,mutation_distr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! History:
!!  - Original Program:   EL Nicolazzi - July      2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculates cumulative poisson distribution,creating a vector of values as output
!! P(x)=((h**x)*(e**-h))/x!, whereas h= expected recombination events (user defined)
!!                                   x= number of effects (used to create a distribution)
!! In this subroutine: one=(h**x),two=(e**-h),three=x! and output[x]=output[x-1]+one*two/three                      
!! Output should be an allocatable (non allocated) vector of real*8, leng should be an integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INTEGER     , PARAMETER    :: i4 = SELECTED_INT_KIND(9)           
!INTEGER     , PARAMETER    :: dp  = SELECTED_REAL_KIND( 15,  307 )
!INTEGER     , PARAMETER    :: largenum=999999 ! This should do for "real" simulations
!REAL    (dp), ALLOCATABLE  :: temp (:), output(:)
!REAL    (sp), INTENT(in)   :: mean
!REAL    (dp)               :: one, two
!INTEGER (i4)               :: i, leng, three
!CHARACTER*1 , INTENT(in)   :: type

!ALLOCATE (temp(1+INT(mean*100)))
!leng=0;temp=0

!DO i = 0,largenum
!   one = mean ** i      !! h**x
!   two = EXP(-mean)     !! e**-h
!   three = i
!   CALL factor(three)   !! x!
!   IF (i == 0) THEN;
!      temp(i+1)=(one*two)/1.
!   ELSE
!      temp(i+1)=temp(i)+(one*two)/real(three)
!   ENDIF
!   IF (temp(i+1) .GE. .99999)EXIT
!   leng = leng + 1
!ENDDO

!IF(type=='R')THEN 
!   ALLOCATE(recomb_distr(leng))
!   DO i = 1,leng
!      recomb_distr(i) = temp(i)
!   ENDDO
!ENDIF

!IF(type=='M')THEN
!   ALLOCATE(mutation_distr(leng))
!   DO i = 1,leng
!      mutation_distr(i) = temp(i)
!   ENDDO
!ENDIF

!--------1---------2---------3---------4---------5---------6---------7---------8

!END SUBROUTINE
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8









!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine factor    ----------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8

SUBROUTINE factor(val)
USE Numeric_Kinds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! History:
!!  - Original Program:   EL Nicolazzi - July      2012
!! Adapted from Prof. Paolo Bison - Corso di Fondamenti di Informatica - Università di Padova.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER(i4) , INTENT(INOUT):: val
INTEGER(i4)                :: a,fat

IF (val<0) THEN;
   PRINT*,"ERROR: ASKING FOR A NEGATIVE FACTOR!!";STOP
ENDIF

fat = 1
DO a = 2,val
  fat = a * fat
ENDDO
val = fat

END SUBROUTINE factor


!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------------------------------- SUBROUTINE intrand -------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8

SUBROUTINE intrand(val,nout,out)
USE Numeric_Kinds
USE RNG
USE Prog_Param, ONLY : idum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! History:
!!  - Original Program:   Gerald Jansen (Date:????)
!!  - Jul 2012 - ELNicolazzi - Adapted to avoid repetition (brute resampling)
!!
!! Variables:
!! val:  space 
!! nout: number of (int) random values to extract
!! out:  variable (vector) of values of the sampling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER(i4), INTENT(IN)       :: val,nout
REAL   (dp), INTENT(OUT)      :: out(nout)
REAL   (dp)                   :: ran_num_1
INTEGER(i4)                   :: i,j,p

IF (nout.GT.val) THEN;
   PRINT*,"ERROR: Impossible to sample w/o replacement from the specified parameters!"
   PRINT*,"VAL(sampling field) should be higher than NOUT (#samples to be extracted)!"
   STOP
ENDIF

DO i=1,nout
   call random_number(ran_num_1)
   out(i) = 1+int(val*ran_num_1)                  ! WARNING: No mutations nor QTL allowed in the first SNP!
   DO j=1,i-1
      DO 
         IF (out(i).NE.out(j)) THEN
            EXIT
         ELSE
            DO p=1,i-1
               IF (out(i).NE.out(p)) THEN
                  CYCLE
               ELSE
                  call random_number(ran_num_1)
                  out(i) = 1+int(val*ran_num_1)   
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO
! Elegant but ineficient for sampling few numbers (ELN)
!j=0;done=0
!DO 
!   ran_num_1=ran1(idum1)
!   temp = 1+INT(val*ran_num_1)       ! WARNING: No mutations/recomb allowed in the first SNP!
!   IF (done(temp) .GT. 0) CYCLE
!   done(temp)=1
!   j=j+1
!   out(j)=temp
!   if(j.eq.nout)EXIT
!ENDDO

END SUBROUTINE intrand
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine initialize_gamma ---------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8

SUBROUTINE initialize_gamma()
USE Numeric_Kinds
USE RNG
USE Prog_param

REAL(sp)      :: init

init=random_gamma(param(1),.TRUE.)

END SUBROUTINE initialize_gamma
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8




!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!---------------------------- Subroutine sort_rec -----------------------------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8

RECURSIVE SUBROUTINE sort_rec(element,first,last)

!--------1---------2---------3---------4---------5---------6---------7---------8
USE Numeric_Kinds

IMPLICIT NONE

INTEGER(i4)             :: i,j
INTEGER(i4), INTENT(IN) :: first,last

REAL   (dp)             :: element(last), median
REAL   (dp)             :: indxval

!--------1---------2---------3---------4---------5---------6---------7---------8

i=first
j=last

median=element(INT((first+last)/2))

DO
   IF( i > j ) EXIT
   DO
      IF( median <=  element(i) ) EXIT
      i=i+1
   END DO
   DO
      IF( element(j) <= median ) EXIT
      j=j-1
   END DO
   IF ( i < j ) THEN

      indxval=element(i)
      element(i)=element(j)
      element(j)=indxval

      i=i+1
      j=j-1
   elseif ( i == j ) THEN
      i=i+1
   END IF
END DO

IF ( first < j ) CALL sort_rec(element,first,j)
IF (  i < last ) CALL sort_rec(element,i,last)

!--------1---------2---------3---------4---------5---------6---------7---------8

END SUBROUTINE sort_rec
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8









!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
! RECURSIVE SUBROUTINE sort_ascending
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8

RECURSIVE SUBROUTINE sort_ascending(element,first,last)

!--------1---------2---------3---------4---------5---------6---------7---------8
USE Numeric_Kinds

IMPLICIT none

INTEGER(i4)             :: i,j
INTEGER(i4), INTENT(in) :: first,last

REAL   (dp)             :: element(last), median
REAL   (dp)             :: indxval

!--------1---------2---------3---------4---------5---------6---------7---------8

i=first
j=last

median=element(int((first+last)/2))

do
   if( i > j ) exit
   do
      if( median     <= element(i) ) exit
      i=i+1
   enddo
   do
      if( element(j) <= median     ) exit
      j=j-1
   enddo
   if ( i < j ) then

      indxval=element(i)
      element(i)=element(j)
      element(j)=indxval

      i=i+1
      j=j-1
   elseif ( i == j ) then
      i=i+1
   endif
enddo

if ( first < j    ) call sort_ascending(element,first,j   )
if ( i     < last ) call sort_ascending(element,i    ,last)

!--------1---------2---------3---------4---------5---------6---------7---------8

END SUBROUTINE sort_ascending
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8










!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
! RECURSIVE SUBROUTINE sort_descending
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8

RECURSIVE SUBROUTINE sort_descending(element,first,last)

!--------1---------2---------3---------4---------5---------6---------7---------8

USE Numeric_Kinds

IMPLICIT none

INTEGER, INTENT(in)     :: first,last
REAL(dp)                :: element(last), median
integer                 :: i,j
REAL(dp)                :: indxval

!--------1---------2---------3---------4---------5---------6---------7---------8

i=first
j=last

median=element(int((first+last)/2))

do
   if( i > j ) exit
   do
      if( median     >= element(i) ) exit
      i=i+1
   enddo
   do
      if( element(j) >= median     ) exit
      j=j-1
   enddo
   if ( i < j ) then

      indxval=element(i)
      element(i)=element(j)
      element(j)=indxval

      i=i+1
      j=j-1
   elseif ( i == j ) then
      i=i+1
   endif
enddo

if ( first < j    ) call sort_descending(element,first,j   )
if ( i     < last ) call sort_descending(element,i    ,last)

!--------1---------2---------3---------4---------5---------6---------7---------8

END SUBROUTINE sort_descending

!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8
!--------1---------2---------3---------4---------5---------6---------7---------8







END MODULE Subs_General
