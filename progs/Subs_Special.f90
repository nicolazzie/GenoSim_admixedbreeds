!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

MODULE Subs_Special
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! HISTORY                                                                                                            
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120 
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Original program: Hossein Jorjani, 2009-2012
! Updates:                                                                    
! Jun 2012 - ELNicolazzi: - Add new variables in param file (previously hard-coded)
!                         - Zeroed "a_phenotype" variable in "generate_phenotypes" subroutine to
!                           allow compilation with gfortran
!                         - Added own sorting subroutine and slightly modified some subroutines using it.
! Jul 2012 - ELN        : - New gen_offspring subroutine with sampling from cumulative Poisson distribution. This 
!                            greatly speeds up things. (suggestion of LVarona)
!                         - Included some checks on user_defined options
!                         - Few bugs were corrected.
! Aug 2012 - ELN        : - Included pedigree file reading and checks (including parts of G.Jansen Simulations)
!                         - Gamma distribution with user defined params for QTL effects
! Sep 2012 - ELN        : - Solved a major OutOfBounds issue in offspring routine
!                         - Included new writing routines.
!                         - Doubled number of phenotypes (fully / half / not correlated)
!                         - Modified choose_n_trait routines
!                         - Solved a phenotype generation problem involving h2.
! Dec 2012 - ELN        : - Fixed a bug in get_Vg subrotuine that cut some QTLs from Vg calculation
! Jan 2013 - ELN        : - Changed the gsen_qtl position. Now Nqtl depend on Me (Ne depending) * coeff_Me (user defined)
!                         - 2 new routines in Subs_General (sort_int_updo & samp_WR) to sample QTL positions
! ...[...]...
! ELN and Luis Varona   : - Many changes.. to many to track!
! ...[...]...
!
! Jan 2014 - ELN        : - Revised mutation and recomb routines, fixing some minor, some major bugs.
!                       : - Included thresholds to output gentoypes and number of animals in pedigree routines.
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
IMPLICIT NONE

CONTAINS

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! subroutine read_user_options
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE read_user_options

USE Numeric_Kinds
USE Prog_Param

IMPLICIT NONE

INTEGER :: open_code,nf=0,nm=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

OPEN(10, FILE=TRIM(path)//'user_options.txt', IOSTAT=open_code)

IF (open_code /= 0) THEN
  DO i=1, 25
    WRITE (*,*) 'FAILED TO OPEN THE USER OPTION FILE'
  ENDDO
ENDIF
writeout=.FALSE.;realped=0

DO i=1,10; READ(10,*);ENDDO
  READ (10,*) species_male
  READ (10,*) species_female
  READ (10,*) breed_male
  READ (10,*) breed_female
  READ (10,*) n_chromosome
  READ (10,*) n_linkage
  READ (10,*) species_gen
  READ (10,*) breed_gen
  READ (10,*) recombination
  READ (10,*) mutation
  READ (10,*) Hcoeff_Me,Lcoeff_Me
  READ (10,*) h2_1,h2_2,h2_3
  READ (10,*) pcorr,mcorr
  READ (10,*) printgeno
  IF(printgeno.GE.1.OR.printgeno.LE.3)writeout=.TRUE.
  READ (10,*) thresh
  READ (10,*) genOUT
  genOUT=breed_gen+1-genOUT
  READ (10,*) pedOUTx
  READ (10,*) Nbreeds
  ALLOCATE(selez(Nbreeds))
  READ (10,*) selez
  READ (10,'(A)')selratio
  READ (10,*),admixprop
  READ (10,*) realped
  ALLOCATE(pedname(realped))
  IF(realped.GT.0)THEN
     READ (10,*) pedname
  ELSE
     READ(10,*)
  ENDIF
  READ (10,*) distr
  IF (distr.EQ.2) THEN
     ALLOCATE(param(1))
     READ (10,*) param(1)
  ENDIF
  

  !! Set selratio for males and females, separately excluding pedigree pops
  IF (Nbreeds-realped>0)THEN
     ALLOCATE(selratio_male(Nbreeds-realped))
     ALLOCATE(selratio_fem(Nbreeds-realped))
     selratio_male=0;selratio_fem=0
     DO 
        pos2 = INDEX(selratio(pos1:),sep)
        IF (pos2 == 0) THEN
           nf = nf+1
           buffer = selratio(pos1:)
           read(buffer,*) selratio_fem(nf)
           EXIT
        ENDIF
        IF (nm.LT.Nbreeds-realped)THEN
           nm = nm+1
           buffer = selratio(pos1:pos1+pos2-2)
           read(buffer,*)selratio_male(nm)
           pos1 = pos2+pos1
        ELSE
           nf = nf+1
           buffer = selratio(pos1:pos1+pos2-2)
           read(buffer,*)selratio_fem(nf)
           pos1 = pos2+pos1
        ENDIF
     ENDDO
  ENDIF
CLOSE(10)
len_genome=n_chromosome*n_linkage*n_locus
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE read_user_options

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Subroutine print_user_options
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE print_user_options
USE Prog_Param
IMPLICIT NONE
CHARACTER*10::dischr,yesno
CHARACTER*50::trait
yesno='No'

WRITE(*,*)'-------------------------------------------------------------------------------------------'
WRITE(*,*)'*******************************************'
WRITE(*,*)'** POPULATION PARAMTERS DEFINED BY USER: **'
WRITE(*,*)'*******************************************';WRITE(*,*)
WRITE(*,'(A46,I7)')' Number of final populations to be simulated:',Nbreeds
IF (realped.GT.0)THEN
WRITE(*,'(A46,I7)')' Number of populations with real pedigree   :',realped
ENDIF
WRITE(*,*)' BASE POPULATION STEP  ' 
WRITE(*,'(A46,I7)')'        male                                :',species_male
WRITE(*,'(A46,I7)')'      female                                :',species_female
WRITE(*,'(A46,I7)')' generations                                :',species_gen

IF (realped.LT.Nbreeds)THEN
WRITE(*,*)' BREED STEP '
WRITE(*,'(A46,I7)')'        male                                :',breed_male
WRITE(*,'(A46,I7)')'      female                                :',breed_female
WRITE(*,'(A46,I7)')' generations (FOR NON PEDIGREE POPs ONLY)   :',breed_gen
WRITE(*,'(A46,I7)')' NEWS**---> ADMIXTURE PROPORTION <---**NEWS :',admixprop
ENDIF
WRITE(*,*);WRITE(*,*)'*******************************************'
WRITE(*,*)'**   GENOMIC PARAMTERS DEFINED BY USER:  **'
WRITE(*,*)'*******************************************';WRITE(*,*)
WRITE(*,'(A46,I7)')' Number of chromosomes                      :', n_chromosome
WRITE(*,'(A46,I7)')' Number of Linkage blocks                   :', n_linkage
WRITE(*,'(A46,I7)')' Total number of SNPs                       :', len_genome
WRITE(*,'(A46,I7)')' -log(Mutation Rate)                        :', mutation
WRITE(*,'(A46,I7)')' Expected recombinations by chromosome      :', recombination
WRITE(*,'(A46,I7)')' Me = 2*SP_Ne*L/log(4*SP_Ne*L) =            :', Me
WRITE(*,'(A46,F7.2)')' HIGH Coeff_Me (Nqtl~ HIGHcoeff*Me)        :', Hcoeff_Me
WRITE(*,'(A46,F7.2)')' LOW  Coeff_Me (Nqtl~ HIGHcoeff*Me)        :', Lcoeff_Me
WRITE(*,'(A46,I7)')' Number of MAX QTLs required (Hcoeff * Me)  :', high_qtl
WRITE(*,'(A46,I7)')' Number of MIN QTLs required (Hcoeff * Me)  :', low_qtl
IF(distr.EQ.1)dischr='Normal'
IF(distr.EQ.2)dischr='Gamma'
IF(distr.EQ.1)WRITE(*,'(A46,1X,A6)')' Distribution of QTL effects                :',dischr
IF(distr.EQ.2)THEN
   WRITE(*,'(A46,2X,A6)')' Distribution of QTL effects                :',dischr
   WRITE(*,'(A46,2X,F5.3)')' Shape parameter of Gamma distribution      :',param(1)
ENDIF
i=0
DO j=1,Nbreeds
   IF(j.LE.realped)THEN
      WRITE(*,'(A10,I2,A,A)')' Breed n.',j,' will use true pedigree file     : ',pedname(j)
   ELSE
      i=i+1
   IF(selez(i).EQ.'H1_A')trait='Nqtl > Me, GROUP A, FIRST h2'
   IF(selez(i).EQ.'H1_B')trait='Nqtl > Me, GROUP B, FIRST h2'
   IF(selez(i).EQ.'H1_C')trait='Nqtl > Me, GROUP C, FIRST h2'
   IF(selez(i).EQ.'H2_A')trait='Nqtl > Me, GROUP A, SECOND h2'
   IF(selez(i).EQ.'H2_B')trait='Nqtl > Me, GROUP B, SECOND h2'
   IF(selez(i).EQ.'H2_C')trait='Nqtl > Me, GROUP C, SECOND h2'
   IF(selez(i).EQ.'H3_A')trait='Nqtl > Me, GROUP A, THIRD h2'
   IF(selez(i).EQ.'H3_B')trait='Nqtl > Me, GROUP B, THIRD h2'
   IF(selez(i).EQ.'H3_C')trait='Nqtl > Me, GROUP C, THIRD h2'
   IF(selez(i).EQ.'L1_A')trait='Nqtl < Me, GROUP A, FIRST h2'
   IF(selez(i).EQ.'L1_B')trait='Nqtl < Me, GROUP B, FIRST h2'
   IF(selez(i).EQ.'L1_C')trait='Nqtl < Me, GROUP C, FIRST h2'
   IF(selez(i).EQ.'L2_A')trait='Nqtl < Me, GROUP A, SECOND h2'
   IF(selez(i).EQ.'L2_B')trait='Nqtl < Me, GROUP B, SECOND h2'
   IF(selez(i).EQ.'L2_C')trait='Nqtl < Me, GROUP C, SECOND h2'
   IF(selez(i).EQ.'L3_A')trait='Nqtl < Me, GROUP A, THIRD h2'
   IF(selez(i).EQ.'L3_B')trait='Nqtl < Me, GROUP B, THIRD h2'
   IF(selez(i).EQ.'L3_C')trait='Nqtl < Me, GROUP C, THIRD h2'
   IF(selez(i).EQ.'NULL') then
        trait='NO SELECTION (random mating)'
        selratio_male(i)=100
        selratio_fem(i)=100
     ENDIF
     WRITE(*,'(A10,I2,A,A)')' Breed n.',j,' will be selected by a trait with    : ',trait
     WRITE(*,'(12X,A,I0)')                ' % of males selected for this breed  : ',selratio_male(i)
     WRITE(*,'(12X,A,I0)')                ' % of females selected for this breed: ',selratio_fem(i)
  ENDIF
ENDDO
IF(printgeno.NE.0)yesno='Yes';WRITE(*,*);
IF(yesno.EQ.'No' )WRITE(*,'(A46,5X,A3)')'  Write genotypes in GTYPES folder           : ',yesno
IF(yesno.EQ.'Yes')WRITE(*,'(A46,4X,A3)')'  Write genotypes in GTYPES folder           : ',yesno
IF(yesno.EQ.'Yes')WRITE(*,'(A46,3X,F4.2)')'  MAF threshold for gentoypes                : ',thresh
IF(yesno.EQ.'Yes')WRITE(*,'(A46,I7)')                 '(last) generations to be written (NO-PED)  :',breed_gen-genOUT+1
IF(yesno.EQ.'Yes'.AND.pedOUTx.EQ.0)WRITE(*,'(A46)')   '--> ALL animals to be written (PED)         '
IF(yesno.EQ.'Yes'.AND.pedOUTx.NE.0)WRITE(*,'(A46,I7)')'(last)animals to be written (PED)          :',pedOUTx
WRITE(*,*)'-------------------------------------------------------------------------------------------'

END SUBROUTINE print_user_options
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Subroutine control_user_options
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE control_user_options

USE Numeric_Kinds
USE Prog_Param

IMPLICIT NONE

INTEGER :: open_code

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

OPEN(10, FILE=TRIM(path)//'user_options.txt', IOSTAT=open_code)

IF (Hcoeff_Me<=0 .OR. Hcoeff_Me>len_genome) THEN
   WRITE(*,*),' ERROR: Invalid HIGH QTL coefficient is specified in param file  '
   write(*,*),' Limits are 0 (excluded) to ',len_genome,'  --> EXIT';STOP
ENDIF

IF (Lcoeff_Me<=0 .OR. Lcoeff_Me>len_genome) THEN
   WRITE(*,*),' ERROR: Invalid LOW QTL coefficient is specified in param file  '
   write(*,*),'        Limits are 0 (excluded) to ',len_genome,'  --> EXIT';STOP
ENDIF

IF (Lcoeff_Me>Hcoeff_Me) THEN
   WRITE(*,*),' ERROR: HIGH QTL_coeff < LOW QTL_coeff. '
   write(*,*),'        You should invert these in user_options.txt  --> EXIT';STOP
ENDIF

IF (genOUT.le.0.and.Nbreeds.GT.realped) THEN
   genOUT=breed_gen
   WRITE(*,*),' WARNING: Number of generations to be printed out is > than actual number of generations for breed step'
   WRITE(*,*),'          The number of generations printed was forced to ',breed_gen
ENDIF

do i=1,Nbreeds-realped
IF (selratio_male(i)<=0 .OR. selratio_male(i)>=101) THEN
   WRITE(*,*),' ERROR: Invalid SELECTION ratio for males is specified in param file  --> EXIT';STOP
ENDIF
enddo

do i=1,Nbreeds-realped
IF (selratio_fem(i)<=0 .OR. selratio_fem(i)>=101) THEN
    WRITE(*,*),' ERROR: Invalid SELECTION ratio for females is specified in param file  --> EXIT';STOP
ENDIF

IF (breed_male*Nbreeds.GT.species_male) THEN
   WRITE(*,*)," ERROR: Number of males in BREED step must be lower than SPECIES males / Number of breeds --> EXIT";STOP
ENDIF
IF (breed_female*Nbreeds.GT.species_female) THEN
   WRITE(*,*)," ERROR: Number of females in BREED step must be lower than SPECIES females / Number of breeds --> EXIT";STOP
ENDIF

IF (admixprop.GT.100.OR.admixprop.LT.0) THEN
   WRITE(*,*)," ERROR: Admixture values are allowed between 0 and 100 --> EXIT";STOP
ENDIF


enddo

END SUBROUTINE control_user_options

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! subroutine readcontrol_pedig
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE readcontrol_pedig(name)
USE Numeric_Kinds
USE Prog_Param

IMPLICIT NONE

!INTEGER(i4), DIMENSION (:,:)   :: temp
INTEGER(i4)                    :: ia,is,id,ig,ix,n
INTEGER(i4)                    :: open_code,ze,preva,prevg
LOGICAL                        :: Xit
CHARACTER*50, INTENT(IN)       :: name               
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
Xit=.FALSE.
pedig_n=0;preva=0;prevg=0;n_gen=0

OPEN(24, FILE=trim(path)//trim(name), IOSTAT=open_code, status='OLD')
IF (open_code /= 0) THEN
    WRITE (*,*) 'FAILED TO OPEN PEDIGREE FILE: ',trim(name);STOP
ENDIF

DO  
   READ(24,*,iostat=open_code)ia,is,id,ix,ig
   IF(open_code/=0)EXIT
   pedig_n=pedig_n+1                       ! COUNT # OF RECORDS AND BASE ANIMALS IN PEDIGFILE

   n_gen(ig+1,1)=n_gen(ig+1,1)+1               ! COUNT # OF RECORDS BY GEN
   IF (is.EQ.0)n_gen(ig+1,2)=n_gen(ig+1,2)+1   ! COUNT # OF MISSING SIRES BY GEN
   IF (id.EQ.0)n_gen(ig+1,3)=n_gen(ig+1,3)+1   ! COUNT # OF MISSING DAMS BY GEN

   IF(ia.LE.preva .and. ig.EQ.prevg)THEN;
      WRITE(*,*);WRITE(*,*),' --- PEDIGREE FILE ',TRIM(name),' MESSAGE -----'
      WRITE(*,*),' PEDIGREE ERROR: PEDIGREE FILE NOT SORTED SEQUENTIALLY IN GENERATION',ig,'(BY ANIMAL)'
      WRITE(*,*),' ANIMAL',ia,'IN ROW',pedig_n,'WHILE ANIMAL',preva,'IN THE PREVIOUS ROW.'
      Xit=.TRUE.
   ENDIF
   IF(ia.LE.is .OR. ia.LE.id) THEN
      WRITE(*,*);WRITE(*,*),' --- PEDIGREE FILE ',TRIM(name),' MESSAGE -----'
      WRITE(*,*),' PEDIGREE ERROR: CONSISTENCY IN PARENTS NUMBER'
      WRITE(*,*),' ANIMALS:',ia,'PARENTS:',is,id
      WRITE(*,*),' PARENTS SHOULD BE NUMBERED BEFORE PROGENY'
      Xit=.TRUE.
   ENDIF
   preva=ia
   prevg=ig
ENDDO
IF(Xit)STOP

REWIND(24)

ALLOCATE(Ped(0:pedig_n),cgens(1000,2),gens(1000,maxval(n_gen(:,1)),2))      ! ALLOCATE PEDIGREE INFORMATION
cgens=0;gens=0

DO ze = 0,pedig_n
   Ped(ze) = ped_type(0,0,0,0,0)
ENDDO
n=0
DO
   READ(24,*,iostat=open_code) ia,is,id,ix,ig
   IF(open_code/=0)EXIT
   ig=ig+1                                  ! AUGMENT Ngen+1 to avoid 0 reference
   n=n+1
   Ped(n) = ped_type(ia,is,id,ix,ig)       ! STORE PEDIGREE

   cgens(ig,ix)=cgens(ig,ix)+1              ! STORE # ANIMALS in GEN divided by sex
   gens(ig,cgens(ig,ix),ix)=ia              ! STORE ANIMALSid (DIVIDED for GEN and SEX)

   IF(ia.GT.pedig_n)THEN;
      WRITE(*,*);WRITE(*,*),' --- PEDIGREE FILE ',TRIM(name),' MESSAGE -----'
      WRITE(*,*),' PEDIGREE ERROR: NON SEQUENTIAL NUMBER OF ANIMALS'
      WRITE(*,*),' ANIMAL #',ia,'IS HIGHER THAN # OF ANIMALS IN PEDIGREE (',pedig_n,')'
      WRITE(*,*),' ANIMALS SHOULD BE RENUMBERED SEQUENTIALLY'
      Xit=.TRUE.
   ENDIF
END DO

!! ALLOW 25 SONS FOR EACH SIRE & 5 FOR EACH DAMset
needSIRE=INT(n_gen(1,1)/25)+50            ! Add 50 random bulls to the dataset
needDAM =INT(n_gen(1,1) /5)+100           ! Add 100 random dams to the data

IF(Xit)STOP
CLOSE(24)
END SUBROUTINE readcontrol_pedig

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! subroutine initialize_h2
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! ELN: Jan 13 - Heavily modified to suit the new distribution of traits
!

SUBROUTINE initialize_h2(h2)
USE Numeric_Kinds
USE Prog_Param, ONLY : h2_1,h2_2,h2_3
IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

REAL, DIMENSION(3*3*2)      :: h2  ! hard coded Only 3 traits and h2 are allowed (3h2*3groups(corr)*2geneticSTRC)
INTEGER(i4)                 :: trait,struct
h2=0.;

DO trait=1,3
   DO struct=0,1
      h2(trait+0+(9*struct))= h2_1      
      h2(trait+3+(9*struct))= h2_2   
      h2(trait+6+(9*struct))= h2_3
   ENDDO
ENDDO

END SUBROUTINE initialize_h2

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! subroutine gen_qtl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE gen_qtl

! Go through all SNPs, assign some of them as QTLs and give an allelic effect to them.
!
! All traits are standardized to have a phenotypic variance of 100 at the base population. Therefore, the additive
! variance would be Vp * h2. For the time being, only additive value is simulated, therefore the additive variance is
! Va = Sigma (2 * p * q) a2 = 2 * n * p * q * a2, where n is number of loci
! a = SQRT (Va/(2 * n * p * q)) Given the fact that in the base population the expected value of p and q are 0.50,
! a = SQRT(Va/(0.5 * n)) = SQRT((2 * Va / n)
!
! Almost re-built and (hopefully) bug free - Jan 2014 (ELN)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE RNG
USE Subs_General
USE Prog_Param, ONLY: idum1, len_genome, n_trait, h2, distr, param, high_qtl, low_qtl, allele_effect, &
                      qtl_position, qtl_size, qtl_type, pcorr, mcorr, HmaxQTL, frequencies
IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ALLpos, qtl_posMeH,qtl_posMeL, corAA, corAB, corAC, corXX, corXY, corXZ 
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: temp_qtl, qtl_temp, tempALL, tempAA, tempAB, tempAC 
INTEGER(i4) :: GrHA=0, GrHB=0, grHC=0, GrLX=0, GrLY=0, GrLZ=0, LmaxQTL
INTEGER(i4) :: NcorAB, NcorAC, NcorXY, NcorXZ, NuncorB, NuncorC, NuncorY, NuncorZ
INTEGER(i4) :: i, j, k, k1, k2, k3, q1, q2, q3, snp, trait, lenMAF
REAL   (dp) :: ran_num_2, Bpct, Ypct, thres
LOGICAL     :: FOUND
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
allele_effect=0.
!! Select SNPs not ~ monomorphic (change to an user_defined min MAF in the future?)
ALLOCATE(ALLpos(len_genome))
thres=.95
k=0
lenMAF=0
DO i=1,len_genome
   IF(frequencies(i,1).GE.thres .OR. frequencies(i,1).LE.(1.-thres)) CYCLE
   IF(frequencies(i,2).GE.thres .OR. frequencies(i,2).LE.(1.-thres)) CYCLE
   lenMAF=lenMAF+1
   ALLpos(lenMAF)=i
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! GROUPS A - B - C and correlations (AB - AC) !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set the (theoretical) total number of SNPs as a function of Me
NuncorB=INT(high_qtl*(1-ABS(pcorr)))
NcorAB =high_qtl-NuncorB
NuncorC=INT(high_qtl*(1-ABS(mcorr)))
NcorAC =high_qtl-NuncorC
HmaxQTL=high_qtl+NuncorB+NuncorC                          ! Nqtl + Nqtl*uncorrQTL(+CORR) + Nqtl*uncorrQTL(-CORR)

ALLOCATE (tempALL( HmaxQTL)                  )
ALLOCATE (tempAA (high_qtl), corAA (high_qtl))
ALLOCATE (tempAB (  NcorAB), corAB (  NcorAB))
ALLOCATE (tempAC (  NcorAC), corAC (  NcorAC))
CALL samp_WR(1, high_qtl, NcorAB, tempAB)                 ! Sample QTLs for group B that correlates with A
CALL samp_WR(1, high_qtl, NcorAC, tempAC)                 ! Sample QTLs for group C that correlates with A

!! Allocate final vectors with new definition of N of QTLs
ALLOCATE (qtl_posMeH  (HmaxQTL), qtl_temp(HmaxQTL))
ALLOCATE (qtl_position(HmaxQTL), qtl_size(HmaxQTL,n_trait), qtl_type(HmaxQTL,6)) 

!! Get positions for ALL QTLs and for Groups: AA, AB and AC
CALL samp_WR(1, lenMAF, HmaxQTL, tempALL)                  ! ALL (temporary)
CALL samp_WR(1, HmaxQTL, high_qtl, tempAA)                 ! Group AA (temporary

qtl_posMeH=ALLpos(tempALL)                                 ! ALL -> Qtl positions on SNPs with < thresh MAF
qtl_position=qtl_posMeH                                    ! ALL (for export)

corAA=qtl_posMeH(tempAA)                                   ! Group AA positions
corAB=corAA(tempAB)                                        ! Group AB positions
corAC=corAA(tempAC)                                        ! Group AC positions
Bpct=REAL(NuncorB)/REAL(high_qtl)                          ! Percentage use to divide BB from CC's

DEALLOCATE(tempAA,tempAB,tempAC)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! GROUPS X - Y - Z and correlations (XY - XZ) !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set the (theoretical) total number of SNPs as a function of Me
NuncorY=INT(low_qtl*(1-ABS(pcorr)))
NcorXY =low_qtl-NuncorY
NuncorZ=INT(low_qtl*(1-ABS(mcorr)))
NcorXZ =low_qtl-NuncorZ
LmaxQTL=low_qtl+NuncorY+NuncorZ                           ! Nqtl + Nqtl*uncorrQTL(+CORR) + Nqtl*uncorrQTL(-CORR)

ALLOCATE (tempAA (low_qtl), corXX (high_qtl))
ALLOCATE (tempAB ( NcorXY), corXY (  NcorXY))
ALLOCATE (tempAC ( NcorXZ), corXZ (  NcorXZ))

CALL samp_WR(1, low_qtl, NcorXY, tempAB)                  ! Sample QTLs for group Y that correlates with X
CALL samp_WR(1, low_qtl, NcorXZ, tempAC)                  ! Sample QTLs for group Z that correlates with X

!! Get positions for ALL QTLs
ALLOCATE (qtl_posMeL (LmaxQTL), temp_qtl (LmaxQTL))
CALL samp_WR(1, HmaxQTL, LmaxQTL, temp_qtl)               ! Get (low)Nqtl positions (sampled form qtl_posMeH)
qtl_posMeL = qtl_posMeH(temp_qtl)                         ! LOW Nqtl positions
 
!! Get positions for Groups: XX, XY and XZ
CALL samp_WR(1, LmaxQTL, low_qtl, tempAA)                 
corXX=qtl_posMeL(tempAA)                                  ! Group XX positions
corXY=     corXX(tempAB)                                  ! Group XY positions
corXZ=     corXX(tempAC)                                  ! Group XZ positions
Ypct=REAL(NuncorY)/REAL(low_qtl)                         ! Percentage use to divide BB from CC's

DEALLOCATE(tempAA,tempAB,tempAC,temp_qtl)

!! Go through all SNPs ... but do something only on those that become QTLs
qtl_size=0.D0
i=1;j=1
k1=1;k2=1;k3=1
q1=1;q2=1;q3=1
qtl_type(:,:)='-'

DO snp=1, len_genome                                     ! Total number of SNPs = n_chromosome*n_linkage*n_locus
   FOUND=.FALSE.
   IF(i .GT. HmaxQTL)CYCLE
   IF (snp .NE. qtl_posMeH(i))CYCLE                          ! Go only for QTLs (avoid SNPs). Nqtl >>> Me step
   IF (snp .EQ. corAB(k1) .AND. k1 .LT. NcorAB) THEN         ! Group AB
      qtl_type(i,1)='A';GrHA=GrHA+1
      qtl_type(i,2)='B';GrHB=GrHB+1
      k1=k1+1
      FOUND=.TRUE.
   ENDIF
   IF (snp .EQ. corAC(k2) .AND. k2 .LT. NcorAC) THEN         ! Group AC (count A only if not accounted for before)
      qtl_type(i,3)='C';GrHC=GrHC+1
      IF (.NOT.FOUND) THEN
         qtl_type(i,1)='A';GrHA=GrHA+1
      ENDIF
      FOUND=.TRUE.
      k2=k2+1
   ENDIF
   IF (snp .EQ. corAA(k3) .AND. k3 .LT. high_qtl) THEN       ! Group AA (count it only if not accounted for)
      IF (.NOT.FOUND) THEN
         qtl_type(i,1)='A';GrHA=GrHA+1
         FOUND=.TRUE.
      ENDIF
      k3=k3+1
   ENDIF
   IF (.NOT.FOUND) THEN                                      ! If not ANY of the above
      ran_num_2=ran1(idum1)                           
      IF (GrHB.LT. high_qtl .AND. GrHC.LT.high_qtl)THEN
         IF (ran_num_2 .LE. Bpct) THEN                       ! Assign BB and CC in order to the correlation with AA
            qtl_type(i,2)='B';GrHB=GrHB+1
         ELSE 
            qtl_type(i,3)='C';GrHC=GrHC+1
         ENDIF
      ELSE
         IF(GrHB.GE.high_qtl.AND.GrHC.LT.high_qtl)THEN       ! If B is full, assign QTL to C
            qtl_type(i,3)='C';GrHC=GrHC+1
         ELSEIF(GrHB.LT.high_qtl.AND.GrHC.GE.high_qtl)THEN   ! If C is full, assign QTL to B
            qtl_type(i,2)='B';GrHB=GrHB+1
         ENDIF
      ENDIF
   ENDIF

   DO trait=1,9,3                                    ! CAREFUL!! THIS SECTION IS HARDCODED!!! (ELN)
      IF (distr.EQ.1)THEN                            ! We want 30 traits with Nqtl= ceff_Me * Me (=2*Ne*L/log(4*Ne*L))
         allele_effect=gasdev(idum1)                 ! Allele effects sampled from a normal distribution
      ELSEIF (distr.EQ.2)THEN
         allele_effect=random_gamma(param(1),.FALSE.)! Allele effects sampled from a gamma  distribution
      ENDIF
      
      IF (qtl_type(i,1)=='A') THEN
         qtl_size(i,trait) = allele_effect * SQRT(2.*h2(trait)/REAL(high_qtl))  ! Assign a size for the QTL effect for A
         IF (qtl_type(i,2)=='B' .AND. pcorr .GT. 0) qtl_size(i,trait+1)=  qtl_size(i,trait)    ! Same effect in A and B
         IF (qtl_type(i,2)=='B' .AND. pcorr .LT. 0) qtl_size(i,trait+1)=-(qtl_size(i,trait))   ! Opposite (if neg corr) effect in A and B
         IF (qtl_type(i,3)=='C' .AND. mcorr .GT. 0) qtl_size(i,trait+2)=  qtl_size(i,trait)    ! Same effect in A and B
         IF (qtl_type(i,3)=='C' .AND. mcorr .LT. 0) qtl_size(i,trait+2)=-(qtl_size(i,trait))   ! Opposite (if neg corr) effect in A and B
      ELSEIF (qtl_type(i,2)=='B') THEN;
         qtl_size(i,trait+1) = allele_effect * SQRT(2.*h2(trait+1)/REAL(high_qtl)) ! Assign a size for the QTL effect for BB
      ELSEIF (qtl_type(i,3)=='C') THEN;
         qtl_size(i,trait+2) = allele_effect * SQRT(2.*h2(trait+2)/REAL(high_qtl)) ! Assign a size for the QTL effect for CC
      ENDIF
   ENDDO
!   print*,i,qtl_posMeH(i),qtl_type(i,1:6),qtl_size(i,1:3),GrHA,GrHB,GrHC, high_qtl
   i=i+1                                               ! Go to the next QTL (Traits 1-60)

   
   IF(j.EQ.LmaxQTL+1)CYCLE
   IF (snp .NE. qtl_posMeL(j))CYCLE                    ! Go only for QTLs (avoid SNPs). Nqtl >>> Me step
   FOUND=.FALSE.
   i=i-1
   IF (snp .EQ. corXY(q1) .AND. q1 .LT. NcorXY) THEN   ! Group XY
      qtl_type(i,4)='X';GrLX=GrLX+1
      qtl_type(i,5)='Y';GrLY=GrLY+1
      q1=q1+1
      FOUND=.TRUE.
   ENDIF
   IF (snp .EQ. corXZ(q2) .AND. q2 .LT. NcorXZ) THEN   ! Group ZX (count X only if not accounted for before) 
      qtl_type(i,6)='Z';GrLZ=GrLZ+1
      q2=q2+1
      IF (.NOT.FOUND) THEN
         qtl_type(i,4)='X';GrLX=GrLX+1
      ENDIF
      FOUND=.TRUE.
   ENDIF
   IF (snp .EQ. corXX(q3) .AND. q3 .LT. low_qtl) THEN  ! Group XX (count X only if not accounted for before) 
      IF (.NOT.FOUND) THEN
         qtl_type(i,4)='X';GrLX=GrLX+1
         FOUND=.TRUE.
      ENDIF
      q3=q3+1
   ENDIF
   IF (.NOT.FOUND) THEN                                ! If not ANY of the above
      ran_num_2=ran1(idum1)                           
      IF (GrLY .LT. low_qtl .AND. GrLZ .LT. low_qtl)THEN
         IF (ran_num_2 .LE. Ypct) THEN                 
            qtl_type(i,5)='Y';GrLY=GrLY+1
         ELSE 
            qtl_type(i,6)='Z';GrLZ=GrLZ+1
         ENDIF
      ELSE
          IF(GrLY.GE.low_qtl.AND.GrLZ.LT.low_qtl)THEN
            qtl_type(i,6)='Z';GrLZ=GrLZ+1
         ELSEIF(GrLY.LT.low_qtl.AND.GrLZ.GE.low_qtl)THEN
            qtl_type(i,5)='Y';GrLY=GrLY+1
         ENDIF
      ENDIF
   ENDIF
   
   DO trait=10,18,3                                  ! CAREFUL!! THIS SECTION IS HARDCODED!!! (ELN)
      IF (distr.EQ.1)THEN                            ! We want 30 traits with Nqtl= ceff_Me * Me (=2*Ne*L/log(4*Ne*L))
         allele_effect=gasdev(idum1)                 ! Allele effects sampled from a normal distribution
      ELSEIF (distr.EQ.2)THEN
         allele_effect=random_gamma(param(1),.FALSE.)! Allele effects sampled from a gamma  distribution
      ENDIF
      
      IF (qtl_type(i,4)=='X') THEN
         qtl_size(i,trait) = allele_effect * SQRT(2.*h2(trait)/REAL(low_qtl))  ! Assign a size for the QTL effect for A
         IF (qtl_type(i,5)=='Y' .AND. pcorr .GT. 0) qtl_size(i,trait+1)=  qtl_size(i,trait)    ! Same effect in A and B
         IF (qtl_type(i,5)=='Y' .AND. pcorr .LT. 0) qtl_size(i,trait+1)=-(qtl_size(i,trait))   ! Opposite (if neg corr) effect in A and B
         IF (qtl_type(i,6)=='Z' .AND. mcorr .GT. 0) qtl_size(i,trait+2)=  qtl_size(i,trait)    ! Same effect in A and B
         IF (qtl_type(i,6)=='Z' .AND. mcorr .LT. 0) qtl_size(i,trait+2)=-(qtl_size(i,trait))   ! Opposite (if neg corr) effect in A and B
      ELSEIF (qtl_type(i,5)=='Y') THEN;
         qtl_size(i,trait+1) = allele_effect * SQRT(2.*h2(trait+1)/REAL(low_qtl))  ! Assign a size for the QTL effect for BB
      ELSEIF (qtl_type(i,6)=='Z') THEN;
         qtl_size(i,trait+2) = allele_effect * SQRT(2.*h2(trait+2)/REAL(low_qtl))  ! Assign a size for the QTL effect for CC
      ENDIF
   ENDDO
!   print*,i,qtl_posMeH(i),qtl_type(i,:),qtl_size(i,10:12),GrLX,GrLY,GrLZ
   j=j+1                                               ! Go to the next QTL (Traits 61-120)
   i=i+1                                              
ENDDO
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
END SUBROUTINE gen_qtl

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! subroutine gen_species_base
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE gen_species_base (n_animal, a_genome)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
!USE Subs_General
USE RNG
USE Prog_Param, ONLY: idum1, n_chromosome, n_linkage, n_locus, n_allele

IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4), INTENT( IN)                                                      :: n_animal
INTEGER(i4), INTENT(OUT), DIMENSION(n_animal,n_chromosome,n_linkage,n_allele) :: a_genome
! For the time being all chromosomes have equal length. However, later the ! chromosomes should have different size. This
! can be accommodated by changing n_linkage to a vector so that each chromosome can have different number of linkage
! groups. Even the introduction of sex chromosome may be done by changing n_linkage from a vector to an array, so that the
! last chromosome has different number ! of linkage groups for the two chromosomes.
INTEGER(i4) :: animal
INTEGER(i4) :: chromosome, linkage, locus, allele, bit_pos
REAL   (dp) :: ran_num_1

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
a_genome=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO animal=1,n_animal
 DO chromosome=1,n_chromosome
  DO linkage=1,n_linkage
   DO locus=1,n_locus
    DO allele=1,n_allele
     bit_pos=locus-1
!     ran_num_1=ran1(idum1)
!     if (ran_num_1.gt..5)then
        a_genome(animal,chromosome,linkage,allele)=ibclr(a_genome(animal,chromosome,linkage,allele),bit_pos)
!     else
!        a_genome(animal,chromosome,linkage,allele)=ibset(a_genome(animal,chromosome,linkage,allele),bit_pos)
!     endif
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE gen_species_base

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE randomize
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE randomize (n_animal,n_selec, sort_animal)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
USE Numeric_Kinds
USE Subs_General
USE RNG
USE Prog_Param, ONLY: idum1

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120                              
INTEGER(i4), INTENT( IN)                       :: n_animal
INTEGER(i4), INTENT( IN)                       :: n_selec
INTEGER(i4), INTENT(OUT), DIMENSION(n_selec)   :: sort_animal
REAL(dp),                 DIMENSION(n_selec)   :: temp
INTEGER(i4)                                    :: i
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
sort_animal=0;temp=0;
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120                               
CALL RANDOM_NUMBER(temp)
DO i= 1, n_selec
   sort_animal(i)=i
ENDDO

call sort_updo(temp,sort_animal,n_selec,'high_to_low')

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE randomize

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE sort_phenotypes
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE sort_phenotypes (n_animal, a_phenotype, sort_animal)
! This subroutine does two things: Sorts the phenotype in ascending order
!                                  "performs" selection
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE RNG
USE Subs_General

IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
INTEGER(i4), INTENT( IN)                       :: n_animal
INTEGER(i4), INTENT(OUT), DIMENSION(n_animal)  :: sort_animal
REAL   (dp), INTENT(IN ), DIMENSION(n_animal)  :: a_phenotype
INTEGER(i4)                       :: i
REAL   (dp), DIMENSION(n_animal)  :: sort_vector

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120                           
sort_animal=0;
sort_vector=0.D0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120        

DO i= 1, n_animal
   sort_vector(i)=a_phenotype(i)
   sort_animal(i)=i
ENDDO

CALL sort_updo (sort_vector, sort_animal, n_animal, 'high_to_low')

END SUBROUTINE sort_phenotypes

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE select_female (n_animal, n_selec, a_phenotype, selec_animal)
! This subroutine does two things: Sorts the phenotype in ascending order
!                                  "performs" selection
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE RNG
USE Subs_General
USE Prog_Param, ONLY: selratioF,idum1

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
INTEGER(i4),INTENT( IN)                     :: n_animal
INTEGER(i4),INTENT( IN)                     :: n_selec
REAL   (dp),INTENT( IN),DIMENSION(n_animal) :: a_phenotype
INTEGER(i4),INTENT(OUT),DIMENSION( n_selec) :: selec_animal
INTEGER(i4)            ,DIMENSION(n_animal) :: sort
INTEGER(i4)                                 :: i
REAL   (dp)            ,DIMENSION(n_animal) :: sort_vector
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

sort=0;
sort_vector=0.D0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO i= 1, n_animal
sort_vector(i)=a_phenotype(i)
sort(i)=i
ENDDO
CALL sort_updo (sort_vector, sort, n_animal, 'high_to_low')
DO i=1,n_selec
   selec_animal(i)=sort(int(ran2(idum1)*n_animal*selratioF/100)+1)
enddo

END SUBROUTINE select_female
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Subroutine SELECT_MALE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
SUBROUTINE select_male (n_animal, n_selec, a_phenotype, selec_animal)
! This subroutine does two things: Sorts the phenotype in ascending order
!                                  "performs" selection
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE RNG
USE Subs_General
USE Prog_Param, ONLY: selratioM,idum1

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
INTEGER(i4),INTENT( IN)                     :: n_animal
INTEGER(i4),INTENT( IN)                     :: n_selec
REAL   (dp),INTENT( IN),DIMENSION(n_animal) :: a_phenotype
INTEGER(i4),INTENT(OUT),DIMENSION( n_selec) :: selec_animal
INTEGER(i4)            ,DIMENSION(n_animal) :: sort
INTEGER(i4)                                 :: i
REAL   (dp)            ,DIMENSION(n_animal) :: sort_vector
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

sort=0;
sort_vector=0.D0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO i= 1, n_animal
sort_vector(i)=a_phenotype(i)
sort(i)=i
ENDDO

CALL sort_updo (sort_vector, sort, n_animal, 'high_to_low')
DO i=1,n_selec
   selec_animal(i)=sort(int(ran2(idum1)*n_animal*selratioM/100.)+1)
enddo

END SUBROUTINE select_male

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE gen_offspring
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE gen_offspring (n_animal, n_male, n_female, s_genome, d_genome, a_genome, myn)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
USE Numeric_Kinds
USE RNG
USE Subs_General
USE Prog_Param,   ONLY: n_chromosome, n_linkage, n_locus, n_allele, sort_sire, sort_dam, gen, &
                        recomb_distr, len_recomb, mutation_distr, len_mutation, idum1, len_genome
IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
LOGICAL    , INTENT( IN)                                                      :: myn
INTEGER(i4), INTENT( IN)                                                      :: n_animal,n_male,n_female
INTEGER(i4), INTENT( IN), DIMENSION(  n_male,n_chromosome,n_linkage,n_allele) :: s_genome
INTEGER(i4), INTENT( IN), DIMENSION(n_female,n_chromosome,n_linkage,n_allele) :: d_genome
INTEGER(i4), INTENT(OUT), DIMENSION(n_animal,n_chromosome,n_linkage,n_allele) :: a_genome
INTEGER(i4)                                                                   :: j, nrec, nmut, animal, exact,c
INTEGER(i4)                                                                   :: chromosome, locus, s_allele, d_allele, bit_pos
INTEGER(i4), ALLOCATABLE, DIMENSION(:)                                        :: s_recomb,d_recomb,s_mutat,d_mutat,order, lgru
REAL   (dp), ALLOCATABLE, DIMENSION(:)                                        :: temp,one
REAL   (dp)                                                                   :: ran_num_1

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
a_genome=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO animal=1,n_animal
   DO chromosome=1,n_chromosome
      !! PATERNAL HAPLOTYPE - START CHROMOSOME SEGMENT
      CALL init_random_seed
      CALL RANDOM_NUMBER(ran_num_1)
      IF(ran_num_1 < 0.5D0) THEN
         s_allele=1                                               ! Choose sire's paternal chromosome
      ELSE
         s_allele=2                                               ! Choose sire's maternal chromosome
      ENDIF
      !! RECOMBINATION for SIRE chromosome
      ran_num_1=ran2(idum1)
      DO nrec = 1, len_recomb                                     ! Define number of recomb in this an/chrom
         IF (ran_num_1 < recomb_distr(nrec))EXIT                  ! Nrec=1, means that NO RECOMBINATION OCCURS
      ENDDO
      IF (ran_num_1 >= recomb_distr(len_recomb)) nrec=nrec+1            ! If rand_num is > to max distr, add 1 extra recomb

      IF (nrec .EQ. 1)THEN;                                       ! NO RECOMBINATION
         DO j=1,n_linkage
            a_genome(animal,chromosome,j,1) = s_genome(sort_sire(animal),chromosome,j,s_allele)
         ENDDO
      ELSE                                                        ! AT LEAST ONE RECOMBINATION
         nrec=nrec-1                                              ! REDUCE THE NUMBER OF NREC (CPoisson distribution starts from p(0))
         ALLOCATE (s_recomb(nrec), order(nrec), temp(nrec),one(1))
         s_recomb=0;order=0;temp=0
         CALL intrand(n_linkage*n_locus, nrec, temp)              ! Create a vector of (random) positions of NREC recombinations
!         print*,"RECOM SIRE",animal,chromosome,temp
         FORALL(j=1:nrec) order(j)=j                              ! This vector is needed by sorting subroutine
         CALL sort_updo (temp, order, nrec, 'low_to_high')        ! Sort positions (these'll be chunk limits for the chromosome)
         ! CREATE CHUNKS OF LINKAGE FOR THIS CHROMOSOME
         s_recomb = int(floor(temp(:)/n_locus))                   ! All blocks where recombination happens along the chromosome
         c=1
         DO j = 1, n_linkage
            IF (s_recomb(c).EQ.0)s_recomb(c)=1
            IF (j .NE. s_recomb(c)) THEN
               a_genome(animal,chromosome,j,1) = s_genome(sort_sire(animal),chromosome,j,s_allele)
            ELSE
               exact=MOD(INT(temp(c)), n_locus)                      ! EXACT POSITION OF SNP RECOMBINING, WITHIN RECOM LINKAGE BLOCK
               IF (s_recomb(c).GE.n_linkage)THEN;
                  CALL intrand(n_locus,1,one)
                  exact=one(1)
               ENDIF
               DO locus = 1, n_locus                                 ! ITERATE FOR SINGLE LOCUS
                  bit_pos=locus-1
                  IF(locus == exact)THEN;                             ! NOT RECOMBINING LOCI
                     IF (s_allele==1) THEN
                        s_allele=2                           ! Shift to the other sire's chromosome (in the next chunk)
                     ELSE
                        s_allele=1                           ! Shift to the other sire's chromosome (in the next chunk)
                     ENDIF
                  ENDIF
                  IF(BTEST(s_genome(sort_sire(animal),chromosome,j,s_allele),bit_pos)) THEN                     ! If bit_pos = 1
                     a_genome(animal,chromosome,j,1)=IBSET(a_genome(animal,chromosome,j,1),bit_pos) ! Keep it as it is
                  ELSE                                                                                                      ! If bit_pos = 0
                     a_genome(animal,chromosome,j,1)=IBCLR(a_genome(animal,chromosome,j,1),bit_pos) ! Keep it as it is
                  ENDIF
               ENDDO
               c=c+1
               IF (c.GE.nrec)c=nrec
            ENDIF
         ENDDO
         DEALLOCATE(s_recomb, order, temp,one)
      ENDIF

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

  !! MATERNAL HAPLOTYPE - START CHROMOSOME SEGMENT
      ran_num_1=ran2(idum1)
      IF(ran_num_1 < 0.5D0) THEN
         d_allele=1                                               ! Choose dam's paternal chromosome
      ELSE
         d_allele=2                                               ! Choose dam's maternal chromosome
      ENDIF

  !! RECOMBINATION for DAM chromosome
      CALL RANDOM_NUMBER(ran_num_1)
      DO nrec = 1, len_recomb                                     ! Define number of recomb in this an/chrom
         IF (ran_num_1 < recomb_distr(nrec))EXIT                  ! Nrec=1, means that NO RECOMBINATION OCCURS
      ENDDO
      IF (ran_num_1 >= recomb_distr(len_recomb)) nrec=nrec+1            ! If rand_num is > to max distr, add 1 extra mutation

      IF (nrec .EQ. 1)THEN;                                       ! NO RECOMBINATION
         do j=1,n_linkage
            a_genome(animal,chromosome,j,2) = d_genome(sort_dam(animal),chromosome,j,d_allele)
         enddo
      ELSE                                                        ! AT LEAST ONE RECOMBINATION
         nrec=nrec-1                                              ! REDUCE THE NUMBER OF NREC (CPoisson distribution starts from p(0))
         ALLOCATE (d_recomb(nrec), order(nrec), temp(nrec),one(1))
         d_recomb=0;order=0;temp=0
         CALL intrand(n_linkage*n_locus, nrec, temp)              ! Create a vector of (random) positions of NREC recombinations
!         print*,"RECOM DAM",animal,chromosome,temp
         FORALL(j=1:nrec) order(j)=j                              ! This vector is needed by sorting subroutine
         CALL sort_updo (temp, order, nrec, 'low_to_high')        ! Sort positions (these'll be chunk limits for the chromosome)

         ! CREATE CHUNKS OF LINKAGE FOR THIS CHROMOSOME
         d_recomb = INT(FLOOR(temp(:)/n_locus))         ! All blocks where recombination happens along the chromosome

         c=1
         DO j = 1, n_linkage
            IF (d_recomb(c).EQ.0)d_recomb(c)=1
            IF (j .NE. d_recomb(c)) THEN
               a_genome(animal,chromosome,j,2) = d_genome(sort_dam(animal),chromosome,j,d_allele)
            ELSE
               exact=MOD(INT(temp(c)), n_locus)                      ! EXACT POSITION OF SNP RECOMBINING, WITHIN RECOM LINKAGE BLOCK
               IF (d_recomb(c).GE.n_linkage)THEN;
                  CALL intrand(n_locus,1,one)
                  exact=one(1)
               ENDIF
               DO locus = 1, n_locus                                 ! ITERATE FOR SINGLE LOCUS
                  bit_pos=locus-1
                  IF(locus == exact)THEN;                             ! NOT RECOMBINING LOCI
                     IF (d_allele==1) THEN
                        d_allele=2                           ! Shift to the other dam's chromosome (in the next chunk)
                     ELSE
                        d_allele=1                           ! Shift to the other dam's chromosome (in the next chunk)
                     ENDIF
                  ENDIF
                  IF(BTEST(d_genome(sort_dam(animal),chromosome,j,d_allele),bit_pos)) THEN                     ! If bit_pos = 1
                     a_genome(animal,chromosome,j,2)=IBSET(a_genome(animal,chromosome,j,2),bit_pos) ! Keep it as it is
                  ELSE                                                                                                      ! If bit_pos = 0
                     a_genome(animal,chromosome,j,2)=IBCLR(a_genome(animal,chromosome,j,2),bit_pos) ! Keep it as it is
                  ENDIF
               ENDDO
               c=c+1
               IF (c .GE. nrec) c=nrec
            ENDIF
         ENDDO
         DEALLOCATE(d_recomb, order, temp,one)
      ENDIF
!               print*
!               print*,"RECOMBINATION!!! LOCUS",locus,"CHUNK",d_recomb(j+1),"DOVE",exact,"STRAND",d_allele
!               print*,"RDAM1:",IBITS(d_genome(sort_dam(animal),chromosome,d_recomb(j+1),1),bit_pos,1)
!               print*,"RDAM1:",IBITS(d_genome(sort_dam(animal),chromosome,d_recomb(j+1),2),bit_pos,1)
!               print*,"RSON1:",IBITS(a_genome(animal,chromosome,d_recomb(j+1),1),bit_pos,1)
!               print*,"RSON2:",IBITS(a_genome(animal,chromosome,d_recomb(j+1),2),bit_pos,1)
!               print*
IF (myn) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! MUTATION for SIRE chromosomes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(len_mutation.EQ.0)CYCLE

      ran_num_1=ran2(idum1)
      DO nmut = 1, len_mutation                                   ! Define number of mutations in this an/chrom
         IF (ran_num_1 < mutation_distr(nmut))EXIT
      ENDDO
      IF (ran_num_1 >= mutation_distr(len_mutation)) nmut=nmut+1          ! If rand_num is > to max distr, add 1 extra mutation

      !!  ## Number of mutations per animal: Position 1 is for 0 mutations (so it is skipped)
      IF (nmut .GT. 1)THEN;                                       ! DO SOMETHING ONLY IF NMUT > 1. Distribution starts from p(0)
         nmut=nmut-1                                              ! Reduce nmut
         ALLOCATE (s_mutat(nmut), order(nmut), temp(nmut),lgru(nmut),one(1))
         s_mutat=0;order=0;temp=0
!         print*,"MUT_SIRE - ANIMAL:",animal,chromosome
         CALL intrand(n_linkage*n_locus, nmut, temp) ! Create a vector of (random) positions of NREC recombinations
!         print*,"MUT SIRE",animal,chromosome,temp
         FORALL(j=1:nmut) order(j)=j                              ! This vector is needed by sorting subroutine
         CALL sort_updo (temp, order, nmut, 'low_to_high')        ! Sort positions
         lgru   = INT(FLOOR(temp(:)/n_locus))       ! All blocks where recombination happens along the chromosome
         s_mutat= MOD(INT(temp(:)), n_locus)
         DO j = 1, nmut                                           ! N_MUTATIONS
            if (lgru(j).eq.0)lgru(j)=1
            IF (lgru(j).GE.n_linkage)THEN;
!                  lgru(j)=n_linkage                 
               CALL intrand(n_locus,1,one)
               s_mutat(j)=one(1)
            ENDIF
!               print*, animal,chromosome,lgru(j),s_mutat(j),temp(j),n_locus,one
            IF(s_mutat(j).EQ.0)CYCLE
            IF(BTEST(a_genome(animal,chromosome,lgru(j),1),s_mutat(j))) THEN                                       ! If bit_pos = 1
               a_genome(animal,chromosome,lgru(j),1)=IBCLR(a_genome(animal,chromosome,lgru(j),1),s_mutat(j))          ! MUTATE
            ELSE                                                                                                ! If bit_pos = 0
               a_genome(animal,chromosome,lgru(j),1)=IBSET(a_genome(animal,chromosome,lgru(j),1),s_mutat(j))         ! MUTATE
            ENDIF
         ENDDO
         DEALLOCATE(lgru, s_mutat, order, temp,one)
      ENDIF ! CLOSING SIRES

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! MUTATION for DAM chromosome
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ran_num_1=ran2(idum1)
      DO nmut = 1, len_mutation                               ! Define number of mutations in this an/chrom
         IF (ran_num_1 < mutation_distr(nmut))EXIT
      ENDDO
      IF (ran_num_1 >= mutation_distr(len_mutation)) nmut=nmut+1          ! If rand_num is > to max distr, add 1 extra mutation

      IF (nmut .GT. 1)THEN;                                   ! DO SOMETHING ONLY IF NMUT > 1. Distribution starts from p(0)
         nmut=nmut-1                                          ! Reduce nmut
         ALLOCATE (d_mutat(nmut), order(nmut), temp(nmut), lgru(nmut),one(1))
         d_mutat=0;order=0;temp=0;lgru=0
         CALL intrand(n_linkage*n_locus, nmut, temp)                 ! Create a vector of (random) positions of NREC recombinations
!         print*,"MUT_DAM",animal,chromosome,temp
         FORALL(j=1:nmut) order(j)=j                          ! This vector is needed by sorting subroutine
         CALL sort_updo (temp, order, nmut, 'low_to_high')    ! Sort positions

         lgru   = INT(FLOOR(temp(:)/n_locus))       ! All blocks where recombination happens along the chromosome
         d_mutat= MOD(INT(temp(:)), n_locus)
         DO j = 1, nmut    ! N_MUTATIONS
            IF(lgru(j).EQ.0) lgru(j)=1
            IF (lgru(j).GE.n_linkage)THEN;
!                  lgru(j)=n_linkage
               CALL intrand(n_locus,1,one)
               d_mutat(j)=one(1)
            ENDIF

            IF(BTEST(a_genome(animal,chromosome,lgru(j),2),d_mutat(j))) THEN                                       ! If bit_pos = 1
               a_genome(animal,chromosome,lgru(j),2)=IBCLR(a_genome(animal,chromosome,lgru(j),2),d_mutat(j))          ! MUTATE
            ELSE                                                                                                ! If bit_pos = 0
               a_genome(animal,chromosome,lgru(j),2)=IBSET(a_genome(animal,chromosome,lgru(j),2),d_mutat(j))         ! MUTATE
            ENDIF
         ENDDO
         DEALLOCATE(lgru, d_mutat, order, temp,one)
      ENDIF ! CLOSING DAMS
   ENDIF ! CLOSING MUTATION
ENDDO ! CLOSING CHROMOSOMES
ENDDO ! CLOSING ANIMALS

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE gen_offspring

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE process_pedigree
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE process_pedigree (n_animal, rp_i_genome, n_ph_sire, n_ph_dam, ph_sire, ph_dam)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE RNG
USE Subs_General
USE Prog_Param,   ONLY: idum1, n_chromosome, n_linkage, n_allele, Ped, sort_sire, sort_dam, cgens, gens

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
INTEGER(i4), INTENT( IN)                                                          :: n_animal, n_ph_sire, n_ph_dam
INTEGER(i4), INTENT( IN), DIMENSION(n_ph_sire,n_chromosome,n_linkage,n_allele)    :: ph_sire
INTEGER(i4), INTENT( IN), DIMENSION(n_ph_dam ,n_chromosome,n_linkage,n_allele)    :: ph_dam
INTEGER(i4), INTENT(OUT), DIMENSION(n_animal ,n_chromosome,n_linkage,n_allele)    :: rp_i_genome
INTEGER(i4),              DIMENSION(1        ,n_chromosome,n_linkage,n_allele)    :: rp_s_genome,rp_d_genome,rp_a_genome
INTEGER(i4),              DIMENSION(n_ph_sire*25)                                 :: phantom_sire ! hard coded?
INTEGER(i4),              DIMENSION(  n_ph_dam*5)                                 :: phantom_dam  ! hard coded?
INTEGER(i4),              DIMENSION(   n_ph_sire)                                 :: conta_sire 
INTEGER(i4),              DIMENSION(    n_ph_dam)                                 :: conta_dam  
INTEGER(i4)                                                                       :: animal, j, y, x, out
INTEGER(i4)                                                                       :: phantom_need
REAL   (dp)                                                                       :: ran_num_1

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
conta_sire=0; conta_dam=0
rp_i_genome=0;
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

! PHANTOM SIRES
phantom_need=n_ph_sire*25

! CREATE A VECTOR OF PHANTOM SIRES FOR ALL PEDIGREE NEEDS, JUST ONCE
j=0
DO 
   ran_num_1=ran1(idum1)
   out = 1 + INT(n_ph_sire*ran_num_1) ! Random phantom sire (animal 0 is not possible)
   IF (conta_sire(out) .EQ. 25) CYCLE
   conta_sire(out)=conta_sire(out)+1
   j=j+1
   phantom_sire(j)=out
   if(j.eq.phantom_need)EXIT
ENDDO

! PHANTOM DAMS
phantom_need=n_ph_dam*5

! CREATE A VECTOR OF PHANTOM DAMS FOR ALL PEDIGREE NEEDS, JUST ONCE

j=0
DO 
   ran_num_1=ran1(idum1)
   out = 1 + INT(n_ph_dam*ran_num_1) ! Random phantom dam (animal 0 is not possible)
   IF (conta_dam(out) .EQ. 5) CYCLE
   conta_dam(out)=conta_dam(out)+1
   j=j+1
   phantom_dam(j)=out
   if(j.eq.phantom_need)EXIT
ENDDO
DEALLOCATE(sort_sire,sort_dam)

!! Erase sort_sire & dam, and create an index pointing to the only animal (useful for gen_offspring)
ALLOCATE(sort_sire(1),sort_dam(1))
sort_sire(1)=1
sort_dam(1) =1

!! PROCESS PEDIGREE
y=0;x=0
DO animal=1,n_animal                                              ! Do it animal by animal. > time, < RAM required
!!!! SIRES
   IF(Ped(animal)%si .EQ. 0)THEN;                                 
      y=y+1
      IF(Ped(animal)%gen .EQ. 1) THEN;                            ! Actually this is GEN0
         rp_s_genome(1,:,:,:)=ph_sire(phantom_sire(y),:,:,:)      ! Pick Phantom sire
         Ped(animal)%si      =n_animal + phantom_sire(y)          ! Keep track of Phantom sire id
      ELSE
         ran_num_1=ran1(idum1)
         out = 1 + INT(cgens((Ped(animal)%gen)-1,1)*ran_num_1)    ! Random new dam  from previous generation
         rp_s_genome(1,:,:,:)=rp_i_genome(gens((Ped(animal)%gen)-1,out,1),:,:,:)           ! Pick new sire from prev gen
         Ped(animal)%si      =gens((Ped(animal)%gen)-1,out,1)     ! Keep track of "new assigned" sire
      ENDIF
   ELSE
      rp_s_genome(1,:,:,:)=rp_i_genome(Ped(animal)%si,:,:,:)      ! Get "Real" sire
   ENDIF

!!!! DAMS
   IF(Ped(animal)%da .EQ. 0)THEN;                                 ! Unknown dam
      x=x+1
      IF(Ped(animal)%gen .EQ. 1) THEN;                            ! Actually this is GEN0
         rp_d_genome(1,:,:,:)=ph_dam(phantom_dam(x),:,:,:)        ! Pick Phantom dam
         Ped(animal)%da      =n_animal + n_ph_sire + phantom_dam(x) ! Keep track of Phantom dam id
      ELSE
         ran_num_1=ran1(idum1)
         out = 1 + INT(cgens((Ped(animal)%gen)-1,2)*ran_num_1)    ! Random new dam  from previous generation
         rp_d_genome(1,:,:,:)=rp_i_genome(gens((Ped(animal)%gen)-1,out,2),:,:,:)           ! Pick new sire from prev gen
         Ped(animal)%da      =gens((Ped(animal)%gen)-1,out,2)     ! Keep track of "new assigned" sire
      ENDIF
   ELSE
      rp_d_genome(1,:,:,:)=rp_i_genome(Ped(animal)%da,:,:,:)      ! Get "Real" dam
   ENDIF


!!!! SONS/DAUGHTERS
   CALL gen_offspring (1, 1, 1, rp_s_genome, rp_d_genome, rp_a_genome,.FALSE.) ! Generate offspring genotype
   rp_i_genome(animal,:,:,:) = rp_a_genome(1,:,:,:)            

ENDDO

END SUBROUTINE process_pedigree

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE calculate_frequencies
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE calculate_frequencies(m_genome,f_genome,o_genotypes)
! This subroutine ONLY counts the number of the two homozygote genotypes, and not the frequency of the two alleles.
! At this stage it also counts one value for the male and female genotypes.

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE Prog_Param, ONLY: len_genome, n_chromosome, n_linkage, n_locus, n_allele, n_male, n_female, homo, ho, cfreq, frequencies

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4), INTENT( IN), DIMENSION(n_male  ,n_chromosome,n_linkage,n_allele) :: m_genome
INTEGER(i4), INTENT( IN), DIMENSION(n_female,n_chromosome,n_linkage,n_allele) :: f_genome
INTEGER(i4),              DIMENSION(len_genome,3)                             :: o_genotypes
INTEGER(i4)                                                                   :: animal, chromosome, linkage, locus
INTEGER(i4)                                                                   :: bit_pos, snp
REAL   (dp)                                                                   :: ffq
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
o_genotypes=0;frequencies=0.D0;homo=0;ho=0;cfreq=0;ffq=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
DO animal=1,n_male
   snp=0
   DO chromosome=1,n_chromosome
      DO linkage=1,n_linkage
         DO locus=1,n_locus
            bit_pos=locus-1
            snp=snp+1
            IF     (      btest(m_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                          btest(m_genome(animal,chromosome,linkage,2),bit_pos))  THEN
                          o_genotypes(snp,1)=o_genotypes(snp,1)+1                        ! One homozygote: AA
            ELSEIF (.NOT. btest(m_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                    .NOT. btest(m_genome(animal,chromosome,linkage,2),bit_pos))  THEN
                          o_genotypes(snp,3)=o_genotypes(snp,3)+1                        ! Other homozygote: BB
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO animal=1,n_female
   snp=0
   DO chromosome=1,n_chromosome
      DO linkage=1,n_linkage
         DO locus=1,n_locus
            bit_pos=locus-1
            snp=snp+1
            IF     (      btest(f_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                          btest(f_genome(animal,chromosome,linkage,2),bit_pos))  THEN
                          o_genotypes(snp,1)=o_genotypes(snp,1)+1                        ! One homozygote: AA
            ELSEIF (.NOT. btest(f_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                    .NOT. btest(f_genome(animal,chromosome,linkage,2),bit_pos))  THEN
                          o_genotypes(snp,3)=o_genotypes(snp,3)+1                        ! Other homozygote: BB
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
DO snp=1, len_genome
   o_genotypes(snp,2) = (n_male + n_female) - (o_genotypes(snp,1) + o_genotypes(snp,3))
   frequencies(snp,1)= REAL(2*(o_genotypes(snp,1)) + o_genotypes(snp,2))/REAL(2*(n_male + n_female))
   frequencies(snp,2)= REAL(2*(o_genotypes(snp,3)) + o_genotypes(snp,2))/REAL(2*(n_male + n_female))

   ho = ho + (REAL(o_genotypes(snp,2))/REAL(n_male + n_female))                     ! Observed heterozygosity 
   IF(frequencies(snp,1) .GT. 0.99 .OR. frequencies(snp,2) .GT. 0.99) homo = homo+1 ! Homozygotes
   ffq=MIN(frequencies(snp,1),frequencies(snp,2))
   IF(ffq .LE. 0.1) THEN
      cfreq(1)=cfreq(1)+1
   ELSEIF(ffq .LE. 0.2) THEN
      cfreq(2)=cfreq(2)+1
   ELSEIF(ffq .LE. 0.3) THEN
      cfreq(3)=cfreq(3)+1
   ELSEIF(ffq .LE. 0.4) THEN
      cfreq(4)=cfreq(4)+1
   ELSE
      cfreq(5)=cfreq(5)+1
   ENDIF
ENDDO

ho = ho / REAL(len_genome)                                                          ! Ho
cfreq = REAL(cfreq) / REAL(len_genome)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE calculate_frequencies

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




SUBROUTINE calculate_frequenciesPED(i_genome,n_animal)

! This subroutine ONLY counts the number of the two homozygote genotypes, and not the frequency of the two alleles.
! At this stage it also counts one value for the male and female genotypes.

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
USE Prog_Param, ONLY: len_genome, n_chromosome, n_linkage, n_locus, n_allele, frequencies, homo, ho, cfreq, pedOUT

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4), INTENT( IN)                                                      :: n_animal
INTEGER(i4), INTENT( IN),DIMENSION( n_animal,n_chromosome,n_linkage,n_allele) :: i_genome
!REAL   (dp), INTENT(OUT),DIMENSION(len_genome,2)                              :: frequencies
INTEGER(i4),             DIMENSION(len_genome,3)                              :: o_genotypes
INTEGER(i4)                                                                   :: animal
INTEGER(i4)                                                                   :: chromosome, linkage, locus, bit_pos, snp
REAL   (dp)                                                                   :: ffq

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
o_genotypes=0;frequencies=0.D0;ho=0;cfreq=0;homo=0;ffq=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO animal=pedOUT+1,n_animal
snp=0
 DO chromosome=1,n_chromosome
  DO linkage=1,n_linkage
   DO locus=1,n_locus
    bit_pos=locus-1
    snp=snp+1
    IF     (      btest(i_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                  btest(i_genome(animal,chromosome,linkage,2),bit_pos))  THEN
                    o_genotypes(snp,1)=o_genotypes(snp,1)+1                                           ! One homozygote: AA
    ELSEIF (.NOT. btest(i_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
            .NOT. btest(i_genome(animal,chromosome,linkage,2),bit_pos))  THEN
                    o_genotypes(snp,3)=o_genotypes(snp,3)+1                                         ! Other homozygote: BB
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO snp=1, len_genome
  o_genotypes(snp,2) = (n_animal-pedOUT) - (o_genotypes(snp,1) + o_genotypes(snp,3))
  frequencies(snp,1)= REAL(2*(o_genotypes(snp,1)) + o_genotypes(snp,2))/REAL(2*(n_animal-pedOUT))
  frequencies(snp,2)= REAL(2*(o_genotypes(snp,3)) + o_genotypes(snp,2))/REAL(2*(n_animal-pedOUT))
  ho = ho + (REAL(o_genotypes(snp,2))/REAL(n_animal-pedOUT))                     ! Observed heterozygosity
  IF(frequencies(snp,1) .GT. 0.999 .OR. frequencies(snp,2) .GT. 0.999) homo = homo+1 ! Homozygotes
  ffq=MIN(frequencies(snp,1),frequencies(snp,2))
  IF(ffq .LE. 0.1) THEN
     cfreq(1)=cfreq(1)+1
  ELSEIF(ffq .LE. 0.2) THEN
     cfreq(2)=cfreq(2)+1
  ELSEIF(ffq .LE. 0.3) THEN
     cfreq(3)=cfreq(3)+1
  ELSEIF(ffq .LE. 0.4) THEN
     cfreq(4)=cfreq(4)+1
  ELSE
     cfreq(5)=cfreq(5)+1
  ENDIF
ENDDO

ho = ho / REAL(len_genome)                                                          ! Ho
cfreq = REAL(cfreq) / REAL(len_genome)


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE calculate_frequenciesPED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120






!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE get_VG
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE get_Vg(Vg)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
USE Numeric_Kinds
USE Prog_Param, ONLY: frequencies, len_genome, n_trait, qtl_position, qtl_size, fixed, HmaxQTL

IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
REAL   (dp), INTENT(OUT), DIMENSION(n_trait)                                  :: Vg
INTEGER(i4)                                                                   :: j, qtp, snp
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

Vg=0.D0;
fixed=0;qtp=1

DO snp=1, len_genome
  IF (snp .EQ. qtl_position(qtp)) THEN
     FORALL(j=1:n_trait) Vg(j)=Vg(j)+2*frequencies(snp,1)*frequencies(snp,2)*(qtl_size(qtp,j)*qtl_size(qtp,j))
     IF ( frequencies(snp,1)==1. .OR. frequencies(snp,2)==1.) fixed(1)=fixed(1)+1
     IF (qtp .LT. HmaxQTL) qtp=qtp+1
  ELSE
     IF ( frequencies(snp,1)==1. .OR. frequencies(snp,2)==1.) fixed(2)=fixed(2)+1
  ENDIF
ENDDO

END SUBROUTINE get_Vg

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE generate_phenotype (n_animal,s_genome,s_phenotype)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE generate_phenotype (n_animal,a_genome,a_phenotype,t_phenotype)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
!USE Subs_General
USE RNG
USE Prog_Param, ONLY: idum1, trait_n, n_chromosome,n_linkage,n_allele, snp, chromosome, linkage, locus, &
                      n_locus, bit_pos, i, animal, frequencies, h2, qtl_size, qtl_position, &
                      trait_n, HmaxQTL

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
INTEGER(i4), INTENT( IN)                                                      :: n_animal
INTEGER(i4), INTENT( IN), DIMENSION(n_animal,n_chromosome,n_linkage,n_allele) :: a_genome
REAL   (dp), INTENT(OUT), DIMENSION(n_animal                                ) :: a_phenotype,t_phenotype
REAL   (dp)                                                                   :: ran_num_1
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
a_phenotype=0.
t_phenotype=0.
snp=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO chromosome=1,n_chromosome
 DO linkage=1,n_linkage
  DO locus=1,n_locus
   bit_pos=locus-1
   snp=snp+1
     DO i=1, HmaxQTL
      IF(qtl_position(i) == snp) THEN
       DO animal=1, n_animal
        IF     (      btest(a_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                      btest(a_genome(animal,chromosome,linkage,2),bit_pos))  THEN
!          a_phenotype(animal) = a_phenotype(animal) + (( 2.D0)*frequencies(snp,2)*qtl_size(i,trait_n))
          a_phenotype(animal) = a_phenotype(animal) +  (( 2.D0)*                   qtl_size(i,trait_n))
    ELSEIF     (.NOT. btest(a_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                .NOT. btest(a_genome(animal,chromosome,linkage,2),bit_pos))  THEN
!          a_phenotype(animal) = a_phenotype(animal) + ((-2.D0)*frequencies(snp,1)*qtl_size(i,trait_n))
          a_phenotype(animal) = a_phenotype(animal)  + ((-2.D0)*                   qtl_size(i,trait_n))
      ELSE
!          a_phenotype(animal) = a_phenotype(animal) + ((frequencies(snp,2)-frequencies(snp,1))*qtl_size(i,trait_n))
          a_phenotype(animal) = a_phenotype(animal) +  qtl_size(i,trait_n)
      ENDIF
     ENDDO
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

t_phenotype=a_phenotype
DO animal=1, n_animal
   ran_num_1=gasdev(idum1)
   a_phenotype(animal) = t_phenotype(animal) + ran_num_1*SQRT((1.-h2(trait_n)))
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE generate_phenotype



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE generate_phenotype (n_animal,s_genome,s_phenotype)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE generate_phenotypeALL (n_animal,a_genome,i_phenotype,it_phenotype,corr)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
!USE Subs_General
USE RNG
USE Prog_Param, ONLY: idum1, trait_n, n_chromosome,n_linkage,n_allele, snp, chromosome, linkage, locus, &
                      n_locus, bit_pos, i, animal, frequencies, h2, n_trait, qtl_size, qtl_position, trait_n, &
                      HmaxQTL

IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
INTEGER(i4), INTENT( IN)                                                      :: n_animal
INTEGER(i4), INTENT( IN), DIMENSION(n_animal,n_chromosome,n_linkage,n_allele) :: a_genome
REAL   (dp), INTENT(OUT), DIMENSION(n_animal,n_trait                        ) :: i_phenotype,it_phenotype
REAL   (dp), INTENT(OUT), DIMENSION(n_trait                                 ) :: corr
REAL   (dp)                                                                   :: ran_num_1
INTEGER(i4)                                                                   :: trait
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
i_phenotype=0.D0
it_phenotype=0.D0
snp=0;corr=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO chromosome=1,n_chromosome
 DO linkage=1,n_linkage
  DO locus=1,n_locus
   bit_pos=locus-1
   snp=snp+1
   DO i=1, HmaxQTL
    IF(qtl_position(i) == snp) THEN
      DO trait=1,n_trait 
       DO animal=1, n_animal
         IF     (      btest(a_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                       btest(a_genome(animal,chromosome,linkage,2),bit_pos))  THEN
!           i_phenotype(animal,trait) = i_phenotype(animal,trait) + (( 2.D0)*frequencies(snp,2)*qtl_size(i,trait))
           i_phenotype(animal,trait) = i_phenotype(animal,trait) + (( 2.D0)*qtl_size(i,trait))
      ELSEIF     (.NOT. btest(a_genome(animal,chromosome,linkage,1),bit_pos) .AND. &
                  .NOT. btest(a_genome(animal,chromosome,linkage,2),bit_pos))  THEN
!           i_phenotype(animal,trait) = i_phenotype(animal,trait)+((-2.D0)*frequencies(snp,1)*qtl_size(i,trait))
           i_phenotype(animal,trait) = i_phenotype(animal,trait)+((-2.D0)*qtl_size(i,trait))
       ELSE
!           i_phenotype(animal,trait) = i_phenotype(animal,trait)+((frequencies(snp,2)-frequencies(snp,1))*qtl_size(i,trait))
           i_phenotype(animal,trait) = i_phenotype(animal,trait)+qtl_size(i,trait)
       ENDIF
      ENDDO
     ENDDO
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

it_phenotype = i_phenotype
DO trait=1,n_trait
   DO animal=1, n_animal
      ran_num_1=gasdev(idum1)
      i_phenotype(animal,trait) = it_phenotype(animal,trait) + ran_num_1*SQRT((1.-h2(trait_n)))
   ENDDO
   CALL get_corr(it_phenotype(:,trait),i_phenotype(:,trait), n_animal, corr(trait))
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE generate_phenotypeALL

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE get_corr
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE get_corr(vect1,vect2,anims,r)
USE Numeric_Kinds

INTEGER ,INTENT(IN)                 ::anims
REAL(dp),INTENT(IN),DIMENSION(anims)::vect1,vect2
REAL(dp)                            ::EX,EXsq,EY,EYsq,EXY,r

EX  =sum(vect1)   
EXsq=sum(vect1**2)
EY  =sum(vect2)
EYsq=sum(vect2**2)
EXY =dot_product(vect1,vect2)
r=((anims*EXY)-(EX*EY)) /  sqrt (((anims*EXsq)-(EX**2)) * ((anims*EYsq)-(EY**2)))

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE get_corr

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



SUBROUTINE diverge_breeds(s_genome,d_genome,breeds_m_genome,breeds_f_genome)!,new_sort_sire,new_sort_dam)
! EXPLAIN SOMETHING HERE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
USE Numeric_Kinds
USE Subs_General
USE Prog_param, ONLY: species_male, species_female, n_chromosome, n_linkage, n_allele, Nbreeds, n_bre_male, &
                      n_bre_female,sort_sire,sort_dam,admixprop
IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4),             DIMENSION(        species_male,n_chromosome,n_linkage,n_allele) :: s_genome
INTEGER(i4),             DIMENSION(      species_female,n_chromosome,n_linkage,n_allele) :: d_genome
INTEGER(i4),             DIMENSION(  Nbreeds,n_bre_male,n_chromosome,n_linkage,n_allele) :: breeds_m_genome
INTEGER(i4),             DIMENSION(Nbreeds,n_bre_female,n_chromosome,n_linkage,n_allele) :: breeds_f_genome
!INTEGER(i4), INTENT(OUT),DIMENSION(                                  Nbreeds,n_bre_male) :: new_sort_sire
!INTEGER(i4), INTENT(OUT),DIMENSION(                                Nbreeds,n_bre_female) :: new_sort_dam

INTEGER(i4) :: i, j, k, ad, mADX
INTEGER(i4), ALLOCATABLE, DIMENSION (:):: sampleposM, sampleposF, skip

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!!! Admixture max (for 2 breeds only!)
mADX = int(n_bre_male*admixprop)/100.
print*, "Max number of admixed animals:",mADX,",from the total:",n_bre_male


allocate(sampleposM(species_male),sampleposF(species_male),skip(species_male))
sampleposM=0;sampleposF=0;skip=0;

call samp_WR(1,species_male  ,mADX,sampleposM)   ! Sample mADX int values, from 1 to species_male, put in sampleposM
call samp_WR(1,species_female,mADX,sampleposF)   ! Sample mADX int values, from 1 to species_female, put in sampleposF

!! Quick and dirty solution for NONE animals admixed
IF (admixprop.EQ.0)THEN;
   sampleposM=99999999
   sampleposF=99999999
ENDIF

k=0
j=1
ad=1
!! FIRST RUN ADMIXTURE ANIMALS
!! Male
open(567,file='RESULTS/admixed_anims.txt')
DO i=1,species_male
   IF (i .EQ. sampleposM(ad)) THEN
      print*, "ADMIXED M ANIM",j," in position:",i,',now position:',ad,',corresp. to ID:',sort_sire(i)
      skip(i)=1
      write(567,*)"Bulls in position:",j,"are admixed ancestors"
      breeds_m_genome(1,j,:,:,:)=s_genome(sort_sire(i),:,:,:)
      breeds_m_genome(2,j,:,:,:)=s_genome(sort_sire(i),:,:,:)
      j = j + 1
      IF (ad .LT. mADX) ad = ad + 1
   ENDIF
ENDDO

!! Quick and dirty solytion for 100% admixed pops
IF (admixprop.LT.100) THEN
   DO i=1,species_male
      if (skip(i).EQ.1)CYCLE
      k=k+1
      print*, "NOT ADMIXED M ANIM",j," in position:",i,'for breed:',k,'corresp. to ID:',sort_sire(i)
      breeds_m_genome(k,j,:,:,:)=s_genome(sort_sire(i),:,:,:)
      IF (k.EQ.Nbreeds)THEN;
         k=0
         j=j+1
      ENDIF
      IF (j.eq.n_bre_male+1)EXIT
   ENDDO
ENDIF
deallocate(skip)

!!!!!!!!!!!!!!! FEMALE POP !!!!!!!!!!!!!!!!!!!!
allocate(skip(species_female))
skip=0
k=0
j=1
ad=1
!! FIRST RUN ADMIXTURE ANIMALS
DO i=1,species_female
   IF (i .EQ. sampleposF(ad)) THEN
      print*, "ADMIXED F ANIM",j," in position:",i,',now position:',ad,'corresp. to ID:',sort_dam(i)
      skip(i)=1
      write(567,*)"Cows in position:",n_bre_male+j,"are admixed ancestors"
      breeds_f_genome(1,j,:,:,:)=d_genome(sort_dam(i),:,:,:)
      breeds_f_genome(2,j,:,:,:)=d_genome(sort_dam(i),:,:,:)
      j = j + 1
      IF (ad .LT. mADX)ad = ad + 1
   ENDIF
ENDDO

!! Quick and dirty solytion for 100% admixed pops
IF(admixprop.LT.100)THEN
   DO i=1,species_female
      if (skip(i).EQ.1)CYCLE
      k=k+1
      print*, "NOT ADMIXED F ANIM",j," in position:",i,'for breed:',k,'corresp. to ID:',sort_dam(i)
      breeds_f_genome(k,j,:,:,:)=d_genome(sort_dam(i),:,:,:)
      IF (k.EQ.Nbreeds)THEN;
         k=0
         j=j+1
      ENDIF
      IF (j.eq.n_bre_female+1)EXIT
   ENDDO
ENDIF

END SUBROUTINE diverge_breeds

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE expand_genome_breed
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE expand_genome_breed (s_genome,d_genome,m_genome,f_genome)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
!USE Subs_General
USE RNG
USE Prog_Param, ONLY: n_chromosome, n_linkage, n_allele, n_bre_male, &
                      n_bre_female,n_male,n_female,sort_sire,sort_dam

IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4), DIMENSION(  n_bre_male,n_chromosome,n_linkage,n_allele) :: s_genome
INTEGER(i4), DIMENSION(n_bre_female,n_chromosome,n_linkage,n_allele) :: d_genome
INTEGER(i4), DIMENSION(      n_male,n_chromosome,n_linkage,n_allele) :: m_genome
INTEGER(i4), DIMENSION(    n_female,n_chromosome,n_linkage,n_allele) :: f_genome

INTEGER  :: i, j, k
REAL     :: m_ratio,f_ratio

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
m_genome=0;f_genome=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

m_ratio=  n_male/n_bre_male
f_ratio=n_female/n_bre_female
   
!! Males
IF (MOD(n_male,n_bre_male)==0)THEN
   k=0
   DO i=1,n_bre_male
      DO j=1,INT(m_ratio)
         k=k+1
         m_genome(k,:,:,:) = s_genome(i,:,:,:)
      ENDDO
   ENDDO
ELSE
   k=0
   DO i=1,n_bre_male
      DO j=1,INT(m_ratio)
         k=k+1
         m_genome(k,:,:,:) = s_genome(i,:,:,:)
      ENDDO
   ENDDO
   DO j=1,(n_male-n_bre_male*INT(m_ratio))
      k=k+1
      !If ratio is not xN, then is more likely that best animals have an extra copy in next gen.
      m_genome(k,:,:,:) = s_genome(j,:,:,:) 
   ENDDO
ENDIF

!! Females
IF (MOD(n_female,n_bre_female)==0)THEN
   K=0
   DO i=1,n_bre_female
      DO j=1,INT(f_ratio)
         k=k+1
         f_genome(k,:,:,:) = d_genome(i,:,:,:)
      ENDDO
   ENDDO
ELSE
   k=0
   DO i=1,n_bre_female
      DO j=1,INT(f_ratio)
         k=k+1
         f_genome(k,:,:,:) = d_genome(i,:,:,:)
      ENDDO
   ENDDO
   DO j=1,(n_female-n_bre_female*INT(f_ratio))
      k=k+1
      !If ratio is not xN, then is more likely that best animals have an extra copy in next gen.
      f_genome(k,:,:,:) = d_genome(j,:,:,:) 
   ENDDO
ENDIF

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
END SUBROUTINE expand_genome_breed

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE expand_genome_country
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE expand_genome_country (s_genome,d_genome,m_genome,f_genome)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE Numeric_Kinds
!USE Subs_General
USE RNG
USE Prog_Param, ONLY: n_chromosome, n_linkage,n_allele, n_male, &
                      n_female,n_bre_male,n_bre_female,sort_sire, sort_dam

IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4), DIMENSION(    n_bre_male,n_chromosome,n_linkage,n_allele) :: s_genome
INTEGER(i4), DIMENSION(  n_bre_female,n_chromosome,n_linkage,n_allele) :: d_genome
INTEGER(i4), DIMENSION(        n_male,n_chromosome,n_linkage,n_allele) :: m_genome
INTEGER(i4), DIMENSION(      n_female,n_chromosome,n_linkage,n_allele) :: f_genome

INTEGER  :: i, j, k
REAL     :: m_ratio,f_ratio
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

m_ratio=  n_male/n_bre_male
f_ratio=n_female/n_bre_female

!! Males
IF (MOD(n_male,n_bre_male)==0)THEN
   k=0
   DO i=1,n_bre_male
      DO j=1,INT(m_ratio)
         k=k+1
         m_genome(k,:,:,:) = s_genome(i,:,:,:)
      ENDDO
   ENDDO
ELSE
   k=0
   DO i=1,n_bre_male
      DO j=1,INT(m_ratio)
         k=k+1
         m_genome(k,:,:,:) = s_genome(i,:,:,:)
      ENDDO
   ENDDO
   DO j=1,(n_male-n_bre_male*INT(m_ratio))
      k=k+1
      !If ratio is not xN, then is more likely that best animals have an extra copy in next gen.
      m_genome(k,:,:,:) = s_genome(sort_sire(j),:,:,:) 
   ENDDO
ENDIF

!! Females
IF (MOD(n_female,n_bre_female)==0)THEN
   K=0
   DO i=1,n_bre_female
      DO j=1,INT(f_ratio)
         k=k+1
         f_genome(k,:,:,:) = d_genome(i,:,:,:)
      ENDDO
   ENDDO
ELSE
   k=0
   DO i=1,n_bre_female
      DO j=1,INT(f_ratio)
         k=k+1
         f_genome(k,:,:,:) = d_genome(i,:,:,:)
      ENDDO
   ENDDO
   DO j=1,(n_female-n_bre_female*INT(f_ratio))
      k=k+1
      !If ratio is not xN, then is more likely that best animals have an extra copy in next gen.
      f_genome(k,:,:,:) = d_genome(sort_dam(j),:,:,:) 
   ENDDO
ENDIF

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE expand_genome_country

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE decompress
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE decompress(in_genome,out_genome)
USE Numeric_Kinds
USE Prog_Param, ONLY: len_genome, n_chromosome, n_linkage, n_locus, n_allele

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! History:
!!  - Original Program:   ELNicolazzi - Keep haplotype info (dunno if used)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER(i4), INTENT(IN),  DIMENSION(1,n_chromosome,n_linkage,n_allele)         ::  in_genome
INTEGER(i1), INTENT(OUT), DIMENSION(len_genome,n_allele)                       :: out_genome
INTEGER(i4):: chromosome,linkage,locus,snp,bit_pos
!--------1---------2---------3---------4---------5---------6---------7---------8
out_genome=0
snp=0
!--------1---------2---------3---------4---------5---------6---------7---------8
DO chromosome=1,n_chromosome
   DO linkage=1,n_linkage
      DO locus = 1, n_locus                                
         bit_pos=locus-1
         snp=snp+1
         IF(BTEST(in_genome(1,chromosome,linkage,1),bit_pos))& !BIT=1 -> SNP=1, otherwise 0 
              out_genome(snp,1)=1                              !Paternal haplotype
         IF(BTEST(in_genome(1,chromosome,linkage,2),bit_pos))&
              out_genome(snp,2)=1                              !Maternal haplotype
      ENDDO
   ENDDO
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE decompress

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE decompress2(in_genome,out_genome)
USE Numeric_Kinds
USE Prog_Param, ONLY: len_genome, n_chromosome, n_linkage, n_locus, n_allele

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! History:
!!  - Original Program:   ELNicolazzi - Keep haplotype info (dunno if used)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER(i4), INTENT(IN),  DIMENSION(1,n_chromosome,n_linkage,n_allele)         ::  in_genome
INTEGER(i1), INTENT(OUT), DIMENSION(len_genome)                                :: out_genome
INTEGER(i1), DIMENSION(len_genome,n_allele)                       :: t_genome
INTEGER(i4):: chromosome,linkage,locus,snp,bit_pos
!--------1---------2---------3---------4---------5---------6---------7---------8
out_genome=9
snp=0
!--------1---------2---------3---------4---------5---------6---------7---------8
DO chromosome=1,n_chromosome
   DO linkage=1,n_linkage
      DO locus = 1, n_locus
         bit_pos=locus-1
         snp=snp+1
         !BIT=1 -> SNP=1, otherwise 0
         IF(BTEST(in_genome(1,chromosome,linkage,1),bit_pos))THEN;
            t_genome(snp,1)=1   !Paternal haplotype
         ELSE
            t_genome(snp,1)=0
         ENDIF
         IF(BTEST(in_genome(1,chromosome,linkage,2),bit_pos))THEN;
            t_genome(snp,2)=1   !Maternal haplotype
         ELSE
            t_genome(snp,2)=0
         ENDIF
         out_genome(snp)=t_genome(snp,1)+t_genome(snp,2)
      ENDDO
   ENDDO
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100----110-------120

END SUBROUTINE decompress2



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE call_breed_trait
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE call_breed_trait(abreed)
USE Prog_Param, ONLY: selez,trait_n

INTEGER, INTENT(IN)   :: abreed
INTEGER               :: i
!--------1---------2---------3---------4---------5---------6---------7---------8                                         

i=abreed
IF(selez(i).EQ.'H1_A') trait_n=1
IF(selez(i).EQ.'H1_B') trait_n=2
IF(selez(i).EQ.'H1_C') trait_n=3
IF(selez(i).EQ.'H2_A') trait_n=4
IF(selez(i).EQ.'H2_B') trait_n=5
IF(selez(i).EQ.'H2_C') trait_n=6
IF(selez(i).EQ.'H3_A') trait_n=7
IF(selez(i).EQ.'H3_B') trait_n=8
IF(selez(i).EQ.'H3_C') trait_n=9
IF(selez(i).EQ.'L1_A') trait_n=10
IF(selez(i).EQ.'L1_B') trait_n=11
IF(selez(i).EQ.'L1_C') trait_n=12
IF(selez(i).EQ.'L2_A') trait_n=13
IF(selez(i).EQ.'L2_B') trait_n=14
IF(selez(i).EQ.'L2_C') trait_n=15
IF(selez(i).EQ.'L3_A') trait_n=16
IF(selez(i).EQ.'L3_B') trait_n=17
IF(selez(i).EQ.'L3_C') trait_n=18

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE call_breed_trait

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE call_country_trait
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE call_country_trait(abreed)
USE Prog_Param, ONLY: selez,trait_n
INTEGER, INTENT(IN)   :: abreed
INTEGER               :: i
!--------1---------2---------3---------4---------5---------6---------7---------8                                         
i=abreed
IF(selez(i).EQ.'H1_A') trait_n=1
IF(selez(i).EQ.'H1_B') trait_n=2
IF(selez(i).EQ.'H1_C') trait_n=3
IF(selez(i).EQ.'H2_A') trait_n=4
IF(selez(i).EQ.'H2_B') trait_n=5
IF(selez(i).EQ.'H2_C') trait_n=6
IF(selez(i).EQ.'H3_A') trait_n=7
IF(selez(i).EQ.'H3_B') trait_n=8
IF(selez(i).EQ.'H3_C') trait_n=9

IF(selez(i).EQ.'L1_A') trait_n=10
IF(selez(i).EQ.'L1_B') trait_n=11
IF(selez(i).EQ.'L1_C') trait_n=12
IF(selez(i).EQ.'L2_A') trait_n=13
IF(selez(i).EQ.'L2_B') trait_n=14
IF(selez(i).EQ.'L2_C') trait_n=15
IF(selez(i).EQ.'L3_A') trait_n=16
IF(selez(i).EQ.'L3_B') trait_n=17
IF(selez(i).EQ.'L3_C') trait_n=18

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE call_country_trait

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_frequency
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_frequency
USE Prog_Param

INTEGER :: x
qtp=1
OPEN(100, FILE=TRIM(path)//outname )

WRITE(100,*)'FIXED SNPs:',fixed(2),'FIXED QTLs:',fixed(1);WRITE(100,*)

x=0
DO snp=1, len_genome
   IF(snp .EQ. qtl_position(qtp))THEN
      x=x+1
      WRITE(100,1000) x, sep, snp, sep, frequencies(snp,1), sep, frequencies(snp,2), sep, 'Q',qtl_type(qtp,:)
      IF (qtp .LT. HmaxQTL) qtp=qtp+1
   ELSE
      IF(frequencies(snp,1) .LT. thresh .OR. frequencies(snp,2) .LT. thresh) CYCLE
      x=x+1
      WRITE(100,1001) x, sep, snp, sep, frequencies(snp,1), sep, frequencies(snp,2), sep, 'M'
   ENDIF
ENDDO
CLOSE(100)

1000 format (I0,A,I0,A,2(F0.3,A),A1,A,6A)
1001 format (I0,A,I0,A,2(F0.3,A),A)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_frequency

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_var
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_var
USE Prog_Param
  
OPEN(112, FILE=TRIM(path)//outname)
WRITE(112,'(A)') 'Trait,SP_first,SP_end,CO_end'
DO i=1, n_trait
   WRITE(112,'(I0,A)',ADVANCE='No') i,sep
   DO j=1,2
      WRITE(112,'(F0.5,A)',ADVANCE='No')g_varALL(i,j),sep
   ENDDO
   WRITE(112,'(F0.5)')g_varALL(i,3)
ENDDO

CLOSE(112)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_var

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_genotypePED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_genotypePED
USE Prog_Param

INTEGER                         ::x
INTEGER,ALLOCATABLE,DIMENSION(:)::outgen,outgen2
CHARACTER*50                    ::outfmt2

ALLOCATE(out_m_genome(len_genome, n_allele))
ALLOCATE(outgen(len_genome),outgen2(len_genome))
outgen=9;outgen2=9
outfmt2='(I0,A1,I0,A1,I0,A1,I0,A1,I0,A1,A2)'

OPEN(101, FILE=TRIM(path)//outname)
OPEN(122, FILE=TRIM(path)//outname2)

DO i=pedOUT+1, pedig_n
   !Do animal by animal to avoid using huge RAM
   CALL decompress(rp_i_genome(i,:,:,:),out_m_genome)

   !! HAPLOTYPES
   IF (printgeno.eq.1) THEN     
      x=0
      qtp=1
      DO j=1, len_genome
         IF(j .EQ. qtl_position(qtp))THEN
            x=x+1
            outgen(x)=out_m_genome(j,1)
            outgen2(x)=out_m_genome(j,2) 
            IF (qtp .LT. HmaxQTL) qtp=qtp+1
         ELSE
            IF(frequencies(j,1) .LT. thresh .OR. frequencies(j,2) .LT. thresh) CYCLE
            x=x+1
            outgen(x)=out_m_genome(j,1)
            outgen2(x)=out_m_genome(j,2) 
         ENDIF
      ENDDO

      outfmt='(I10,1X,                   I1)'
      WRITE(outfmt(10:20),'(I0)') x
      WRITE(101,outfmt)  i,outgen(1:x)
      WRITE(101,outfmt)  i,outgen2(1:x)
      WRITE(122,outfmt2) i,sep,PED(i)%an,sep,Ped(i)%si,sep,Ped(i)%da,sep,Ped(i)%sex,sep,'HT'

   !!! GENOTYPES
   ELSE
      outgen=9
      x=0
      qtp=1
      DO j=1, len_genome
         IF(j .EQ. qtl_position(qtp))THEN
            x=x+1
            outgen(x)=out_m_genome(j,1)+out_m_genome(j,2)
            IF (qtp .LT. HmaxQTL) qtp=qtp+1
         ELSE
            IF(frequencies(j,1) .LT. thresh .OR. frequencies(j,2) .LT. thresh) CYCLE
            x=x+1
            outgen(x)=out_m_genome(j,1)+out_m_genome(j,2)
         ENDIF
      ENDDO
      outfmt='(I10,1X,                   I1)'
      WRITE(outfmt(10:20),'(I0)') x
      WRITE(101,outfmt)  i,outgen(1:x)
      WRITE(122,outfmt2) i,sep,Ped(i)%an,sep,Ped(i)%si,sep,Ped(i)%da,sep,Ped(i)%sex,sep,'GT'
   ENDIF
ENDDO
DEALLOCATE(out_m_genome)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_genotypePED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_phenotypePED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_phenotypePED
USE Prog_Param

!! Open and create header
OPEN(102, FILE=TRIM(path)//outname)
WRITE(102,'(A12)',ADVANCE='No') 'INDIVIDUAL,T'
DO i=1,n_trait-1
   WRITE(102,'(I0,A4,F0.2,A3)',ADVANCE='No') i,'(h2=' , h2(i) , '),T'
ENDDO
WRITE(102,'(I0,A4,F0.2,A1)') n_trait,'(h2=' , h2(n_trait) , ')'

DO j=1,pedig_n
   WRITE(102,'(I10,A1)',ADVANCE='No')j,sep
   DO i=1,n_trait-1
      WRITE(102,'(F0.6,A1)',ADVANCE='No')i_phenotype(j,i),sep
   ENDDO
   WRITE(102,'(F0.6)')i_phenotype(j,n_trait)
  
ENDDO
CLOSE(102)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_phenotypePED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_tbvPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_tbvPED
USE Prog_Param

OPEN(103, FILE=TRIM(path)//outname)

!! Open and create header
OPEN(103, FILE=TRIM(path)//outname)
WRITE(103,'(A12)',ADVANCE='No') 'INDIVIDUAL,T'
DO i=1,n_trait-1
   WRITE(103,'(I0,A4,F0.2,A3)',ADVANCE='No') i,'(h2=' , h2(i) , '),T'
ENDDO
WRITE(103,'(I0,A4,F0.2,A1)') n_trait,'(h2=' , h2(n_trait) , ')'

DO j=1,pedig_n
   WRITE(103,'(I0,A1)',ADVANCE='No')j,sep
   DO i=1,n_trait-1
      WRITE(103,'(F0.6,A1)',ADVANCE='No')it_phenotype(j,i),sep
   ENDDO
   WRITE(103,'(F0.6)')it_phenotype(j,n_trait)
ENDDO

CLOSE(103)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_tbvPED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120





!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_corrPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_corrPED
USE Prog_Param

OPEN(104, FILE=TRIM(path)//outname)
WRITE(104,'(A)') 'Trait_#,h2,r(tbv,phen)'
DO i=1, n_trait
   WRITE(104,'(I0,A,F0.2,A,F0.4)') i,sep,h2(i),sep,corr(i)
ENDDO
CLOSE(104)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_corrPED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_genotypeNOPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

SUBROUTINE write_genotypeNOPED
USE Prog_Param
INTEGER                           ::x
INTEGER, ALLOCATABLE, DIMENSION(:)::outgen,outgen2
CHARACTER*50                      ::outfmt2,outfmt3

ALLOCATE(out_m_genome(len_genome, n_allele))
ALLOCATE(out_f_genome(len_genome, n_allele))
ALLOCATE(outgen(len_genome),outgen2(len_genome))

outgen=9;outgen2=9
outfmt2='(I0,1A,I0,1A,I0,1A,I0,1A,I1,1A,A2)'
outfmt3='(I0,1A,I0,5A1,I1,1A,A2)'

DO i=1, n_male
   !Do animal by animal to avoid using huge RAM
   CALL decompress(m_genome(i,:,:,:),out_m_genome)
!! HAPLOTYPES
   IF (printgeno.eq.1 .or. printgeno.eq.3) THEN
      x=0
      qtp=1
      DO j=1, len_genome
         IF(j .EQ. qtl_position(qtp))THEN
            x=x+1
            outgen(x)=out_m_genome(j,1)
            outgen2(x)=out_m_genome(j,2)
            IF (qtp .LT. HmaxQTL) qtp=qtp+1
         ELSE
            IF(frequencies(j,1) .LT. thresh .OR. frequencies(j,2) .LT. thresh) CYCLE
            x=x+1
            outgen(x)=out_m_genome(j,1)
            outgen2(x)=out_m_genome(j,2)
         ENDIF
      ENDDO
      outfmt='(I10,1X,                   I1)'
      WRITE(outfmt(10:20),'(I0)') x
      if (printgeno.eq.3)then
         IF (breed.eq.1) then;
            WRITE(1011,outfmt)  seq+i*2-1,outgen(1:x)
            WRITE(1011,outfmt)  seq+i*2-1,outgen2(1:x)
         ELSE
            WRITE(1011,outfmt)  seq+i*2,outgen(1:x)
            WRITE(1011,outfmt)  seq+i*2,outgen2(1:x)
         ENDIF
      else
         IF (breed.eq.1) then;
            WRITE(101,outfmt)  seq+i*2-1,outgen(1:x)
            WRITE(101,outfmt)  seq+i*2-1,outgen2(1:x)
         ELSE
            WRITE(101,outfmt)  seq+i*2,outgen(1:x)
            WRITE(101,outfmt)  seq+i*2,outgen2(1:x)
         ENDIF
         IF(gen.GT.genOUT-1)THEN
            IF (breed.eq.1) then;
               WRITE(122,outfmt2) seq+i*2-1,',',gen,sep,(seq-n_male*4)+(sort_sire1(i)*2-1),sep,&
                    (seq-n_male*4)+(n_male*2+(sort_dam1(i)*2-1)),sep,1,sep,'HT'
            ELSE
               WRITE(122,outfmt2) seq+i*2,',',gen,sep,(seq-n_male*4)+(sort_sire1(i)*2),sep,&
                    (seq-n_male*4)+(n_male*2+(sort_dam1(i)*2)),sep,1,sep,'HT'
            ENDIF
            
         ELSE
            IF (breed.eq.1) then;
               WRITE(122,outfmt3) i*2-1,',',gen,sep,'0',sep,'0',sep,1,sep,'HT'
            else
               WRITE(122,outfmt3) i*2,',',gen,sep,'0',sep,'0',sep,1,sep,'HT'
            endif
         endif
      ENDIF
!! GENOTYPES
   ENDIF
   IF (printgeno.ge.2) THEN
      outgen=9
      x=0
      qtp=1
      DO j=1, len_genome
         IF(j .EQ. qtl_position(qtp))THEN
            x=x+1
            outgen(x)=out_m_genome(j,1)+out_m_genome(j,2)
            IF (qtp .LT. HmaxQTL) qtp=qtp+1
         ELSE
            IF(frequencies(j,1) .LT. thresh .OR. frequencies(j,2) .LT. thresh) CYCLE
            x=x+1
            outgen(x)=out_m_genome(j,1)+out_m_genome(j,2)
         ENDIF
      ENDDO
      outfmt='(I10,1X,                   I1)'
      WRITE(outfmt(10:20),'(I0)') x
      IF (breed.eq.1) then;
         WRITE(101,outfmt)  seq+i*2-1,outgen(1:x)
      ELSE
         WRITE(101,outfmt)  seq+i*2,outgen(1:x)
      ENDIF
      IF(gen.GT.genOUT-1)THEN
         IF (breed.eq.1) then;
            WRITE(122,outfmt2) seq+i*2-1,',',gen,sep,(seq-n_male*4)+(sort_sire1(i)*2-1),sep,&
                 (seq-n_male*4)+(n_male*2+(sort_dam1(i)*2-1)),sep,1,sep,'GT'
         ELSE
            WRITE(122,outfmt2) seq+i*2,',',gen,sep,(seq-n_male*4)+(sort_sire1(i)*2),sep,&
                 (seq-n_male*4)+(n_male*2+(sort_dam1(i)*2)),sep,1,sep,'GT'
         ENDIF
      ELSE
         IF (breed.eq.1) then;
            WRITE(122,outfmt3) i*2-1,',',gen,sep,'0',sep,'0',sep,1,sep,'GT'
         else
            WRITE(122,outfmt3) i*2,',',gen,sep,'0',sep,'0',sep,1,sep,'GT'
         endif
      ENDIF
   ENDIF
ENDDO

DO i=1, n_female
   ! Do animal by animal to aviod using huge RAM
   CALL decompress(f_genome(i,:,:,:),out_f_genome)
   !! HAPLOTYPES
   outgen(:)=5
   outgen2(:)=5
   x=0
   qtp=1
   IF (printgeno.eq.1.or.printgeno.eq.3) THEN
      DO j=1, len_genome
         IF(j .EQ. qtl_position(qtp))THEN
            x=x+1
            outgen(x)=out_f_genome(j,1)
            outgen2(x)=out_f_genome(j,2)
            IF (qtp .LT. HmaxQTL) qtp=qtp+1
         ELSE
            IF(frequencies(j,1) .LT. thresh .OR. frequencies(j,2) .LT. thresh) CYCLE
            x=x+1
            outgen(x)=out_f_genome(j,1)
            outgen2(x)=out_f_genome(j,2)
         ENDIF
      ENDDO
      outfmt='(I10,1X,                   I1)'
      WRITE(outfmt(10:20),'(I0)') x
      if (printgeno.eq.3) then
         IF (breed.eq.1) then;
            WRITE(1011,outfmt)  seq+(n_male+i)*2-1,outgen(1:x)
            WRITE(1011,outfmt)  seq+(n_male+i)*2-1,outgen2(1:x)
         ELSE
            WRITE(1011,outfmt)  seq+(n_male+i)*2,outgen(1:x)
            WRITE(1011,outfmt)  seq+(n_male+i)*2,outgen2(1:x)
         ENDIF
      else
         IF (breed.eq.1) then;
            WRITE(101,outfmt)  seq+(n_male+i)*2-1,outgen(1:x)
            WRITE(101,outfmt)  seq+(n_male+i)*2-1,outgen2(1:x)
         ELSE
            WRITE(101,outfmt)  seq+(n_male+i)*2,outgen(1:x)
            WRITE(101,outfmt)  seq+(n_male+i)*2,outgen2(1:x)
         ENDIF
         IF(gen.GT.genOUT-1)THEN
            IF (breed.eq.1) then;
               WRITE(122,outfmt2) seq+(n_male+i)*2-1,',',gen,sep,(seq-n_male*4)+(sort_sire(i)*2-1),sep,&
                    (seq-n_male*4)+(n_male*2+(sort_dam(i)*2-1)),sep,2,sep,'HT'
            ELSE
               WRITE(122,outfmt2) seq+(n_male+i)*2,',',gen,sep,(seq-n_male*4)+(sort_sire(i)*2),sep,&
                    (seq-n_male*4)+(n_male*2+(sort_dam(i)*2)),sep,2,sep,'HT'
            ENDIF
         ELSE
            IF (breed.eq.1) then;        
               WRITE(122,outfmt3) (n_male+i)*2-1,',',gen,sep,'0',sep,'0',sep,2,sep,'HT'
            ELSE
               WRITE(122,outfmt3) (n_male+i)*2,',',gen,sep,'0',sep,'0',sep,2,sep,'HT'
            ENDIF
         ENDIF
      endif
   endif
!! GENOTYPES
   IF (printgeno.ge.2)then
   outgen(:)=5
   x=0
   qtp=1
      DO j=1, len_genome
         IF(j .EQ. qtl_position(qtp))THEN
            x=x+1
            outgen(x)=out_f_genome(j,1)+out_f_genome(j,2)
            IF (qtp .LT. HmaxQTL) qtp=qtp+1
         ELSE
            IF(frequencies(j,1) .LT. thresh .OR. frequencies(j,2) .LT. thresh) CYCLE
            x=x+1
            outgen(x)=out_f_genome(j,1)+out_f_genome(j,2)
         ENDIF
      ENDDO
      outfmt='(I10,1X,                   I1)'
      WRITE(outfmt(10:20),'(I0)') x
      IF (breed.eq.1) then;
         WRITE(101,outfmt)  seq+(n_male+i)*2-1,outgen(1:x)
      ELSE
         WRITE(101,outfmt)  seq+(n_male+i)*2,outgen(1:x)
      ENDIF
      IF(gen.GT.genOUT-1)THEN
         IF (breed.eq.1) then;
            WRITE(122,outfmt2) seq+(n_male+i)*2-1,',',gen,sep,(seq-n_male*4)+(sort_sire(i)*2-1),sep,&
                 (seq-n_male*4)+(n_male*2+(sort_dam(i)*2-1)),sep,2,sep,'GT'
         ELSE
            WRITE(122,outfmt2) seq+(n_male+i)*2,',',gen,sep,(seq-n_male*4)+(sort_sire(i)*2),sep,&
                 (seq-n_male*4)+(n_male*2+(sort_dam(i)*2)),sep,2,sep,'GT'
         ENDIF
      ELSE
         IF (breed.eq.1) then;
            WRITE(122,outfmt3) (n_male+i)*2-1,',',gen,sep,'0',sep,'0',sep,2,sep,'GT'
         ELSE
            WRITE(122,outfmt3) (n_male+i)*2,',',gen,sep,'0',sep,'0',sep,2,sep,'GT'
         ENDIF
      ENDIF
   ENDIF
ENDDO
seq=seq+((n_male+n_female)*2)

DEALLOCATE(out_m_genome,out_f_genome,outgen)

END SUBROUTINE write_genotypeNOPED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_phenotypeNOPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_phenotypeNOPED
USE Prog_Param


!! Open and create header
DO j=1,n_male
   WRITE(102,'(I0,A1,I0,A1)',ADVANCE='No')j,'_',gen,sep
   DO i=1,n_trait-1
      IF (i.EQ.trait_n)THEN;
         WRITE(102,'(A,F0.6,A2)',ADVANCE='No')'*',sALL_phenotype(j,i),'*,'
      ELSE
         WRITE(102,'(F0.6,A1)',ADVANCE='No')sALL_phenotype(j,i),sep         
      ENDIF
   ENDDO
   IF(n_trait.EQ.trait_n)THEN
      WRITE(102,'(A,F0.6,A)')'*',sALL_phenotype(j,n_trait),'*'
   ELSE
      WRITE(102,'(F0.6)')sALL_phenotype(j,n_trait)
   ENDIF
ENDDO

DO j=1,n_female
   WRITE(102,'(I0,A1,I0,A1)',ADVANCE='No')n_male+j,'_',gen,sep
   DO i=1,n_trait-1
      IF (i.EQ.trait_n)THEN;
         WRITE(102,'(A,F0.6,A2)',ADVANCE='No')'*',dALL_phenotype(j,i),'*,'
      ELSE
         WRITE(102,'(F0.6,A1)',ADVANCE='No')dALL_phenotype(j,i),sep
      ENDIF
   ENDDO
   IF(n_trait.EQ.trait_n)THEN
      WRITE(102,'(A,F0.6,A)')'*',dALL_phenotype(j,n_trait),'*'
   ELSE
      WRITE(102,'(F0.6)')dALL_phenotype(j,n_trait)
   ENDIF
ENDDO

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_phenotypeNOPED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_tbvNOPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_tbvNOPED
USE Prog_Param


DO j=1,n_male
   WRITE(103,'(I0,A1,I0,A1)',ADVANCE='No')j,'_',gen,sep
   DO i=1,n_trait-1
      IF (i.EQ.trait_n)THEN;
         WRITE(103,'(A,F0.6,A2)',ADVANCE='No')'*',tsALL_phenotype(j,i),'*,'
      ELSE
         WRITE(103,'(F0.6,A1)',ADVANCE='No')tsALL_phenotype(j,i),sep
      ENDIF
   ENDDO
   IF(n_trait.EQ.trait_n)THEN
      WRITE(103,'(A,F0.6,A)')'*',tsALL_phenotype(j,n_trait),'*'
   ELSE
      WRITE(103,'(F0.6)')tsALL_phenotype(j,n_trait)
   ENDIF
ENDDO

DO j=1,n_female
   WRITE(103,'(I0,A1,I0,A1)',ADVANCE='No')n_male+j,'_',gen,sep
   DO i=1,n_trait-1
      IF (i.EQ.trait_n)THEN;
         WRITE(103,'(A,F0.6,A2)',ADVANCE='No')'*',tdALL_phenotype(j,i),'*,'
      ELSE
         WRITE(103,'(F0.6,A1)',ADVANCE='No')tdALL_phenotype(j,i),sep
      ENDIF
   ENDDO
   IF(n_trait.EQ.trait_n)THEN
      WRITE(103,'(A,F0.6,A)')'*',tdALL_phenotype(j,n_trait),'*'
   ELSE
      WRITE(103,'(F0.6)')tdALL_phenotype(j,n_trait)
   ENDIF
ENDDO
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_tbvNOPED

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_corrPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
SUBROUTINE write_corrNOPED
USE Prog_Param

DO i=1, n_trait
      WRITE(104,'(I0,A,I0,A,F0.2,A,F0.4,A,F0.4)') gen,sep,i,sep,h2(i),sep,corr(i),sep,corrB(i)
ENDDO
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END SUBROUTINE write_corrNOPED
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE gen_breed_base
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!SUBROUTINE gen_breed_base (n_animal, s_genome, d_genome, a_genome)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!USE Numeric_Kinds
!USE Subs_General
!USE RNG
!USE Prog_Param, ONLY: idum1, n_chromosome, n_linkage, n_locus, n_allele, &
!                      sort_sire, sort_dam, recombination_rate, mutation_rate

!IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!INTEGER(i4), INTENT( IN)                                                      :: n_animal
!INTEGER(i4), INTENT( IN), DIMENSION(n_animal,n_chromosome,n_linkage,n_allele) :: s_genome, d_genome
!INTEGER(i4), INTENT(OUT), DIMENSION(n_animal,n_chromosome,n_linkage,n_allele) :: a_genome

!INTEGER(i4) :: animal
!INTEGER(i4) :: chromosome, linkage, locus, s_allele, d_allele, bit_pos, mutation
!REAL   (dp) :: ran_num_1

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!a_genome=0;


!DO animal=1,n_animal


! DO chromosome=1,n_chromosome
!  ran_num_1=ran1(idum1)
!  IF(ran_num_1 < 0.5D0) THEN
!   s_allele=1                                                                          ! Choose sire's paternal chromosome
!  ELSE
!   s_allele=2                                                                          ! Choose sire's maternal chromosome
!  ENDIF
!  ran_num_1=ran1(idum1)
!  IF(ran_num_1 < 0.5D0) THEN
!   d_allele=1                                                                           ! Choose dam's paternal chromosome
!  ELSE
!   d_allele=2                                                                           ! Choose dam's maternal chromosome
!  ENDIF
!  DO linkage=1,n_linkage
!   DO locus=1,n_locus
!     bit_pos=locus-1
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!     ran_num_1=ran1(idum1)
!     IF(ran_num_1 < recombination_rate) THEN
!      IF (s_allele==1) s_allele=2                                                   ! Shift to the other sire's chromosome
!      IF (s_allele==2) s_allele=1                                                   ! Shift to the other sire's chromosome
!     ENDIF
!     ran_num_1=ran1(idum1)
!     IF(BTEST(s_genome(sort_sire(animal),chromosome,linkage,s_allele),bit_pos)) THEN                      ! If bit_pos = 1
!      IF(ran_num_1 > mutation_rate) THEN
!       a_genome(animal,chromosome,linkage,1)=IBSET(a_genome(animal,chromosome,linkage,1),bit_pos)       ! Keep it as it is
!      ELSE
!       a_genome(animal,chromosome,linkage,1)=IBCLR(a_genome(animal,chromosome,linkage,1),bit_pos); mutation=mutation+1
!      ENDIF
!     ELSE                                                                                                 ! If bit_pos = 0
!      IF(ran_num_1 > mutation_rate) THEN
!       a_genome(animal,chromosome,linkage,1)=IBCLR(a_genome(animal,chromosome,linkage,1),bit_pos)       ! Keep it as it is
!      ELSE
!       a_genome(animal,chromosome,linkage,1)=IBSET(a_genome(animal,chromosome,linkage,1),bit_pos); mutation=mutation+1
!      ENDIF
!     ENDIF
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!     ran_num_1=ran1(idum1)
!     IF(ran_num_1 < recombination_rate) THEN
!      IF (d_allele==1) d_allele=2                                                    ! Shift to the other dam's chromosome
!      IF (d_allele==2) d_allele=1                                                    ! Shift to the other dam's chromosome
!     ENDIF
!     ran_num_1=ran1(idum1)
!     IF(BTEST(d_genome(sort_dam(animal),chromosome,linkage,d_allele),bit_pos)) THEN                       ! If bit_pos = 1
!      IF(ran_num_1 > mutation_rate) THEN
!       a_genome(animal,chromosome,linkage,2)=IBSET(a_genome(animal,chromosome,linkage,2),bit_pos)       ! Keep it as it is
!      ELSE
!       a_genome(animal,chromosome,linkage,2)=IBCLR(a_genome(animal,chromosome,linkage,2),bit_pos); mutation=mutation+1
!      ENDIF
!     ELSE                                                                                               ! Set bit_pos to 0
!      IF(ran_num_1 > mutation_rate) THEN
!       a_genome(animal,chromosome,linkage,2)=IBCLR(a_genome(animal,chromosome,linkage,2),bit_pos)       ! Keep it as it is
!      ELSE
!       a_genome(animal,chromosome,linkage,2)=IBSET(a_genome(animal,chromosome,linkage,2),bit_pos); mutation=mutation+1
!      ENDIF
!     ENDIF
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!END SUBROUTINE gen_breed_base

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! subroutine sel_snp
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!SUBROUTINE sel_snp(snp_array)
!
! Go through all SNPs, choose NUMSNP(e.g. 54001) of them as SNPs on a SNP chip
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!USE Numeric_Kinds
!USE Subs_General
!USE RNG
!USE Prog_Param, ONLY: idum1, chromosome, len_genome, n_chromosome, linkage, n_linkage, locus, n_locus, numsnp
!
!IMPLICIT NONE
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!INTEGER(i4), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: snp_array
!INTEGER(i4) :: i=0, position=0
!REAL   (dp) :: ran_num_1
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!ALLOCATE (snp_array(numsnp,2))
!snp_array=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!DO chromosome=1,n_chromosome
!  DO linkage=1,n_linkage
!    DO locus=1,n_locus
!      position=position+1
!      ran_num_1=ran1(idum1)
!      IF (ran_num_1 < (REAL(numsnp)/REAL(len_genome))) THEN
!        i=i+1
!        IF (i > numsnp) EXIT
!        snp_array(i,1)=chromosome
!        snp_array(i,2)=position
!      ENDIF
!    ENDDO
!  ENDDO
!ENDDO
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!END SUBROUTINE sel_snp
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!
!

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE choose_a_fitness_trait
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!SUBROUTINE choose_a_fitness_trait (trait_n, group)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!USE Numeric_Kinds
!USE RNG
!USE Prog_Param, ONLY: idum1, h2

!IMPLICIT NONE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!CHARACTER*1, INTENT(IN )   :: group
!!INTEGER(i4), INTENT(OUT)   :: trait_n
!REAL   (dp)                :: ran_num_1

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!trait_n=0
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!ran_num_1=ran1(idum1)
!IF (group .EQ. 'A') trait_n =     INT (ran_num_1*10)+1 !range h2=1-11 (included)
!IF (group .EQ. 'B') trait_n= 30 + INT (ran_num_1*10)+1
!IF (group .EQ. 'C') trait_n= 60 + INT (ran_num_1*10)+1
!IF (group .EQ. 'D') trait_n= 90 + INT (ran_num_1*10)+1

!IF (trait_n .EQ. 11) trait_n=10                        !reset range h2=1-10
!IF (trait_n .EQ. 41) trait_n=40                        !reset range h2=1-10
!IF (trait_n .EQ. 71) trait_n=70                        !reset range h2=1-10
!IF (trait_n .EQ. 101) trait_n=100                      !reset range h2=1-10

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!END SUBROUTINE choose_a_fitness_trait

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE choose_a_medium_h2_trait
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!SUBROUTINE choose_a_medium_h2_trait (trait_n, group)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!USE Numeric_Kinds
!USE RNG
!USE Prog_Param, ONLY: idum1, h2

!IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!CHARACTER*1, INTENT(IN )   :: group
!INTEGER(i4), INTENT(OUT)   :: trait_n
!REAL   (dp)                :: ran_num_1
!trait_n=0

!ran_num_1=ran1(idum1)

!IF (group .EQ. 'A') trait_n =10  + INT (ran_num_1*10)+1 !range h2=12-31 (included)
!IF (group .EQ. 'B') trait_n= 40  + INT (ran_num_1*10)+1
!IF (group .EQ. 'C') trait_n= 70  + INT (ran_num_1*10)+1
!IF (group .EQ. 'D') trait_n= 100 + INT (ran_num_1*10)+1

!IF (trait_n .EQ. 21)  trait_n=20                        !reset range h2=12-30
!IF (trait_n .EQ. 51)  trait_n=50                        !reset range h2=12-30
!IF (trait_n .EQ. 81)  trait_n=80                        !reset range h2=12-30
!IF (trait_n .EQ. 111) trait_n=110                       !reset range h2=12-30

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!END SUBROUTINE choose_a_medium_h2_trait

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE choose_a_high_h2_trait
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!SUBROUTINE choose_a_high_h2_trait (trait_n, group)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!USE Numeric_Kinds
!USE RNG
!USE Prog_Param, ONLY: idum1, h2

!IMPLICIT NONE

!CHARACTER*1, INTENT(IN )   :: group
!INTEGER(i4), INTENT(OUT)   :: trait_n
!REAL   (dp)                :: ran_num_1
!trait_n=0

!ran_num_1=ran1(idum1)

!IF (group .EQ. 'A') trait_n =20  + INT (ran_num_1*10)+1 !range h2=35-81 (included)
!IF (group .EQ. 'B') trait_n= 50  + INT (ran_num_1*10)+1
!IF (group .EQ. 'C') trait_n= 80  + INT (ran_num_1*10)+1
!IF (group .EQ. 'D') trait_n= 110 + INT (ran_num_1*10)+1

!IF (trait_n .EQ. 31 ) trait_n=30                        !reset range h2=35-80
!IF (trait_n .EQ. 61 ) trait_n=60                        !reset range h2=35-80
!IF (trait_n .EQ. 91 ) trait_n=90                        !reset range h2=35-80
!IF (trait_n .EQ. 121) trait_n=120                       !reset range h2=35-80

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!END SUBROUTINE choose_a_high_h2_trait

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120








!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! SUBROUTINE write_genotype
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!SUBROUTINE write_genotype
!USE Prog_Param
!
!INTEGER, ALLOCATABLE, DIMENSION(:)::outgen
!
!ALLOCATE(out_m_genome(len_genome, n_allele))
!ALLOCATE(out_f_genome(len_genome, n_allele))
!ALLOCATE(outgen(len_genome))
!
!outfmt='(I0,A,I,A,A,           I1)'
!WRITE(outfmt(13:22),'(I0)')len_genome
!
!OPEN(101, FILE=outname)
!DO i=1, n_male
!   !Do animal by animal to avoid using huge RAM
!   CALL decompress(m_genome(i,:,:,:),out_m_genome)
!   IF (printgeno.eq.1) THEN
!      WRITE(101,outfmt) i,sep,1,sep,'ph',sep,out_m_genome(:,1)
!      WRITE(101,outfmt) i,sep,1,sep,'mh',sep,out_m_genome(:,2)
!   ELSEIF (printgeno.eq.2) THEN
!      outgen(:)=9
!      DO j=1,len_genome
!         outgen(j)=out_m_genome(j,1)+out_m_genome(j,2)
!      ENDDO
!      WRITE(101,outfmt) i,sep,1,sep,'gt',sep,outgen
!   ENDIF
!ENDDO

!DO i=1, n_female
!   CALL decompress(f_genome(i,:,:,:),out_f_genome)
!   IF (printgeno.eq.1) THEN
!      WRITE(101,outfmt) n_male+i,sep,2,sep,'ph',sep,out_f_genome(:,1)
!      WRITE(101,outfmt) n_male+i,sep,2,sep,'mh',sep,out_f_genome(:,2)
!   ELSEIF (printgeno.eq.2) THEN
!      outgen(:)=9
!      DO j=1,len_genome
!         outgen(j)=out_f_genome(j,1)+out_f_genome(j,2)
!      ENDDO
!      WRITE(101,outfmt) n_male+i,sep,2,sep,'gt',sep,outgen
!   ENDIF
!ENDDO
!CLOSE(101)
!
!DEALLOCATE(out_m_genome,out_f_genome,outgen)
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
!END SUBROUTINE write_genotype
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120






!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

END MODULE Subs_Special
