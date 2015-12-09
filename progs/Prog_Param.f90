MODULE prog_param
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! REVISIONS
! 2012-08-30: Many things had been hard coded in the programs. These shouls be moved out to here and read from the file
!             "user_options.txt". The first one idetified by Ezequiel is the number of SNPs on the chip.
!             The variable "numsnp" is for the number of SNPs.


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

USE numeric_kinds

IMPLICIT NONE

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
                                   !      VALID EXAMPLES
INTEGER(i4)  :: species_male  ,  & !  =   50
                species_female,  & !  =   50
                breed_male    ,  & !  =   50
                breed_female  ,  & !  =   50
                country_male  ,  & !  =   50
                country_female,  & !  =   50
                n_chromosome  ,  & !  =   30
                n_linkage     ,  & !  = 1000
                species_gen   ,  & !  =   10
                breed_gen     ,  & !  =    5
                country_gen   ,  & !  =    5
                recombination ,  & !  =    1
                mutation      ,  & !  =    4
                numsnp        ,  & !  =54001
                Nbreeds       ,  & !  =    2
                printgeno     ,  & !  =    1
                genOUT        ,  & !  =    1
                pedOUTx       ,  & !  =  100
                realped       ,  & !  =    0
                distr         ,  &
                admixprop


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! THERE   IS   NO   NEED   TO   CHANGE   ANYTHING   BEYOND   THIS   POINT
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

INTEGER(i4)                                  :: i          =      0, & ! Counters
                                                j          =      0, &
                                                k          =      0, &
                                                l          =      0, &
                                                m          =      0, &
                                                bit_pos    =      0, &
                                                snp        =      0, &
                                                n_male     =      0, &
                                                n_female   =      0, &
                                                n_max      =      0, &
                                                chromosome =      0, &
                                                linkage    =      0, &
                                                locus      =      0, &
                                                allele     =      0, &
                                                s_allele   =      0, &
                                                d_allele   =      0, &
                                                gen        =      0, &
                                                trait_n    =      0, &
                                                animal     =      0, &
                                                origbm     =      0, &
                                                origbf     =      0, &
                                                qtp        =      0, &
                                                pos1       =      1

INTEGER(i4)                                  :: n_locus    =     32, &
                                                n_allele   =      2, &
                                                n_trait    =     18, &
                                                idum1      = -99991, &
                                                idum2      = -99989, &
                                                idum3      = -99961

INTEGER(i4)                                  :: n_bre_male         , &
                                                n_bre_female       , &
                                                breed              , &
                                                seq                , &
                                                pos2


INTEGER(i4)                                  :: pedig_n            , &
                                                baseSIRE           , &
                                                baseDAM            , &
                                                needSIRE           , &
                                                needDAM            , &
                                                pedOUT

INTEGER(i4)                                  :: len_recomb         , &
                                                len_mutation       , &
                                                len_genome         , &
                                                Me                 , &
                                                HmaxQTL            , &
                                                high_qtl           , &
                                                low_qtl            , &
                                                homo               , &
                                                selratioM          , &
                                                selratioF

INTEGER(i4), ALLOCATABLE, DIMENSION(:)       :: qtl_position       , &
                                                fixed              


INTEGER(i4), ALLOCATABLE, DIMENSION(:)       :: anim               , &
                                                sire               , &
                                                dam                , &
                                                sort_sire          , &
                                                sort_dam           , &
                                                sort_sire1         , &
                                                sort_dam1          , &
                                                selratio_male      , &
                                                selratio_fem

INTEGER(i1), ALLOCATABLE, DIMENSION(:,:)     :: out_m_genome       , &
                                                out_f_genome

INTEGER(i4), ALLOCATABLE, DIMENSION(:,:)     :: o_genotypes        , &
                                                cgens              , &
                                                n_gen              , &
                                                newsort_sire       , &
                                                newsort_dam       

INTEGER(i4), ALLOCATABLE, DIMENSION(:,:,:)   :: gens

INTEGER(i4), ALLOCATABLE, DIMENSION(:,:,:,:) :: f_genome           , &
                                                d_genome           , &
                                                m_genome           , &
                                                s_genome           , &
                                                rp_i_genome

INTEGER(i4), ALLOCATABLE, DIMENSION(:,:,:,:,:)::breeds_m_genome    , &
                                                breeds_f_genome

REAL   (sp)                                  :: t1                 , &
                                                delta              , &
                                                exp_mutation       , &
                                                Hcoeff_Me          , &
                                                Lcoeff_Me          , &
                                                h2_1               , &
                                                h2_2               , &
                                                h2_3               , &
                                                pcorr              , &
                                                mcorr              , &
                                                ho                 , &
                                                thresh

REAL   (sp)                                  :: mutation_rate      , &
                                                recombination_rate , &  !!! TO ELIMINATE!!! 
                                                allele_effect

REAL   (dp)                                  :: ran_num_1          , &
                                                ran_num_2

REAL   (sp), EXTERNAL                        :: ran12              , &
                                                ran2               , &
                                                gasdev

REAL   (dp), ALLOCATABLE, DIMENSION(:)       :: s_phenotype        , &
                                                d_phenotype        , &
                                                ts_phenotype       , &
                                                td_phenotype       , &
                                                t_phenotype        , &
                                                corr               , &
                                                corrB              , &
                                                cfreq

REAL   (dp), ALLOCATABLE, DIMENSION(:)       :: recomb_distr       , &
                                                mutation_distr     , &
                                                g_variance     

REAL   (dp), ALLOCATABLE, DIMENSION(:,:)     :: frequencies        , &
                                                g_varALL           , &
                                                i_phenotype        , &
                                                it_phenotype       , &
                                                sALL_phenotype     , &
                                                dALL_phenotype     , &
                                                tsALL_phenotype    , &
                                                tdALL_phenotype   
     
REAL   (sp), ALLOCATABLE, DIMENSION(:,:)     :: qtl_size          

REAL   (sp), ALLOCATABLE, DIMENSION(:)       :: param              , &
                                                h2                 

INTEGER(i4), ALLOCATABLE, DIMENSION(:,:)     :: snp_array     

CHARACTER*1                                  :: sep=','            ! Separator for output files
CHARACTER*3                                  :: buffer             

CHARACTER*1, ALLOCATABLE, DIMENSION(:,:)     :: qtl_type          

CHARACTER*50                                 :: outname            , &
                                                outname2           , &
                                                outfmt            

CHARACTER*200                                :: path               , &
                                                selratio

CHARACTER*4, ALLOCATABLE, DIMENSION(:)       :: selez

CHARACTER*50, ALLOCATABLE, DIMENSION(:)      :: pedname

TYPE :: ped_type                                     ! pedigree structure
     INTEGER(i4)                             :: an   ! animal sequetial id number
     INTEGER(i4)                             :: si   ! sire sequential id number
     INTEGER(i4)                             :: da   ! dam sequential id number 
     INTEGER(i4)                             :: sex  ! sex number (1male 2female)
     INTEGER(i4)                             :: gen  ! generation number
END TYPE ped_type

TYPE (ped_type),ALLOCATABLE,DIMENSION (:)    :: Ped

LOGICAL                                      :: writeout






!TYPE genetic_effect
! INTEGER(i4), DIMENSION((n_chromosome*n_linkage*n_locus)/20)         :: qtl_position
! REAL   (dp), DIMENSION((n_chromosome*n_linkage*n_locus)/20,n_trait) :: qtl_size
!END TYPE genetic_effect
!TYPE (genetic_effect)                                                :: qtl_effect


!TYPE fulle_record
! INTEGER                                                             :: aid
! INTEGER                                                             :: sid
! INTEGER                                                             :: did
! INTEGER                                                             :: order
! INTEGER, DIMENSION(n_chromosome,n_linkage,n_allele)                 :: genome
! INTEGER, DIMENSION(n_trait)                                         :: geneotypic_value
! REAL(dp), DIMENSION(n_trait)                                        :: pheneotypic_value
!END TYPE full_record
!TYPE(full_record), DIMENSION((country_male+country_female)*country_gen) :: animal_record


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! C O N S T A N T S

INTEGER(i4), PARAMETER               :: int_0         =   0
INTEGER(i4), PARAMETER               :: int_1         =   1
INTEGER(i4), PARAMETER               :: int_2         =   2

REAL(dp)   , PARAMETER               :: ddeci         =   0.1D0
REAL(dp)   , PARAMETER               :: dcenti        =   0.01D0
REAL(dp)   , PARAMETER               :: dmilli        =   0.001D0
REAL(dp)   , PARAMETER               :: dmicro        =   0.000001D0
REAL(dp)   , PARAMETER               :: dnano         =   0.000000001D0
REAL(dp)   , PARAMETER               :: dpico         =   0.000000000001D0

REAL(dp)   , PARAMETER               :: dp_0          =   0.D0
REAL(dp)   , PARAMETER               :: dp_1          =   1.D0
REAL(dp)   , PARAMETER               :: dp_2          =   2.D0
REAL(dp)   , PARAMETER               :: dp_4          =   4.D0
REAL(dp)   , PARAMETER               :: dp_10         =  10.D0
REAL(dp)   , PARAMETER               :: dp_100        = 100.D0
REAL(dp)   , PARAMETER               :: POIS(6) =(/.3579,.7358,.9197,.9810,.9963,.9999/) ! To be updated with LV's funcion 

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120




!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
END MODULE prog_param
