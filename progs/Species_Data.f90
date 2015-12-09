!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

PROGRAM species_data
! This program is supposed to simulate data for genomic studies.
! It should start from (a) the base population and go through (2) species data, and (c) breed data, and finally to the
! country data.
! Therefore the name should be changed to something more appropriate.
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! HISTORY                                                                                                            
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Original program: Hossein Jorjani, 2009-2011
! Updates:
! June 2010 - Ezequiel Nicolazzi: - Add new variables in param file (previously hard-coded)
!         July,  1st 2012 - ELN : - New varibables(breeds_f/m_genome) to allow more breeds to be simulated               
!         July,23-24 2012 - ELN : - Cumulative Poisson to determine recomb/mutation in simulation (suggested by LVarona)
!         Aug,    06 2012 - ELN : - Stable version of multiple breeds with/without REAL pedigree
!                                 - Some bugs fixed
!                                 - Frequencies consider fixed SNPs & QTLs
!         Sep,    17 2012 - ELN : - Put all writing in Subs_Special
!                                 - Solved major issue in breed step (which gave segmentation fault)
!                                 - Checked consistency across platforms
!                 27 2012 - ELN : - Solved a phenotype generation error.
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
USE Prog_Param
USE Subs_Special
USE Subs_General
USE RNG 

IMPLICIT NONE
integer::run,ii,tim(3),DIVERGEITER=100,wil
real::rnd
logical::start=.TRUE.,DIVERGEMORE=.TRUE.
INTEGER(i1),ALLOCATABLE,DIMENSION(:)::outgen
call init_random_seed
call GET_ENVIRONMENT_VARIABLE(name="WORKPATH",value=path) ! Useful if run in different folders (linux/unix/mac only)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
CALL read_user_options

!! Me calculation, following Goddard (2009) in Genetica.      
Me=int((2*(Species_male+Species_female)*n_chromosome)/log(real(4*(Species_male+Species_female)*n_chromosome)))
high_qtl=int(Hcoeff_Me*Me)
low_qtl =int(Lcoeff_Me*Me)

CALL print_user_options

!--------0--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! Few checks before starting
WRITE(*,*);WRITE(*,*) "CHECKING USER_OPTIONS AND PEDIGREES (IF ANY)"
CALL control_user_options 

ALLOCATE(n_gen(1000,3))         

DO run=1,realped      ! This checks pedigrees and prints some simple stats 
   WRITE(*,*);WRITE(*,*) "CHECHING PEDIGREE FILE: ",TRIM(pedname(run))
   CALL readcontrol_pedig(pedname(run))    !!! PEDIGREE CONTROL ONLY (DISCARD ALL THE INFORMATION SEIZED)
   WRITE(*,*)"DISTRIBUTION OF ANIMALS WITHIN GENERATIONS IN PEDIGREE";WRITE(*,*)
   DO ii=1,1000
      IF (n_gen(ii,1).GT.0)THEN
         WRITE(*,'(2X,A,I3,A,I6,A,I6,A,I6,A,I6)') "Gen:",ii-1," # animals:",n_gen(ii,1),&
                                               " UNK sires",n_gen(ii,2)," UNKN dams",n_gen(ii,3)
      ENDIF
   END DO
   DEALLOCATE(ped,cgens,gens)
ENDDO

DEALLOCATE(n_gen)
WRITE(*,*) "CHECKS COMPLETED SUCCESSFULLY";WRITE(*,*)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

!!! Define mutation and recombination rates/distributions
mutation_rate        = 2.5/(10.**(mutation))           
exp_mutation         = (n_linkage*n_locus)*mutation_rate               !! Exp(mutations)=Chrom_length*mutation_rat

!!! Call Cumulative Poisson distribution to get the thresholds for recombination and mutation
CALL Cpoisson(REAL(recombination),  recomb_distr,  len_recomb) ! C-Poisson distribution for recombination
CALL Cpoisson(       exp_mutation,mutation_distr,len_mutation) ! C-Poisson distribution for mutation

n_max=MAX(species_male,species_female)
n_male  =species_male
n_female=species_female
ALLOCATE(outgen      ( len_genome))
ALLOCATE(o_genotypes ( len_genome,3 ))
ALLOCATE(frequencies ( len_genome,2 ))
ALLOCATE(h2          (   n_trait    ))
ALLOCATE(g_variance  (   n_trait    ))
ALLOCATE(g_varALL    (   n_trait ,3 ))
ALLOCATE(fixed       (   2          ))
ALLOCATE(cfreq       (   5          ))

ALLOCATE(f_genome(n_female,n_chromosome,n_linkage,n_allele))
ALLOCATE(d_genome(n_female,n_chromosome,n_linkage,n_allele))
ALLOCATE(m_genome(  n_male,n_chromosome,n_linkage,n_allele))
ALLOCATE(s_genome(  n_male,n_chromosome,n_linkage,n_allele))

ALLOCATE(d_phenotype(n_female))
ALLOCATE(s_phenotype(n_male  ))
ALLOCATE(t_phenotype(n_max   )) ! Used later on

ALLOCATE(sort_sire(n_max     ))
ALLOCATE(sort_dam (n_max     ))

IF(distr.EQ.2) CALL initialize_gamma         ! Initialize gamma distribution (if desired)
CALL initialize_h2(h2)                       ! Initialize h2's
CALL gen_species_base (n_male,   m_genome)   ! Generate first set of males
CALL gen_species_base (n_female, f_genome)   ! Generate first set of females

!CALL calculate_frequencies(m_genome,f_genome,o_genotypes)

                                                                   OPEN(243, FILE=TRIM(path)//'RESULTS/runtime_stats.txt')

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!
! ****************************************************** SPECIES *********************************************************
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

DO gen=1,species_gen
   call init_random_seed
   t1 = SECNDS(0.0)
   ! Previous generations of males and females are now sires and dams.
   s_genome=m_genome; m_genome=0
   d_genome=f_genome; f_genome=0

   IF(MOD(gen,1000)==0) THEN
      CALL calculate_frequencies(s_genome,d_genome,o_genotypes)
      WRITE(243,*),"Gen:",gen,"Homoz:",homo,"OBShetero:",ho,"FreqDistr:",cfreq
      WRITE(*,*),"Gen:",gen,"Homoz:",homo,"OBShetero:",ho,"FreqDistr:",cfreq
   ENDIF

   ! This generates males for next Gen
   CALL randomize (n_male  ,n_male, sort_sire)
   CALL randomize (n_female,n_male, sort_dam)

   CALL gen_offspring (n_male, n_male, n_female , s_genome, d_genome, m_genome, .TRUE.)
   !!!
!   CALL decompress2(s_genome(1,:,:,:),outgen)
!   print*,gen,"  END:",outgen(len_genome-15:len_genome)
!   CALL decompress2(s_genome(2,:,:,:),outgen)
!   print*,gen,"  END:",outgen(len_genome-15:len_genome)
   
   ! This generates females for next Gen
   CALL randomize (n_male  ,n_female, sort_sire)
   CALL randomize (n_female,n_female, sort_dam)
   CALL gen_offspring (n_female, n_male, n_female , s_genome, d_genome, f_genome, .TRUE.)

                                                                                                        delta = SECNDS(t1)
                                                                                                 IF(MOD(gen,1000)==0) THEN
                                  WRITE (*,'(A,I0,A,F8.2,A)') ' SPECIES GENERATION ',gen,' FINISHED IN: ',delta,' SECONDS'
                                                                                                                     ENDIF
END DO

!!! Assign QTLs to non-homozygous SNPs after the species generation.
CALL calculate_frequencies(m_genome,f_genome,o_genotypes)
CALL gen_qtl
CALL get_Vg(g_variance);g_varALL(:,1)=g_variance

                                                                   OPEN(110, FILE=TRIM(path)//'RESULTS/gtl_positions.txt')
                                                                                 outfmt='(I8,1x,          (F0.4,1x),A,6A)'
                                                                                         WRITE(outfmt(8:15),'(I0)')n_trait
                                                                                                            DO i=1,HmaxQTL
                                             WRITE(110,outfmt) qtl_position(i), (qtl_size(i,j), j=1,n_trait),qtl_type(i,:)
                                                                                                                     ENDDO
                                                                                                                CLOSE(110)

                                                                                        outname="RESULTS/freq_basePOP.end"
                                                                                                     CALL write_frequency
!                                                                                                         IF(writeout)THEN
!                                                                                   WRITE(*,*),"WRITING OUT GENOTYPES... "
!                                                                                          outname="GTYPES/allspecies.end"
!                                                                                             call write_genotype(outname)
!                                                                                                                    ENDIF
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

s_genome=m_genome
d_genome=f_genome

DEALLOCATE(m_genome, f_genome)
DEALLOCATE(d_phenotype,s_phenotype,t_phenotype)


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! ******************************************************* DIVERGE BREEDS *************************************************
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

n_bre_male  =  breed_male   
n_bre_female=breed_female   

ALLOCATE(breeds_m_genome   (  Nbreeds,n_bre_male,n_chromosome,n_linkage,n_allele))
ALLOCATE(breeds_f_genome   (Nbreeds,n_bre_female,n_chromosome,n_linkage,n_allele))
!ALLOCATE(newsort_sire                                      (  Nbreeds,n_bre_male))
!ALLOCATE(newsort_dam                                       (Nbreeds,n_bre_female))

DO breed=1,Nbreeds
   call init_random_seed
   IF(start)THEN
      origbm=  breed_male
      origbf=breed_female
   ENDIF

   breed_male  =origbm   
   breed_female=origbf 
   g_varALL(:,2:3)=0.D0

   IF (breed .GT. realped)THEN
      WRITE(*,*);WRITE(*,*),"-------------------------------------------------------------------------"
      WRITE(*,'(1X,A22,I3,A)'),"** PROCESSING BREED N.",breed," - NO REAL PEDIGREE FILE FOR THIS BREED "
      WRITE(*,*),"------------------------------------------------------------------------";WRITE(*,*)
   ELSE
      WRITE(*,*);WRITE(*,*),"--------------------------------------------------------------------------"
      WRITE(*,'(1X,A22,I3,A23,1X,A)'),"** PROCESSING BREED N.",breed," - REAL PEDIGREE FILE -> ",TRIM(pedname(breed))
      WRITE(*,*),"--------------------------------------------------------------------------";WRITE(*,*)
   ENDIF
   

   ! Divide N "SPECIES" animals into smaller breeds, following sort_sire/dam
   ! Do this only one time
   ! Get animals divided equally into Nbreeds,and with new (reduced) sort variable.  
   IF(start)THEN
      ALLOCATE(d_phenotype(species_female))
      ALLOCATE(s_phenotype(species_male  ))
      ALLOCATE(t_phenotype(n_max         ))

      CALL ITIME(tim)
      DO i=1,tim(3)
         CALL random_number(rnd)
      ENDDO
      trait_n=1+int(17*rnd)              !! Imperfect random sampler of trait to diverge breeds
      CALL calculate_frequencies    (s_genome,d_genome,o_genotypes)
      write(243,*);write(243,*)," GENERATION BREEDS ";write(243,*);
      write(243,*),"Gen:",gen,"Homoz:",homo,"OBShetero:",ho,"FreqDistr:",cfreq
      write(243,*)
      CALL generate_phenotype (species_male  ,s_genome,s_phenotype,t_phenotype)
      CALL generate_phenotype (species_female,d_genome,d_phenotype,t_phenotype)
      ! Use the same "sort_sire" & "sort_dam" for both male and female offspring
      CALL sort_phenotypes (species_male  ,s_phenotype,sort_sire)
      CALL sort_phenotypes (species_female,d_phenotype,sort_dam)
      CALL diverge_breeds(s_genome,d_genome,breeds_m_genome,breeds_f_genome) 
      DEALLOCATE(d_phenotype,s_phenotype,t_phenotype)
   ENDIF
 
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120
! *****************************RE-BUILD VECTORS FOR 1 BREED TO FOLLOW NEXT STEPS   ***************************************
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120

   DEALLOCATE(s_genome,  d_genome)
   DEALLOCATE(sort_sire, sort_dam)

   ALLOCATE(s_genome   (  n_bre_male,n_chromosome,n_linkage,n_allele))
   ALLOCATE(d_genome   (n_bre_female,n_chromosome,n_linkage,n_allele))
!   ALLOCATE(sort_sire  (  n_bre_male                                ))
!   ALLOCATE(sort_dam   (n_bre_female                                ))
   
   ! Assign variables to a (reduced, without migration) subset of SPECIES DATA
   s_genome =breeds_m_genome(breed,:,:,:,:)
   d_genome =breeds_f_genome(breed,:,:,:,:)
!   sort_sire=   newsort_sire(breed,:      )
!   sort_dam =    newsort_dam(breed,:      )

   !! NOTE: This should probably be changed if multiple countries are considered! (ELN)
   gen=0
   IF (breed .LE. realped) THEN
      ALLOCATE(n_gen(1000,3))          ! Coming from Prog_Param
      CALL readcontrol_pedig(pedname(breed)) 
      n_max=MAX(needSIRE,needDAM)
      IF (n_bre_male.LT.n_max.OR.n_bre_female.LT.n_max)THEN
         WRITE(*,*)   
         WRITE(*,*) 'WARNING:   Number of BASE (fe)males was ',  n_bre_male,'. Set to',n_max
         WRITE(*,*) 'WARNING:    (this to fit the need of phantom animals in the REAL pedigree).'
         WRITE(*,*)      
         n_male  =n_max 
         n_female=n_max 
      ELSE
         n_max=MAX(n_bre_male,n_bre_female)
         IF(n_bre_male.NE.n_bre_female)THEN
            WRITE(*,*),'WARNING: breed male and female MUST be = (for now). Normalizing it to:',n_max
            n_male=n_max
            n_female=n_max
         ELSE
            n_male  =n_bre_male
            n_female=n_bre_female
         ENDIF
      ENDIF
   ELSE
      ! n_bre_male cannot be != n_bre_female because they come from the same animal#(species) / n_breeds
      n_male  =n_bre_male
      n_female=n_bre_female
   ENDIF

   ALLOCATE(s_phenotype(  n_male                                ))
   ALLOCATE(d_phenotype(n_female                                ))
   ALLOCATE(t_phenotype(   n_max                                ))
   ALLOCATE(m_genome   (  n_male,n_chromosome,n_linkage,n_allele))
   ALLOCATE(f_genome   (n_female,n_chromosome,n_linkage,n_allele))
   
   ! Expand the size of m_genome and f_genome in the species populations to
   !        the size of s_genome anf d_genome in the breed   populations
   CALL expand_genome_breed(s_genome,d_genome,m_genome,f_genome) 

   ! Calculate initial frequencies for each breed (name it as number of breed
   CALL calculate_frequencies(m_genome,f_genome,o_genotypes)
   write(*,*),"--->PREDIVERGENCE Gen:",gen,"Homoz:",homo,"OBShetero:",ho,"FreqDistr:",cfreq
   CALL get_Vg(g_variance);g_varALL(:,2)=g_variance

   DEALLOCATE(s_genome, d_genome)
!   DEALLOCATE(sort_sire, sort_dam)

   ALLOCATE(s_genome   (  n_male,n_chromosome,n_linkage,n_allele))
   ALLOCATE(d_genome   (n_female,n_chromosome,n_linkage,n_allele))
   ALLOCATE(sort_sire  (n_max                                   ))
   ALLOCATE(sort_dam   (n_max                                   ))


   ! Diverve breeds more if wanted
   IF (DIVERGEMORE) THEN
      DO wil=1,DIVERGEITER
         call init_random_seed
         print*,'EXTRA DIVERGENCE ITER:',wil
         s_genome=m_genome; m_genome=0;sort_sire=0
         d_genome=f_genome; f_genome=0;sort_dam=0
         
         CALL randomize (n_male, n_male, sort_sire)
         CALL randomize (n_female, n_male, sort_dam)
         CALL gen_offspring (n_male, n_male, n_female, s_genome, d_genome, m_genome,.FALSE.)
         
         CALL randomize (n_male, n_female, sort_sire)
         CALL randomize (n_female, n_female, sort_dam)
         CALL gen_offspring (n_female, n_male, n_female, s_genome, d_genome, f_genome,.FALSE.)
         
         CALL calculate_frequencies(m_genome,f_genome,o_genotypes)
         write(243,*),"Gen:",gen,"Homoz:",homo,"OBShetero:",ho,"FreqDistr:",cfreq
      ENDDO
   ENDIF
!   s_genome=m_genome
!   d_genome=f_genome

   DEALLOCATE(s_genome,d_genome,sort_sire,sort_dam)
   DEALLOCATE(d_phenotype,s_phenotype,t_phenotype)

   gen=0  
   n_max=MAX(breed_male,breed_female)
   n_male  =breed_male
   n_female=breed_female
   ALLOCATE(s_phenotype(  n_male                                ))
   ALLOCATE(d_phenotype(n_female                                ))
   ALLOCATE(ts_phenotype(  n_male                               ))
   ALLOCATE(td_phenotype(n_female                               ))
!     ALLOCATE(m_genome   (  n_male,n_chromosome,n_linkage,n_allele))
!     ALLOCATE(f_genome   (n_female,n_chromosome,n_linkage,n_allele))
     
     ! Expand the size of m_genome and f_genome in the species populations to
     !        the size of s_genome anf d_genome in the breed   populations
     
!     CALL expand_genome_country(s_genome,d_genome,m_genome,f_genome)
     
!     DEALLOCATE(s_genome, d_genome)
     ALLOCATE(s_genome   (  n_male,n_chromosome,n_linkage,n_allele))
     ALLOCATE(d_genome   (n_female,n_chromosome,n_linkage,n_allele))
     ALLOCATE(sort_sire  (n_max))
     ALLOCATE(sort_dam   (n_max))
     ALLOCATE(sort_sire1  (n_max))
     ALLOCATE(sort_dam1   (n_max))
     
!     s_genome=m_genome; m_genome=0
!     d_genome=f_genome; f_genome=0
     
     t1 = SECNDS(0.0)

     delta = SECNDS(t1)
     WRITE (*,'(A,I0,A,F8.2,A)') ' BREED GENERATION ',gen,' FINISHED IN: ',delta,' SECONDS'

                                                                                                          IF(writeout)THEN
                                                                                                              !! GENOTYPES
                                                                                                                     seq=0
                                                            WRITE (outname, '("GTYPES/NOPED_final_", I0, ".end")' )  breed
                                                         WRITE (outname2,'("GTYPES/ID_GenealNOPED_", I0, ".end")' )  breed
                                                                                                    if(printgeno.eq.3)then
                                                                      OPEN(1011, FILE=TRIM(path)//TRIM(outname)//"_HAPLO")
                                                                        OPEN(101, FILE=TRIM(path)//TRIM(outname)//"_GENO")
                                                                                                                      else
                                                                                 OPEN(101, FILE=TRIM(path)//TRIM(outname))
                                                                                                                     endif
                                                                                OPEN(122, FILE=TRIM(path)//TRIM(outname2))
                                                                                                  CALL write_genotypeNOPED
                                                                                                                     ENDIF

                                                                                                            !! PHENOTYPES
                                                            WRITE (outname, '("RESULTS/phen_finalNOPED_",I0,".txt")')breed
                                                                                 OPEN(102, FILE=TRIM(path)//TRIM(outname))
                                                                            WRITE(102,'(A12)',ADVANCE='No') 'INDIVIDUAL,T'
                                                                                                          DO i=1,n_trait-1
                                                        WRITE(102,'(I0,A4,F0.2,A3)',ADVANCE='No') i,'(h2=' , h2(i) , '),T'
                                                                                                                     ENDDO
                                                           WRITE(102,'(I0,A4,F0.2,A1)') n_trait,'(h2=' , h2(n_trait) , ')'


                                                                                                                   !! TBVs
                                                             WRITE (outname, '("RESULTS/tbv_finalNOPED_",I0,".txt")')breed
                                                                                 OPEN(103, FILE=TRIM(path)//TRIM(outname))
                                                                            WRITE(103,'(A12)',ADVANCE='No') 'INDIVIDUAL,T'
                                                                                                          DO i=1,n_trait-1
                                                         WRITE(103,'(I0,A4,F0.2,A3)',ADVANCE='No') i,'(h2=' , h2(i), '),T'
                                                                                                                     ENDDO
                                                           WRITE(103,'(I0,A4,F0.2,A1)') n_trait,'(h2=' , h2(n_trait) , ')'

                                                                                                           !!CORRELATIONS
                                                   WRITE (outname, '("RESULTS/corr_phen-tbv_finalNOPED_",I0,".txt")')breed
                                                                                 OPEN(104, FILE=TRIM(path)//TRIM(outname))
                                              WRITE(104,'(A)') 'Generation,Trait_#,h2,MALE_r(tbv,phen),FEMALE_r(tbv,phen)'

                                                                                 
     ALLOCATE(corr (n_trait),  corrB(n_trait))
     ALLOCATE(sALL_phenotype (  n_male,n_trait))
     ALLOCATE(tsALL_phenotype(  n_male,n_trait))
     ALLOCATE(dALL_phenotype (n_female,n_trait))
     ALLOCATE(tdALL_phenotype(n_female,n_trait))

     CALL call_country_trait(breed-realped)

     DO gen=1,breed_gen
        call init_random_seed
        t1 = SECNDS(0.0)
        ! Previous generations of males and females are now sires and dams.
        s_genome=m_genome; m_genome=0
        d_genome=f_genome; f_genome=0
        
        ! Selection for a randomly chosen medium h2 trait every generation
        CALL calculate_frequencies(s_genome,d_genome,o_genotypes)
        write(243,*);write(243,*)," GENERATION NO_PEDIGREE ";write(243,*);
        write(243,*),"Gen:",gen,"Homoz:",homo,"OBShetero:",ho,"FreqDistr:",cfreq
        write(243,*)
        
        IF(writeout .AND. gen.GE.genOUT)THEN
           CALL generate_phenotypeALL(n_male  ,s_genome,sALL_phenotype,tsALL_phenotype,corr)
           CALL generate_phenotypeALL(n_female,d_genome,dALL_phenotype,tdALL_phenotype,corrB)
           s_phenotype = sALL_phenotype(:,trait_n) 
           ts_phenotype=tsALL_phenotype(:,trait_n)
           d_phenotype = dALL_phenotype(:,trait_n) 
           td_phenotype=tdALL_phenotype(:,trait_n)
        ELSE
           CALL generate_phenotype  (n_male  ,s_genome,s_phenotype,ts_phenotype)
           CALL generate_phenotype  (n_female,d_genome,d_phenotype,td_phenotype)
        ENDIF
 
        selratioM=selratio_male(breed-realped)
        selratioF=selratio_fem(breed-realped)

        call select_male(n_male, n_male, s_phenotype,sort_sire1)
        call select_female(n_female, n_male, d_phenotype, sort_dam1)
        sort_sire=sort_sire1
        sort_dam=sort_dam1
        CALL gen_offspring (n_male, n_male, n_female, s_genome, d_genome, m_genome,.FALSE.)
        
        call select_male(n_male, n_female, s_phenotype,sort_sire)
        call select_female(n_female, n_female, d_phenotype, sort_dam)
        CALL gen_offspring (n_female, n_male, n_female, s_genome, d_genome, f_genome,.FALSE.)
        delta = SECNDS(t1)

                                                                                                    IF( gen.GE.genOUT)THEN
                                                                                                 CALL write_phenotypeNOPED
                                                                                                 CALL       write_tbvNOPED
                                                                                                 CALL      write_corrNOPED
                                                                                     IF(writeout)CALL  write_genotypeNOPED
                                                                                                                     ENDIF

                                    WRITE (*,'(A,I0,A,F8.2,A)') ' BREED GENERATION ',gen,' FINISHED IN: ',delta,' SECONDS'
     END DO
                                                                                                    IF(writeout)CLOSE(101)
                                                                                                    IF(writeout)CLOSE(122)
                                                                                                                CLOSE(102)
                                                                                                                CLOSE(103)
                                                                                                                CLOSE(104)
     CALL calculate_frequencies(m_genome,f_genome,o_genotypes)
     CALL get_Vg(g_variance);g_varALL(:,3)=g_variance


                                                                   WRITE (outname, '("RESULTS/g_var_pop",I0,".txt")')breed
                                                                                                            CALL write_var
                                                        WRITE (outname, '("RESULTS/freq_finalNOPED_",I0,"_end.txt")')breed
                                                                                                      CALL write_frequency

      DEALLOCATE(corr,corrB)
      DEALLOCATE(sALL_phenotype,dALL_phenotype,tsALL_phenotype,tdALL_phenotype)
      DEALLOCATE(s_phenotype,d_phenotype,m_genome,f_genome,ts_phenotype,td_phenotype)
      DEALLOCATE(sort_sire1,sort_dam1)
   start=.FALSE.
                                                                                       END DO !! BREED-COUNTRY COMBINATION

WRITE(*,*);WRITE(*,*),"--------------------------------------------------------------------------"
WRITE(*,*),"                   PROGRAM FINISHED SUCCESSFULLY!!!!                      "
WRITE(*,*),"--------------------------------------------------------------------------";WRITE(*,*)

END PROGRAM species_data
