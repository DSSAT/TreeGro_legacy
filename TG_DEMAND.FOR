!=======================================================================
!  TG_DEMAND, Subroutine, J.W. Jones and G. Hoogenboom.
!-----------------------------------------------------------------------
!  Calculates potential demand for C and N based upon new growth and
!  existing N deficiency in old tissue.
!-----------------------------------------------------------------------
!  REVISION       HISTORY
!  01/01/1990 JWJ Written.
!  02/01/1993 GH  Revised.
!  04/24/1994 NBP Changed TAIRHR to TGRO.
!  08/22/1995 GH  Added seed composition routine from KJB & ELPiper
!  04/02/1996 JWJ Modified partitioning during early growth
!  01/10/1997 GH  Added TURFAC effect on seed growth and pod addition
!  09/15/1998 CHP Modified for modular format
!  05/10/1999 GH  Incorporated in CROPGRO
!-----------------------------------------------------------------------
!  Called by:  PLANT
!  Calls:      SDCOMP, TREE_IPDMND
!=======================================================================

      SUBROUTINE TG_DEMAND(DYNAMIC, CONTROL,
     &  AGRLF, AGRRT, AGRSH2, AGRSTM, CROP, DRPP, DXR57,  !Input
     &  FILECC, FILEGC, FILEIO, FNINSH, FRACDN, LAGSD,    !Input
     &  LNGPEG, NDLEAF, NSTRES, PAR, PCNL, PCNRT, PCNST,  !Input
     &  PGAVL, PUNCSD, PUNCTR, PLTPOP, RPROAV, RTWT,      !Input
     &  SDDES, SDNO, SDVAR, SHELN, SHVAR, STMWT, SWFAC,   !Input
     &  TAVG, TDUMX, TDUMX2, TGRO, TURFAC, VSTAGE, WCRLF, !Input
     &  WCRRT, WCRST, WNRLF, WNRRT, WNRSH, WNRST, WTLF,   !Input
     &  WTSD, WTSHE, XPOD, NVEG0, NR1, NR2, NR5, NR7,     !Input

     &  AGRSD1, AGRSD2, AGRVG, AGRVG2, CDMREP, F, FNINL,  !Output
     &  FNINR, FNINS, FNINSD, FRLF, FRRT, FRSTM, GDMSD,   !Output
     &  GRRAT1, NDMNEW,  NDMOLD, NDMREP, NDMSDR, NDMTOT,  !Output
     &  NDMVEG, NMINEP, NMOBR, PHTIM, PNTIM, POTCAR,      !Output
     &  POTLIP, SDGR, TURADD, XFRT, YREND)                !Output

!-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData
      IMPLICIT NONE
      SAVE

      CHARACTER*2 CROP
      CHARACTER*3 TYPSDT
      CHARACTER*6   ERRKEY
      PARAMETER (ERRKEY = 'TG_DEMAND')
      CHARACTER*30 FILEIO
      CHARACTER*78 MSG(2)
      CHARACTER*92 FILECC, FILEGC

      INTEGER DYNAMIC   !, TIMDIF
      INTEGER NPP, I, NAGE, DAS
      INTEGER NDLEAF, NR1, NR2, NR5, NR7, NVEG0, YREND

      REAL FRLFM, FRSTMM, YY, XX, TMPFAC
      REAL REDPUN,TMPFCS,PAGE,REDSHL,SDMAX,CDMSH,GDMSH,ADDSHL
      REAL TEMXFR,CAVTOT,GDMSDO,CNOLD
      REAL NVSTL,NVSTS,NVSTR,FRNLFT
      REAL TABEX,CURV
      REAL POTLIP, POTCAR

      REAL TPHFAC,PARSLA,FFVEG
      REAL GROYES,GAINNW,GAINWT
      REAL SLAVAR, SLAREF, FINREF
      REAL SLAMAX, SLAMIN, THRESH
      REAL AGRSH2, AGRRT, AGRSTM, FRLFF, FRSTMF
      REAL CARMIN, LIPOPT, LIPTB, SLOSUM

      REAL AGRLF, AGRSD1, AGRSD2, AGRVG, AGRVG2,
     &  CDMREP, CDMSD, CDMSDR, CDMTOT,
     &  CDMVEG, DRPP, DUMFAC, DXR57, F,
     &  FNINL, FNINR, FNINS, FNINSD, FNINSH,
     &  FRACDN, FRLF, FRLFMX,
     &  FRRT, FRSTM, FVEG,
     &  GDMSD, GDMSDR,
     &  GROMAX, GRRAT1, LAGSD, LNGPEG, LNGSH,
     &  NDMNEW, NDMOLD, NDMREP,
     &  NDMSD, NDMSDR, NDMSH, NDMTOT, NDMVEG,
     &  NMINEP, NMOBMX, NMOBR, NRCVR, NSTRES,
     &  NVSMOB,
     &  PAR, PCNL, PCNRT, PCNST,
     &  PGAVL, PLIGSD, PLTPOP, PMINSD, POASD,
     &  PROLFF, PROLFI,
     &  PRORTF, PRORTI, PROSTF, PROSTI, RCH2O,
     &  RLIG, RLIP, RMIN, RNO3C,
     &  ROA, RPRO, RPROAV, RTWT, SDGR,
     &  SDLIP, SDPRO, SDVAR, SHLAG, SHVAR,
     &  SIZELF, SIZREF, SLAMN, SLAMX, SLAPAR,
     &  SRMAX, STMWT, SWFAC, TAVG, TDUMX,
     &  SIZRAT, TDUMX2,
     &  TURADD, TURFAC, TURSLA, TURXFR,
     &  VSSINK, VSTAGE, WCRLF, WCRRT, WCRST, WNRLF,
     &  WNRRT, WNRSH, WNRST, WTLF, XFRMAX,
     &  XFRT, XFRUIT, XPOD

      REAL FNSDT(4)
      REAL XVGROW(6), YVGROW(6), YVREF(6)
      REAL XSLATM(10), YSLATM(10), XTRFAC(10), YTRFAC(10),
     &         XXFTEM(10), YXFTEM(10)
      REAL XLEAF(25), YLEAF(25), YSTEM(25)
      REAL TGRO(TS)
      REAL SDDES(NCOHORTS), SDNO(NCOHORTS), SHELN(NCOHORTS) 
      REAL WTSD(NCOHORTS), WTSHE(NCOHORTS)
      REAL PHTIM(NCOHORTS), PNTIM(NCOHORTS)

      REAL TURFSL

!CHP - puncture variables, not functional
      REAL PUNCSD, PUNCTR, RPRPUN

      TYPE (ControlType) CONTROL

!***********************************************************************
!***********************************************************************
!     Run Initialization - Called once per simulation
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------

! Add the following if statement by JZW
      ! IF (CONTROL % RUN == 1 .OR. INDEX('FQ',CONTROL%RNMODE) < 0) THEN
      CALL IPDMND(
     &  FILECC, FILEGC, FILEIO,                           !Input
     &  CARMIN, FINREF, FNSDT, FRLFF, FRLFMX,             !Output
     &  FRSTMF, LIPOPT, LIPTB, LNGSH, NMOBMX,             !Output
     &  NRCVR, NVSMOB, PLIGSD, PMINSD, POASD,             !Output
     &  PROLFF, PROLFI, PRORTF, PRORTI, PROSTF, PROSTI,   !Output
     &  RCH2O, RLIG, RLIP, RMIN, RNO3C, ROA,              !Output
     &  RPRO, SDLIP, SDPRO, SHLAG, SLAMAX, SLAMIN,        !Output
     &  SLAPAR, SLAREF, SLAVAR, SLOSUM, SIZELF, SIZREF,   !Output
     &  SRMAX, THRESH, TURSLA, TYPSDT, VSSINK, XFRMAX,    !Output
     &  XFRUIT, XLEAF, XSLATM, XTRFAC, XVGROW, XXFTEM,    !Output
     &  YLEAF, YSLATM, YSTEM, YTRFAC, YVREF, YXFTEM)      !Output
      ! EndIf

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      !If (CONTROL % MULTI <2) then
      IF (CONTROL%RUN .EQ. 1 .OR. INDEX('QPF',CONTROL%RNMODE).LE.0) THEN
        CDMSDR = 0.0  !CH2O required for seed growth from mobilized N
        GDMSDR = 0.0  !Potential seed growth from NDMSDR 
        FNINSD = 0.0  !Maximum fraction of N for growing seed tissue based on temperature
        NDMNEW = 0.0  !Total N demand for new growth
        NDMSD  = 0.0  !Total N demand to grow seed demand (GDMSD)
        NDMSDR = 0.0  !Amount of Mobilized N which can be used for seed growth
        NDMVEG = 0.0  !N required for vegetative growth if all PGAVL is used as computed
        NMOBR  = 0.0  !Stage-dependent potential N mining rate expressed as a fraction of the maximum rate 
        SDGR   = 0.0  !Potential growth rate per seed
        FNINL  = 0.0  !Maximum fraction of N for growing leaf tissue
        FNINR  = 0.0  !Maximum fraction of N for growing root tissue
        NMINEP = 0.0  !Potential N mobilization from storage
        RPRPUN = 1.0  !Puncture damage reduction variable
        TMPFAC = 1.0  !Modifies maximum growth rate for seed and shells depending on temperature       
      endif
        NDMREP = 0.0  !Total N needed for potential reproductive growth
        NDMSH  = 0.0  !Total N demand to grow shell demand (GDMSH)
        FNINS  = 0.0  !Maximum fraction of N for growing stem tissue
!-----------------------------------------------------------------------
!     SET VARIETY SPECIFIC LEAF PARAMETERS
!-----------------------------------------------------------------------
      IF (CROP .NE. 'FA') THEN
        DUMFAC = SLAVAR / SLAREF
        F      = DUMFAC * FINREF
        FVEG   = DUMFAC * SLAMAX
        SLAMN  = DUMFAC * SLAMIN
        SLAMX  = DUMFAC * SLAMAX
        GROMAX = 0.0
        SIZRAT = SIZELF / SIZREF

        DO I = 1,6
          YVGROW(I) = SIZRAT * YVREF(I)
        ENDDO

!-----------------------------------------------------------------------
!     INITIALIZE PARTITIONING PARAMETERS
!-----------------------------------------------------------------------
        FRLF = TABEX(YLEAF,XLEAF,0.0,8)
        FRSTM = TABEX(YSTEM,XLEAF,0.0,8)
        FRRT = 1.0 - FRLF - FRSTM

      ENDIF

!***********************************************************************
!***********************************************************************
!     EMERGENCE CALCULATIONS - Performed once per season upon emergence
!         or transplanting of plants
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. EMERG) THEN
!-----------------------------------------------------------------------
      IF (CONTROL % RUN == 1 .OR. INDEX('FPQ',CONTROL%RNMODE) < 0) THEN
        XFRT   = XFRUIT ! XFRT: Current day's partitioning to reproductive growth (0-1)
                        ! XFRUIT: Maximum fraction of daily available gross photosynthate (PG) 
        ADDSHL = 0.0    ! Today's growth demand for shells of age NPP (g[shell] / m2 / d)
        TURXFR = 0.0    ! Turgor water stress factor used to modify partitioning to 
        GDMSD  = 0.0    ! Seed growth demand based on temperature and photoperiod
        CDMSD  = 0.0    ! Total CH2O demand to grow seed demand (GDMSD)
        NDMSD  = 0.0    ! Total N demand to grow seed demand (GDMSD) (g[N] / m2 / d) 
        GDMSDR = 0.0    ! Potential seed growth from NDMSDR (amount of Mobilized N which 
        CDMSDR = 0.0    ! CH2O required for seed growth from mobilized N
        NDMSDR = 0.0    ! Amount of Mobilized N which can be used for seed growth
        CDMREP = 0.0    ! Total CH2O needed for potential reproductive growth
        NAGE   = 0      ! Age of cohort (d)
        DO NPP = 1,NCOHORTS
          PHTIM(NPP) = 0. !Cumulative photothermal time ages of seeds and shells 
          PNTIM(NPP) = 0. ! Photothermal days from first flower when flowers in age group I 
        END DO
        FNINSD = SDPRO * 0.16   ! FNINSD: Maximum fraction of N for growing seed tissue based on temperature (g[N] / g[seed])
        ! SDPROR    Ratio to adjust lipid and carbohydrate proportions when seed 
        !           protein differs from protein composition of standard cultivar (SDPROS) 
        FNINL  = PROLFI * 0.16  ! FNINL: Maximum fraction of N for growing leaf tissue (g[N] / g[leaf])
        ! PROLFI    Maximum protein composition in leaves during growth with 
        !           luxurious supply of N (g[protein] / g[leaf tissue])
        FNINS  = PROSTI * 0.16  ! FNINS: Maximum fraction of N for growing stem tissue (g[N] / g[stem])
        ! PROSTI    Maximum protein composition in stems during growth with 
        !           luxurious supply of N (g[protein] / g[stem])
        FNINR  = PRORTI * 0.16  ! FNINR: Maximum fraction of N for growing root tissue (g[N] / g[root])
        ! PRORTI    Maximum protein composition in roots during growth with 
        !             luxurious supply of N (g[protein] / g[root])
      ELSE
        ADDSHL = 0.0    ! Today's growth demand for shells of age NPP (g[shell] / m2 / d)
        TURXFR = 0.0    ! Turgor water stress factor used to modify partitioning to
        GDMSD  = 0.0    ! Seed growth demand based on temperature and photoperiod
        CDMSD  = 0.0    ! Total CH2O demand to grow seed demand (GDMSD)
        NDMSD  = 0.0    ! Total N demand to grow seed demand (GDMSD) (g[N] / m2 / d)
        GDMSDR = 0.0    ! Potential seed growth from NDMSDR (amount of Mobilized N which
        CDMSDR = 0.0    ! CH2O required for seed growth from mobilized N
        CDMREP = 0.0    ! Total CH2O needed for potential reproductive growth
        NAGE   = 0      ! Age of cohort (d)
        DO NPP = 1,NCOHORTS
          PHTIM(NPP) = 0. !Cumulative photothermal time ages of seeds and shells
          PNTIM(NPP) = 0. ! Photothermal days from first flower when flowers in age group I
        END DO
!        FNINSD = SDPRO * 0.16   ! FNINSD: Maximum fraction of N for growing seed tissue based on temperature (g[N] / g[seed])
!        ! SDPROR    Ratio to adjust lipid and carbohydrate proportions when seed
!        !           protein differs from protein composition of standard cultivar (SDPROS)
!        FNINL  = PROLFI * 0.16  ! FNINL: Maximum fraction of N for growing leaf tissue (g[N] / g[leaf])
!        ! PROLFI    Maximum protein composition in leaves during growth with
!        !           luxurious supply of N (g[protein] / g[leaf tissue])
!        FNINS  = PROSTI * 0.16  ! FNINS: Maximum fraction of N for growing stem tissue (g[N] / g[stem])
!        ! PROSTI    Maximum protein composition in stems during growth with
!        !           luxurious supply of N (g[protein] / g[stem])
!        FNINR  = PRORTI * 0.16  ! FNINR: Maximum fraction of N for growing root tissue (g[N] / g[root])
!        ! PRORTI    Maximum protein composition in roots during growth with
!        !             luxurious supply of N (g[protein] / g[root])
!
      ENDIF
!***********************************************************************
!***********************************************************************
!     DAILY RATE/INTEGRATION
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
!     DAS = MAX(0,TIMDIF(YRSIM,YRDOY))
      CALL GET(CONTROL)
      DAS = CONTROL % DAS
!-----------------------------------------------------------------------
!     Compute max N mining, NMINEP, based on stage-dependent mining
!     rate, NMOBR
!-----------------------------------------------------------------------
!     Assume that a Maximum Fraction (NMOBMX) of N can be Mobilized per Day
!     NVSMOB is the relative N mobil rate in veg stage, rel to reprod. stage
!-----------------------------------------------------------------------
!     9/27/95 ACCELERATE N MOBILIZATION AFTER R5, FUNCTION OF (1-SWFAC)
!     ALLOWS ACCELERATING BY 50% IF MAX DEFICIT.
!     2/6/96 SOMETIMES SEEDS FILL, XPOD IS LOW, THEN N MOBILIZATION SLOWS
!     I DON'T REALLY WANT THAT, LATE IN CYCLE.  KJB
!     NOW, DXR57 HITS CLOSE TO 1 AT MATURITY AND PREVENTS THAT
!-----------------------------------------------------------------------
      NMOBR  = NVSMOB * NMOBMX * TDUMX
      IF (DAS .GT. NR5) THEN
        NMOBR = NMOBMX * TDUMX2 * (1.0 + 0.5*(1.0 - SWFAC))
     &      * (1.0 + 0.3*(1.0 - NSTRES)) * (NVSMOB + (1. - NVSMOB)
     &      * MAX(XPOD,DXR57**2.))
      ENDIF
      NMINEP = NMOBR * (WNRLF + WNRST + WNRRT + WNRSH)

!-----------------------------------------------------------------------
      IF (DAS .GE. NR1) THEN
!-----------------------------------------------------------------------
!     Accumulate physiological age of flower (PNTIM) and pod (PHTIM) cohorts
!-----------------------------------------------------------------------
        IF (DAS - NR1 + 1 > NCOHORTS) THEN
          WRITE(MSG(1),'(A,I5)')
     &      'Number of flower cohorts exceeds maximum limit of',NCOHORTS
          CALL WARNING(1,ERRKEY,MSG)
          CALL ErrorCode(CONTROL, 100, ERRKEY, YREND)
          RETURN
        ENDIF

        IF (DAS .EQ. NR1) THEN
          PNTIM(1) = 0.0
        ELSE
          PNTIM(DAS - NR1 + 1) = PNTIM(DAS - NR1) + TDUMX
        ENDIF

        IF (DAS .LE. NR2) THEN
          PHTIM(1) = 0.0
        ELSE
          PHTIM(DAS - NR2 + 1) = PHTIM(DAS - NR2) + TDUMX
        ENDIF

!-----------------------------------------------------------------------
!     Calculate function for modifying seed growth rate with temperature
!-----------------------------------------------------------------------
        TMPFAC = 0.
        TMPFCS = 0.
        DO I = 1,24
          TMPFAC =CURV(TYPSDT,FNSDT(1),FNSDT(2),FNSDT(3),FNSDT(4),
     &                  TGRO(I))
        TMPFCS = TMPFCS + TMPFAC
        ENDDO
        TMPFAC = TMPFCS / 24.
!-----------------------------------------------------------------------
!       Calculate reduction in seed growth due to insect punctures
!-----------------------------------------------------------------------
        IF (PUNCSD .GT. 0.001) THEN
          REDPUN = 1.0 - (PUNCTR/PUNCSD) * RPRPUN
          REDPUN = MAX(0.0,REDPUN)
        ELSE
          REDPUN = 1.0
        ENDIF
!-----------------------------------------------------------------------
!       Water stress factor (TURADD) effect on reproductive growth and
!       pod addition.  Stress is defined to INCREASE growth and addition.
!-----------------------------------------------------------------------
        TURADD = TABEX (YTRFAC,XTRFAC,TURFAC,4)
!-----------------------------------------------------------------------
!     Calculate maximum growth per seed based on temp and seed punctures
!-----------------------------------------------------------------------
        SDGR = SDVAR * TMPFAC * REDPUN * (1.-(1.-DRPP)*SRMAX) *
     &       (1. + TURADD)
!-----------------------------------------------------------------------
!     Initialize Seed Growth Demands and CH2O and N required for seed
!       growth
!-----------------------------------------------------------------------
        GDMSD  = 0.0
        CDMSD  = 0.0
        NDMSD  = 0.0
        GDMSDR = 0.0
        CDMSDR = 0.0
        NDMSDR = 0.0
!-----------------------------------------------------------------------
        IF (DAS .GT. NR2) THEN
          DO NPP = 1, DAS - NR2
!-----------------------------------------------------------------------
!     Calculate physiol age of seed cohort.  Do not allow seed to grow
!     until shells are greater than LAGSD physiol age.
!-----------------------------------------------------------------------
            PAGE = PHTIM(DAS - NR2 + 1) - PHTIM(NPP)
            IF (PAGE .GE. LAGSD) THEN
!-----------------------------------------------------------------------
!     Allow cohort growth until threshing limit (seed wt./pod wt) occurs
!     taking into account damage by pests to seed and shells
!-----------------------------------------------------------------------
              REDSHL = 0
              IF (SDDES(NPP).GT.0) THEN
                REDSHL = WTSHE(NPP)*SDDES(NPP)/(SDDES(NPP)+SDNO(NPP))
              ENDIF
              SDMAX = (WTSHE(NPP)-REDSHL)*THRESH/(100.-THRESH)-WTSD(NPP)
              SDMAX = MAX(0.0,SDMAX)
!-----------------------------------------------------------------------
!     Compute Seed Growth Demand, GDMSD, and N required for seed, NDMSD
!-----------------------------------------------------------------------
              GDMSD  = GDMSD  + MIN(SDGR*SDNO(NPP)*REDPUN, SDMAX)
            ENDIF
          ENDDO
!-----------------------------------------------------------------------
!     Call seed composition routine
!-----------------------------------------------------------------------
          CALL SDCOMP(
     &      CARMIN, LIPOPT, LIPTB, PLIGSD, PMINSD, POASD, !Input
     &      RCH2O, RLIG, RLIP, RMIN, RNO3C, ROA, SDLIP,   !Input
     &      SDPRO, SLOSUM, TAVG,                          !Input
     &      AGRSD1, AGRSD2, FNINSD, POTCAR, POTLIP)       !Output

          NDMSD  = FNINSD * GDMSD
!-----------------------------------------------------------------------
!     Calculate Amount of Mobilized N Which Can be Used for Seed Growth,
!     NDMSDR, potential seed growth from this source of N, GDMSDR,
!     and CH2O required for this seed growth from mobilized N, CDMSDR
!-----------------------------------------------------------------------
          IF (NDMSD .GT. NMINEP) THEN
            NDMSDR = NMINEP
          ELSE
            NDMSDR = NDMSD
          ENDIF
          GDMSDR = NDMSDR/FNINSD
          CDMSDR = GDMSDR * (AGRSD1 + FNINSD*6.25 * RPRO)
!-----------------------------------------------------------------------
!    Compute Total CH2O Demand to Grow GDMSD g Tissue
!-----------------------------------------------------------------------
          CDMSD = (MAX(0.0,(GDMSD - GDMSDR))) * AGRSD2 + CDMSDR
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Compute max growth per shell, depending on temp, daylength
!-----------------------------------------------------------------------
      GRRAT1 = SHVAR * TMPFAC * (1.- (1.-DRPP) * SRMAX)
     & * (1.0 + TURADD)
!-----------------------------------------------------------------------
!     Initialize Shell Growth Demand, N (NDMSH) and C (CDMSH) needed for growth
!-----------------------------------------------------------------------
      GDMSH = 0.0
      NDMSH = 0.0
      CDMSH = 0.0
!-----------------------------------------------------------------------
!     Compute growth demand for shells, GDMSH, allowing slow growth
!     until LNGPEG age, then potential growth until LNGSH
!-----------------------------------------------------------------------
      IF (DAS .GT. NR2) THEN
        DO NPP = 1,DAS - NR2
          NAGE = DAS - NR2 + 1 - NPP  !NAGE not used - chp
          PAGE = PHTIM(DAS - NR2 + 1) - PHTIM(NPP)
          IF (PAGE .LE. LNGSH .AND. SHELN(NPP) .GE. 0.001 .AND.
     &       GRRAT1 .GE. 0.001) THEN
            IF (PAGE .GE. LNGPEG) THEN
              !Shells between LNGPEG and LNGSH
              ADDSHL = GRRAT1 * SHELN(NPP)
            ELSE
              !Shells < LNGPEG
              ADDSHL = GRRAT1 * SHELN(NPP) * SHLAG
            ENDIF
          ENDIF
          GDMSH  = GDMSH + ADDSHL
        ENDDO
!-----------------------------------------------------------------------
!     Compute CH2O required for the potential shell growth
!-----------------------------------------------------------------------
        CDMSH = GDMSH * AGRSH2
      ENDIF
!-----------------------------------------------------------------------
!     Compute TEMXFR, the temp effect on partitioning to pods
!     High temp would increase fraction growth to vegetative tissue
!-----------------------------------------------------------------------
      TEMXFR = 0.
      DO I = 1,24
        TEMXFR = TEMXFR + TABEX(YXFTEM,XXFTEM,TGRO(I),6)
      ENDDO
      TEMXFR = TEMXFR/24.
!-----------------------------------------------------------------------
!     Partitioning to pods is increased under drought stress conditions
!        depending on XFRMAX, an input parameter
!-----------------------------------------------------------------------
      TURXFR = XFRMAX * (1. - TURFAC)
      TURXFR = MIN(TURXFR,1.0)
      TURXFR = MAX(TURXFR,0.0)
!-----------------------------------------------------------------------
!     Night length and temperature are multiplicative
!     but turgor effect adds to the partitioning
!-----------------------------------------------------------------------
      XFRT = XFRUIT * TEMXFR + XFRUIT * TURXFR
!     XFRT = XFRUIT * RNIT * TEMXFR   !NEED TO FIX FOR DAYLENGTH EFFECT
      XFRT = MIN(XFRT,1.0)
      XFRT = MAX(XFRT,0.0)
!-----------------------------------------------------------------------
!    Total Potential Available CH2O for Reprod Growth (CAVTOT)
!    and total CH2O needed for potential reproductive growth (CDMREP)
!-----------------------------------------------------------------------
      CAVTOT = PGAVL * XFRT
      CDMREP = CDMSH + CDMSD
!-----------------------------------------------------------------------
!    Adjust C-Demand for New Growth if C-Available is Less than C Demand
!    Also adjust tissue growth demand for seeds and shells
!-----------------------------------------------------------------------
      GDMSDO = GDMSD
      IF (CDMREP .GT. CAVTOT) THEN
        IF (CDMSD .GT. CAVTOT) THEN
          CDMSH = 0.0
          GDMSH = 0.0
          CDMSD = CAVTOT
          IF (CDMSDR .GT. CAVTOT) THEN
            CDMSDR = CAVTOT
          ENDIF
          GDMSD = (MAX(0.0,(CDMSD-CDMSDR)))/AGRSD2 +
     &    CDMSDR/(AGRSD1+FNINSD*6.25*RPRO)
          NDMSDR = GDMSDR * FNINSD
        ELSE
          CDMSH = CAVTOT - CDMSD
          GDMSH = CDMSH/AGRSH2
        ENDIF
        CDMREP = CDMSD + CDMSH
      ENDIF
!-----------------------------------------------------------------------
!     Compute N demand for seed, shell, and total reproductive growth
!-----------------------------------------------------------------------
      NDMSD  = GDMSD * FNINSD
      NDMSH  = GDMSH * FNINSH
      NDMREP = NDMSD + NDMSH

!-----------------------------------------------------------------------
!     Vegetative partitioning factors and demand for C and N for new
!     growth before VSSINK, assume leaf expansion is fixed, compute
!     SLA based on function of light, temp, etc, then compute
!     FRLF (leaf partitioning), then FRRT, FRSTM
!-----------------------------------------------------------------------
!     Check to See if New Vegetative Tissue Can Be Grown, Using PGAVL
!-----------------------------------------------------------------------
      CDMVEG = MAX(0.0,(1.-XFRT)*PGAVL)
      NDMVEG = 0.0
      CDMVEG = (PGAVL * XFRT - CDMREP) + CDMVEG

!-----------------------------------------------------------------------
!       This is from documentation:  check no longer needed?? chp
!-----------------------------------------------------------------------
!      CDMVEG = MAX(0.0,(1.-XFRT)*PGAVL)
!      IF (PGAVL * XFRT .GT. CDMREP) THEN
!        IF(N .LE. NDLEAF) CDMVEG = (PGAVL * XFRT - CDMREP) + CDMVEG
!      ENDIF
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      IF (DAS .EQ. NR1) THEN
!-----------------------------------------------------------------------
!     Fraction of growth going to leaves and roots decreases
!     linearly between R1 and NDLEAF.
!-----------------------------------------------------------------------
        FRLFM  = TABEX (YLEAF, XLEAF, VSTAGE, 8)
        FRSTMM = TABEX (YSTEM, XLEAF, VSTAGE, 8)
        YY = FRLFM - FRLFF 
        XX = FRSTMM - FRSTMF
      ENDIF
!-----------------------------------------------------------------------
      IF (DAS .LT. NR1) THEN
!-----------------------------------------------------------------------
!     Calculate Pattern of Vegetative Partitioning, a function of V-STAGE
!-----------------------------------------------------------------------
        FRLF  = TABEX(YLEAF,XLEAF,VSTAGE,8)
        FRSTM = TABEX(YSTEM,XLEAF,VSTAGE,8)
      ELSE
!-----------------------------------------------------------------------
!     Partitioning between vegetative tissues depends on development
!     as expressed by FRACDN, the relative development between R1 and NDLEAF
!-----------------------------------------------------------------------
        FRLF = FRLFM - YY * FRACDN
        FRSTM = FRSTMM - XX * FRACDN
        IF ( DAS .GE. NDLEAF) THEN
          FRLF = FRLFF
          FRSTM = FRSTMF
        ENDIF
      ENDIF

!     This is where to modify partitioning for extra root growth:
!     check units!!! fraction vs percentage
!     FRLF = FRLF - FRLF/(FRLF+FRSTM) * (extra root value)
!     FRSTM= FRSTM - FRSTM/(FRLF+FRSTM) * (extra root value)
      FRRT = 1. - FRLF - FRSTM

!-----------------------------------------------------------------------
!     Compute F, specific leaf area for new leaf weight
!-----------------------------------------------------------------------
      TPHFAC = 0.
      DO I = 1,24
        TPHFAC = TPHFAC + TABEX (YSLATM,XSLATM,TGRO(I),5)
      ENDDO
      TPHFAC = TPHFAC/24.
!-----------------------------------------------------------------------
      PARSLA = (SLAMN+(SLAMX-SLAMN)*EXP(SLAPAR*PAR))/SLAMX
      TURFSL = MAX(0.1, (1.0 - (1.0 - TURFAC)*TURSLA))
!-----------------------------------------------------------------------
!     Compute overall effect of TMP, PAR, water stress on SLA (F), first
!     for veg stages, then transition to rep stage from R1 to end leaf
!     effect of PAR on SLA, COX PEANUT SCI. 5:27, 1978
!-----------------------------------------------------------------------
      FFVEG = FVEG * TPHFAC * PARSLA * TURFSL

      F = FFVEG
      IF (XFRT*FRACDN .GE. 0.05) F = FFVEG * (1.0 - XFRT * FRACDN)
!-----------------------------------------------------------------------
!     For determinate plants (XFRUIT=1.) leaf expansion stops at NDLEAF
!-----------------------------------------------------------------------
      IF (XFRUIT .GT. 0.9999 .AND. DAS .GE. NDLEAF) F = 0.0

!-----------------------------------------------------------------------
!     During early vegetative growth, leaf area expansion depends on
!     VSTAGE (Prior to VSSINK).  This sets FRLF, partitioning of d.m.
!     to leaves.  FRRT and FRSTM are then computed by left over C.  When
!     an upper limit of d.m. goes to leaves, leaf area expansion is
!     restricted so that F is maintained as computed and minimal amounts
!     of C is partitioned to FRSTM and FRRT  (JWJ 4/1/96)
!-----------------------------------------------------------------------
      IF (VSTAGE .LT. VSSINK) THEN
        GROYES = GROMAX
        GROMAX = TABEX(YVGROW,XVGROW,VSTAGE,6) * SIZELF/SIZREF
        GAINNW = (GROMAX - GROYES) * PLTPOP
!-----------------------------------------------------------------------
!     CALCULATE MINIMUM WEIGHT NEEDED TO ADD GAINNW LEAF AREA/M2,
!     AND AMOUNT OF LEAF WEIGHT WHICH CAN BE GROWN WITH PG AVAILABLE
!-----------------------------------------------------------------------
        IF (F .GT. 1.E-5) THEN
          GAINWT = GAINNW/F
        ELSE
          GAINWT = 0.0
        ENDIF
!-----------------------------------------------------------------------
!     Compute fraction of C partitioned to leaves, based on F, VSSINK
!     Limit leaf pertitioning to FRLFMX (i.e., FRLFMX = 0.7)
!-----------------------------------------------------------------------
        FRLF = (AGRLF*GAINWT)/(CDMVEG + 0.0001)
        IF (FRLF .GT. FRLFMX) THEN
          GAINWT = (CDMVEG/AGRLF) * FRLFMX
          GAINNW = GAINWT * F
          FRLF = FRLFMX
        ENDIF
!-----------------------------------------------------------------------
!     Recompute FRSTM and FRRT based on FRLF
!-----------------------------------------------------------------------
        FRSTM = (1. - FRLF) * FRSTM / (FRSTM + FRRT)
        FRRT  = 1. - FRLF - FRSTM
!-----------------------------------------------------------------------
      ENDIF
!-----------------------------------------------------------------------
!     Compute CH2O cost per g of tissue, excluding cost for protein (AGRVG)
!     and total CH2O cost per g of veg tissue (AGRVG2)
!-----------------------------------------------------------------------
      AGRVG = AGRLF * FRLF + AGRRT * FRRT + AGRSTM * FRSTM
      AGRVG2 = AGRVG + (FRLF*PROLFI+FRRT*PRORTI+FRSTM*PROSTI)*RPROAV
!-----------------------------------------------------------------------
!    Compute N Demand for New Tissue, including reproductive and vegetative
!-----------------------------------------------------------------------
      NDMVEG = (CDMVEG/AGRVG2) * (FRLF*FNINL+FRSTM*FNINS+
     &   FRRT*FNINR)
      NDMNEW = NDMREP + NDMVEG
!-----------------------------------------------------------------------
!    Check to See if Any C is Left After Reproductive Growth for
!    Reducing N to Re-Fill Old Tissue, if N Can Be Taken up by Roots
!-----------------------------------------------------------------------
      CNOLD = MAX(0.0,PGAVL-CDMREP)
      NDMOLD = 0.0
!-----------------------------------------------------------------------
!    Nitrogen Demand for Old Tissue
!-----------------------------------------------------------------------
      IF (DAS .GT. NVEG0 .AND. DAS .LT. NR7 .AND.
     &          CNOLD .GT. 0.0) THEN
        NVSTL = FNINL
        NVSTS = FNINS
        NVSTR = FNINR
        IF (DXR57 .GT.0.0) THEN
           FRNLFT = (NRCVR + (1. - NRCVR) * (1. - DXR57**2))
           NVSTL = PROLFF*0.16 + (FNINL-PROLFF*0.16) * FRNLFT
           NVSTS = PROSTF*0.16 + (FNINS-PROSTF*0.16) * FRNLFT
           NVSTR = PRORTF*0.16 + (FNINR-PRORTF*0.16) * FRNLFT
        ENDIF
        NDMOLD = (WTLF  - WCRLF) * MAX(0.0,(NVSTL - PCNL /100.))
     &         + (STMWT - WCRST) * MAX(0.0,(NVSTS - PCNST/100.))
     &         + (RTWT  - WCRRT) * MAX(0.0,(NVSTR - PCNRT/100.))
        IF (NDMOLD .GT. (CNOLD/RNO3C*0.16)) THEN
          NDMOLD = CNOLD/RNO3C*0.16
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!    Total N Demand
!-----------------------------------------------------------------------
      NDMTOT = NDMREP + NDMVEG + NDMOLD
!-----------------------------------------------------------------------
!    Compute Total Demand for C, and Max. C that Could be Mined
!     CDMTOT not used - chp
!-----------------------------------------------------------------------
      CDMTOT = CDMREP + CDMVEG + NDMOLD*RNO3C/0.16 
      GDMSD = GDMSDO
!-----------------------------------------------------------------------
!    At this point, PGAVL will be used entirely, assuming that N can be
!    made available in the ratio described.
!     Growth Demands : GDMSD, GDMSH
!     N-Demands      : NDMREP, NDMVEG, NDMOLD, NDMTOT, NDMNEW
!     C-Demands      : CDMREP, CDMVEG, CDMTOT, CNOLD

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!-----------------------------------------------------------------------

      RETURN
      END SUBROUTINE TG_DEMAND
!=======================================================================

