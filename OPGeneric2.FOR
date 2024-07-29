C=======================================================================
C  OPGENERIC2, Subroutine
C
C  Generates output for simulated data, up to 10 variables
C-----------------------------------------------------------------------
C  Revision history
C
C  09/22/2008 CHP Written
C=======================================================================

      SUBROUTINE OPGENERIC2 (CONTROL, 
     &    NVars, Width, HeaderTxt, FormatTxt,  
     &    VAR1, VAR2, VAR3, VAR4, VAR5, VAR6, VAR7, VAR8, VAR9, VAR10)

C-----------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT  NONE
      SAVE

      TYPE (ControlType), INTENT(IN) :: CONTROL
      CHARACTER*80      , INTENT(IN) :: FormatTxt
      CHARACTER*120     , INTENT(IN) :: HeaderTxt
      INTEGER           , INTENT(IN) :: NVars
      REAL              , INTENT(IN) :: Var1, Var2, Var3, Var4, Var5
      REAL              , INTENT(IN) :: Var6, Var7, Var8, Var9, Var10

      CHARACTER*12, PARAMETER :: OUTG = 'Generic2.OUT'
      CHARACTER*95 FMT_STRING
      INTEGER DAS, DOY, DYNAMIC, ERRNUM, INCDAT, LUN
      INTEGER RUN, Width, YEAR, YRDOY
      LOGICAL FEXIST

      DAS     = CONTROL % DAS
      DYNAMIC = CONTROL % DYNAMIC
      RUN     = CONTROL % RUN
      YRDOY   = CONTROL % YRDOY

!***********************************************************************
!***********************************************************************
!     Seasonal Initialization - Called once per season
!***********************************************************************
      IF (DYNAMIC == SEASINIT) THEN
!-----------------------------------------------------------------------
        CALL GETLUN('Generic2', LUN)
      
        INQUIRE (FILE = OUTG, EXIST = FEXIST)
        IF (FEXIST) THEN
          OPEN (UNIT = LUN, FILE = OUTG, STATUS = 'OLD',
     &      IOSTAT = ERRNUM, POSITION = 'APPEND')
        ELSE
          OPEN (UNIT = LUN, FILE = OUTG, STATUS = 'NEW',
     &      IOSTAT = ERRNUM)
          WRITE(LUN,'("*Generic2 daily output")')
        ENDIF

!       Write headers
        CALL HEADER(SEASINIT, LUN, RUN)
        WRITE (LUN,'(A,A)') "@YEAR DOY   DAS", HeaderTxt(1:Width)

        FMT_STRING = "(1X,I4,1X,I3.3,1X,I5)"  

!***********************************************************************
!***********************************************************************
!     Daily OUTPUT
!***********************************************************************
      ELSEIF (DYNAMIC == OUTPUT .OR. 
     &        DYNAMIC == INTEGR .OR. 
     &        DYNAMIC == RATE) THEN

        CALL YR_DOY(YRDOY, YEAR, DOY)
        CALL YR_DOY(INCDAT(YRDOY,-1), YEAR, DOY)
        WRITE (LUN,FMT_STRING,ADVANCE='NO') YEAR, DOY, DAS
        SELECT CASE (NVars)
        CASE(1); WRITE(LUN,FormatTxt) VAR1
        CASE(2); WRITE(LUN,FormatTxt) VAR1,VAR2
        CASE(3); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3
        CASE(4); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4
        CASE(5); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4,VAR5
        CASE(6); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4,VAR5,VAR6
        CASE(7); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4,VAR5,VAR6,VAR7
        CASE(8); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4,VAR5,VAR6,
     &                                VAR7,VAR8
        CASE(9); WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4,VAR5,VAR6,
     &                                VAR7,VAR8,VAR9
        CASE(10);WRITE(LUN,FormatTxt) VAR1,VAR2,VAR3,VAR4,VAR5,VAR6,
     &                                VAR7,VAR8,VAR9,VAR10
        END SELECT

!***********************************************************************
!***********************************************************************
!     Daily OUTPUT
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASEND) THEN
!***********************************************************************
        CLOSE (LUN)

!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE OPGENERIC2
