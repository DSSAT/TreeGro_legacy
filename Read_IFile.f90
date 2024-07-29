          SUBROUTINE Read_IFile (CONTROL,                &    
          LAII, LWAI, SWAI, RWAI, PWAI, RDPI, CHTI, CWIDI, GSTI ) 

      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      IMPLICIT     NONE

      CHARACTER*1   PLDS
      CHARACTER*2   CROP
      CHARACTER*6   VARTY
     ! CHARACTER*10  STNAME(20) 
      CHARACTER*12  FILEIB, FILEX   
      CHARACTER*20  VARNAME
      CHARACTER*30  FILEIO 
      CHARACTER*78  MSG(2)
      CHARACTER*80  PATHEX
      CHARACTER*92  PATHIB
      CHARACTER*128 CHAR

      CHARACTER*6  ERRKEY, SECTION
      PARAMETER (ERRKEY = 'READ_I')

      INTEGER ISECT, LINC, LNUM, LUNIO, ERR, FOUND, STGDOY(20), TRT
      REAL, INTENT(OUT), OPTIONAL :: LAII, LWAI, SWAI, RWAI, PWAI
      REAL, INTENT(OUT), OPTIONAL :: RDPI, GSTI, CHTI, CWIDI 
      REAL LAII0, GSTI0, LWAI0,SWAI0,RWAI0,PWAI0,RDPI0,CHTI0,CWIDI0
      LOGICAL GoodFile

!     The variable "CONTROL" is of type "ControlType".
      TYPE (ControlType) CONTROL
      
      FILEIO  = CONTROL % FILEIO
      LUNIO   = CONTROL % LUNIO
      !FILEX = CONTROL % FILEX 
      GoodFile = .TRUE.
      LNUM = 0
!-----------------------------------------------------------------------
!       Read data from FILEIO for use in PLANT module
!-----------------------------------------------------------------------
      OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,0)
     ! LNUM = LNUM + 3  
!-----------------------------------------------------------------------
!    Read name of FILEX - to get name of intial biomass file (if present)
!-----------------------------------------------------------------------
      READ (LUNIO,'(3(/),15X,A12,1X,A80)',IOSTAT=ERR) FILEX, PATHEX
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)
     
!-----------------------------------------------------------------------
!   Find and read initial biomass file
    FILEIB = FILEX(1:11) // "I"
    PathIB = TRIM(PATHEX) // FILEIB
    INQUIRE (FILE = PATHIB, EXIST = GoodFile)

    OPEN (LUNIO, FILE = PathIB, STATUS = 'OLD', IOSTAT=ERR)
    IF (ERR .NE. 0) GoodFile = .FALSE.

    IF (GoodFile) THEN
      SECTION = '*EXP. DATA (I)'
      CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
      IF (FOUND .EQ. 0) GoodFile = .FALSE.
    ENDIF

    DO WHILE (GoodFile)
      CALL IGNORE (LUNIO,LNUM,ISECT,CHAR)
      IF (ISECT .NE. 1) THEN
        GoodFile = .FALSE.; EXIT
      ENDIF

      READ(CHAR,'(I2,9F8.0)', IOSTAT=ERR) TRT,     &
        LAII0, GSTI0, LWAI0, SWAI0, RWAI0, PWAI0, RDPI0, CHTI0, CWIDI0
      IF (ERR .NE. 0) THEN
        GoodFile = .FALSE.; EXIT
      ENDIF

      IF (TRT == CONTROL % TRTNUM) then
        IF (PRESENT(LAII)) LAII = LAII0
        IF (PRESENT(GSTI)) GSTI = GSTI0
        IF (PRESENT(LWAI)) LWAI = LWAI0
        IF (PRESENT(SWAI)) SWAI = SWAI0
        IF (PRESENT(RWAI)) RWAI = RWAI0
        IF (PRESENT(PWAI)) PWAI = PWAI0
        IF (PRESENT(RDPI)) RDPI = RDPI0
        IF (PRESENT(CHTI)) CHTI = CHTI0
        IF (PRESENT(CWIDI)) CWIDI = CWIDI0
        EXIT
      Endif
    ENDDO

!   No valid initial biomass file - use default values depending
!       on planting method
    IF (.NOT. GoodFile) THEN
      CALL ERROR(ERRKEY,0,"I file is missing",0)
      !CALL ERROR(ERRKEY,ErrCode,' ',0)
      WRITE(MSG(1),100) 
      WRITE(MSG(2),101) 
  100 FORMAT('Error in Read_IFile, I file is missing ')
  101 FORMAT('Program will stop.')
      CALL WARNING(2, 'Read_IFile', MSG)
      WRITE(*,'(A78)') MSG(1)
      STOP
    ENDIF

       RETURN
       END SUBROUTINE Read_IFile
!=====================================================================
  MODULE Interface_Read_IFile
!     Interface needed for optional arguments with IPSOIL
      INTERFACE
        SUBROUTINE Read_IFile (CONTROL,                &    
          LAII, LWAI, SWAI, RWAI, PWAI, RDPI, CHTI, CWIDI, GSTI )   !Output
          USE ModuleDefs
          TYPE (ControlType), INTENT(IN) :: CONTROL
          REAL, INTENT(OUT), OPTIONAL :: LAII, LWAI, SWAI, RWAI, PWAI
          REAL, INTENT(OUT), OPTIONAL :: RDPI, CHTI, CWIDI, GSTI 
        END SUBROUTINE
      END INTERFACE
      END MODULE Interface_Read_IFile
!========================================================================
!     Read_IFile VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! CHTI = Initial canopy height m
! CWIDI= Initial canopy width m
! GSTI   Initial development stage. It is VSTAGE (i.e. DVSI)
! LAII   Initial leaf area per plant (m2 plant-1) (LAPE)
! LWAI  Initial leaf weight" (kg/ha) (WLVGI)
! RWAI   Initial root weight (kg/ha) (WRTI)
! PWAI   Initial weight storage organs (kg/ha) (WSOI) (pod weight)
! SWAI   Initial stem weight (kg/ha) (i.e. WSTI)
! RDPI   Initial root depth (m) (ZRTI)
!-----------------------------------------------------------------------
!     END SUBROUTINE Read_IFile
!=======================================================================  
