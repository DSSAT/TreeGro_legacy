!  CapillaryRitchie.for 
!
!  SUBROUTINE:
!  Capillary - Entry point of console application.
! The upflow will go to the line 475 of WATBAL.for
! Based on the following paper
! Ayars, J.E., Christen, E.W., Soppe, R.W., Meyer, W.S., 2006. The resource potential 
! of in-situ shallow ground water use in irrigated agriculture: a review. 
! Irrigation Sci. 24, 147-160.

!****************************************************************************
!
!  PROGRAM: Capillary Raise
!****************************************************************************

      SUBROUTINE CapillaryRitchie_2D(DYNAMIC, MgmtWTD, SOILPROP, SW,  !Input 
     &      Cells, FLOWUPd, CapiCell)                                  !Output
      USE Cells_2D
      USE ModuleDefs
      USE ModuleData
      implicit none
      SAVE

      Type(CellType), DIMENSION(MaxRows,MaxCols) :: CELLS
      INTEGER, INTENT(IN) :: DYNAMIC
      REAL, INTENT(IN) :: MgmtWTD
      TYPE (SoilType), INTENT(IN) :: SOILPROP
      REAL, DIMENSION(NL), INTENT(IN) :: SW
      REAL, INTENT(OUT) :: FLOWUPd
      Integer I, J, i1, i2, NLAYR, LIMIT_2D, FurCol1, FurRow1
      Real DUL(NL), LL(NL), DS(NL) 
      REAL DLAYR(NL),  SAT(NL), DULavg, SWavg, SATavg, LLavg, ET
      REAL, DIMENSION(MaxRows,MaxCols) :: SWV,  CapiCell 
      REAL, DIMENSION(MaxCols) :: SWBedCOLavg, SWFurCOLavg 
      REAL, DIMENSION(MaxCols) :: SWBedCOLavgOrd, SWFurCOLavgOrd
      REAL, DIMENSION(MaxCols) :: CapiBedCol, CapiFurCol
      REAL ActWTD, a, b, c, Zmax, Zr, RTDEP, x
      REAL CapiBed, CapiFur, BEDWD, BEDHT, ROWSPC, temp
      REAL, DIMENSION(NL) :: ThetaCap
      REAL, DIMENSION(MaxRows,MaxCols) :: Width
!***********************************************************************
!***********************************************************************
! Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      
      NLAYR  = SOILPROP % NLAYR
      DLAYR  = SOILPROP % DLAYR  
      LL     = SOILPROP % LL
      SAT    = SOILPROP % SAT 
      DUL    = SOILPROP % DUL    
      DS    = SOILPROP % DS  
      BEDWD = BedDimension % BEDWD 
      BEDHT = BedDimension % BEDHT 
      ROWSPC = BedDimension%ROWSPC_cm
      FurCol1 = BedDimension % FurCol1
      FurRow1 = BedDimension % FurRow1
      Width = CELLS % Struc % Width
      
	a = 3.9
      b = 3.8
      c = 0.5 
      
      FLOWUPd = 0.0

!***********************************************************************
! DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
      FLOWUPd = 0.0
      DULavg = 0.0
      LLavg = 0.0
      SWavg = 0.0
      SATavg = 0.0
      SWV = CELLS%STATE%SWV
      LIMIT_2D =  BedDimension % LIMIT_2D ! Jin Wu is it changed daily?
      CALL WTDEPT(NLAYR, DLAYR, DS, DUL, SAT, SW,                !Input
     &    ActWTD) 
      IF (MgmtWTD > 9999.) RETURN

      DO I = 1, BedDimension % LIMIT_2D !NLAYR
        ! if (DLAYR(I) . GT. MgmtWTD) exit
         DULavg = DULavg + DUL(I) * DLAYR(I)
         LLavg = LLavg + LL(I) * DLAYR(I)    
         SWavg = SWavg + SW(I) * DLAYR(I)
         SATavg = SATavg + SAT(I) * DLAYR(I)
      END DO
      DULavg = DULavg/DS(NLAYR)
      LLavg = LLavg/DS(NLAYR)
      SWavg = SWavg/DS(NLAYR)
      SATavg = SATavg/DS(NLAYR)
    !  LIMIT_2D
       ! ThetaCap_top
      ! ThetaCap_botum
      ! ThetaCap_top(L+1) not equal to ThetaCap_botum(L)
      CALL GET('SPAM', 'ET',  ET)  ! GET('SPAM', 'EP',  EP) ! ET is refer to crop ET  
      ! Make sure to set CALL PUT('PLANT', 'RTDEP',  RTDEP) in CropGro.for   
      CALL GET('PLANT', 'RTDEP',  RTDEP)
        x = (SATavg - SWavg)/(SATavg - LLavg)
        Zmax=100*(-278.2*DULavg**3 + 33.6 * DULavg **2 + 37.53 * DULavg)
        Zr = ActWTD -  RTDEP/3
    !  write(*,*)"Zmax=",Zmax, " Zr=",Zr," ActWD=",ActWTD," RTDEP=",RTDEP
       !denominator = EXP(b * (ZR/Zmax) * (1+ EXP(c/(x+0.01))
        FLOWUPd=ET* a /(exp( b * Zr / Zmax )*( 1 + exp( c / (x + .01))))
        !  mm/d=mm/d
       Write(*,*)"FLOWUPd is", FLOWUPd, " ET=", ET, " x=", x ! upflow is about 10-9 if water table depth is 80, ET, a,b,Zr,Zmax,c,x !x=0, ET=0?
    !    Write(*,*) "DULavg=",DULavg, " SATavg=",SATavg,"LLavg=", LLavg !(ETPHOT put ET & spam.for )
        !CellFLOWUPmin =CELLS % Struc % Width * FLOWUP/24/60/10 //? /(bedwidth+furrowwidth)
        !WidthFrac(Row,Col)
        ! cm2/min = cm * cm/min
        ! DeltaT is in min
        !================================================================
        ! Calculatuion the distribution to each column(the model of distribution used here should be improved later)
        ! The reference paper for FLOWUPd is not for bed system. We propose some simple modify here 
        ! Capillary rise distribution between bed and furrow
        ! Jin Wu, assume that the water table surface is below bed
        ! CapiBed * BEDWD  + CapiFur * (ROWSPC - BEDWD) = FLOWUPd * ROWSPC
        ! CapiBed / CapiFur = MgmtWTD/(MgmtWTD -  BEDHT)
         temp = BEDWD + (ROWSPC- BEDWD) *(MgmtWTD -  BEDHT)/MgmtWTD
         CapiBed = FLOWUPd * ROWSPC /temp
         CapiFur = CapiBed * (MgmtWTD -  BEDHT) /MgmtWTD 
         !  mm/d = mm/d
         ! CapiFur = FLOWUPd * ROWSPC* (MgmtWTD -  BEDHT) /(temp*MgmtWTD)
        ! CellStruc(FurRow1,J) % Width 
         ! Capillary rise distribution between column
         
        SWBedCOLavg = 0.
        SWFurCOLavg = 0.
        DO J = 1, NColsTot
        ! Jin Wu, if LIMIT_2D is above water table layer, then we should not use LIMIT_2D
          DO I = 1, LIMIT_2D !NLAYR
            If (J. LT. FurCol1) then
              SWBedCOLavg(J) = SWBedCOLavg(J)  +  SWV(I, J)
            else
              SWFurCOLavg(J) = SWFurCOLavg(J)  +  SWV(I, J)
            endif
          EndDo
          if (J. LT. FurCol1) then
            SWBedCOLavg(J) =  SWBedCOLavg(J)/(FurCol1-1)
          else
            SWFurCOLavg(J) =  SWFurCOLavg(J)/(NColsTot - FurCol1 + 1)
          endif
        EndDo
      
        call MySort(SWBedCOLavg, SWBedCOLavgOrd, FurCol1-1 )
        call MySort(SWFurCOLavg, SWFurCOLavgOrd, NColsTot-FurCol1 )
        ! Make order from small to large of column soil water average
        ! Distribute more capillary rise for drier column whil keep the total capillary rised water unchanged
        CapiBedCol = 0.
        CapiFurCol = 0.
        ! Jin Wu: need debug here. when/why got SWBedCOLavgOrd(J)=0
        DO J = 1, FurCol1-1
            CapiBedCol(SWBedCOLavgOrd(J)) = CapiBed *
     &                SWBedCOLavg(SWBedCOLavgOrd(FurCol1-1-J))/SWavg
            !        mm/d = mm/d    
            ! here we assume that both bed and Furrow has same SATavg
         EndDo
         DO J = FurCol1, NColsTot
            !        mm/d = mm/d    
            CapiFurCol(SWFurCOLavgOrd(J)) = CapiFur *
     &          SWFurCOLavg(SWFurCOLavgOrd(NColsTot-FurCol1-J))/SWavg
            ! here we assume that both bed and Furrow has same SATavg
         EndDo
        !===============================================================
        ! Calculatuion the distribution to each layer (the model of distribution used here should be improved later)
         !   Calculate water content within capillary fringe (for water table)
         CALL CapFringe(            !& 
     &        MgmtWTD,  SOILPROP,   !&    !Input
     &        ThetaCap)   
         DO J = 1, NColsTot  
            CapiCell(LIMIT_2D, J) = CapiBedCol(J)
            if (J. LT. FurCol1) then
               i1 = 1
            else  
               i1 = FurRow1
            Endif
            DO I = LIMIT_2D, i1
    !            CapiCell(I-1, J) = CapiBedCell(I,J)-(SAT(I)-SWV(I,J))
    ! &                                / CELLS % Struc % Width
                CapiCell(I-1, J) = 
     &              CapiCell(I,J)-(ThetaCap(I)-SWV(I,J)) / Width(I,J)
            Enddo
!If the capillary rise water can not be used up until come to the soil surface
! We remove these water from the FLOWUPd
! Jin Wu question, if divided the capillary into time step, these unused water may be able to use
            if (CapiCell(1, J). GT. 0.) then
              DO I = LIMIT_2D, 1
                 CapiCell(I, J) = CapiCell(I,J)- CapiCell(1, J)
              enddo
              if (J. LT. FurCol1) then
                 CapiBedCol(J) = CapiBedCol(J) - CapiCell(1, J)
              else  
                 CapiFurCol(J) = CapiFurCol(J) - CapiCell(1, J)
              Endif
              
              FLOWUPd = FLOWUPd - CapiCell(1, J)
            endif
         EndDo
       Write(*,*) "Cell Capi",CapiBed, CapiBedCol(J),CapiCell(2,1)
         
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE CapillaryRitchie_2D
C=====================================================================
!     Capillary_2D VARIABLE DEFINITIONS: (updated Dec 2009)
!-----------------------------------------------------------------------
! ActWTD    Depth to water table (cm)
! a          Regression coefficients  
! b          Regression coefficients 
! c          Regression coefficients  
! DLAYR(L)   Soil thickness in layer L (cm)
! DS(L)      Depth to bottom of soil layer L	cm
! ET         Total daily Crop ET(mm/d)
! FLOWUP     Upflux (mm/day)  
! FLOWUPd    Daily Upflux (mm/day)  
! NL         Defines the total number of soil layers (puddled and nonpuddled) and has a maximum value of 10 (including NLPUD). 
!                For example, NL and NLPUD can be set at 8 and 3, respectively, that is, a soil profile with three 
!                puddled soil layers (of which the third represents the plow sole) and five layers in the nonpuddled subsoil.
! RTDEP     Root depth (cm)
! x          The relative water content described by the relation
! ZR         Depth from 1/3rd of the depth of the root zone to the ground water level (cm), 
! Zmax       The threshold water table depth below which upflow would be less than 1 mm/day as defined by Talsma (1963) (cm)
!-----------------------------------------------------------------------
!     END SUBROUTINE Capillary_2D
!=======================================================================
    