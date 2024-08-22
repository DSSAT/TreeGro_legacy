##########################################################################################
# CMake project file for DSSAT-CSM
##########################################################################################
# Define the project and the dependencies that it has
##########################################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
   message(FATAL_ERROR "DSSAT in-source builds are not permitted. \nMake a separate folder"
   						"for building:\nmkdir build; cd build; cmake ..\nBefore that,"
   						"remove the files already created:\n"
   						"CMakeCache.txt and CMakeFiles\n")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

PROJECT(DSSAT-CSM Fortran)

##########################################################################################
# DSSAT_VERSION
##########################################################################################

SET(MAJOR 4)
SET(MINOR 5)
SET(MODEL 0)
SET(BUILD 0)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%H
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH_ABBR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

MESSAGE(STATUS "MAJOR: ${MAJOR} MINOR: ${MINOR} MODEL: ${MODEL} BUILD: ${BUILD} COMMIT: ${GIT_COMMIT_HASH_ABBR}")
MESSAGE(STATUS "BRANCH: ${BRANCH}")
MESSAGE(STATUS "COMMIT: ${GIT_COMMIT_HASH}")
# Set the version
SET(VERSION ${MAJOR}.${MINOR}.${MODEL})

IF(APPLE)
    SET(ENV{MACOSX_DEPLOYMENT_TARGET} 10.10)
ENDIF(APPLE)

##########################################################################################

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Organize output files.  On Windows this also keeps .dll files next
# to the .exe files that need them, making tests easy to run.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Keep Uncomment if Fortran 90 support is required
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

##########################################################################################
# Define the actual files and folders that make the build
##########################################################################################

# Define the executable name and OSDefinitions file
SET(EXECUTABLE_NAME dscsm0${MAJOR}${MINOR})


# Add the source files
set( 
  source_list 
 "ASMDM.FOR"
 "AUTHAR.FOR"
 "AUTPLT.FOR"
 "CANOPY.FOR"
 "Capillary.for"
 "CE_RATIO_C.FOR"
 "CENTURY.for"
 "CHEMICAL.FOR"
 "CO2VAL.FOR"
 "CROPGRO.for"
 "CRSIMDEF.FOR"
 "CSCER040.FOR"
 "CSCERES_Interface.for"
 "CSCRP_Interface.for"
 "CSCRP045.FOR"
 "CSDIS045.FOR"
 "CSM.for"
 "CSP_CANOPY.FOR"
 "CSP_CASUPRO.FOR"
 "CSP_DEMAND.FOR"
 "CSP_GROW.FOR"
 "CSP_GROW_CANE.FOR"
 "CSP_HRes.for"
 "CSP_INCOMP.FOR"
 "CSP_INCOMP_OUT.FOR"
 "CSP_INPHENOL.FOR"
 "CSP_IPDMND_OUT.FOR"
 "CSP_IPPHENOL.FOR"
 "CSP_IPPLNT.FOR"
 "CSP_MOBIL.FOR"
 "CSP_NUPTAK.FOR"
 "CSP_OPGROW.FOR"
 "CSP_OPHARV.FOR"
 "CSP_PHENOL.FOR"
 "CSP_PHOTO.for"
 "CSP_RESPIR.FOR"
 "CSP_ROOTS.FOR"
 "CSP_SENES.FOR"
 "CSP_VEGGR.FOR"
 "CSREA045.FOR"
 "CSUTS045.FOR"
 "DATES.FOR"
 "DECRAT_C.FOR"
 "DEMAND.FOR"
 "Drainage_2D.FOR"
 "Drip_Irrig.FOR"
 "EFLOW_C.FOR"
 "EQUIL2.FOR"
 "ERROR.FOR"
 "ESR_SoilEvap.for"
 "ESR_SoilEvap_2D.for"
 "ETPHOT.for"
 "ETPHR.for"
 "FCHEM.for"
 "Fert_Place.for"
 "Flood_Chem.for"
 "Flood_Irrig.FOR"
 "FLOODI.FOR"
 "FREEZE.FOR"
 "FreshWt.FOR"
 "GROW.FOR"
 "HMET.for"
 "HRes_CGRO.for"
 "HResCeres.for"
 "IMMOBLIMIT_C.for"
 "INCOMP.FOR"
 "INCORPOR_C.for"
 "INFIL.for"
 "Info.FOR"
 "input_sub.for"
 "INSOIL.FOR"
 "INSTGE.FOR"
 "INTRO.FOR"
 "INVAR.FOR"
 "IPCHEM.FOR"
 "IPECO.FOR"
 "IPENV.FOR"
 "ipexp.for"
 "IPHedley_C.for"
 "IPHedley_inorg.for"
 "IPIBS.FOR"
 "IPMAN.FOR"
 "IPPARM.FOR"
 "IPPEST.FOR"
 "Ipphenol.for"
 "IPPLNT.FOR"
 "IPPROG.FOR"
 "ipsim.for"
 "IPSLIN.FOR"
 "IPSOIL.FOR"
 "IPSOIL_Inp.FOR"
 "IPTILL.FOR"
 "IPVAR.FOR"
 "IPWTH_alt.for"
 "IRRIG.FOR"
 "LAND.FOR"
 "LINDM.FOR"
 "LITDEC_C.FOR"
 "LMATCH.FOR"
 "MgmtOps.for"
 "ML_CERES.for"
 "ML_GROSUB.FOR"
 "ML_NFACT.FOR"
 "ml_NUPTAK.FOR"
 "ML_opharv.FOR"
 "ML_PHASEI.FOR"
 "ML_PHENOL.FOR"
 "ML_rootgr.for"
 "ML_TILLSUB.FOR"
 "MOBIL.FOR"
 "ModuleDefs.for"
 "MULCHEVAP.FOR"
 "MULCHLAYER.FOR"
 "MULCHWAT.FOR"
 "MZ_CERES.FOR"
 "MZ_GROSUB.for"
 "MZ_IX_GROSUB.for"
 "MZ_IX_KNUMBER.for"
 "MZ_IX_LEAFAREA.for"
 "MZ_IX_NUPTAK.FOR"
 "MZ_IX_PHENOL.for"
 "MZ_IX_PHOTSYNT.for"
 "MZ_IX_PLANTG.for"
 "MZ_IX_RADABS.for"
 "MZ_IX_RESPIR.FOR"
 "MZ_KUPTAK.FOR"
 "MZ_NFACTO.for"
 "MZ_NUPTAK.FOR"
 "MZ_OPGROW.FOR"
 "MZ_OPHARV.FOR"
 "MZ_OPNIT.FOR"
 "MZ_PHENOL.for"
 "MZ_ROOTS.FOR"
 "NCHECK_C.FOR"
 "NCHECK_inorg.FOR"
 "NCHECK_inorg_2D.FOR"
 "NCHECK_organic.FOR"
 "NFIX.FOR"
 "NFLUX.FOR"
 "NFLUX_2D.FOR"
 "NUPTAK.FOR"
 "NUPTAK_2D.FOR"
 "OM_Place.for"
 "OPETPHOT.FOR"
 "OpFlood.for"
 "OPFLOODN.FOR"
 "OPGEN.FOR"
 "OPGeneric.FOR"
 "OPGeneric2.FOR"
 "Opgrow.for"
 "OPHARV.FOR"
 "OPHEAD.FOR"
 "OPMULCH.FOR"
 "OPPEST.FOR"
 "OpPlantP.for"
 "OpSoilKi.for"
 "OPSOILNI.FOR"
 "OPSOILNI_2D.FOR"
 "OpSoilOrg.for"
 "OpSoilPi.for"
 "OPSOMLIT_C.FOR"
 "OpStemp.for"
 "OPSTRESS.FOR"
 "OPSUM.for"
 "OPTEMPXY2K.FOR"
 "optempy2k.for"
 "OPVIEW.FOR"
 "OPWBAL.for"
 "OPWBAL_2D.for"
 "OPWEATH.FOR"
 "OXLAYER.FOR"
 "P_CASUPRO.FOR"
 "P_CERES.FOR"
 "P_CGRO.FOR"
 "P_IPPLNT.FOR"
 "p_plant.for"
 "P_Uptake.for"
 "Paddy_Mgmt.for"
 "PARTIT_C.FOR"
 "PATH.FOR"
 "PEST.FOR"
 "PESTCP.FOR"
 "PET.FOR"
 "PHENOL.FOR"
 "PHOTO.FOR"
 "plant.for"
 "PlantNBal.FOR"
 "PODDET.FOR"
 "PODS.FOR"
 "PPlantSubs.for"
 "PT_GROSUB.for"
 "PT_NFACTO.for"
 "PT_NUPTAK.for"
 "PT_OPGROW.for"
 "PT_OPHARV.for"
 "PT_PHASEI.for"
 "PT_PHENOL.for"
 "PT_ROOTGR.for"
 "PT_ROOTGR_2D.for"
 "PT_SUBSTOR.FOR"
 "PT_THTIME.for"
 "Read_IFile.f90"
 "READS.FOR"
 "RESPIR.FOR"
 "RETC_VG.for"
 "RI_Calcshk.for"
 "RI_GNURSE.FOR"
 "RI_Grosub.for"
 "RI_Ipcrop.for"
 "RI_KUPTAK.FOR"
 "RI_Nfacto.for"
 "RI_Nuptak.for"
 "RI_Opgrow.for"
 "RI_Opharv.for"
 "RI_Phenol.for"
 "RI_Rootgr.for"
 "RI_Tillsub.for"
 "RI_Transpl_g.for"
 "RI_Transpl_p.for"
 "RICE.for"
 "RNOFF.for"
 "ROOTDM.FOR"
 "ROOTS.for"
 "ROOTS_2D.for"
 "RootSoilVol.for"
 "ROOTWU.FOR"
 "ROOTWU_2D.FOR"
 "ROOTWU_HR.FOR"
 "RPLACE_C.FOR"
 "RStages.for"
 "RunList.for"
 "SAL_Stemp.for"
 "SATFLO.for"
 "SC_CanesimCanopy.for"
 "SC_Canop3.for"
 "SC_CNG_mods.for"
 "SC_CNGRO.FOR"
 "SC_COEFFS.FOR"
 "SC_OPHARV.FOR"
 "SC_OUTPUT.for"
 "SC_PARTIT.for"
 "SC_PHENOL.for"
 "SC_PHOTOS.for"
 "SC_Poplt3.for"
 "SC_ROOTG.for"
 "SDCOMP.FOR"
 "SECLI.FOR"
 "SECROP.FOR"
 "SEEDDM.FOR"
 "SEFERT.FOR"
 "SEFLD.FOR"
 "SEFREQ.FOR"
 "SEHARV.FOR"
 "SEINIT.FOR"
 "SEIRR.FOR"
 "SENES.FOR"
 "SENESADD_C.for"
 "SENS.FOR"
 "SEPEST.FOR"
 "SEPLT.FOR"
 "SERES.FOR"
 "SESIM.FOR"
 "SESOIL.FOR"
 "SETIME.FOR"
 "SEVAR.FOR"
 "SEWTH.FOR"
 "SG_CERES.for"
 "SG_GROSUB.FOR"
 "SG_NFACT.FOR"
 "sg_NUPTAK.FOR"
 "SG_OPHARV.FOR"
 "SG_PHASEI.FOR"
 "SG_PHENOL.FOR"
 "SG_ROOTGR.for"
 "SLigCeres.for"
 "SOIL.FOR"
 "SoilCBal_C.for"
 "SoilCellInit_2D.F90"
 "SoilCellUtils_2D.F90"
 "SoilCNPinit_C.FOR"
 "SOILDYN.FOR"
 "SOILEV.for"
 "SoilK_init.FOR"
 "SoilKi.for"
 "SoilMixing.for"
 "SoilNBalSum.for"
 "SOILNI.for"
 "SOILNI_2D.for"
 "SoilNi_init.for"
 "SoilNi_init_2D.for"
 "SoilNiBal.for"
 "SoilNiBal_2D.for"
 "SoilNoBal.FOR"
 "SoilNoBal_C.FOR"
 "SoilNoPoBal.FOR"
 "SoilOrg.FOR"
 "SoilOrg_init.for"
 "SoilPBalSum.for"
 "SoilPi.for"
 "SoilPi_init.FOR"
 "soilpibal.for"
 "SoilPoBal.FOR"
 "SoilPoBal_C.FOR"
 "SOLAR.FOR"
 "SOMDEC_C.FOR"
 "SOMFIX_C.FOR"
 "SOMINIT_c.for"
 "SOMLITPRINT_C.for"
 "SPAM.for"
 "SPAM_2D.for"
 "SPSUBS.for"
 "STEMP.for"
 "SUBSL2.for"
 "SUWCMS2.for"
 "SW_FreshWt.for"
 "sw_GROSUB.for"
 "SW_Sensor.for"
 "TextureClass.for"
 "TG_CANOPY.FOR"
 "TG_DEMAND.FOR"
 "TG_GROW.FOR"
 "TG_Opgrow.for"
 "TG_OPHARV.FOR"
 "TG_PHENOL.FOR"
 "TG_PODDET.FOR"
 "TG_PODS.FOR"
 "TG_ROOTS.for"
 "TG_ROOTS_2D.for"
 "TG_RStages.for"
 "TG_SENES.FOR"
 "TG_TREEGRO.for"
 "TG_VEGGR.FOR"
 "TILEDRAIN.FOR"
 "Tillage.for"
 "TillEvent.for"
 "TR_Calcshk.for"
 "TR_Grosub.for"
 "TR_Ipcrop.for"
 "TR_Nfacto.for"
 "TR_Nuptak.for"
 "TR_OPGROW.FOR"
 "TR_Opharv.for"
 "TR_Phenol.for"
 "TR_Rootgr.for"
 "TR_SUBSTOR.FOR"
 "TR_Tillsub.for"
 "TR_Transpl_g.for"
 "TR_Transpl_p.for"
 "trans.for"
 "TSOMLIT_C.for"
 "UTILS.FOR"
 "VEGDM.FOR"
 "VEGGR.FOR"
 "Warning.FOR"
 "WATBAL.FOR"
 "WatBal2D.for"
 "WaterTable.F90"
 "WaterTable_2D.F90"
 "WBAL.for"
 "WBAL_2D.for"
 "WBAL_2D_ts.for"
 "WBSUBS.for"
 "WEATHR.for"
 "WEATHR_Inp.FOR"
 "WGEN.for"
 "WTHMOD.for"
 "WTHSET.FOR"
)
ADD_EXECUTABLE(${EXECUTABLE_NAME} ${source_list})

foreach(source IN LISTS source_list)
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" source_path_msvc "${source_path}")
    source_group("${source_path_msvc}" FILES "${source}")
endforeach()


install(TARGETS ${EXECUTABLE_NAME} DESTINATION
	RUNTIME DESTINATION .)

install(FILES Utilities/run_dssat DESTINATION .
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	)

install(DIRECTORY Data/ DESTINATION .
	PATTERN "*.in" EXCLUDE)

##################################################
# Resume information about flags and option used #
##################################################

message( "-- Flags" )
message( "   FFLAGS       ${CMAKE_Fortran_FLAGS}" )
message( "   RELEASE      ${CMAKE_Fortran_FLAGS_RELEASE}" )
message( "   DEBUG        ${CMAKE_Fortran_FLAGS_DEBUG}")
message( "   LINKER       ${CMAKE_EXE_LINKER_FLAGS}")
message( "-- Build Info" )
message( "   BUILD TYPE   ${CMAKE_BUILD_TYPE}" )
message( "   VERSION      ${MAJOR}.${MINOR}.${MODEL}.${BUILD}" )
message( "   I. PREFIX    ${CMAKE_INSTALL_PREFIX}" )
message( "   Executable   ${EXECUTABLE_NAME}" )
