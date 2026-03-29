dset ^fm_grav2d_50m.dat
options  byteswapped
undef 1.e30
title  OUTPUT FROM WRF V4.2 MODEL
xdef  1024 linear 0 0.0001
ydef     1 linear 0 0.0001
zdef  128 linear 1 1  
tdef    7 linear 00Z01JAN2000       1MN      
VARS  104
XLAT           1  0  LATITUDE, SOUTH IS NEGATIVE (degree_north)
XLONG          1  0  LONGITUDE, WEST IS NEGATIVE (degree_east)
LU_INDEX       1  0  LAND USE CATEGORY (-)
VAR_SSO        1  0  variance of subgrid-scale orography (m2)
PFM          128  0  pressure (Pa)
PRESSURE_B   128  0  BASE pressure (Pa)
RHOFM        128  0  DENSITY (Kg m-3)
RHO_BASE     128  0  BASE DENSITY (Kg m-3)
ZZ           128  0  height values on full (w) levels (-)
U            128  0  x-wind component (m s-1)
V            128  0  y-wind component (m s-1)
W            128  0  z-wind component (m s-1)
PH           128  0  perturbation geopotential (m2 s-2)
PHB          128  0  base-state geopotential (m2 s-2)
T            128  0  perturbation potential temperature theta-t0 (K)
THM          128  0  either 1) pert moist pot temp=(1+Rv/Rd Qv)*(theta)-T0, or 2) pert dry pot temp=t (K)
MU             1  0  perturbation dry air mass in column (Pa)
MUB            1  0  base state dry air mass in column (Pa)
NEST_POS       1  0  - (-)
P            128  0  perturbation pressure (Pa)
PB           128  0  BASE STATE PRESSURE (Pa)
ZHALF        128  0  height values on half (prho) levels (-)
ZFULL        128  0  height values on full (w) levels (-)
P_HYD        128  0  hydrostatic pressure (Pa)
Q2             1  0  QV at 2 M (kg kg-1)
T2             1  0  TEMP at 2 M (K)
TH2            1  0  POT TEMP at 2 M (K)
PSFC           1  0  SFC PRESSURE (Pa)
U10            1  0  U at 10 M (m s-1)
V10            1  0  V at 10 M (m s-1)
AREA2D         1  0  Horizontal grid cell area, using dx, dy, and map factors (m2)
DX2D           1  0  Horizontal grid distance: sqrt(area2d) (m)
QVAPOR       128  0  Water vapor mixing ratio (kg kg-1)
SHDMAX         1  0  ANNUAL MAX VEG FRACTION (-)
SHDMIN         1  0  ANNUAL MIN VEG FRACTION (-)
SNOALB         1  0  ANNUAL MAX SNOW ALBEDO IN FRACTION (-)
TSLB           5  0  SOIL TEMPERATURE (K)
SMOIS          5  0  SOIL MOISTURE (m3 m-3)
SH2O           5  0  SOIL LIQUID WATER (m3 m-3)
SEAICE         1  0  SEA ICE FLAG (-)
XICEM          1  0  SEA ICE FLAG (PREVIOUS STEP) (-)
SFROFF         1  0  SURFACE RUNOFF (mm)
UDROFF         1  0  UNDERGROUND RUNOFF (mm)
IVGTYP         1  0  DOMINANT VEGETATION CATEGORY (-)
ISLTYP         1  0  DOMINANT SOIL CATEGORY (-)
VEGFRA         1  0  VEGETATION FRACTION (-)
GRDFLX         1  0  GROUND HEAT FLUX (W m-2)
ACGRDFLX       1  0  ACCUMULATED GROUND HEAT FLUX (J m-2)
ACSNOM         1  0  ACCUMULATED MELTED SNOW (kg m-2)
SNOW           1  0  SNOW WATER EQUIVALENT (kg m-2)
SNOWH          1  0  PHYSICAL SNOW DEPTH (m)
CANWAT         1  0  CANOPY WATER (kg m-2)
SSTSK          1  0  SKIN SEA SURFACE TEMPERATURE (K)
COSZEN         1  0  COS of SOLAR ZENITH ANGLE (dimensionless)
LAI            1  0  LEAF AREA INDEX (m-2/m-2)
VAR            1  0  OROGRAPHIC VARIANCE (-)
MAPFAC_M       1  0  Map scale factor on mass grid (-)
MAPFAC_MX      1  0  Map scale factor on mass grid, x direction (-)
MAPFAC_MY      1  0  Map scale factor on mass grid, y direction (-)
MF_VX_INV      1  0  Inverse map scale factor on v-grid, x direction (-)
F              1  0  Coriolis sine latitude term (s-1)
E              1  0  Coriolis cosine latitude term (s-1)
SINALPHA       1  0  Local sine of map rotation (-)
COSALPHA       1  0  Local cosine of map rotation (-)
HGT            1  0  Terrain Height (m)
TSK            1  0  SURFACE SKIN TEMPERATURE (K)
RAINC          1  0  ACCUMULATED TOTAL CUMULUS PRECIPITATION (mm)
RAINSH         1  0  ACCUMULATED SHALLOW CUMULUS PRECIPITATION (mm)
RAINNC         1  0  ACCUMULATED TOTAL GRID SCALE PRECIPITATION (mm)
SNOWNC         1  0  ACCUMULATED TOTAL GRID SCALE SNOW AND ICE (mm)
GRAUPELNC      1  0  ACCUMULATED TOTAL GRID SCALE GRAUPEL (mm)
HAILNC         1  0  ACCUMULATED TOTAL GRID SCALE HAIL (mm)
CLDFRA       128  0  CLOUD FRACTION (-)
SWDOWN         1  0  DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE (W m-2)
GLW            1  0  DOWNWARD LONG WAVE FLUX AT GROUND SURFACE (W m-2)
SWNORM         1  0  NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT) (W m-2)
OLR            1  0  TOA OUTGOING LONG WAVE (W m-2)
ALBEDO         1  0  ALBEDO (-)
CLAT           1  0  COMPUTATIONAL GRID LATITUDE, SOUTH IS NEGATIVE (degree_north)
ALBBCK         1  0  BACKGROUND ALBEDO (-)
EMISS          1  0  SURFACE EMISSIVITY (-)
NOAHRES        1  0  RESIDUAL OF THE NOAH SURFACE ENERGY BUDGET (W m{-2})
TMN            1  0  SOIL TEMPERATURE AT LOWER BOUNDARY (K)
XLAND          1  0  LAND MASK (1 FOR LAND, 2 FOR WATER) (-)
UST            1  0  U* IN SIMILARITY THEORY (m s-1)
PBLH           1  0  PBL HEIGHT (m)
HFX            1  0  UPWARD HEAT FLUX AT THE SURFACE (W m-2)
QFX            1  0  UPWARD MOISTURE FLUX AT THE SURFACE (kg m-2 s-1)
LH             1  0  LATENT HEAT FLUX AT THE SURFACE (W m-2)
ACHFX          1  0  ACCUMULATED UPWARD HEAT FLUX AT THE SURFACE (J m-2)
ACLHF          1  0  ACCUMULATED UPWARD LATENT HEAT FLUX AT THE SURFACE (J m-2)
SNOWC          1  0  FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER) (-)
SR             1  0  fraction of frozen precipitation (-)
PCB            1  0  base state dry air mass in column (Pa)
PC             1  0  perturbation dry air mass in column (Pa)
LANDMASK       1  0  LAND MASK (1 FOR LAND, 0 FOR WATER) (-)
LAKEMASK       1  0  LAKE MASK (1 FOR LAKE, 0 FOR NON-LAKE) (-)
SST            1  0  SEA SURFACE TEMPERATURE (K)
SST_INPUT      1  0  SEA SURFACE TEMPERATURE FROM WRFLOWINPUT FILE (K)
pressure     128  0  Model pressure (hPa)
height       128  0  Model height (km)
tk           128  0  Temperature (K)
tc           128  0  Temperature (C)
theta        128  0  Potential Temperature (K)
ENDVARS
@ global String comment TITLE =  OUTPUT FROM WRF V4.2 MODEL
@ global String comment START_DATE = 0001-01-01_00:00:00
@ global String comment SIMULATION_START_DATE = 0001-01-01_00:00:00
@ global String comment WEST-EAST_GRID_DIMENSION =  1025
@ global String comment SOUTH-NORTH_GRID_DIMENSION =     3
@ global String comment BOTTOM-TOP_GRID_DIMENSION =   129
@ global String comment DX =        50.00
@ global String comment DY =        50.00
@ global String comment AERCU_OPT =     0
@ global String comment AERCU_FCT =         1.00
@ global String comment IDEAL_CASE =     6
@ global String comment DIFF_6TH_SLOPEOPT =     0
@ global String comment AUTO_LEVELS_OPT =     2
@ global String comment DIFF_6TH_THRESH =         0.10
@ global String comment DZBOT =        50.00
@ global String comment DZSTRETCH_S =         1.30
@ global String comment DZSTRETCH_U =         1.10
@ global String comment SKEBS_ON =     0
@ global String comment SPEC_BDY_FINAL_MU =     1
@ global String comment USE_Q_DIABATIC =     0
@ global String comment GRIDTYPE = C
@ global String comment DIFF_OPT =     1
@ global String comment KM_OPT =     1
@ global String comment DAMP_OPT =     0
@ global String comment DAMPCOEF =         0.03
@ global String comment KHDIF =        75.00
@ global String comment KVDIF =        75.00
@ global String comment MP_PHYSICS =     0
@ global String comment RA_LW_PHYSICS =     0
@ global String comment RA_SW_PHYSICS =     0
@ global String comment SF_SFCLAY_PHYSICS =     0
@ global String comment SF_SURFACE_PHYSICS =     0
@ global String comment BL_PBL_PHYSICS =     0
@ global String comment CU_PHYSICS =     0
@ global String comment SF_LAKE_PHYSICS =     0
@ global String comment SURFACE_INPUT_SOURCE =     3
@ global String comment SST_UPDATE =     0
@ global String comment GRID_FDDA =     0
@ global String comment GFDDA_INTERVAL_M =     0
@ global String comment GFDDA_END_H =     0
@ global String comment GRID_SFDDA =     0
@ global String comment SGFDDA_INTERVAL_M =     0
@ global String comment SGFDDA_END_H =     0
@ global String comment HYPSOMETRIC_OPT =     1
@ global String comment USE_THETA_M =     1
@ global String comment GWD_OPT =     0
@ global String comment SF_URBAN_PHYSICS =     0
@ global String comment SF_SURFACE_MOSAIC =     0
@ global String comment SF_OCEAN_PHYSICS =     0
@ global String comment SHCU_PHYSICS =     0
@ global String comment MFSHCONV =     0
@ global String comment FEEDBACK =     1
@ global String comment SMOOTH_OPTION =     2
@ global String comment SWRAD_SCAT =         1.00
@ global String comment W_DAMPING =     0
@ global String comment DT =         0.50
@ global String comment RADT =         0.00
@ global String comment BLDT =         0.00
@ global String comment CUDT =         0.00
@ global String comment AER_OPT =     0
@ global String comment SWINT_OPT =     0
@ global String comment AER_TYPE =     1
@ global String comment AER_AOD550_OPT =     1
@ global String comment AER_ANGEXP_OPT =     1
@ global String comment AER_SSA_OPT =     1
@ global String comment AER_ASY_OPT =     1
@ global String comment AER_AOD550_VAL =         0.12
@ global String comment AER_ANGEXP_VAL =         1.30
@ global String comment AER_SSA_VAL =         0.85
@ global String comment AER_ASY_VAL =         0.90
@ global String comment MOIST_ADV_OPT =     1
@ global String comment SCALAR_ADV_OPT =     1
@ global String comment TKE_ADV_OPT =     1
@ global String comment DIFF_6TH_OPT =     0
@ global String comment DIFF_6TH_FACTOR =         0.12
@ global String comment OBS_NUDGE_OPT =     0
@ global String comment BUCKET_MM =        -1.00
@ global String comment BUCKET_J =        -1.00
@ global String comment PREC_ACC_DT =         0.00
@ global String comment ISFTCFLX =     0
@ global String comment ISHALLOW =     0
@ global String comment ISFFLX =     1
@ global String comment ICLOUD =     1
@ global String comment ICLOUD_CU =     0
@ global String comment TRACER_PBLMIX =     1
@ global String comment SCALAR_PBLMIX =     0
@ global String comment YSU_TOPDOWN_PBLMIX =     0
@ global String comment GRAV_SETTLING =     0
@ global String comment DFI_OPT =     0
@ global String comment NTASKS_X =    -1
@ global String comment NTASKS_Y =    -1
@ global String comment NTASKS_TOTAL =     1
@ global String comment SIMULATION_INITIALIZATION_TYPE = IDEALIZED DATA
@ global String comment WEST-EAST_PATCH_START_UNSTAG =     1
@ global String comment WEST-EAST_PATCH_END_UNSTAG =  1024
@ global String comment WEST-EAST_PATCH_START_STAG =     1
@ global String comment WEST-EAST_PATCH_END_STAG =  1025
@ global String comment SOUTH-NORTH_PATCH_START_UNSTAG =     1
@ global String comment SOUTH-NORTH_PATCH_END_UNSTAG =     2
@ global String comment SOUTH-NORTH_PATCH_START_STAG =     1
@ global String comment SOUTH-NORTH_PATCH_END_STAG =     3
@ global String comment BOTTOM-TOP_PATCH_START_UNSTAG =     1
@ global String comment BOTTOM-TOP_PATCH_END_UNSTAG =   128
@ global String comment BOTTOM-TOP_PATCH_START_STAG =     1
@ global String comment BOTTOM-TOP_PATCH_END_STAG =   129
@ global String comment GRID_ID =     1
@ global String comment PARENT_ID =     0
@ global String comment I_PARENT_START =     1
@ global String comment J_PARENT_START =     1
@ global String comment PARENT_GRID_RATIO =     1
@ global String comment CEN_LAT =         0.00
@ global String comment CEN_LON =         0.00
@ global String comment TRUELAT1 =         0.00
@ global String comment TRUELAT2 =         0.00
@ global String comment MOAD_CEN_LAT =         0.00
@ global String comment STAND_LON =         0.00
@ global String comment POLE_LAT =         0.00
@ global String comment POLE_LON =         0.00
@ global String comment GMT =         0.00
@ global String comment JULYR =     0
@ global String comment JULDAY =     1
@ global String comment MAP_PROJ =     0
@ global String comment MAP_PROJ_CHAR = Cartesian
@ global String comment MMINLU =
@ global String comment NUM_LAND_CAT =    21
@ global String comment ISWATER =     0
@ global String comment ISLAKE =     0
@ global String comment ISICE =     0
@ global String comment ISURBAN =     0
@ global String comment ISOILWATER =     0
@ global String comment HYBRID_OPT =     0
@ global String comment ETAC =         0.00
