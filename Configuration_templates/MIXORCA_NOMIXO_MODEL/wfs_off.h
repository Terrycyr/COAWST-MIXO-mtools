/*
 * ** svn $Id: ias.h 830 2017-01-24 21:21:11Z arango $
 * *******************************************************************************
 * ** Copyright (c) 2002-2018 The ROMS/TOMS Group                               **
 * **   Licensed under a MIT/X style license                                    **
 * **   See License_ROMS.txt                                                    **
 * *******************************************************************************
 * **
 * ** Options for NGOM 1/25 DEG resolution.
 * **
 * ** Application flag:   oil_03
 * ** Input script:       ocean_oil_03.in
 * */

#undef WRF_MODEL
#undef SWAN_MODEL
#define ROMS_MODEL

# define BOUNDARY_ALLREDUCE /* use mpi_allreduce in mp_boundary */
# undef COLLECT_ALLGATHER  /* use mpi_allgather in mp_collect  */
# define COLLECT_ALLREDUCE  /* use mpi_allreduce in mp_collect  */
# undef REDUCE_ALLGATHER   /* use mpi_allgather in mp_reduce   */
# define REDUCE_ALLREDUCE   /* use mpi_allreduce in mp_reduce   */

/*
 *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *   * **  Offline settings
 *    *    *    *    * **-----------------------------------------------------------------------------
 *     *     *     *     * */
#define OFFLINE
#ifdef OFFLINE
#define OUT_DOUBLE
#define OCLIMATOLOGY
#define ATCLIMATOLOGY 
/* best results with MPDATA offline */
!#define TS_U3HADVECTION
!#define TS_C4VADVECTION
#define TS_MPDATA

/* For inclusion of mixing variables from history file, 
 *  *  * you need to either define the following individual mixing variables 
 *   *   * or define AKXCLIMATOLOGY and/or MIXCLIMATOLOGY
 *    *    *
 *     *     * AKTCLIMATOLOGY and AKSCLIMATOLOGY and AKVCLIMATLOGY define AKXCLIMATOLOGY
 *      *      * TKECLIMATOLOGY and GLSCLIMATOLOGY define MIXCLIMATOLOGY
 *       *      */

#undef MIXCLIMATOLOGY
#undef AKXCLIMATOLOGY
#undef AKVCLIMATLOGY
#undef AKTCLIMATOLOGY
#define AKSCLIMATOLOGY    /* get best performance forcing realistic Aks */


/* for offline, turn off forcings (bulk forcing undefined elsewhere in file)
 * All forcing is coming in through climatology from online case. */

/*
 *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   * **  Adding offline passive tracers
 *    *    *    * **-----------------------------------------------------------------------------
 *     *     *     * */
#undef OFFLINE_TPASSIVE
#ifdef OFFLINE_TPASSIVE
# define T_PASSIVE
# define ANA_BPFLUX
# define ANA_SPFLUX
#endif

/*
 *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *   * **  Adding offline floats
 *    *    *    *    * **-----------------------------------------------------------------------------
 *     *     *     *     * */
#undef OFFLINE_FLOATS
#undef FLOAT_VWALK
#ifdef OFFLINE_FLOATS
# define FLOATS
# undef FLOAT_OIL
# undef WOIL_INTEGRATED
# undef OIL_DEBUG
# undef OIL_EULR
#endif

#define OFFLINE_BIOLOGY
#ifdef OFFLINE_BIOLOGY
# define BIOLOGY
# define BIO_RIVER
#endif

#endif


#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#undef  NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */
#undef  VERIFICATION            /* NL Observation Verification Driver */
#undef  NORMALIZATION           /* Background error Covariance Normalization */
#undef  AD_SENSITIVITY          /* Adjoint Sensitivity Driver */

/*
 * **-----------------------------------------------------------------------------
 * **  Nonlinear basic state settings.
 * **-----------------------------------------------------------------------------
 * */
#undef  AVERAGES               /*Write out time-averaged data*/
#undef  AVERAGES_FLUXES

#define DIAGNOSTICS_BIO
#define DIAGNOSTICS_TS
#define DIAGNOSTICS

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define TS_DIF2
#define SOLVE3D
#define SALINITY
#define TEMPERATURE
#define NONLIN_EOS
#define UV_LOGDRAG        
#define DJ_GRADPS
#define MIX_GEO_UV
#define MIX_GEO_TS
#define CURVGRID
#define SPHERICAL
#define MASKING

#undef LIMIT_VDIFF

#undef LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#define  MY25_MIXING
#ifdef MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES
#endif

/*
 *  *  * **-----------------------------------------------------------------------------
 *   *   * **  MIXO-RCA model settings.
 *    *    * **-----------------------------------------------------------------------------
 *     *     * */
#undef OPT_IF
#define  BIO_MIXORCA
#ifdef BIO_MIXORCA

/*
 *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *  * **    Basic setting
 *    *    *    *  * **-----------------------------------------------------------------------------
 *     *     *     *  * */
# define ANA_SPFLUX
# define ANA_BPFLUX
# undef ANA_SPONGE
# undef VERTICAL_MIGRATION
#  ifdef VERTICAL_MIGRATION
#  define VERTICAL_MIGRATION_MIXO
#  endif
# define SAL_CDOM
# undef REAERATION_COSINE
# define SSC_OFFLINE
# undef COLD_MORTALITY
# define NUTR_MORTALITY


# define PHYT_G1
#  ifdef PHYT_G1
#  define QUADRATIC_GRAZ_P1
#  endif

# undef PHYT_G2
#  ifdef PHYT_G2
#  define QUADRATIC_GRAZ_P2
#  endif

# undef PHYT_G3
#  ifdef PHYT_G3
#  define QUADRATIC_GRAZ_P3
#  endif

/*
 *  *  *  *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *   *   *   *  * **    Sediment model setting
 *    *    *    *    *    *    *  * **-----------------------------------------------------------------------------
 *     *     *     *     *     *     *  * */
# define MIXO_SED
# undef ZERO_SED_FLUX
# undef NO_SED_NO23

/*
 *  *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *   *  * **    Kelvin Flynn's model setting
 *    *    *    *    *  * **-----------------------------------------------------------------------------
 *     *     *     *     *  * */
# define MIXOTROPHY
# undef MIXO_TURB
# if defined MIXO_TURB
#  define TKECLIMATOLOGY
# endif
# undef SBARRIER_WEAK_MIXO
# define SBARRIER_STRONG_MIXO
# define QUADRATIC_GRAZ_PREY
# define QUADRATIC_GRAZ_MIXO
# undef PAT_KB_INGEST
# undef FAN_LIGHT_MODEL
# define SIMPLE_MX_GRAZ
/*
 *  *  *  *  *  *  *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *   *   *   *   *   *   *  * **    Zooplankton model setting
 *    *    *    *    *    *    *    *    *    *  * **-----------------------------------------------------------------------------
 *     *     *     *     *     *     *     *     *     *  * */
# undef ZOOPLANKTON


/*
 *  *  *  *  *  *  *  *  * **-----------------------------------------------------------------------------
 *   *   *   *   *   *   *   *  * **    Carbon model setting
 *    *    *    *    *    *    *    *  * **-----------------------------------------------------------------------------
 *     *     *     *     *     *     *     *  * */
# undef MIXO_CARBON

#endif

#undef  BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE_OUT
# define ANA_RAIN
# define ANA_CLOUD
#else
# undef  QCORRECTION
# undef  SOLAR_SOURCE
# undef  DIURNAL_SRFLUX
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
#endif

#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

#undef FORWARD_WRITE
#undef FORWARD_READ
#undef FORWARD_MIXING

/*
 * **-----------------------------------------------------------------------------
 * **  Variational Data Assimilation.
 * **-----------------------------------------------------------------------------
 * */

#ifdef NORMALIZATION
# undef  MULTIPLE_TLM
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  W4DVAR
# undef  R_SYMMETRY
# define CORRELATION
# undef  CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  TLM_CHECK
# undef  BALANCE_OPERATOR
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#if defined IS4DVAR || defined IS4DVAR_OLD
# undef  MULTIPLE_TLM
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  W4DVAR
# undef  R_SYMMETRY
# undef  CORRELATION
# undef  CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  TLM_CHECK
# undef  BALANCE_OPERATOR
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#ifdef W4DVAR
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  IS4DVAR
# undef  R_SYMMETRY
# undef  CORRELATION
# define CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# define RPM_RELAXATION
# undef  TLM_CHECK
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#ifdef W4DPSAS
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  IS4DVAR
# undef  R_SYMMETRY
# undef  CORRELATION
# define CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  TLM_CHECK
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#ifdef SANITY_CHECK
# define FULL_GRID
# define FORWARD_READ
# define FORWARD_WRITE
# define FORWARD_MIXING
# define OUT_DOUBLE
# define ANA_PERTURB
# define ANA_INITIAL
#endif
