/*
** Options for  tokin  Grid.
**
** Application flag:   tokin
** Input script:       ocean_tokin.in
*/
#define ROMS_MODEL
#define WF_SHELF

#define BOUNDARY_ALLREDUCE /* use mpi_allreduce in mp_boundary */
#undef  COLLECT_ALLGATHER  /* use mpi_allgather in mp_collect  */
#define COLLECT_ALLREDUCE  /* use mpi_allreduce in mp_collect  */
#define REDUCE_ALLGATHER   /* use mpi_allgather in mp_reduce   */
#undef  REDUCE_ALLREDUCE   /* use mpi_allreduce in mp_reduce   */

#define WET_DRY
#undef STATIONS
#define	UV_ADV
#define	UV_COR
#define	UV_VIS2
#define UV_SMAGORINSKY
#define	TS_DIF2
#define TS_SMAGORINSKY
#define SOLVE3D
#define	SALINITY
#define TEMPERATURE
#define	NONLIN_EOS
#define UV_LOGDRAG 	  
#undef FLOATS
!#define FLOAT_STICKY
#undef FLOAT_VWALK
#undef VWALK_FORWARD

/*#define TS_U3HADVECTION And #define TS_SVADVECTION replaced by  TS_MPDATA  */
#define DJ_GRADPS
#define MIX_GEO_UV
#define MIX_GEO_TS
#define CURVGRID
#define SPHERICAL
# define SPLINES_VVISC
# define SPLINES_VDIFF
#define MASKING
!#define LIMIT_BSTRESS

#define T_PASSIVE

#undef SEDIMENT
#ifdef SEDIMENT
!# define BEDLOAD_MPM
# define ANA_BPFLUX
# define ANA_SPFLUX
# undef BEDLOAD_SOULSBY
!# define RIVER_SEDIMENT /* Process river sediment point-sources */
!# define SED_DENS /* Activate sediment to affect equation of state */
# define SED_MORPH /* allow bottom model elevation to evolve */
# define SUSPLOAD /* Activate suspended load transport */
# define SED_TAU_CD_CONST   /* use constant critical stress for deposition */
#endif

#define BULK_FLUXES
#ifdef   BULK_FLUXES
# define LONGWAVE_OUT
# define ATM_PRESS
# undef EMINUSP
# define ANA_RAIN    /* could undef this and use rain data from moorings */
# define ANA_CLOUD
# define WTYPE_GRID
# define ANA_WTYPE 
# define SOLAR_SOURCE  /* add shortwave with exp(z) profile */
# define WIND_MINUS_CURRENT
!# define QCORRECTION         /* use if net heat flux correction */
!# define SCORRECTION         /* use if freshwater flux correction */ 
#else
# define ANA_STFLUX
!# define ANA_SMFLUX
#endif

#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX

!#define MY25_MIXING
!# define LMD_MIXING
#define  GLS_MIXING
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
!#  define CRAIG_BANNER
#  define CHARNOK
!#define LIMIT_VDIFF
!#define LIMIT_VVISC
# endif

# ifdef LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_DDMIX
#  undef  LMD_NONLOCAL
# endif

#define RADIATION_2D
 
#undef Add_Tide

#ifdef Add_Tide

!#define RAMP_TIDES
#define SSH_TIDES        
#define UV_TIDES        

#if defined SSH_TIDES || defined UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
#endif

#endif
 
#define	AVERAGES         
#define AVERAGES_FLUXES    /* for full diagnostics of mean and eddy heat fluxes */
!#define DIAGNOSTICS_UV
!#define DIAGNOSTICS_TS

