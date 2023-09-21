      PROGRAM create_SSC
      USE netcdf
!
!
!  Imported and exported variable declarations. 
!
      INTEGER ::   ncID, varID, DAY, YEAR, DAYS, REC_INT
      INTEGER ::   ncid1,I,J,K
      INTEGER ::   Dimid1, Dimid2,Dimid3,Dimid4,Dimid5
      INTEGER ::   Dimid6, Dimid7,Dimid8,Dimid9,Dimid10
      INTEGER ::   Varid1, Varid2, Varid3, Varid4, Varid5
      INTEGER, PARAMETER :: NX=456,NY=184,NZ=21,NT=48
      CHARACTER(LEN=100) :: SUBNA, HISNA, AVGNA, filename
      CHARACTER(LEN=5) :: str, str2, str3
!
!  Variables in ROMS NetCDF file
!
      REAL*8     ::  pm(NX,NY)
      REAL*8     ::  pn(NX,NY)
      REAL*8     ::  h_roms(NX,NY)
      REAL*8     ::  masku(NX,NY)
      REAL*8     ::  maskv(NX,NY)
      REAL*8     ::  maskr(NX,NY)

      REAL*8     ::  ocean_time(NT),ocean_time_ref
      REAL*8     ::  mud_01(NX,NY,NZ,NT)
      REAL*8     ::  mud_02(NX,NY,NZ,NT)
      REAL*8     ::  SSC(NX,NY,NZ,NT)
    
      REAL*8     :: dzero

      dzero = 243.
      YEAR = 2022
      DAYS = 365-dzero
      REC_INT = 1
!-----------------------------------------------------------------------
! Define output files
!-----------------------------------------------------------------------
      SUBNA = 'CREATE'

!     Define Dimensions
      WRITE(str2,"(I4.4)") YEAR

      filename =  './WFS_'//TRIM(ADJUSTL(str2))//'_SSC.nc'

      status = nf90_create(trim(adjustl(filename))
     . , nf90_clobber,ncid1)
      CALL nccheck_status(status,'CREATE',SUBNA)

      status = nf90_def_dim(ncid = ncid1, name = "xi_rho",len = NX,
     .dimid = DimID1)
      CALL nccheck_status(status,'DIM1',SUBNA)
      status = nf90_def_dim(ncid = ncid1, name = "xi_u",len = NX-1,
     .dimid = DimID2)
      CALL nccheck_status(status,'DIM2',SUBNA)
      status = nf90_def_dim(ncid = ncid1, name = "xi_v",len = NX,
     .dimid = DimID3)
      CALL nccheck_status(status,'DIM3',SUBNA)

      status = nf90_def_dim(ncid = ncid1, name = "eta_rho",len = NY,
     .dimid = DimID4)
      CALL nccheck_status(status,'DIM4',SUBNA)
      status = nf90_def_dim(ncid = ncid1, name = "eta_u",len = NY,
     .dimid = DimID5)
      CALL nccheck_status(status,'DIM5',SUBNA)
      status = nf90_def_dim(ncid = ncid1, name = "eta_v",len = NY-1,
     .dimid = DimID6)
      CALL nccheck_status(status,'DIM6',SUBNA)

      status = nf90_def_dim(ncid = ncid1, name = "s_rho",len = NZ,
     .dimid = DimID7)
      CALL nccheck_status(status,'DIM7',SUBNA)
      status = nf90_def_dim(ncid = ncid1, name = "s_w",len = NZ+1,
     .dimid = DimID8)
      CALL nccheck_status(status,'DIM8',SUBNA)

      status = nf90_def_dim(ncid = ncid1, name = "ocean_time",
     .len = NT*DAYS,dimid = DimID9)
      CALL nccheck_status(status,'DIM9',SUBNA)

      status = nf90_def_dim(ncid = ncid1, name = "SSC_time",
     .len = NT*DAYS,dimid = DimID10)
      CALL nccheck_status(status,'DIM10',SUBNA)

!     Define variables
      status = nf90_def_var(ncid = ncid1, name = "ocean_time",
     .xtype = nf90_double,
     .dimids = (/DimID9/),varID = VarID1)
      CALL nccheck_status(status,'VAR1',SUBNA)

      status = nf90_put_att(ncid1,VarID1,
     . 'long_name','time since initialization')
      status = nf90_put_att(ncid1,VarID1,
     . 'units','seconds since 2022-01-01 00:00:00')
      status = nf90_put_att(ncid1,VarID1,
     . 'calendar','gregorian')
      status = nf90_put_att(ncid1,VarID1,
     . 'field','time, scalar, series')

      status = nf90_def_var(ncid = ncid1, name = "SSC_time",
     .xtype = nf90_double,
     .dimids = (/DimID10/),varID = VarID2)
      CALL nccheck_status(status,'VAR2',SUBNA)

      status = nf90_put_att(ncid1,VarID2,
     . 'long_name','time')
      status = nf90_put_att(ncid1,VarID2,
     . 'units','seconds')
      status = nf90_put_att(ncid1,VarID2,
     . 'field','time, scalar, series')

      status = nf90_def_var(ncid = ncid1, name = "SSC",
     .xtype = nf90_double,
     .dimids = (/ DimID1, DimID4, DimID7, DimID10/),varID = VarID3)
      CALL nccheck_status(status,'VAR3',SUBNA)

      status = nf90_enddef(ncid1)
      CALL nccheck_status(status,'ENDDEF',SUBNA)

      k=0
      DO DAY= 1+dzero, DAYS+dzero
      k=k+1
!-----------------------------------------------------------------------
! name of this subroutine!-----------------------------------------------------------------------
!
      SUBNA= 'EXTRACT'
      WRITE(str,"(I5.5)") DAY
      HISNA = 
     .'../WFS_'//TRIM(ADJUSTL(str2))//'_avg_'//TRIM(ADJUSTL(str))//'.nc'
      PRINT*, HISNA
!-----------------------------------------------------------------------
! Get ROMS input
!-----------------------------------------------------------------------
!
! ROMS grid file
      status=nf90_open(TRIM(ADJUSTL(HISNA)),nf90_nowrite,ncID)
      status=nf90_inq_varid(ncID,'ocean_time',varID)
      status=nf90_get_var(ncID,varID,ocean_time,
     .                    start=(/1/),
     .                    count=(/NT/),
     .                    stride=(/REC_INT/))
      CALL nccheck_status(status,'ocean_time',SUBNA)
      
      status=nf90_inq_varid(ncID,'mud_01',varID)      ! mud_01
      status=nf90_get_var(ncID, varID, mud_01,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX,NY,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'mud_01',SUBNA)

      status=nf90_inq_varid(ncID,'mud_02',varID)      ! mud_02
      status=nf90_get_var(ncID, varID, mud_02,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX,NY,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'mud_02',SUBNA)

      status = nf90_close(ncID)
      CALL nccheck_status(status,'END READING',SUBNA)

!-----------------------------------------------------------------------
! Calculate SSC
!-----------------------------------------------------------------------
      SSC = mud_01+mud_02

!-----------------------------------------------------------------------
! Ouput ROMS variable into single file
!----------------------------------------------------------------------- 

      IF(k==1) THEN
         ocean_time_ref = ocean_time(1)
      ENDIF
      DO I = 1,NT
         ocean_time(I) = ocean_time(I)-ocean_time_ref
      ENDDO

      status = nf90_put_var(ncid1,VarId1,ocean_time,
     .start = (/NT*(k-1)+1/),
     .count = (/NT/))
      CALL nccheck_status(status,'OUTPUT1',SUBNA)

      status = nf90_put_var(ncid1,VarId2,ocean_time,
     .start = (/NT*(k-1)+1/),
     .count = (/NT/))
      CALL nccheck_status(status,'OUTPUT2',SUBNA)

      status = nf90_put_var(ncid1,VarId3,SSC,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT3',SUBNA)

      ENDDO

      status = nf90_close(ncid1)
      CALL nccheck_status(status,'END',SUBNA)

      END
