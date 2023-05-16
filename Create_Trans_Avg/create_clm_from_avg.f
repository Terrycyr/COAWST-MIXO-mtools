      PROGRAM create_trans
      USE netcdf
!
!
!  Imported and exported variable declarations. 
!
      INTEGER ::   ncID, varID, DAY, YEAR, DAYS, REC_INT
      INTEGER ::   ncid1(20),I,J,K
      INTEGER ::   Dimid1, Dimid2,Dimid3,Dimid4,Dimid5
      INTEGER ::   Dimid6, Dimid7,Dimid8,Dimid9,Dimid10
      INTEGER ::   Dimid11, Dimid12,Dimid13,Dimid14,Dimid15
      INTEGER ::   Dimid16
      INTEGER ::   Varid1, Varid2, Varid3, Varid4, Varid5
      INTEGER ::   Varid6, Varid7, Varid8, Varid9, Varid10
      INTEGER ::   Varid11, Varid12, Varid13, Varid14, Varid15
      INTEGER ::   Varid16, Varid17, Varid18, Varid19, Varid20
      INTEGER, PARAMETER :: NX=456,NY=184,NZ=21,NT=24
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
      REAL*8     ::  zeta(NX,NY,NT)
      REAL*8     ::  salt(NX,NY,NZ,NT)
      REAL*8     ::  temp(NX,NY,NZ,NT)
      REAL*8     ::  u(NX-1,NY,NZ,NT)
      REAL*8     ::  v(NX,NY-1,NZ,NT)
      REAL*8     ::  Huon(NX-1,NY,NZ,NT)
      REAL*8     ::  Hvom(NX,NY-1,NZ,NT)
      REAL*8     ::  ubar(NX-1,NY,NT)
      REAL*8     ::  vbar(NX,NY-1,NT)
      REAL*8     ::  w(NX,NY,NZ+1,NT)
      REAL*8     ::  omega(NX,NY,NZ+1,NT)
      REAL*8     ::  AKs(NX,NY,NZ+1,NT)
      REAL*8     ::  swrad(NX,NY,NT)
      REAL*8     ::  shflux(NX,NY,NT)
    
      REAL*8     :: dzero

      dzero = 88.
      YEAR = 2021
      DAYS = 365-dzero
      REC_INT = 2
!-----------------------------------------------------------------------
! Define output files
!-----------------------------------------------------------------------
      SUBNA = 'CREATE'

!     Define Dimensions
      WRITE(str2,"(I4.4)") YEAR

      DO I=1,11

      WRITE(str3,"(I2.2)") I

      filename =  './WFS_'//TRIM(ADJUSTL(str2))//'_TRANS_'
     .//TRIM(ADJUSTL(str3))//'.nc'

      status = nf90_create(trim(adjustl(filename))
     . , nf90_clobber,ncid1(i))
      CALL nccheck_status(status,'CREATE',SUBNA)

      status = nf90_def_dim(ncid = ncid1(i), name = "xi_rho",len = NX,
     .dimid = DimID1)
      CALL nccheck_status(status,'DIM1',SUBNA)
      status = nf90_def_dim(ncid = ncid1(i), name = "xi_u",len = NX-1,
     .dimid = DimID2)
      CALL nccheck_status(status,'DIM2',SUBNA)
      status = nf90_def_dim(ncid = ncid1(i), name = "xi_v",len = NX,
     .dimid = DimID3)
      CALL nccheck_status(status,'DIM3',SUBNA)

      status = nf90_def_dim(ncid = ncid1(i), name = "eta_rho",len = NY,
     .dimid = DimID4)
      CALL nccheck_status(status,'DIM4',SUBNA)
      status = nf90_def_dim(ncid = ncid1(i), name = "eta_u",len = NY,
     .dimid = DimID5)
      CALL nccheck_status(status,'DIM5',SUBNA)
      status = nf90_def_dim(ncid = ncid1(i), name = "eta_v",len = NY-1,
     .dimid = DimID6)
      CALL nccheck_status(status,'DIM6',SUBNA)

      status = nf90_def_dim(ncid = ncid1(i), name = "s_rho",len = NZ,
     .dimid = DimID7)
      CALL nccheck_status(status,'DIM7',SUBNA)
      status = nf90_def_dim(ncid = ncid1(i), name = "s_w",len = NZ+1,
     .dimid = DimID8)
      CALL nccheck_status(status,'DIM8',SUBNA)

      status = nf90_def_dim(ncid = ncid1(i), name = "ocean_time",
     .len = NT*DAYS,dimid = DimID9)
      CALL nccheck_status(status,'DIM9',SUBNA)

      status = nf90_def_dim(ncid = ncid1(i), name = "salt_time",
     .len = NT*DAYS,dimid = DimID10)
      CALL nccheck_status(status,'DIM10',SUBNA)

      status = nf90_def_dim(ncid = ncid1(i), name = "temp_time",
     .len = NT*DAYS,dimid = DimID11)
      CALL nccheck_status(status,'DIM11',SUBNA)

!     Define variables
      status = nf90_def_var(ncid = ncid1(i), name = "ocean_time",
     .xtype = nf90_double,
     .dimids = (/DimID9/),varID = VarID1)
      CALL nccheck_status(status,'VAR1',SUBNA)

      status = nf90_put_att(ncid1(i),VarID1,
     . 'long_name','time since initialization')
      status = nf90_put_att(ncid1(i),VarID1,
     . 'units','seconds since 2021-03-30 00:00:00')
      status = nf90_put_att(ncid1(i),VarID1,
     . 'calendar','gregorian')
      status = nf90_put_att(ncid1(i),VarID1,
     . 'field','time, scalar, series')

      status = nf90_def_var(ncid = ncid1(i), name = "temp_time",
     .xtype = nf90_double,
     .dimids = (/DimID11/),varID = VarID11)
      CALL nccheck_status(status,'VAR11',SUBNA)

      status = nf90_put_att(ncid1(i),VarID11,
     . 'long_name','time')
      status = nf90_put_att(ncid1(i),VarID11,
     . 'units','seconds')
      status = nf90_put_att(ncid1(i),VarID11,
     . 'field','time, scalar, series')

      status = nf90_def_var(ncid = ncid1(i), name = "salt_time",
     .xtype = nf90_double,
     .dimids = (/DimID10/),varID = VarID12)
      CALL nccheck_status(status,'VAR12',SUBNA)

      status = nf90_put_att(ncid1(i),VarID12,
     . 'long_name','time')
      status = nf90_put_att(ncid1(i),VarID12,
     . 'units','seconds')
      status = nf90_put_att(ncid1(i),VarID12,
     . 'field','time, scalar, series')


      SELECT CASE (I)
      CASE(1)
      status = nf90_def_var(ncid = ncid1(i), name = "zeta",
     .xtype = nf90_double,
     .dimids = (/ DimID1, DimID4, DimID9 /),varID = VarID2)
      CALL nccheck_status(status,'VAR2',SUBNA)

      CASE(2)
      status = nf90_def_var(ncid = ncid1(i), name = "ubar",
     .xtype = nf90_double,
     .dimids = (/ DimID2, DimID5, DimID9/),varID = VarID3)
      CALL nccheck_status(status,'VAR3',SUBNA)

      CASE(3)
      status = nf90_def_var(ncid = ncid1(i), name = "vbar",
     .xtype = nf90_double,
     .dimids = (/ DimID3, DimID6, DimID9/),varID = VarID4)
      CALL nccheck_status(status,'VAR4',SUBNA)

      CASE(4)
      status = nf90_def_var(ncid = ncid1(i), name = "u",
     .xtype = nf90_double,
     .dimids = (/ DimID2, DimID5, DimID7, DimID9/),varID = VarID5)
      CALL nccheck_status(status,'VAR5',SUBNA)

      CASE(5)
      status = nf90_def_var(ncid = ncid1(i), name = "v",
     .xtype = nf90_double,
     .dimids = (/ DimID3, DimID6, DimID7, DimID9/),varID = VarID6)
      CALL nccheck_status(status,'VAR6',SUBNA)

      CASE(6)
      status = nf90_def_var(ncid = ncid1(i), name = "omega",
     .xtype = nf90_double,
     .dimids = (/ DimID1, DimID4, DimID8, DimID9/),varID = VarID7)
      CALL nccheck_status(status,'VAR7',SUBNA)

      CASE(7)
      status = nf90_def_var(ncid = ncid1(i), name = "temp",
     .xtype = nf90_double,
     .dimids = (/ DimID1, DimID4, DimID7, DimID11/),varID = VarID8)
      CALL nccheck_status(status,'VAR8',SUBNA)

      CASE(8)
      status = nf90_def_var(ncid = ncid1(i), name = "salt",
     .xtype = nf90_double,
     .dimids = (/ DimID1, DimID4, DimID7, DimID10/),varID = VarID9)
      CALL nccheck_status(status,'VAR9',SUBNA)

      CASE(9)
      status = nf90_def_var(ncid = ncid1(i), name = "AKs",
     .xtype = nf90_double,
     .dimids = (/ DimID1, DimID4, DimID8, DimID9/),varID = VarID10)
      CALL nccheck_status(status,'VAR10',SUBNA)

      CASE(10)
      status = nf90_def_var(ncid = ncid1(i), name = "Huon",
     .xtype = nf90_double,
     .dimids = (/ DimID2, DimID5, DimID7, DimID9/),varID = VarID13)
      CALL nccheck_status(status,'VAR13',SUBNA)

      CASE(11)
      status = nf90_def_var(ncid = ncid1(i), name = "Hvom",
     .xtype = nf90_double,
     .dimids = (/ DimID3, DimID6, DimID7, DimID9/),varID = VarID14)
      CALL nccheck_status(status,'VAR14',SUBNA)

      END SELECT
      status = nf90_enddef(ncid1(i))
      CALL nccheck_status(status,'ENDDEF',SUBNA)

      ENDDO

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
      
      status=nf90_inq_varid(ncID,'zeta',varID)     ! free surface
      status=nf90_get_var(ncID, varID, zeta,
     .                    start=(/1,1,1/),
     .                    count=(/NX,NY,NT/),
     .                    stride=(/1,1,REC_INT/))
      CALL nccheck_status(status,'zeta',SUBNA)
      
      status=nf90_inq_varid(ncID,'salt',varID)     ! salinity
      status=nf90_get_var(ncID, varID, salt,
     .                    start=(/1,1,1,1/),
     .                    count=(/NX,NY,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'salt',SUBNA)
      
      status=nf90_inq_varid(ncID,'temp',varID)     ! temperature
      status=nf90_get_var(ncID,varID, temp,
     .                    start=(/1,1,1,1/),
     .                    count=(/NX,NY,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'temp',SUBNA)
      
      status=nf90_inq_varid(ncID,'ubar',varID)        ! ubar-momentum
      status=nf90_get_var(ncID,varID,ubar,
     .                    start = (/1,1,1/),
     .                    count = (/NX-1,NY,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'ubar',SUBNA)

      status=nf90_inq_varid(ncID,'vbar',varID)        ! vbar-momentum
      status=nf90_get_var(ncID,varID,vbar,
     .                    start = (/1,1,1/),
     .                    count = (/NX,NY-1,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'vbar',SUBNA)

      status=nf90_inq_varid(ncID,'u',varID)        ! u-momentum
      status=nf90_get_var(ncID,varID,u, 
     .                    start = (/1,1,1,1/),
     .                    count = (/NX-1,NY,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'u',SUBNA)
      
      status=nf90_inq_varid(ncID,'v',varID)        ! v-momentum
      status=nf90_get_var(ncID,varID,v,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX,NY-1,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'v',SUBNA)

      status=nf90_inq_varid(ncID,'Huon',varID)        ! u-mass transport
      status=nf90_get_var(ncID,varID,Huon,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX-1,NY,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'Huon',SUBNA)

      status=nf90_inq_varid(ncID,'Hvom',varID)        ! v-mass transport
      status=nf90_get_var(ncID,varID,Hvom,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX,NY-1,NZ,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'Hvom',SUBNA)

!      status=nf90_inq_varid(ncID,'w',varID)    ! w-momentum
!      status=nf90_get_var(ncID, varID, w,
!     .                    start = (/1,1,1,1/),
!     .                    count = (/NX,NY,NZ+1,NT/),
!     .                    stride=(/1,1,1,REC_INT/))
!      CALL nccheck_status(status,'w',SUBNA)


      status=nf90_inq_varid(ncID,'omega',varID)    ! omega-momentum
      status=nf90_get_var(ncID, varID, omega,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX,NY,NZ+1,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'omega',SUBNA)
      
      status=nf90_inq_varid(ncID,'AKs',varID)      ! diffusivity
      status=nf90_get_var(ncID, varID, AKs,
     .                    start = (/1,1,1,1/),
     .                    count = (/NX,NY,NZ+1,NT/),
     .                    stride=(/1,1,1,REC_INT/))
      CALL nccheck_status(status,'AKs',SUBNA)

!      status=nf90_inq_varid(ncID,'swrad',varID)      ! short wave radiation
!      status=nf90_get_var(ncID, varID, swrad,
!     .                    start = (/1,1,1/),
!     .                    count = (/NX,NY,NT/),
!     .                    stride=(/1,1,REC_INT/))
!      CALL nccheck_status(status,'swrad',SUBNA)
!
!      status=nf90_inq_varid(ncID,'shflux',varID)      ! net heat flux
!      status=nf90_get_var(ncID, varID, shflux,
!     .                    start = (/1,1,1/),
!     .                    count = (/NX,NY,NT/),
!     .                    stride=(/1,1,REC_INT/))
!      CALL nccheck_status(status,'shflux',SUBNA)

!-----------------------------------------------------------------------
! Ouput ROMS variable into single file
!----------------------------------------------------------------------- 

      IF(k==1) THEN
         ocean_time_ref = ocean_time(1)
      ENDIF
      DO I = 1,NT
         ocean_time(I) = ocean_time(I)-ocean_time_ref
      ENDDO

      DO I=1,11
      status = nf90_put_var(ncid1(i),VarId1,ocean_time,
     .start = (/NT*(k-1)+1/),
     .count = (/NT/))
      CALL nccheck_status(status,'OUTPUT1',SUBNA)

      status = nf90_put_var(ncid1(i),VarId11,ocean_time,
     .start = (/NT*(k-1)+1/),
     .count = (/NT/))
      CALL nccheck_status(status,'OUTPUT11',SUBNA)

      status = nf90_put_var(ncid1(i),VarId12,ocean_time,
     .start = (/NT*(k-1)+1/),
     .count = (/NT/))
      CALL nccheck_status(status,'OUTPUT12',SUBNA)
      ENDDO

      status = nf90_put_var(ncid1(1),VarId2,zeta,
     .start = (/1,1,NT*(k-1)+1/),
     .count = (/NX,NY,NT/))
      CALL nccheck_status(status,'OUTPUT2',SUBNA)

      status = nf90_put_var(ncid1(2),VarId3,ubar,
     .start = (/1,1,NT*(k-1)+1/),
     .count = (/NX-1,NY,NT/))
      CALL nccheck_status(status,'OUTPUT3',SUBNA)

      status = nf90_put_var(ncid1(3),VarId4,vbar,
     .start = (/1,1,NT*(k-1)+1/),
     .count = (/NX,NY-1,NT/))
      CALL nccheck_status(status,'OUTPUT4',SUBNA)

      status = nf90_put_var(ncid1(4),VarId5,u,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX-1,NY,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT5',SUBNA)

      status = nf90_put_var(ncid1(5),VarId6,v,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY-1,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT6',SUBNA)

      status = nf90_put_var(ncid1(6),VarId7,omega,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY,NZ+1,NT/))
      CALL nccheck_status(status,'OUTPUT7',SUBNA)
      
      status = nf90_put_var(ncid1(7),VarId8,temp,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT8',SUBNA)

      status = nf90_put_var(ncid1(8),VarId9,salt,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT9',SUBNA)

      status = nf90_put_var(ncid1(9),VarId10,AKs,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY,NZ+1,NT/))
      CALL nccheck_status(status,'OUTPUT10',SUBNA)

      status = nf90_put_var(ncid1(10),VarId13,Huon,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX-1,NY,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT13',SUBNA)

      status = nf90_put_var(ncid1(11),VarId14,Hvom,
     .start = (/1,1,1,NT*(k-1)+1/),
     .count = (/NX,NY-1,NZ,NT/))
      CALL nccheck_status(status,'OUTPUT14',SUBNA)

      ENDDO

      DO I=1,11
        status = nf90_close(ncid1(I))
        CALL nccheck_status(status,'END',SUBNA)
      ENDDO

      END
