function create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg,NBT,NDYE,bio_sed, dataset_source)

%create init file
nc_init=netcdf.create(init_file,'clobber');
 
%% Global attributes:
myName= mfilename;
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by ',myName,' on ' datestr(now)]);
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'type', dataset_source);


%% Dimensions:

disp(' ## Defining Dimensions...')
 
%get some grid info
  [LP,MP]=size(gn.lon_rho);
  L=LP-1;
  Lm=L-1;
  M=MP-1;
  Mm=M-1;
  L  = Lm+1;
  M  = Mm+1;
  xpsi  = L;
  xrho  = LP;
  xu    = L;
  xv    = LP;
  epsi = M;
  erho = MP;
  eu   = MP;
  ev   = M;
  N       = gn.N;
  
psidimID = netcdf.defDim(nc_init,'xpsi',L);
xrhodimID = netcdf.defDim(nc_init,'xrho',LP);
xudimID = netcdf.defDim(nc_init,'xu',L);
xvdimID = netcdf.defDim(nc_init,'xv',LP);

epsidimID = netcdf.defDim(nc_init,'epsi',M);
erhodimID = netcdf.defDim(nc_init,'erho',MP);
eudimID = netcdf.defDim(nc_init,'eu',MP);
evdimID = netcdf.defDim(nc_init,'ev',M);

s_rhodimID = netcdf.defDim(nc_init,'s_rho',N);
s_wdimID = netcdf.defDim(nc_init,'s_w',N+1);
if(Nbed>0)
    NbeddimID = netcdf.defDim(nc_init,'Nbed',Nbed);
end
timedimID = netcdf.defDim(nc_init,'ocean_time',1);
if(Nveg>0)
    NvegdimID = netcdf.defDim(nc_init,'Nveg',Nveg);
end

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

sphericalID = netcdf.defVar(nc_init,'spherical','short',timedimID);
netcdf.putAtt(nc_init,sphericalID,'long_name','grid type logical switch');
netcdf.putAtt(nc_init,sphericalID,'flag_meanings','spherical Cartesian');
netcdf.putAtt(nc_init,sphericalID,'flag_values','1, 0');

VtransformID = netcdf.defVar(nc_init,'Vtransform','long',timedimID);
netcdf.putAtt(nc_init,VtransformID,'long_name','vertical terrain-following transformation equation');

VstretchingID = netcdf.defVar(nc_init,'Vstretching','long',timedimID);
netcdf.putAtt(nc_init,VstretchingID,'long_name','vertical terrain-following stretching function');
 
theta_bID = netcdf.defVar(nc_init,'theta_b','double',timedimID);
netcdf.putAtt(nc_init,theta_bID,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(nc_init,theta_bID,'units','1');

theta_sID = netcdf.defVar(nc_init,'theta_s','double',timedimID);
netcdf.putAtt(nc_init,theta_sID,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(nc_init,theta_sID,'units','1');

tcline_ID = netcdf.defVar(nc_init,'Tcline','double',timedimID);
netcdf.putAtt(nc_init,tcline_ID,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(nc_init,tcline_ID,'units','meter');

hc_ID = netcdf.defVar(nc_init,'hc','double',timedimID);
netcdf.putAtt(nc_init,hc_ID,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(nc_init,hc_ID,'units','meter');

Cs_rID = netcdf.defVar(nc_init,'Cs_r','double',s_rhodimID);
netcdf.putAtt(nc_init,Cs_rID,'long_name','S-coordinate stretching curves at RHO-points');
netcdf.putAtt(nc_init,Cs_rID,'units','1');
netcdf.putAtt(nc_init,Cs_rID,'valid_min',-1);
netcdf.putAtt(nc_init,Cs_rID,'valid_max',0);
netcdf.putAtt(nc_init,Cs_rID,'field','Cs_r, scalar');

Cs_wID = netcdf.defVar(nc_init,'Cs_w','double',s_wdimID);
netcdf.putAtt(nc_init,Cs_wID,'long_name','S-coordinate stretching curves at W-points');
netcdf.putAtt(nc_init,Cs_wID,'units','1');
netcdf.putAtt(nc_init,Cs_wID,'valid_min',-1);
netcdf.putAtt(nc_init,Cs_wID,'valid_max',0);
netcdf.putAtt(nc_init,Cs_wID,'field','Cs_w, scalar');

s_rhoID = netcdf.defVar(nc_init,'s_rho','double',s_rhodimID);
netcdf.putAtt(nc_init,s_rhoID,'long_name','S-coordinate at RHO-points');
netcdf.putAtt(nc_init,s_rhoID,'units','1');
netcdf.putAtt(nc_init,s_rhoID,'valid_min',-1);
netcdf.putAtt(nc_init,s_rhoID,'valid_max',0);
netcdf.putAtt(nc_init,s_rhoID,'field','s_rho, scalar');

s_wID = netcdf.defVar(nc_init,'s_w','double',s_wdimID);
netcdf.putAtt(nc_init,s_wID,'long_name','S-coordinate at W-points');
netcdf.putAtt(nc_init,s_wID,'units','1');
netcdf.putAtt(nc_init,s_wID,'valid_min',-1);
netcdf.putAtt(nc_init,s_wID,'valid_max',0);
netcdf.putAtt(nc_init,s_wID,'field','s_w, scalar');

ocean_timeID = netcdf.defVar(nc_init,'ocean_time','double',timedimID);
netcdf.putAtt(nc_init,ocean_timeID,'long_name','time since initialization');
netcdf.putAtt(nc_init,ocean_timeID,'units','days');
netcdf.putAtt(nc_init,ocean_timeID,'field','ocean_time, scalar, series');

saltID = netcdf.defVar(nc_init,'salt','float',[xrhodimID erhodimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,saltID,'long_name','salinity');
netcdf.putAtt(nc_init,saltID,'units','PSU');
netcdf.putAtt(nc_init,saltID,'field','salinity, scalar, series');

tempID = netcdf.defVar(nc_init,'temp','float',[xrhodimID erhodimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,tempID,'long_name','temperature');
netcdf.putAtt(nc_init,tempID,'units','C');
netcdf.putAtt(nc_init,tempID,'field','temperature, scalar, series');

if(NDYE>0)
    dyeID = netcdf.defVar(nc_init,'dye_01','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,dyeID,'long_name','tracer 01');
    netcdf.putAtt(nc_init,dyeID,'units','mg/L');
    netcdf.putAtt(nc_init,dyeID,'field','tracer 01, scalar, series');
end

uID = netcdf.defVar(nc_init,'u','float',[xudimID eudimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,uID,'long_name','u-momentum component');
netcdf.putAtt(nc_init,uID,'units','meter second-1');
netcdf.putAtt(nc_init,uID,'field','u-velocity, scalar, series');

ubarID = netcdf.defVar(nc_init,'ubar','float',[xudimID eudimID timedimID]);
netcdf.putAtt(nc_init,ubarID,'long_name','vertically integrated u-momentum component');
netcdf.putAtt(nc_init,ubarID,'units','meter second-1');
netcdf.putAtt(nc_init,ubarID,'field','ubar-velocity, scalar, series');

vID = netcdf.defVar(nc_init,'v','float',[xvdimID evdimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,vID,'long_name','v-momentum component');
netcdf.putAtt(nc_init,vID,'units','meter second-1');
netcdf.putAtt(nc_init,vID,'field','v-velocity, scalar, series');

vbarID = netcdf.defVar(nc_init,'vbar','float',[xvdimID evdimID timedimID]);
netcdf.putAtt(nc_init,vbarID,'long_name','vertically integrated v-momentum component');
netcdf.putAtt(nc_init,vbarID,'units','meter second-1');
netcdf.putAtt(nc_init,vbarID,'field','vbar-velocity, scalar, series');
 
zetaID = netcdf.defVar(nc_init,'zeta','float',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,zetaID,'long_name','free-surface');
netcdf.putAtt(nc_init,zetaID,'units','meter');
netcdf.putAtt(nc_init,zetaID,'field','free-surface, scalar, series');
 
for mm=1:NCS
    count=['00',num2str(mm)];
    count=count(end-1:end);

    eval(['mud_',count,'ID = netcdf.defVar(nc_init,''mud_',count,''',''double'',[xrhodimID erhodimID s_rhodimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''long_name'',''suspended cohesive sediment, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''units'',''kilogram meter-3'');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''field'',''mud_',count,', scalar, series'');'])

    eval(['mudfrac_',count,'ID = netcdf.defVar(nc_init,''mudfrac_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''long_name'',''cohesive sediment fraction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''units'',''nondimensional'');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''field'',''mudfrac_',count,', scalar, series'');'])

     eval(['mudmass_',count,'ID = netcdf.defVar(nc_init,''mudmass_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''long_name'',''cohesive sediment mass, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''units'',''kilogram meter-2'');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''field'',''mudmass_',count,', scalar, series'');'])

end
for mm=1:NNS
    count=['00',num2str(mm)];
    count=count(end-1:end);

    eval(['sand_',count,'ID = netcdf.defVar(nc_init,''sand_',count,''',''double'',[xrhodimID erhodimID s_rhodimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''long_name'',''suspended noncohesive sediment, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''units'',''kilogram meter-3'');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''field'',''sand_',count,', scalar, series'');'])

    eval(['sandfrac_',count,'ID = netcdf.defVar(nc_init,''sandfrac_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''long_name'',''noncohesive sediment fraction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''units'',''nondimensional'');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''field'',''sandfrac_',count,', scalar, series'');'])

    eval(['sandmass_',count,'ID = netcdf.defVar(nc_init,''sandmass_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''long_name'',''noncohesive sediment mass, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''units'',''kilogram meter-2'');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''field'',''sandmass_',count,', scalar, series'');'])

    eval(['bedload_Usand_',count,'ID = netcdf.defVar(nc_init,''bedload_Usand_',count,''',''double'',[xudimID eudimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''long_name'',''bed load flux of sand in U-direction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''units'',''kilogram meter-1 s-1'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''field'',''bedload_Usand_',count,', scalar, series'');'])

    eval(['bedload_Vsand_',count,'ID = netcdf.defVar(nc_init,''bedload_Vsand_',count,''',''double'',[xvdimID evdimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''long_name'',''bed load flux of sand in V-direction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''units'',''kilogram meter-1 s-1'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''field'',''bedload_Vsand_',count,', scalar, series'');'])
end

if(Nbed>0)
    bed_thicknessID = netcdf.defVar(nc_init,'bed_thickness','double',[xrhodimID erhodimID NbeddimID timedimID]);
    netcdf.putAtt(nc_init,bed_thicknessID,'long_name','sediment layer thickness');
    netcdf.putAtt(nc_init,bed_thicknessID,'units','meter');
    netcdf.putAtt(nc_init,bed_thicknessID,'time','ocean_time');
    netcdf.putAtt(nc_init,bed_thicknessID,'field','bed thickness, scalar, series');

    bed_ageID = netcdf.defVar(nc_init,'bed_age','double',[xrhodimID erhodimID NbeddimID timedimID]);
    netcdf.putAtt(nc_init,bed_ageID,'long_name','sediment layer age');
    netcdf.putAtt(nc_init,bed_ageID,'units','day');
    netcdf.putAtt(nc_init,bed_ageID,'time','ocean_time');
    netcdf.putAtt(nc_init,bed_ageID,'field','bed age, scalar, series');

    bed_porosityID = netcdf.defVar(nc_init,'bed_porosity','double',[xrhodimID erhodimID NbeddimID timedimID]);
    netcdf.putAtt(nc_init,bed_porosityID,'long_name','sediment layer porosity');
    netcdf.putAtt(nc_init,bed_porosityID,'units','nondimensional');
    netcdf.putAtt(nc_init,bed_porosityID,'time','ocean_time');
    netcdf.putAtt(nc_init,bed_porosityID,'field','bed porosity, scalar, series');

    bed_biodiffID = netcdf.defVar(nc_init,'bed_biodiff','double',[xrhodimID erhodimID NbeddimID timedimID]);
    netcdf.putAtt(nc_init,bed_biodiffID,'long_name','biodiffusivity at bottom of each layer');
    netcdf.putAtt(nc_init,bed_biodiffID,'units','meter2 second-1');
    netcdf.putAtt(nc_init,bed_biodiffID,'time','ocean_time');
    netcdf.putAtt(nc_init,bed_biodiffID,'field','bed biodiffusivity, scalar, series');

    bed_taucritID = netcdf.defVar(nc_init,'bed_tau_crit','double',[xrhodimID erhodimID NbeddimID timedimID]);
    netcdf.putAtt(nc_init,bed_taucritID,'long_name','tau critical in each layer');
    netcdf.putAtt(nc_init,bed_taucritID,'units','Newtons meter-2');
    netcdf.putAtt(nc_init,bed_taucritID,'time','ocean_time');
    netcdf.putAtt(nc_init,bed_taucritID,'field','bed tau critical, scalar, series');
end

grain_diameterID = netcdf.defVar(nc_init,'grain_diameter','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,grain_diameterID,'long_name','sediment median grain diameter size');
netcdf.putAtt(nc_init,grain_diameterID,'units','meter');
netcdf.putAtt(nc_init,grain_diameterID,'time','ocean_time');
netcdf.putAtt(nc_init,grain_diameterID,'field','grain diameter, scalar, series');

grain_densityID = netcdf.defVar(nc_init,'grain_density','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,grain_densityID,'long_name','sediment median grain density');
netcdf.putAtt(nc_init,grain_densityID,'units','kilogram meter-3');
netcdf.putAtt(nc_init,grain_densityID,'time','ocean_time');
netcdf.putAtt(nc_init,grain_densityID,'field','grain density, scalar, series');

settling_velID = netcdf.defVar(nc_init,'settling_vel','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,settling_velID,'long_name','sediment median grain settling velocity');
netcdf.putAtt(nc_init,settling_velID,'units','meter second-1');
netcdf.putAtt(nc_init,settling_velID,'time','ocean_time');
netcdf.putAtt(nc_init,settling_velID,'field','settling vel, scalar, series');

erosion_stressID = netcdf.defVar(nc_init,'erosion_stress','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,erosion_stressID,'long_name','sediment median critical erosion stress');
netcdf.putAtt(nc_init,erosion_stressID,'units','meter2 second-2');
netcdf.putAtt(nc_init,erosion_stressID,'time','ocean_time');
netcdf.putAtt(nc_init,erosion_stressID,'field','erosion stress, scalar, series');

if(Nbed>0)
    ripple_lengthID = netcdf.defVar(nc_init,'ripple_length','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,ripple_lengthID,'long_name','bottom ripple length');
    netcdf.putAtt(nc_init,ripple_lengthID,'units','meter');
    netcdf.putAtt(nc_init,ripple_lengthID,'time','ocean_time');
    netcdf.putAtt(nc_init,ripple_lengthID,'field','ripple length, scalar, series');

    ripple_heightID = netcdf.defVar(nc_init,'ripple_height','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,ripple_heightID,'long_name','bottom ripple height');
    netcdf.putAtt(nc_init,ripple_heightID,'units','meter');
    netcdf.putAtt(nc_init,ripple_heightID,'time','ocean_time');
    netcdf.putAtt(nc_init,ripple_heightID,'field','ripple height, scalar, series');

    dmix_offsetID = netcdf.defVar(nc_init,'dmix_offset','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,dmix_offsetID,'long_name','dmix erodibility profile offset');
    netcdf.putAtt(nc_init,dmix_offsetID,'units','meter');
    netcdf.putAtt(nc_init,dmix_offsetID,'time','ocean_time');
    netcdf.putAtt(nc_init,dmix_offsetID,'field','dmix_offset, scalar, series');

    dmix_slopeID = netcdf.defVar(nc_init,'dmix_slope','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,dmix_slopeID,'long_name','dmix erodibility profile slope');
    netcdf.putAtt(nc_init,dmix_slopeID,'units','_');
    netcdf.putAtt(nc_init,dmix_slopeID,'time','ocean_time');
    netcdf.putAtt(nc_init,dmix_slopeID,'field','dmix_slope, scalar, series');

    dmix_timeID = netcdf.defVar(nc_init,'dmix_time','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,dmix_timeID,'long_name','dmix erodibility profile time scale');
    netcdf.putAtt(nc_init,dmix_timeID,'units','seconds');
    netcdf.putAtt(nc_init,dmix_timeID,'time','ocean_time');
    netcdf.putAtt(nc_init,dmix_timeID,'field','dmix_time, scalar, series');
end

if(Nveg>0)
    vegID = netcdf.defVar(nc_init,'plant_height','double',[xrhodimID erhodimID NvegdimID timedimID]);
    netcdf.putAtt(nc_init,vegID,'long_name','plant height');
    netcdf.putAtt(nc_init,vegID,'units','meter');
    netcdf.putAtt(nc_init,vegID,'time','ocean_time');
    netcdf.putAtt(nc_init,vegID,'field','plant_height, scalar, series');

    vegID = netcdf.defVar(nc_init,'plant_density','double',[xrhodimID erhodimID NvegdimID timedimID]);
    netcdf.putAtt(nc_init,vegID,'long_name','plant density');
    netcdf.putAtt(nc_init,vegID,'units','plant-meter2');
    netcdf.putAtt(nc_init,vegID,'time','ocean_time');
    netcdf.putAtt(nc_init,vegID,'field','plant_density, scalar, series');

    vegID = netcdf.defVar(nc_init,'plant_diameter','double',[xrhodimID erhodimID NvegdimID timedimID]);
    netcdf.putAtt(nc_init,vegID,'long_name','plant diameter');
    netcdf.putAtt(nc_init,vegID,'units','meter');
    netcdf.putAtt(nc_init,vegID,'time','ocean_time');
    netcdf.putAtt(nc_init,vegID,'field','plant_diameter, scalar, series');

    vegID = netcdf.defVar(nc_init,'plant_thickness','double',[xrhodimID erhodimID NvegdimID timedimID]);
    netcdf.putAtt(nc_init,vegID,'long_name','plant thickness');
    netcdf.putAtt(nc_init,vegID,'units','meter');
    netcdf.putAtt(nc_init,vegID,'time','ocean_time');
    netcdf.putAtt(nc_init,vegID,'field','plant_thickness, scalar, series');

    vegID = netcdf.defVar(nc_init,'marsh_mask','double',[xrhodimID erhodimID NvegdimID timedimID]);
    netcdf.putAtt(nc_init,vegID,'long_name','marsh mask');
    netcdf.putAtt(nc_init,vegID,'units','nondimensional');
    netcdf.putAtt(nc_init,vegID,'time','ocean_time');
    netcdf.putAtt(nc_init,vegID,'field','marsh_mask, scalar, series');
end

if(NBT>0)
    acID = netcdf.defVar(nc_init,'A_C','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,acID,'long_name','MX model prey');
    netcdf.putAtt(nc_init,acID,'units','mg C L-1');
    netcdf.putAtt(nc_init,acID,'field','A_C, scalar, series');

    achlcID = netcdf.defVar(nc_init,'A_CHLC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,achlcID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_init,achlcID,'units','gChl/gC');
    netcdf.putAtt(nc_init,achlcID,'field','A_CHLC, scalar, series');

    ancID = netcdf.defVar(nc_init,'A_NC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,ancID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_init,ancID,'units','gN/gC');
    netcdf.putAtt(nc_init,ancID,'field','A_NC, scalar, series');

    apcID = netcdf.defVar(nc_init,'A_PC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,apcID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_init,apcID,'units','gP/gC');
    netcdf.putAtt(nc_init,apcID,'field','A_PC, scalar, series');

    bsiID = netcdf.defVar(nc_init,'BSI','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,bsiID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_init,bsiID,'units','mg SI L-1');
    netcdf.putAtt(nc_init,bsiID,'field','BSI, scalar, series');

    doID = netcdf.defVar(nc_init,'DO','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,doID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_init,doID,'units','mg O2 L-1');
    netcdf.putAtt(nc_init,doID,'field','DO, scalar, series');

    exdocID = netcdf.defVar(nc_init,'EXDOC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,exdocID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_init,exdocID,'units','mg C L-1');
    netcdf.putAtt(nc_init,exdocID,'field','exdoc, scalar, series');

    ldocID = netcdf.defVar(nc_init,'LDOC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,ldocID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_init,ldocID,'units','mg C L-1');
    netcdf.putAtt(nc_init,ldocID,'field','ldoc, scalar, series');

    ldonID = netcdf.defVar(nc_init,'LDON','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,ldonID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_init,ldonID,'units','mg N L-1');
    netcdf.putAtt(nc_init,ldonID,'field','ldon, scalar, series');

    ldopID = netcdf.defVar(nc_init,'LDOP','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,ldopID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_init,ldopID,'units','mg P L-1');
    netcdf.putAtt(nc_init,ldopID,'field','ldop, scalar, series');


    lpocID = netcdf.defVar(nc_init,'LPOC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,lpocID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_init,lpocID,'units','mg C L-1');
    netcdf.putAtt(nc_init,lpocID,'field','lpoc, scalar, series');

    lponID = netcdf.defVar(nc_init,'LPON','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,lponID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_init,lponID,'units','mg N L-1');
    netcdf.putAtt(nc_init,lponID,'field','lpon, scalar, series');

    lpopID = netcdf.defVar(nc_init,'LPOP','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,lpopID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_init,lpopID,'units','mg P L-1');
    netcdf.putAtt(nc_init,lpopID,'field','lpop, scalar, series');

    mavguID = netcdf.defVar(nc_init,'M_AVGCU','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mavguID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_init,mavguID,'units','C/C/d');
    netcdf.putAtt(nc_init,mavguID,'field','mavgu, scalar, series');

    mcID = netcdf.defVar(nc_init,'M_C','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mcID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_init,mcID,'units','mg C L-1');
    netcdf.putAtt(nc_init,mcID,'field','mc, scalar, series');

    mchlcID = netcdf.defVar(nc_init,'M_CHLC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mchlcID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_init,mchlcID,'units','gChl/gC');
    netcdf.putAtt(nc_init,mchlcID,'field','mchlc, scalar, series');

    mfcID = netcdf.defVar(nc_init,'M_FC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mfcID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_init,mfcID,'units','gC/gC');
    netcdf.putAtt(nc_init,mfcID,'field','mfc, scalar, series');

    mfchlcID = netcdf.defVar(nc_init,'M_FCHLC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mfchlcID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_init,mfchlcID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_init,mfchlcID,'field','mfchlc, scalar, series');

    mncID = netcdf.defVar(nc_init,'M_NC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mncID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_init,mncID,'units','gN/gC');
    netcdf.putAtt(nc_init,mncID,'field','mnc, scalar, series');

    mpcID = netcdf.defVar(nc_init,'M_PC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mpcID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_init,mpcID,'units','gP/gC');
    netcdf.putAtt(nc_init,mpcID,'field','mpc, scalar, series');

    mfncID = netcdf.defVar(nc_init,'MFNC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mfncID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_init,mfncID,'units','gN/gC');
    netcdf.putAtt(nc_init,mfncID,'field','mfnc, scalar, series');

    mfpcID = netcdf.defVar(nc_init,'MFPC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,mfpcID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_init,mfpcID,'units','gP/gC');
    netcdf.putAtt(nc_init,mfpcID,'field','mfpc, scalar, series');

    nh4tID = netcdf.defVar(nc_init,'NH4T','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,nh4tID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_init,nh4tID,'units','mg N L-1');
    netcdf.putAtt(nc_init,nh4tID,'field','nh4t, scalar, series');

    no23ID = netcdf.defVar(nc_init,'NO23','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,no23ID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_init,no23ID,'units','mg N L-1');
    netcdf.putAtt(nc_init,no23ID,'field','no23, scalar, series')

    o2eqID = netcdf.defVar(nc_init,'O2EQ','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,o2eqID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_init,o2eqID,'units','mg O2 L-1');
    netcdf.putAtt(nc_init,o2eqID,'field','o2eq, scalar, series')

    phyt1ID = netcdf.defVar(nc_init,'PHYT1','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,phyt1ID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_init,phyt1ID,'units','mg C L-1');
    netcdf.putAtt(nc_init,phyt1ID,'field','phyt1, scalar, series')

    phyt2ID = netcdf.defVar(nc_init,'PHYT2','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,phyt2ID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_init,phyt2ID,'units','mg C L-1');
    netcdf.putAtt(nc_init,phyt2ID,'field','phyt2, scalar, series')

    phyt3ID = netcdf.defVar(nc_init,'PHYT3','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,phyt3ID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_init,phyt3ID,'units','mg C L-1');
    netcdf.putAtt(nc_init,phyt3ID,'field','phyt3, scalar, series')

    zoo1ID = netcdf.defVar(nc_init,'ZOO1','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,zoo1ID,'long_name','SMALL ZOOPLANKTON (ZOO1)');
    netcdf.putAtt(nc_init,zoo1ID,'units','mg C L-1');
    netcdf.putAtt(nc_init,zoo1ID,'field','zoo1, scalar, series')

    zoo2ID = netcdf.defVar(nc_init,'ZOO2','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,zoo2ID,'long_name','LARGE ZOOPLANKTON (ZOO2)');
    netcdf.putAtt(nc_init,zoo2ID,'units','mg C L-1');
    netcdf.putAtt(nc_init,zoo2ID,'field','zoo2, scalar, series')

    po4tID = netcdf.defVar(nc_init,'PO4T','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,po4tID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_init,po4tID,'units','mg P L-1');
    netcdf.putAtt(nc_init,po4tID,'field','po4t, scalar, series')

    rdocID = netcdf.defVar(nc_init,'RDOC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,rdocID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_init,rdocID,'units','mg C L-1');
    netcdf.putAtt(nc_init,rdocID,'field','rdoc, scalar, series')

    rdonID = netcdf.defVar(nc_init,'RDON','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,rdonID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_init,rdonID,'units','mg N L-1');
    netcdf.putAtt(nc_init,rdonID,'field','rdon, scalar, series')

    rdopID = netcdf.defVar(nc_init,'RDOP','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,rdopID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_init,rdopID,'units','mg P L-1');
    netcdf.putAtt(nc_init,rdopID,'field','rdop, scalar, series')

    TAID = netcdf.defVar(nc_init,'TA','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,TAID,'long_name','TOTAL ALKALINITY');
    netcdf.putAtt(nc_init,TAID,'units','umol L-1');
    netcdf.putAtt(nc_init,TAID,'field','TA, scalar, series')

    DICID = netcdf.defVar(nc_init,'DIC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,DICID,'long_name','DISSOLVED INORGANIC CARBON');
    netcdf.putAtt(nc_init,DICID,'units','umol L-1');
    netcdf.putAtt(nc_init,DICID,'field','DIC, scalar, series')

    CACO3ID = netcdf.defVar(nc_init,'CACO3','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,CACO3ID,'long_name','CALCIUM CARBONATE');
    netcdf.putAtt(nc_init,CACO3ID,'units','mmol L-1');
    netcdf.putAtt(nc_init,CACO3ID,'field','CACO3, scalar, series')

    rpocID = netcdf.defVar(nc_init,'RPOC','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,rpocID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_init,rpocID,'units','mg C L-1');
    netcdf.putAtt(nc_init,rpocID,'field','rpoc, scalar, series')

    rponID = netcdf.defVar(nc_init,'RPON','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,rponID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_init,rponID,'units','mg N L-1');
    netcdf.putAtt(nc_init,rponID,'field','rpon, scalar, series')

    rpopID = netcdf.defVar(nc_init,'RPOP','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,rpopID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_init,rpopID,'units','mg P L-1');
    netcdf.putAtt(nc_init,rpopID,'field','rpop, scalar, series')

    salID = netcdf.defVar(nc_init,'SAL','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,salID,'long_name','SALINITY');
    netcdf.putAtt(nc_init,salID,'units','psu');
    netcdf.putAtt(nc_init,salID,'field','sal, scalar, series')

    sitID = netcdf.defVar(nc_init,'SIT','float',[xrhodimID erhodimID s_rhodimID timedimID]);
    netcdf.putAtt(nc_init,sitID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_init,sitID,'units','mg SI L-1');
    netcdf.putAtt(nc_init,sitID,'field','sit, scalar, series')
end

if(bio_sed>0)
    CTEMPID = netcdf.defVar(nc_init,'CTEMP','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CTEMPID,'long_name','temperature');
    netcdf.putAtt(nc_init,CTEMPID,'units','degree');
    netcdf.putAtt(nc_init,CTEMPID,'time','ocean_time');
    netcdf.putAtt(nc_init,CTEMPID,'field','CTEMP, scalar, series');

    CPOPG1ID = netcdf.defVar(nc_init,'CPOPG1','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOPG1ID,'long_name','Labile particulate organic P');
    netcdf.putAtt(nc_init,CPOPG1ID,'units','MG P/M^3');
    netcdf.putAtt(nc_init,CPOPG1ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOPG1ID,'field','CPOPG1, scalar, series');

    CPOPG2ID = netcdf.defVar(nc_init,'CPOPG2','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOPG2ID,'long_name','Refractory particulate organic P');
    netcdf.putAtt(nc_init,CPOPG2ID,'units','MG P/M^3');
    netcdf.putAtt(nc_init,CPOPG2ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOPG2ID,'field','CPOPG2, scalar, series');

    CPOPG3ID = netcdf.defVar(nc_init,'CPOPG3','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOPG3ID,'long_name','Conservative particulate organic P');
    netcdf.putAtt(nc_init,CPOPG3ID,'units','MG P/M^3');
    netcdf.putAtt(nc_init,CPOPG3ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOPG3ID,'field','CPOPG3, scalar, series');

    CPONG1ID = netcdf.defVar(nc_init,'CPONG1','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPONG1ID,'long_name','Labile particulate organic N');
    netcdf.putAtt(nc_init,CPONG1ID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,CPONG1ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPONG1ID,'field','CPONG1, scalar, series');

    CPONG2ID = netcdf.defVar(nc_init,'CPONG2','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPONG2ID,'long_name','Refractory particulate organic N');
    netcdf.putAtt(nc_init,CPONG2ID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,CPONG2ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPONG2ID,'field','CPONG2, scalar, series');

    CPONG3ID = netcdf.defVar(nc_init,'CPONG3','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPONG3ID,'long_name','Conservative particulate organic N');
    netcdf.putAtt(nc_init,CPONG3ID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,CPONG3ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPONG3ID,'field','CPONG3, scalar, series');

    CPOCG1ID = netcdf.defVar(nc_init,'CPOCG1','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOCG1ID,'long_name','Labile particulate organic C');
    netcdf.putAtt(nc_init,CPOCG1ID,'units','MG C/M^3');
    netcdf.putAtt(nc_init,CPOCG1ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOCG1ID,'field','CPOCG1, scalar, series');

    CPOCG2ID = netcdf.defVar(nc_init,'CPOCG2','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOCG2ID,'long_name','Refractory particulate organic C');
    netcdf.putAtt(nc_init,CPOCG2ID,'units','MG C/M^3');
    netcdf.putAtt(nc_init,CPOCG2ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOCG2ID,'field','CPOCG2, scalar, series');

    CPOCG3ID = netcdf.defVar(nc_init,'CPOCG3','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOCG3ID,'long_name','Conservative particulate organic C');
    netcdf.putAtt(nc_init,CPOCG3ID,'units','MG C/M^3');
    netcdf.putAtt(nc_init,CPOCG3ID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOCG3ID,'field','CPOCG3, scalar, series');

    CPOSID = netcdf.defVar(nc_init,'CPOS','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CPOSID,'long_name','biogenic Si');
    netcdf.putAtt(nc_init,CPOSID,'units','MG Si/M^3');
    netcdf.putAtt(nc_init,CPOSID,'time','ocean_time');
    netcdf.putAtt(nc_init,CPOSID,'field','CPOS, scalar, series');

    PO4T2TM1SID = netcdf.defVar(nc_init,'PO4T2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,PO4T2TM1SID,'long_name','dissolved PO4');
    netcdf.putAtt(nc_init,PO4T2TM1SID,'units','MG P/M^3');
    netcdf.putAtt(nc_init,PO4T2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,PO4T2TM1SID,'field','PO4T2TM1S, scalar, series');

    NH4T2TM1SID = netcdf.defVar(nc_init,'NH4T2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,NH4T2TM1SID,'long_name','dissolved NH4');
    netcdf.putAtt(nc_init,NH4T2TM1SID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,NH4T2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,NH4T2TM1SID,'field','NH4T2TM1S, scalar, series');

    NO3T2TM1SID = netcdf.defVar(nc_init,'NO3T2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,NO3T2TM1SID,'long_name','dissolved NO3');
    netcdf.putAtt(nc_init,NO3T2TM1SID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,NO3T2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,NO3T2TM1SID,'field','NO3T2TM1S, scalar, series');

    HST2TM1SID = netcdf.defVar(nc_init,'HST2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,HST2TM1SID,'long_name','dissolved HS');
    netcdf.putAtt(nc_init,HST2TM1SID,'units','MG O2/M^3');
    netcdf.putAtt(nc_init,HST2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,HST2TM1SID,'field','HST2TM1S, scalar, series');

    SIT2TM1SID = netcdf.defVar(nc_init,'SIT2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,SIT2TM1SID,'long_name','dissolved SI');
    netcdf.putAtt(nc_init,SIT2TM1SID,'units','MG Si/M^3');
    netcdf.putAtt(nc_init,SIT2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,SIT2TM1SID,'field','SIT2TM1S, scalar, series');
    
    BNTHSTR1SID = netcdf.defVar(nc_init,'BNTHSTR1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,BNTHSTR1SID,'long_name','benthic stress');
    netcdf.putAtt(nc_init,BNTHSTR1SID,'units','null');
    netcdf.putAtt(nc_init,BNTHSTR1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,BNTHSTR1SID,'field','BNTHSTR1S, scalar, series');

    PO4T1TM1SID = netcdf.defVar(nc_init,'PO4T1TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,PO4T1TM1SID,'long_name','dissolved PO4');
    netcdf.putAtt(nc_init,PO4T1TM1SID,'units','MG P/M^3');
    netcdf.putAtt(nc_init,PO4T1TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,PO4T1TM1SID,'field','PO4T1TM1S, scalar, series');

    NH4T1TM1SID = netcdf.defVar(nc_init,'NH4T1TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,NH4T1TM1SID,'long_name','dissolved NH4');
    netcdf.putAtt(nc_init,NH4T1TM1SID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,NH4T1TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,NH4T1TM1SID,'field','NH4T1TM1S, scalar, series');

    NO3T1TM1SID = netcdf.defVar(nc_init,'NO3T1TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,NO3T1TM1SID,'long_name','dissolved NO3');
    netcdf.putAtt(nc_init,NO3T1TM1SID,'units','MG N/M^3');
    netcdf.putAtt(nc_init,NO3T1TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,NO3T1TM1SID,'field','NO3T1TM1S, scalar, series');

    HST1TM1SID = netcdf.defVar(nc_init,'HST1TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,HST1TM1SID,'long_name','dissolved HS');
    netcdf.putAtt(nc_init,HST1TM1SID,'units','MG O1/M^3');
    netcdf.putAtt(nc_init,HST1TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,HST1TM1SID,'field','HST1TM1S, scalar, series');

    SIT1TM1SID = netcdf.defVar(nc_init,'SIT1TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,SIT1TM1SID,'long_name','dissolved SI');
    netcdf.putAtt(nc_init,SIT1TM1SID,'units','MG Si/M^3');
    netcdf.putAtt(nc_init,SIT1TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,SIT1TM1SID,'field','SIT1TM1S, scalar, series');

    CH4T1TM1SID = netcdf.defVar(nc_init,'CH4T1TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CH4T1TM1SID,'long_name','methane');
    netcdf.putAtt(nc_init,CH4T1TM1SID,'units','MG O2/M^3');
    netcdf.putAtt(nc_init,CH4T1TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,CH4T1TM1SID,'field','CH4T1TM1S, scalar, series');

    CH4T2TM1SID = netcdf.defVar(nc_init,'CH4T2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,CH4T2TM1SID,'long_name','methane');
    netcdf.putAtt(nc_init,CH4T2TM1SID,'units','MG O2/M^3');
    netcdf.putAtt(nc_init,CH4T2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,CH4T2TM1SID,'field','CH4T2TM1S, scalar, series');

    SO4T2TM1SID = netcdf.defVar(nc_init,'SO4T2TM1S','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,SO4T2TM1SID,'long_name','sulfate');
    netcdf.putAtt(nc_init,SO4T2TM1SID,'units','MG O2/M^3');
    netcdf.putAtt(nc_init,SO4T2TM1SID,'time','ocean_time');
    netcdf.putAtt(nc_init,SO4T2TM1SID,'field','SO4T2TM1S, scalar, series');

    HSEDID = netcdf.defVar(nc_init,'HSED','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,HSEDID,'long_name','sediment depth');
    netcdf.putAtt(nc_init,HSEDID,'units','M');
    netcdf.putAtt(nc_init,HSEDID,'time','ocean_time');
    netcdf.putAtt(nc_init,HSEDID,'field','HSED, scalar, series');

    VSEDID = netcdf.defVar(nc_init,'VSED','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,VSEDID,'long_name','sedimentation velocity');
    netcdf.putAtt(nc_init,VSEDID,'units','CM/YEAR');
    netcdf.putAtt(nc_init,VSEDID,'time','ocean_time');
    netcdf.putAtt(nc_init,VSEDID,'field','VSED, scalar, series');

    VPMIXID = netcdf.defVar(nc_init,'VPMIX','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,VPMIXID,'long_name','particulate mixing velocity');
    netcdf.putAtt(nc_init,VPMIXID,'units','M^2/DAY');
    netcdf.putAtt(nc_init,VPMIXID,'time','ocean_time');
    netcdf.putAtt(nc_init,VPMIXID,'field','VPMIX, scalar, series');

    VDMIXID = netcdf.defVar(nc_init,'VDMIX','double',[xrhodimID erhodimID timedimID]);
    netcdf.putAtt(nc_init,VDMIXID,'long_name','dissloved mixing velocity');
    netcdf.putAtt(nc_init,VDMIXID,'units','M^2/DAY');
    netcdf.putAtt(nc_init,VDMIXID,'time','ocean_time');
    netcdf.putAtt(nc_init,VDMIXID,'field','VDMIX, scalar, series');
end

netcdf.close(nc_init)