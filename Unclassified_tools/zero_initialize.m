function zero_initialize(fn)
%   Detailed explanation goes here
ncid = netcdf.open(fn,'WRITE');
flag = 0;
i=0;
while(flag==0)
    try
        varname = netcdf.inqVar(ncid,i);
        varid = netcdf.inqVarID(ncid,varname);
        data = netcdf.getVar(ncid,varid);
        data = zeros(size(data));
        netcdf.putVar(ncid,varid,data);
        i=i+1;
    catch
        flag=1;
    end
end
netcdf.close(ncid);
end

