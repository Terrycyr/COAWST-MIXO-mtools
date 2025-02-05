clear all;close all;

[~,~,loc] = xlsread('./all locations data in surfer.xls');
[~,~,IN] = xlsread('./Inorganicnutrientscomposite.xls');
[~,~,DON] = xlsread("FINAL EB DON DATA.xls");

[r1,c1] = size(loc);

for i=2:r1
   if(~ischar(loc{i,1}))
       loc{i,1} = num2str(loc{i,1});
   end
end

[r2,c2] = size(IN);

%Inorganic nitrogen
k=0;
k2=0;
for i = 1:r2
    if(strcmp(IN{i,1},'Date'))
        continue;
    end 
    tmp = strsplit(IN{i,1},'/');
    try
        year = str2num(tmp{3});
        mon = str2num(tmp{1});
        day = str2num(tmp{2});
    catch
        continue;
    end
    k=k+1;
    IN_tmp{k,1} = datenum(year,mon,day);
    if(isnan(IN{i,5}))
        if(ischar(IN{i,2}))
            IN_tmp{k,2} = IN{i,2};
        else
            IN_tmp{k,2} = num2str(IN{i,2});
        end
    else
        if(ischar(IN{i,5}))
            IN_tmp{k,2} = IN{i,5};
        else
            IN_tmp{k,2} = num2str(IN{i,5});
        end
    end
    pos = find(strcmp(loc(:,1),IN_tmp{k,2}));
    if(~isempty(pos))
        k2=k2+1;
        IN_out(k2,1:2) = IN_tmp(k,1:2);
        IN_out(k2,3) = loc(pos,2);
        IN_out(k2,4) = loc(pos,3);
        if(ischar(IN{i,3}))
            IN_out{k2,5} = str2num(IN{i,3});
        else
            IN_out{k2,5} = IN{i,3};
        end
        for col = [6 7 10]
            if(isnumeric(IN{i,col}))
                if(IN{i,col}<0.0001)
                    IN{i,col} = NaN;
                end
            else
                try
                    IN{i,col} = str2num(IN{i,col});
                    if(IN{i,col}<0.0001||isempty(IN{i,col}))
                        IN{i,col} = NaN;
                    end
                catch
                    IN{i,col} = NaN;
                end

            end
        end
        IN_out(k2,6) = IN(i,6);
        IN_out(k2,7) = IN(i,7);
        IN_out(k2,8) = IN(i,10);
    end      
end

IN_num = cell2mat(IN_out(:,[1,3:8]));
save('Inorganicnutrients.mat','IN_out','IN_num');

figure;
scatter(IN_num(:,3),IN_num(:,2));

%DON
[r3,c3] = size(DON);
k=0;
for i=3:r3
    if(~isnumeric(DON{i,5}(1))||isnan(DON{i,5}(1)))
        continue;
    end
    year = str2num(DON{i,1}(end-1:end));
    if(year>90)
        year=year+1900;
    else
        year=year+2000;
    end
    mon = str2num(DON{i,1}(end-3:end-2));
    day = 1;
    pos = find(strcmp(loc(:,1),num2str(DON{i,2})));
    if(~isempty(pos))
        k=k+1;
        DON_out{k,1} = datenum(year,mon,day);
        if(isnumeric(DON{i,2}))
            DON{i,2} = num2str(DON{i,2});
        end
        DON_out{k,2} = DON{i,2};
        DON_out{k,3} = loc{pos,2};
        DON_out{k,4} = loc{pos,3};
        DON_out{k,5} = DON{i,3};
        tmp1 = DON{i,5};
        tmp2 = DON{i,6};
        if(isnan(tmp2))
            DON_out{k,6} = tmp1;
        else
            DON_out{k,6} = tmp2;
        end
    else
        continue;
    end
end

DON_num = cell2mat(DON_out(:,[1 3:6]));
save('Organicnutrients.mat','DON_out','DON_num');

figure;
scatter(DON_num(:,3),DON_num(:,2));




