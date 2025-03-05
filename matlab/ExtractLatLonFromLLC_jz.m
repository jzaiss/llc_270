% Example extraction of latlon region from LLC domain

addpath('/nobackup/jzaissbo/matlabcode/lib');

% grid dimension
nx=270;

% latlon region
minlat=-70;
maxlat=57;


% grid directory locations
gin='/nobackup/jzaissbo/llc_270_global/rawgrid/';
gout= '/nobackup/jzaissbo/llc_270_global/grid/';

% extract indices for latlon region
prec='real*4';
fnam=[gin 'YC.data'];
tmp=read_llc_fkij(fnam,nx,1,1,1);
tmp(~tmp)=nan;
fc=[1 2 4 5];
ix{1}=1:nx;
ix{2}=1:nx;
ix{3}=[];
ix{4}=1:nx;
ix{5}=1:nx;
jx{1}=find(tmp>=minlat&tmp<=maxlat);
jx{2}=find(tmp>=minlat&tmp<=maxlat);
jx{3}=[];
jx{4}=find(tmp>=minlat&tmp<=maxlat);
jx{5}=find(tmp>=minlat&tmp<=maxlat);
m(1)=0;
for f=1:length(fc)
    m(f+1)=length(ix{fc(f)});
end
n=length(jx{fc(1)});

% get and save grid information
eval(['cd ' gout])
for fnm={'XC','XG','YC','YG','hFacC'}
    fin=[gin fnm{1} '.data'];
    fout=[fnm{1} '_' int2str(sum(m)) 'x' int2str(n)];
    for f=1:length(fc)
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
    switch fnm{1}
      case {'AngleCS','U2zonDir'}, fld(:)=1;
      case {'AngleSN','V2zonDir'}, fld(:)=0;
      case 'RAC'
        fld2=fld;
        for j=1:(nx*2)
            tmp=double(fld(:,j));
            fld2(:,j)=mean(tmp(find(tmp)));
        end
        fld=fld2;
        case {'XC','XG'}, fld=fld(:,500)*ones(1,2*nx);
        case {'YC','YG'}, fld=ones(4*nx,1)*fld(500,:);
    end
    writebin(fout,fld);
end


% example extraction of a tracer field
%fnam='/nobackup/dcarrol2/LLC_540_africa/MITgcm/run/diags/hourly_SST.0001209756.data';
%fld=zeros(sum(m),n);
%for f=1:length(fc)
%    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
%        read_llc_fkij(fnam,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
%end
%figure(1), clf, quikpcolor(fld'), colorbar

% example extraction of a vector field
% note that zonal velocity is U in faces 1/2 and V in faces 4/5
% and meridional velocity is V in faces 1/2 and -U in faces 4/5

%%

%Convert 3D velocity binary files to 3 1D binary files: [270, 3150, 50, 3] --> 3x [270, 3150, 50]

dataDir = '/nobackup/hzhang1/ITER/iter42_1992_2023/for_Zaiss/daily/'; %change to hong directory

saveDir1 = '/nobackup/jzaissbo/llc_270_global/uvel/';
saveDir2 = '/nobackup/jzaissbo/llc_270_global/vvel/';
saveDir3 = '/nobackup/jzaissbo/llc_270_global/wvel/';

VelFiles = dir([dataDir 'daily_velmass_3d*data']);

for i = 1:length(VelFiles)

	Vel = readbin([dataDir VelFiles(i).name], [270, 3150, 50, 3]);

    U=Vel(:,:,:,1);
    V=Vel(:,:,:,2);
    W=Vel(:,:,:,3);

    str=VelFiles(i).name;
    newStr=extractBetween(str,".","data",'Boundaries','inclusive');

    fnu=['uvel' newStr{1}];
    fnv=['vvel' newStr{1}];
    fnw=['wvel' newStr{1}];

    writebin([saveDir1 fnu],U);
    writebin([saveDir2 fnv],V);
    writebin([saveDir3 fnw],W);

    disp(VelFiles(i).name);

end


dataDir1 = '/nobackup/jzaissbo/llc_270_global/uvel/raw';
dataDir2 = '/nobackup/jzaissbo/llc_270_global/vvel/raw';
dataDir3 = '/nobackup/hzhang1/ITER/iter42_1992_2023/for_Zaiss/daily'; 

saveDir1 = '/nobackup/jzaissbo/llc_270_global/uvel/';
saveDir2 = '/nobackup/jzaissbo/llc_270_global/vvel/';
saveDir3 = '/nobackup/jzaissbo/llc_270_global/no3/';

uVelFiles = dir([dataDir1 'uvel*data']);
vVelFiles = dir([dataDir2 'vvel*data']);
NO3Files = dir([dataDir3 'daily_NO3*data']);


for i = 1:length(uVelFiles)

finu = [dataDir1 uVelFiles(i).name];
finv = [dataDir2 vVelFiles(i).name];
finn = [dataDir3 NO3Files(i).name];

%Velo files
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finu,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)});
        fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finv,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finv,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)});
        fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = - ...               % <<<<<<<<
            read_llc_fkij(finu,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
    
end

%Nitrate Files
fld=zeros(sum(m),n);
for f=1:length(fc)
   fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
       read_llc_fkij(finn,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
end

str=NO3Files(i).name;
newStr=extractBetween(str,".","data",'Boundaries','inclusive');

fnn = ['nitrate' newStr{1}];

writebin([saveDir1 uVelFiles(i).name],fldu);
writebin([saveDir2 vVelFiles(i).name],fldv);
writebin([saveDir3 fnn],fld);

disp(uVelFiles(i).name);

end

%figure(2), clf, quikpcolor(fldu'), colorbar
%figure(3), clf, quikpcolor(fldv'), colorbar

