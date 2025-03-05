% Example extraction of latlon region from LLC domain

%addpath(genpath('/nobackup/dcarrol2/MATLAB'));

% grid dimension
nx=270;
ny=3150;
nz=50;

% latlon region
minlat=-70;
maxlat=57;

% grid directory locations
gin='C:\Users\jessicaz\Documents\research\llc_270\rawgrid\';
gout= 'C:\Users\jessicaz\Documents\research\llc_270\grid\';

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

% Southwest center points, no direction, fill blank tiles
fnx='DXC';
fny='DYC';
finx=[gin fnx '.data'];
finy=[gin fny '.data'];
foutx=[fnx '_' int2str(sum(m)) 'x' int2str(n)];
fouty=[fny '_' int2str(sum(m)) 'x' int2str(n)];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
end
fldx=ones(4*nx,1)*fldx(1,:);
fldy=ones(4*nx,1)*fldy(1,:);
writebin(foutx,fldx);
writebin(fouty,fldy);

% Southwest center points, no direction, fill blank tiles
fnx='RAW';
fny='RAS';
finx=[gin fnx '.data'];
finy=[gin fny '.data'];
foutx=[fnx '_' int2str(sum(m)) 'x' int2str(n)];
fouty=[fny '_' int2str(sum(m)) 'x' int2str(n)];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
end
fldx2=fldx;
for j=1:(nx*2)
    tmp=double(fldx(:,j));
    fldx2(:,j)=mean(tmp(find(tmp)));
end
fldx=fldx2;
fldy2=fldy;
for j=1:(nx*2)
    tmp=double(fldy(:,j));
    fldy2(:,j)=mean(tmp(find(tmp)));
end
fldy=fldy2;
writebin(foutx,fldx);
writebin(fouty,fldy);

% Southwest center points, no direction
fnx='hFacW';
fny='hFacS';
finx=[gin fnx '.data'];
finy=[gin fny '.data'];
foutx=[fnx '_' int2str(sum(m)) 'x' int2str(n)];
fouty=[fny '_' int2str(sum(m)) 'x' int2str(n)];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);

% Southwest edges, no direction, fill blank tiles
fnx='DXG';
fny='DYG';
finx=[gin fnx '.data'];
finy=[gin fny '.data'];
foutx=[fnx '_' int2str(sum(m)) 'x' int2str(n)];
fouty=[fny '_' int2str(sum(m)) 'x' int2str(n)];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
end
fldx=ones(4*nx,1)*fldx(1,:);
fldy=ones(4*nx,1)*fldy(1,:);
writebin(foutx,fldx);
writebin(fouty,fldy);

% Southwest corner (vorticity) points, no direction, fill blank tiles
fnm='RAZ';
fin=[gin fnm '.data'];
finy=[gin fny '.data'];
fout=[fnm '_' int2str(sum(m)) 'x' int2str(n)];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
end
fld=ones(4*nx,1)*fld(1,:);
writebin(fout,fld);

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

dataDir = 'C:\Users\jessicaz\Documents\research\llc_270\daily3d\'; %change to hong directory

saveDir1 = 'C:\Users\jessicaz\Documents\research\llc_270\uvel\';
saveDir2 = 'C:\Users\jessicaz\Documents\research\llc_270\vvel\';
saveDir3 = 'C:\Users\jessicaz\Documents\research\llc_270\wvel\';

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


dataDir1 = 'C:\Users\jessicaz\Documents\research\llc_270\uvel\';
dataDir2 = 'C:\Users\jessicaz\Documents\research\llc_270\vvel\';
dataDir3 = 'C:\Users\jessicaz\Documents\research\llc_270\daily3d\'; %change to hong directory

saveDir1 = 'C:\Users\jessicaz\Documents\research\llc_270\uvel\';
saveDir2 = 'C:\Users\jessicaz\Documents\research\llc_270\vvel\';
saveDir3 = 'C:\Users\jessicaz\Documents\research\llc_270\no3\';

uVelFiles = dir([dataDir1 'uvel*data']);
vVelFiles = dir([dataDir2 'vvel*data']);
NO3Files = dir([dataDir3 'daily_NO3*data']); %----> change to nitrate name


for i = 1:length(uVelFiles)

finu = [dataDir1 uVelFiles(i).name];
finv = [dataDir2 vVelFiles(i).name];
finn = [dataDir3 NO3Files(i).name];

%Velo files
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finu,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)});
        fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finv,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finv,nx,fc(f),1:10,ix{fc(f)},jx{fc(f)});
        fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:) = - ...               % <<<<<<<<
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

disp(uVelFiles(i).name);

end

%figure(2), clf, quikpcolor(fldu'), colorbar
%figure(3), clf, quikpcolor(fldv'), colorbar

