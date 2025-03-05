addpath('C:\Users\jessicaz\Documents\MATLAB\lib')
cd('C:\Users\jessicaz\Documents\research\llc_270')

dataDir = 'C:\Users\jessicaz\Documents\research\llc_270\daily3d\';

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

    disp(VelFiles(i).name)

end




test=quikread_llc([saveDir1 VelFiles(1).name],270);

quikplot_llc(test)
c=colorbar();
caxis([-0.5 0.5])