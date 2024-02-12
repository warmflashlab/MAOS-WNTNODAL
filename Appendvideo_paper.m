%% append files to their corresponding tif (only for videos)

clear; clc; clear memory;

pathtoread = '//Volumes//';%where are the files 
EdataDir = '//Volumes/MAOS 4/LR-BcatNodalhethom/';

AppendDir=[EdataDir 'Appended/'];
if ~exist(AppendDir,'dir')
    mkdir(AppendDir);
end

addpath(genpath('~/Documents/GitHub/stemcells')); 

uppath=strfind(pathtoread,'/');



% YOU HAVE TU INPUT THE TIMES IN ONE POSITION AND THE NUMBER OF CHANNELS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

timeinonexy = 2; %how many time files of one xy position are there.
channels = 2; 

%  YOU HAVE TO INPUT THE TIMES IN ONE POSITION AND THE NUMBER OF CHANNELS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

files= dir(fullfile(pathtoread,'*.tif'));
nfiles=length(files);
filesname=cell(1,nfiles);

for ii=1:nfiles
    filesname{ii}=files(ii).name;
end

fullDF={};
for ii=1:nfiles
fullDF{ii}=fullfile(pathtoread,files(ii).name);
end

nT=23; %just for creating the variable 
nZ=10; %just for making a variable

numstosum=0:timeinonexy-1;

divisor=1:timeinonexy:nfiles;
list=zeros(length(divisor),timeinonexy);

for ii=1:length(divisor)
list(ii,:)=divisor(ii)+numstosum;
end



 %% The main appender with multiple Z and T 
tic
%erase the last frame?
eraselastframe=1;


for hh = 1:length(divisor)%numero de files de 3 en 3
    
    hh
    
    disp('Extracting data');
   
for gg = 1:timeinonexy
  
   reader{gg}=bfGetReader(fullDF{list(hh,gg)});
    nT(gg)=reader{gg}.getSizeT;
    nZ(gg)=reader{gg}.getSizeZ;
    nX(gg)=reader{gg}.getSizeX;
    nY(gg)=reader{gg}.getSizeY;
   
   
   imgholder{gg}=bfopen(fullDF{list(hh,gg)});
    
   
    
    list(hh,gg)

end

catholder={};


for gg=1:timeinonexy
    
    if gg==1
    %catholder=cat(1,imgholder{1}{1},imgholder{gg}{1});
    catholder=imgholder{1}{1};
    else %timeinonexy>2
    catholder=cat(1,catholder,imgholder{gg}{1});
    end

end

disp('Temporary concatenated image created');


chlistidxs=[]; %images idxs divided by channels
for cc=1:channels
chlistidxs(cc,:)=cc:channels:length(catholder);
end
chlistidxs=chlistidxs';

disp('Channels indexes created');


ztlistidxs=[]; 
 itemplist=1:nZ(1):length(chlistidxs);
 ftemplist=nZ(1):nZ(1):length(chlistidxs);

for tt=1:sum(nT)
    counter=1;
for zzss=itemplist(tt):ftemplist(tt)

ztlistidxs(counter,:,tt)=chlistidxs(zzss,:);
counter=counter+1;

end
end

disp('Z and T indexes created');
disp('re-arranging bf file');

 if eraselastframe==1
 bigimg=uint16(zeros(1024,1024,channels,nZ(1),sum(nT)-1));
 else
 bigimg=uint16(zeros(1024,1024,channels,nZ(1),sum(nT)));
 end

 
for cc=1:channels
for zz=1:nZ(1)
    
    if eraselastframe==1    
        for tt=1:sum(nT)-1
 bigimg(:,:,cc,zz,tt)=catholder{ztlistidxs(zz,cc,tt),1};
 [cc zz tt]
        end      
    else
        for tt=1:sum(nT)
 bigimg(:,:,cc,zz,tt)=catholder{ztlistidxs(zz,cc,tt),1};
 [cc zz tt]
        end
    end
 
end
end

 
 filename=strcat('a-h-',filesname{divisor(hh)}(1:end-10),'.tif');

hh
disp('Saving file...');
bfsave(bigimg,[AppendDir filename],'Compression','Uncompressed','dimensionOrder','XYCZT'); % it has to be order even if it doesn't make sense to the matrix! 

end

disp('Files Appended');
imshow('done meme.png')
   timetoc= toc;
disp('Elapsed time is');
    timetoc/60
    disp(' minutes')
    
   