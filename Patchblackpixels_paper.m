%% Patch black pixels 

%this is for patching the black pixels with a background value to aid in
%the lookup tables finding the mins and maxs

clear; close all; clc; clear memory;
addpath(genpath('~/Documents/GitHub/stemcells')); 

readpath =['/Volumes//'];

patchDir = fullfile(readpath,'Patched/');
if ~exist(patchDir,'dir')
    mkdir(patchDir);
end


files= dir(fullfile(readpath,'*.tif'));
nfiles=length(files);


% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡

bgvalue= 550; % this is manually input, a value of background to patch zeros

% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡


%% Read files / locate black pixels / Patch them with background value


for ii=1:nfiles
clear reader; clear original; clear img; clear slice; clear Yidx; clear Xidx;
disp('Reading original file...');
    reader=bfGetReader([files(ii).folder '/' files(ii).name]);
    
nC = reader.getSizeC;
nT = reader.getSizeT;
nZ = reader.getSizeZ;
sX = reader.getSizeX;
sY = reader.getSizeY;

disp('file - Y - X - Channels - Time - Z');
[ii sY sX nC nT nZ]

original=bfGetPlane(reader,reader.getIndex(0,0,0)+1);
[Yidx, Xidx]=find(original<1);

if sum(any(original<1))>0
    
img=uint16(zeros(sY,sX,nC,nZ,nT));
                %for reader : z,c,t 
                %for saving correctly make a matrix (y,x,c,z,t)
                disp('Patching zeros...');
                
                clear Yidx; clear Xidx;
for cc=1:nC
 for zz=1:nZ
     for tt=1:nT

         slice=bfGetPlane(reader,reader.getIndex(zz-1,cc-1,tt-1)+1);
         aux=find(slice<1);
         
         slice(aux)=bgvalue;
         img(:,:,cc,zz,tt)=slice;
         
         [cc zz tt]
         
     end
 end
end
disp('Saving...');
disp(ii);
bfsave(img,[patchDir 'patched-' files(ii).name(1:end-4) '.tif'],'Compression','Uncompressed','dimensionOrder','XYCZT');
 

   
disp('Moving to next file');
else
cprintf('*blue','this file does not have black pixels');
    disp(files(ii).name)
cprintf('*blue','this file does not have black pixels');
end
end
disp('All files patched');
