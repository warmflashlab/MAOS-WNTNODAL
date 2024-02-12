
%% Make Max Intensity of a video 


tic
clear; close all; clc; clear memory;
addpath(genpath('~/Documents/GitHub/stemcells')); 

readpath = ['//Volumes/'];


maxDir=readpath;


% !!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!

takeoutfirstplane = 0; % do you want to delete the first Z?


% !!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!


files=dir(fullfile(readpath, '*.tif'));
nfiles=length(files);


maxDir = fullfile(maxDir,'Max/');
if ~exist(maxDir,'dir')
    mkdir(maxDir);
end



            for kk=1:nfiles
                clear img;
            fileIn = [files(kk).name];
            fileOut = ['MAXProj_',fileIn(1:end-4),'.tif'];
            fileDir= [[files(kk).folder],'/'];

            reader = bfGetReader([fileDir,fileIn]);

            nT = reader.getSizeT;
            nC = reader.getSizeC;
            nZ = reader.getSizeZ;

                sX = reader.getSizeX;
                sY = reader.getSizeY;
             
kk
disp('of');
nfiles

           
            img=uint16(zeros(sY,sX,1,nC,nT));
           

           for tt=1:nT
               for cc=1:nC
            
                   if takeoutfirstplane > 0
                   
            img(:,:,1,cc,tt) =  bfMaxIntensity_nofirst(reader,tt,cc);
            fileOut = ['nf_MAXProj_',fileIn(1:end-4),'.tif'];
                       
                   else
                   
            img(:,:,1,cc,tt) =  bfMaxIntensity(reader,tt,cc); % así se tienen que guardar en la variable para que respete el dimension order (x,y,z,c,t)
                   
                   end
               end
           end
 bfsave(img,[maxDir fileOut],'Compression','Uncompressed','dimensionOrder','XYZCT'); %XYZTC actually works but it believes z=t and ch=t
 
            end


   timetoc= toc;
   
disp('Max Intensity done');
disp('Elapsed time is');
    timetoc/60
    disp(' minutes')

    