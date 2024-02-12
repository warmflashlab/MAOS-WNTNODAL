%% for Nodal movies, read all pixels
%First, locate files

tic
clear; close all; clc; clear memory;
addpath(genpath('~/Documents/GitHub/stemcells')); 

readpath ='//Volumes//'; 
maxpath = '//Volumes/MAOS 4/LR-BcatNodalhethom/Appended/Max/';

experiment = 'nodals-';

files=dir(fullfile(readpath,'*.tif'));
maxfiles=dir(fullfile(maxpath,'*.tif'));
nfiles=length(files); nmaxfiles=length(maxfiles);

if nfiles == nmaxfiles
    disp('same number of files for max and stacks');
else
    error('not same number of max and 3d files');
end

resultsDir = fullfile(readpath, [experiment '-DFCMImat/']);
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
else
    
load([resultsDir experiment '-DFCMI.mat']);
load([resultsDir experiment '-DFCMI_max.mat']);
load([resultsDir experiment '-DFCMI_files.mat']);
load([resultsDir experiment '-bins.mat']);
load([resultsDir experiment '-spacebinned.mat']);

disp('previous matrix loaded');
end


%% Read the pixel intensities of movies if you do not have read them
% SKIP if you have loaded mats

channeltoread = 2;
objfactor=0.625;
bins=20;

spacebinned=0:bins:360;

DFCMI={};
DFCMI_file=cell(1,nfiles);

            for ii=1:nfiles
                clear reader; clear iPlane; clear frame; clear dfcenter; clear MI;
                
            fileIn = [files(ii).name];
            fileDir= [[files(ii).folder],'/'];

            reader = bfGetReader([fileDir,fileIn]);
            nT = reader.getSizeT;
            %nC = reader.getSizeC;
            nZ = reader.getSizeZ;
                sX = reader.getSizeX;
                sY = reader.getSizeY;
             
ii
disp('of');
nfiles

center = [sY/2, sX/2];
 idxs = 1:(sY*sX);
[y1,x1] = ind2sub([sY,sX], idxs);


dfcenter = sqrt((y1-center(1)).^2 + (x1-center(2)).^2).*objfactor;

           for tt=1:nT
               for zz=1:nZ

            
                   iPlane=reader.getIndex(zz-1,channeltoread-1,tt-1)+1;
                   frame = bfGetPlane(reader,iPlane);
               
                        temp.DFCMI = [dfcenter(:), frame(:)];
                        temp.DFCMI(temp.DFCMI(:,2)==0,:)=[]; % to delete zero values for quantif
                        % for binning
                        
                        temp.binMI=zeros(length(spacebinned)-1, 2);
                        
                            for kk=1:length(spacebinned)-1
                            
                                temp.idx=find(temp.DFCMI(:,1)>spacebinned(kk) & temp.DFCMI(:,1)<=spacebinned(kk+1));
                                temp.list=temp.DFCMI(temp.idx,2);
                                
                                temp.binMI(kk,:)=[spacebinned(kk+1) mean(temp.list)+std(single(temp.list))];

                            end
                           
                            temp.fDFCMI{tt}{zz}=temp.binMI ;           
               end
           end

 DFCMI_file{ii}=[files(ii).name];
 DFCMI{ii}=temp.fDFCMI;
 temp.fDFCMI={};
           
            end
            
save([resultsDir experiment '-DFCMI.mat'],'DFCMI');
save([resultsDir experiment '-DFCMI_files.mat'],'DFCMI_file');
save([resultsDir experiment '-bins.mat'],'bins');
save([resultsDir experiment '-spacebinned.mat'],'spacebinned');
                   
   timetoc= toc;
   
disp('Single Voxel quantification done');
disp('Elapsed time is');
    timetoc/60
    disp(' minutes')
    
%% Quantification for max projections
% SKIP if you have loaded mats

channeltoread = 2;
objfactor=0.625;
bins=20;

spacebinned=0:bins:360;
DFCMI_max={};



            for ii=1:nfiles
                clear reader; clear iPlane; clear frame; clear dfcenter; clear MI;
                
            fileIn = [maxfiles(ii).name];
            reader = bfGetReader([maxpath,fileIn]); %make sure there's the /
            nT = reader.getSizeT;
            %nC = reader.getSizeC;
            nZ = reader.getSizeZ;
                sX = reader.getSizeX;
                sY = reader.getSizeY;
             
ii
disp('of');
nfiles

center = [sY/2, sX/2];
 idxs = 1:(sY*sX);
[y1,x1] = ind2sub([sY,sX], idxs);

dfcenter = sqrt((y1-center(1)).^2 + (x1-center(2)).^2).*objfactor;

           for tt=1:nT
               for zz=1:nZ


                   iPlane=reader.getIndex(zz-1,channeltoread-1,tt-1)+1;
                   frame = bfGetPlane(reader,iPlane);
                   
        
                        temp.DFCMI = [dfcenter(:), frame(:)];
                        temp.DFCMI(temp.DFCMI(:,2)==0,:)=[]; % to delete zero values for quantif
                        % for binning
                        
                        temp.binMI=zeros(length(spacebinned)-1, 2);
                        
                            for kk=1:length(spacebinned)-1
                            
                                temp.idx=find(temp.DFCMI(:,1)>spacebinned(kk) & temp.DFCMI(:,1)<=spacebinned(kk+1));
                                temp.list=temp.DFCMI(temp.idx,2);
                                
                                temp.binMI(kk,:)=[spacebinned(kk+1) mean(temp.list)+std(single(temp.list))];

                            end
                           
                            temp.fDFCMI{tt}{zz}=temp.binMI ;           
               end
           end

 DFCMI_max{ii}=temp.fDFCMI;
 temp.fDFCMI={};
           
            end
   
save([resultsDir experiment '-DFCMI_max.mat'],'DFCMI_max');

                   
   timetoc= toc;
   
disp('Max projections quantifications done');
disp('Elapsed time is');
    timetoc/(60)
    disp(' minutes')
    
    
%% Group

channeltoread = 2;
objfactor=0.625; %objective factor um/pixel


%groups
 gd.description = {'hetNODAL-ctrl','hetNODAL-ACT','hetNODAL-BMP'};
% 
 gd.days = [1,2,3,4,5,6,7,8];
 gd.conditions = [1,2,3,4,5,6,7,8];
 
 d1time = 30;
 d2time = 30;
 d3time = 30;


if numel(gd.days) ~= length(gd.description)
    error('the agroupation of days doesnt match the description');
else
     cprintf('blue','days match description');
     disp('days match description');
end

if numel(gd.conditions) ~= length(gd.description)
    error('the agroupation of conditions doesnt match the description');
else
    disp(' ');
    cprintf('blue','conditions match description');
     disp('days match description');
end

%indicate which files are grouped
 groups{1}= [];  groups{2}= [];  groups{3}= [];  groups{4}= [];
  groups{5}= [];  groups{6}= [];  groups{7}= [];  groups{8}= []; 

check=length([groups{:}]);

if check==nfiles

    disp('all files are within a group');
else
disp(' WARNING : NOT ALL FILES ARE IN GROUPS !!! ');
end


%% calculations 

gd.DFCMI={}; gd.DFCMI_max={};

for gg = 1:length(groups)
for ssgg=1:length(groups{gg})    
    
    for tt=1:length(DFCMI{groups{gg}(ssgg)})
        for zz=1:length(DFCMI{groups{gg}(ssgg)}{tt})
        
    if ssgg==1
gd.DFCMI{gg}{tt}{zz}(:,ssgg)=DFCMI{groups{gg}(ssgg)}{tt}{zz}(:,2);          
    else
gd.DFCMI{gg}{tt}{zz}(:,ssgg)=DFCMI{groups{gg}(ssgg)}{tt}{zz}(:,2);      
    end
        end
    
    end
    
end
end


for gg = 1:length(groups)
for ssgg=1:length(groups{gg})    
    
    for tt=1:length(DFCMI_max{groups{gg}(ssgg)})
        
        
    if ssgg==1
gd.DFCMI_max{gg}{tt}{1}(:,ssgg)=DFCMI_max{groups{gg}(ssgg)}{tt}{1}(:,2);          
    else
gd.DFCMI_max{gg}{tt}{1}(:,ssgg)=DFCMI_max{groups{gg}(ssgg)}{tt}{1}(:,2);      
    end
       
    
    end
    
end
end


temp.size=length(spacebinned)-1;

for gg=1:length(groups)
    for tt = 1:length(gd.DFCMI{gg})
        for zz = 1:length(gd.DFCMI{gg}{tt})
   
            if sum(sum(gd.DFCMI{gg}{tt}{zz}==0))>1
                temp.idx=find(gd.DFCMI{gg}{tt}{zz}==0);
                gd.DFCMI{gg}{tt}{zz}(temp.idx)=[];
                temp.size2=length(gd.DFCMI{gg}{tt}{zz})/temp.size;
                gd.DFCMI{gg}{tt}{zz}=reshape(gd.DFCMI{gg}{tt}{zz},[temp.size, temp.size2]);
               
                disp('some data has extra t points');
            end
            
            gd.DFCMI{gg}{tt}{zz} = [spacebinned(2:end)' mean(gd.DFCMI{gg}{tt}{zz},2)+std(gd.DFCMI{gg}{tt}{zz},0,2)];
            º
        end
    end
   
end
%for maX:
for gg=1:length(groups)
    for tt = 1:length(gd.DFCMI_max{gg})
   
            if sum(sum(gd.DFCMI_max{gg}{tt}{1}==0))>1
                temp.idx=find(gd.DFCMI_max{gg}{tt}{1}==0);
                gd.DFCMI_max{gg}{tt}{1}(temp.idx)=[];
                temp.size2=length(gd.DFCMI_max{gg}{tt}{1})/temp.size;
                gd.DFCMI_max{gg}{tt}{1}=reshape(gd.DFCMI_max{gg}{tt}{1},[temp.size, temp.size2]);
               
                disp('some data has extra t points');
            end
            
            gd.DFCMI_max{gg}{tt}{1} = [spacebinned(2:end)' mean(gd.DFCMI_max{gg}{tt}{1},2)+std(gd.DFCMI_max{gg}{tt}{1},0,2)];

    end
   
end

%you dont need the next one because you have extracted the max data already

donotcountfirstZ = 0; %in zero the max proj will be used, in 1 the 3D quant will be used
gd.maxZ_DFCMI={};

for gg=1:length(groups)
    for tt = 1:length(gd.DFCMI{gg})
    nT(gg)=tt;
        
        
        
        if donotcountfirstZ >0
             allZextracted = cell2mat(gd.DFCMI{gg}{tt});
           allZextracted = allZextracted(:,3:end); 
           gd.maxZ_DFCMI{gg}{tt} = [allZextracted(:,1) max(allZextracted,[],2)];
        else
            allZextracted = cell2mat(gd.DFCMI_max{gg}{tt});
        gd.maxZ_DFCMI{gg}{tt} = [allZextracted(:,1) max(allZextracted,[],2)];
        end

    end
end

%assemble grid
grid={};

for gg=1:length(groups)
    grid{gg} = zeros(length(spacebinned)-2,nT(gg));
    for tt = 1:length(gd.DFCMI{gg})
      for bb = 1:length(spacebinned)-2
   
            grid{gg}(bb,tt)=gd.maxZ_DFCMI{gg}{tt}(bb,2);

      end
    end

end


for gg=1:length(groups)
    allmaxes(gg)=max(max(cell2mat(gd.maxZ_DFCMI{gg})));
end

bridge = 5; %how many x points you want it to be bridged between days

%now put this for all grids
for dd=1:length(days)-1
for ii=1:length(gd.conditions)
    temp_size=size(grid{gd.conditions(ii,dd)});
    leftchunk=zeros(bridge,temp_size(2));
    rightchunk=zeros(bridge,temp_size(2));
    
    leftchunk=grid{gd.conditions(ii,dd)}(:,end-(bridge-1):end);
    rightchunk=grid{gd.conditions(ii,dd)+1}(:,1:bridge);
    
    temp_chunk=(leftchunk+rightchunk)/2;
   % [dd ii]
    
    grid{gd.conditions(ii,dd)}(:,end-(bridge-1):end) = temp_chunk;
    grid{gd.conditions(ii,dd)+1}(:,1:bridge) = temp_chunk;
end
end


smoothlevel = 15;
smoothgrid={};
for gg=1:length(groups)
smoothgrid{gg}=movmean(grid{gg},smoothlevel);
smoothgrid{gg}=movmean(smoothgrid{gg},smoothlevel,2);

end

%% make the graph

%close all;
makeblack=0;

screensize = get(groot, 'ScreenSize');

cp=colormap(hot);

colormap(cp);


tiles=reshape(1:length(groups),length(days),length(conditions))';
posforbar = round(length(conditions)/2); %position of the tile for the color bar
if length(days)<3
    posford23 = tiles(:,end);
else
posford23 = tiles(:,2:3); %positions for all in day 2 and 3
end
posford23 = reshape(posford23,1,[]);
posforconditions = tiles(:,1)';
posford1 = tiles(:,1)';

% 

if makeblack==0
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4))]);
lettercolor='k'; oppositeLC='w'; oppositetriplet=[1 1 1]; colortriplet = [0 0 0];


colormap(cp);
else
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4))],'color','black');
lettercolor='w'; oppositeLC='k'; oppositetriplet=[0 0 0]; colortriplet = [1 1 1];
colormap(heatforblack);
colormap(cp);
end
%for black
t=tiledlayout(length(conditions),length(days),'TileSpacing','none','Padding','loose');
sgtitle('Kymograph Distance from Center','FontSize',20,'FontWeight','bold','color',lettercolor);
ylabel(t,'Distance from center (µm)','FontWeight','bold','color',lettercolor);
xlabel(t,'Hours','FontWeight','bold','FontSize',15,'color',lettercolor);

for gg=1:length(groups)
    nexttile

 
pplot(gg)=imagesc(smoothgrid{gg},'Interpolation','bilinear');

pplot(gg).YData = [0:bins:370]';

caxis([mean(allmaxes(1:3)) mean(allmaxes(4:end))]);


if makeblack==1
axis=gca
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
set(gca,'color','black')
else
    axis=gca
    set(gca,'YDir','normal');
end


if any(gg==posford23)
set(gca,'Ytick',[],'YColor','none','box','off','YTickLabel',[],'layer','top');
pplot(gg).Parent.YColor=oppositetriplet;

end

if gg==tiles(posforbar,end)
cb(gg)=colorbar;
ylabel(cb(gg),'Nodal Intensity (A.U.)','FontWeight','bold','color',lettercolor);
cb(gg).Color = colortriplet;
cb(gg).Box='off'
end

if any(gg==posforconditions)
    tempidx=find(gg==posforconditions);
    ylabel(conditions{tempidx},'FontSize',15);
end

if any(gg==posford1)
pplot(gg).XData=pplot(gg).XData*(d1time/60);
xlim([0.5 22]);
axis.YAxis.FontSize=15;
set(gca,'box','off','Ycolor',lettercolor,'layer','top');
%
axis.XTickLabel=[];
axis.YTick=[0:50:370];
axis.YTickLabelMode='auto';
end

if length(days)>1
if any(gg == gd.days(2,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d2time/60);
     xlim([0.5 22]);   
axis.XTickLabel=[];
end

if length(days)>2
if any(gg == gd.days(3,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d3time/60);
     xlim([0.5 22]);   
axis.XTickLabel=[];
end
end
end

if length(days)<2
    xlim([0.5 22]);
end

 ylim([0 351]);
axis.XAxis.FontSize=12;
axis.LineWidth=1;
end


if makeblack==0
filename=strcat(experiment,'-','Kymograph of DFC-Time-Nuc');
 saveas(fig, fullfile(resultsDir, filename), 'png');
   saveas(fig, fullfile(resultsDir, filename), 'svg'); 
 savefig(fig,fullfile(resultsDir, [filename '.fig']),'compact');
else
    fig.InvertHardcopy = 'off';
   filename=strcat(experiment,'-','Kymograph of DFC-Time-Nuc in black');
 saveas(fig,fullfile(resultsDir, filename), 'png');
 savefig(fig,fullfile(resultsDir, [filename '.fig']),'compact');
end

%% now the other graph... ZvT

gd.fZ_DFCMI={};

for gg=1:length(groups)
  
    for tt = 1:length(gd.DFCMI{gg})
      for zz = 1:length(gd.DFCMI{gg}{tt})
   
            gd.fZ_DFCMI{gg}{tt}(zz,:)=mean(gd.DFCMI{gg}{tt}{zz}(:,2));

      end
    end

end


Zgrid={};
for gg=1:length(groups)
Zgrid{gg}=cell2mat(gd.fZ_DFCMI{gg});
end


%if there's a Z-drift, you can correct it here:
shifted = gd.days(1,:);
Zdrift = 1; %is the first plane drifted?
for gg=1:length(shifted)
Zgrid{shifted(gg)}(Zdrift,:)=[]; 
end


for gg=1:length(groups)
    stupidsize=size(Zgrid{gg});
    
bigZgrid{gg}=bgvalue*ones(actualmaxZ,nT(gg));
bigZgrid{gg}(1:stupidsize(1),1:stupidsize(2))=Zgrid{gg};
end

smoothbigZgrid={};
for gg=1:length(groups)
smoothbigZgrid{gg}=movmean(bigZgrid{gg},2);
smoothbigZgrid{gg}=movmean(smoothbigZgrid{gg},2,2);

end


%close all;

tiles=reshape(1:length(groups),length(days),length(conditions))'; % it has to be like this for it to be in the order of the graph tiles

makeblack=0;

if makeblack==0
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4)/1.1)]);
lettercolor='k'; oppositeLC='w'; oppositetriplet=[1 1 1]; colortriplet = [0 0 0];

cp=colormap(hot);
colormap(cp);

else
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4)/1.1)],'color','black');
lettercolor='w'; oppositeLC='k'; oppositetriplet=[0 0 0]; colortriplet = [1 1 1];

cp=colormap(hot);
cp(1,:)=[0,0,0]; 

colormap(cp);

end
%for black
t=tiledlayout(length(conditions),length(days),'TileSpacing','tight','Padding','loose');

sgtitle('Kymographs','FontSize',20,'FontWeight','bold','color',lettercolor);
ylabel(t,'Z (µm)','FontWeight','bold','color',lettercolor);
xlabel(t,'Hours','FontWeight','bold','FontSize',15,'color',lettercolor);

for gg=1:length(groups)
    nexttile

pplot(gg)=imagesc(smoothbigZgrid{gg});
pplot(gg).YData = pplot(gg).YData.*3;

caxis([round(min(nuc.mins)) round(max(nuc.maxes)/100)*100]);

if makeblack==1
axis=gca
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
set(gca,'color','black')
else
    axis=gca
    axis.XColor='k';
    axis.YColor='k';
    set(gca,'YDir','normal');
end

%if gg==2 | gg==3 | gg==5 | gg==6 | gg==8 | gg==9
if any(gg==posford23)
set(gca,'Ytick',[],'YColor','none','box','off','YTickLabel',[],'layer','top');
pplot(gg).Parent.YColor=oppositetriplet;
 %set(gca,'Visible','off');
end

if gg==tiles(posforbar,end)
cb(gg)=colorbar;
ylabel(cb(gg),'Nodal Mean Intensity (A.U.)','FontWeight','bold','color',lettercolor);
cb(gg).Color = colortriplet;
cb(gg).Box='off'
%cb(gg).Ticks = [600 850 1100];
%cb(gg).TickLabels = {'0','0.5','1'};
end

if any(gg==posforconditions)
    tempidx=find(gg==posforconditions);
    ylabel(conditions{tempidx},'FontSize',15);
end

if any(gg==posford1)
pplot(gg).XData=pplot(gg).XData*(d1time/60);
xlim([0.5 22]);
axis.YAxis.FontSize=15;
axis.XTickLabel = [];
set(gca,'box','off','Ycolor',lettercolor,'layer','top');
axis.YTick = [0:10:60];
end

if any(gg == gd.days(2,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d2time/60);
     xlim([0.5 22]);
     axis.XTickLabel = [];
end

if length(days)>2
if any(gg == gd.days(3,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d3time/60);
     xlim([.5 21]);
     axis.XTickLabel = [];
     
end
end


  ylim([2 60]);
 axis.XAxis.FontSize=12;
end


if makeblack==0
filename=strcat(experiment,'-','Kymograph of Z-Time-Nuc');
 saveas(fig, fullfile(resultsDir, filename), 'png');
  saveas(fig, fullfile(resultsDir, filename), 'svg'); 
 savefig(fig,fullfile(resultsDir, [filename '.fig']),'compact');
else
    fig.InvertHardcopy = 'off';
   filename=strcat(experiment,'-','Kymograph of Z-Time-Nuc in black');
 saveas(fig,fullfile(resultsDir, filename), 'png');
 savefig(fig,fullfile(resultsDir, [filename '.fig']),'compact');
end



