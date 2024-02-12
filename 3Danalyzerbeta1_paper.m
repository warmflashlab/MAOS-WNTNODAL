%% 3D analyzer for beta cat movies

clc; close all; clear; clear memory;


dataDir = '/Volumes/'; %this is the path where the images are
h5Dir= '/Volumes';
EdataDir = '/Volumes//'; %where do you want the results


addpath(genpath('~/Documents/GitHub/stemcells')); 
pathtofile =[dataDir];
experiment = '-';

nucfiles=dir(fullfile(pathtofile,'*.csv'));

memfiles=nucfiles; %in case you do not need membrane quantifications

 conditions={'1','2','3','4','5','6'};

objfactor=0.625; % objective magnitifaction factor micrometer to pixels for 20x SD

screensize=get(0,'ScreenSize');
spacing=2; %how much spacing between each Z


% THIS IS SUPER IMPORTANT !!!!!!

loaddata = 0;

% THIS IS SUPER IMPORTANT !!!!!!


if loaddata<1
h5files=dir(fullfile(h5Dir,'*.h5*'));
end



nfiles=length(nucfiles);
mem.nfiles=length(memfiles);

if nfiles~=mem.nfiles
    disp('nuclear tables do not match the number of membrane tables');
end

days={'day 1','day 2','day 3'};

resultsDir = fullfile(EdataDir,experiment);
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

graphsRD=[resultsDir '/graphs/']; 

if ~exist(graphsRD,'dir')
    mkdir(graphsRD);
end

if loaddata<1
for ii=1:nfiles
idxh5=strfind(h5files(ii).name,'_');
idxnuc=strfind(nucfiles(ii).name,'_');

compare=strcmp(h5files(ii).name(1:idxh5(1)),nucfiles(ii).name(1:idxnuc(1)));

if compare <1
    error('ERROR : a h5 file is not matching its table : ERROR');

end
end
end

disp('all files have h5 and table correspondant, you can proceed');
%%  read the tables and extract data

nuctable={};
memtable={};
h5s=[];

chantoread=1; %which channel do you want to read? this is the default
 
    
for ii=1:nfiles
   nuctable{ii}=readtable([nucfiles(ii).folder '/' nucfiles(ii).name]); 
   nT(ii)=max(table2array(nuctable{ii}(:,2)))+1; %this will look on the timestep specifically because of the column number
   memtable{ii}=readtable([memfiles(ii).folder '/' memfiles(ii).name]); 
   
   if loaddata<1
   h5s(ii,:)=size(h5read([h5files(ii).folder '/' h5files(ii).name], '/exported_data'));
   end
end

%h5 files sizes is [x y c z t]

if loaddata<1
save([dataDir 'h5s.mat'],'h5s');
else
    load([dataDir 'h5s.mat']);
end

%% arrange extracted data from tables
% 
 %ilabel={}; xpos={}; ypos={}; zpos={}; linid={}; trackid={}; meanI={}; Alljuntits={}; Alldata={};
 nuc.ilabel={}; nuc.xpos={}; nuc.ypos={}; nuc.zpos={}; nuc.meanI={}; nuc.class={}; nuc.Alljuntits={}; nuc.Alldata={}; nuc.frame={}; nuc.area={};
 mem.ilabel={}; mem.xpos={}; mem.ypos={}; mem.zpos={}; mem.meanI={}; mem.class={}; mem.Alljuntits={}; mem.Alldata={}; mem.frame={}; mem.area={};
% 

for ii=1:nfiles
   
for jj=1:nT(ii)
nuc.ilabel{ii}{jj}=nuctable{ii}.labelimage_oid(find(nuctable{ii}.timestep==jj-1)); mem.ilabel{ii}{jj}=memtable{ii}.object_id(find(memtable{ii}.timestep==jj-1));

nuc.xpos{ii}{jj}=nuctable{ii}.ObjectCenter_0(find(nuctable{ii}.timestep==jj-1)); mem.xpos{ii}{jj}=memtable{ii}.ObjectCenter_0(find(memtable{ii}.timestep==jj-1));
nuc.ypos{ii}{jj}=nuctable{ii}.ObjectCenter_1(find(nuctable{ii}.timestep==jj-1)); mem.ypos{ii}{jj}=memtable{ii}.ObjectCenter_1(find(memtable{ii}.timestep==jj-1));
nuc.zpos{ii}{jj}=nuctable{ii}.ObjectCenter_2(find(nuctable{ii}.timestep==jj-1)); mem.zpos{ii}{jj}=memtable{ii}.ObjectCenter_2(find(memtable{ii}.timestep==jj-1));

nuc.class{ii}{jj}=nuctable{ii}.PredictedClass(find(nuctable{ii}.timestep==jj-1)); mem.Pclass{ii}{jj}=memtable{ii}.PredictedClass(find(memtable{ii}.timestep==jj-1));
nuc.frame{ii}{jj}=nuctable{ii}.timestep(find(nuctable{ii}.timestep==jj-1)); mem.frame{ii}{jj}=memtable{ii}.timestep(find(memtable{ii}.timestep==jj-1));

nuc.meanI{ii}{jj}=nuctable{ii}.MeanIntensity_0(find(nuctable{ii}.timestep==jj-1)); %0 is the first channel, 1 is the second
mem.meanI{ii}{jj}=memtable{ii}.MeanIntensity_0(find(memtable{ii}.timestep==jj-1)); %of the first channel

nuc.dfc{ii}{jj}= sqrt((nuc.xpos{ii}{jj}-h5s(ii,1)/2).^2 + (nuc.ypos{ii}{jj}-h5s(ii,2)/2).^2).*objfactor;
mem.dfc{ii}{jj}= sqrt((mem.xpos{ii}{jj}-h5s(ii,1)/2).^2 + (mem.ypos{ii}{jj}-h5s(ii,2)/2).^2).*objfactor;

nuc.area{ii}{jj}=nuctable{ii}.ObjectArea(find(nuctable{ii}.timestep==jj-1));
mem.area{ii}{jj}=memtable{ii}.ObjectArea(find(memtable{ii}.timestep==jj-1));

nuc.Alljuntits{ii}{jj}=[nuc.frame{ii}{jj}+1 nuc.ilabel{ii}{jj} repmat(ii,length(nuc.ilabel{ii}{jj}),1) nuc.area{ii}{jj} nuc.xpos{ii}{jj} nuc.ypos{ii}{jj} nuc.dfc{ii}{jj} nuc.zpos{ii}{jj} nuc.meanI{ii}{jj} ];
mem.Alljuntits{ii}{jj}=[mem.frame{ii}{jj}+1 mem.ilabel{ii}{jj} repmat(ii,length(mem.ilabel{ii}{jj}),1) mem.area{ii}{jj} mem.xpos{ii}{jj} mem.ypos{ii}{jj} mem.dfc{ii}{jj} mem.zpos{ii}{jj} mem.meanI{ii}{jj} ];

%that +1 is so the time does not start with zero


[ii jj]
end
end


description='frame(1) - ilabel(2) - file number(3) - area(4) - xpos(5) - ypos(6) - distance from center(7)- zpos(8) - mean intensity(9) ';

unfiltered.nuc.Alljuntits=nuc.Alljuntits;
unfiltered.mem.Alljuntits=mem.Alljuntits;

%% Filters and re-grouping
 
% % AREA FILTER
% % DFC FILTER
% % Z FILTER
% 
 filterarea=500; %any area below this will be eliminated;
 filterdfc=450; %any distance from the center greater than this will be deleted
 filterzc=1; %any z-centroid below this will be eliminated
 spacing = 2;
 repetitions = 1; %if there are less than this numbers of cells will be ignored

% % AREA FILTER
% % DFC FILTER
% % Z FILTER
% 
% 

for ii=1:nfiles
    for jj=1:nT(ii)

       nuc.Alljuntits{ii}{jj}(nuc.Alljuntits{ii}{jj}(:,4)<filterarea,:)=[];
       nuc.Alljuntits{ii}{jj}(nuc.Alljuntits{ii}{jj}(:,end-2)>filterdfc,:)=[];
       nuc.Alljuntits{ii}{jj}(nuc.Alljuntits{ii}{jj}(:,end-1)<filterzc,:)=[];
       
       mem.Alljuntits{ii}{jj}(mem.Alljuntits{ii}{jj}(:,4)<filterarea,:)=[];
      % mem.Alljuntits{ii}{jj}(mem.Alljuntits{ii}{jj}(:,4)>filterdfc,:)=[];
        
    end
end

% start arranging data
nuc.meanAlljuntits={}; nuc.stdAlljuntits={}; mem.meanAlljuntits={}; mem.stdAlljuntits={};
for ii=1:nfiles
    for jj=1:nT(ii)

nuc.meanAlljuntits{ii}{jj}=mean(nuc.Alljuntits{ii}{jj}(:,end)); %last row because last row is mean intensity
nuc.stdAlljuntits{ii}{jj}=std(nuc.Alljuntits{ii}{jj}(:,end))./2;

mem.meanAlljuntits{ii}{jj}=mean(mem.Alljuntits{ii}{jj}(:,end));
mem.stdAlljuntits{ii}{jj}=std(mem.Alljuntits{ii}{jj}(:,end))./2;

if height(mem.stdAlljuntits{ii}{jj})<1
mem.stdAlljuntits{ii}{jj}=0;
end

    end
end

nuc.meanAlldata={}; nuc.stdAlldata={}; mem.meanAlldata={}; mem.stdAlldata={};
for ii = 1:nfiles
for jj = 1:nT(ii)

if jj==1
nuc.Alldata{ii}=nuc.Alljuntits{ii}{jj};
nuc.meanAlldata{ii}=[jj nuc.meanAlljuntits{ii}{jj}];
nuc.stdAlldata{ii}=[jj nuc.stdAlljuntits{ii}{jj}];

mem.Alldata{ii}=mem.Alljuntits{ii}{jj};
mem.meanAlldata{ii}=[jj mem.meanAlljuntits{ii}{jj}];
mem.stdAlldata{ii}=[jj mem.stdAlljuntits{ii}{jj}];


else
    nuc.Alldata{ii}=cat(1,nuc.Alldata{ii},nuc.Alljuntits{ii}{jj});
    nuc.meanAlldata{ii}=cat(1,nuc.meanAlldata{ii},[jj nuc.meanAlljuntits{ii}{jj}]);
    nuc.stdAlldata{ii}=cat(1,nuc.stdAlldata{ii},[jj nuc.stdAlljuntits{ii}{jj}]);
    
      mem.Alldata{ii}=cat(1,mem.Alldata{ii},mem.Alljuntits{ii}{jj});
    mem.meanAlldata{ii}=cat(1,mem.meanAlldata{ii},[jj mem.meanAlljuntits{ii}{jj}]);
    mem.stdAlldata{ii}=cat(1,mem.stdAlldata{ii},[jj mem.stdAlljuntits{ii}{jj}]);
end


end
end

divisor= nfiles/length(conditions);

%curate all the scatter data

curate= 99; % what percentage of all data want to extract the maximum from?
nuc.checkformax={}; nuc.maxes=[]; nuc.mins=[];
mem.checkformax={}; mem.maxes=[]; mem.mins=[];

%for mean intensity
for ii=1:nfiles
nuc.checkformax{ii}=sort(nuc.Alldata{ii}(:,end));
mem.checkformax{ii}=sort(mem.Alldata{ii}(:,end));

nuc.maxes(ii)=nuc.checkformax{ii}(round(length(nuc.checkformax{ii})*curate/100));
mem.maxes(ii)=mem.checkformax{ii}(round(length(mem.checkformax{ii})*curate/100));
nuc.mins(ii)=min(nuc.Alldata{ii}(:,end));
mem.mins(ii)=min(mem.Alldata{ii}(:,end));

nuc.Alldata{ii}(nuc.Alldata{ii}(:,end)>max(nuc.maxes),:)=[];

end


%% group files by conditions

gd.description = {'day1-ctrl','day2-ctrl', 'day3-ctrl', 'day1-10x', 'day2-10x', 'day3-10x','day1-SB10x', 'day2-SB10x', 'day3-SB10x','day1-1x-2dSB10x', 'day2-1x-2dSB10x', 'day3-1x-2dSB10x',...
    'day1-CHIR3', 'day2-CHIR3', 'day3-CHIR3','day1-CHIR10', 'day2-CHIR10', 'day3-CHIR10'};

gd.days = [1,4,7,10,13,16;2,5,8,11,14,17;3,6,9,12,15,18]; %separate by days
gd.conditions = [1,2,3;4,5,6;7,8,9;10,11,12;13,14,15;16,17,18]; %separate by conditions

if numel(gd.days) ~= length(gd.description)
    error('the agroupation of days doesnt match the description');
else
     cprintf('blue','days match description');
end

if numel(gd.conditions) ~= length(gd.description)
    error('the agroupation of conditions doesnt match the description');
else
    disp(' ');
    cprintf('blue','conditions match description');
end


%indicate which are the files into each group. You need to look at the
%order that comes out of files(ii).name

%control          10x               SB            %2dSB    CHIR  3        CHIR10
groups{1}= [];    groups{4}= [];           groups{7}=[];     groups{10}= [];          groups{13}= [];   groups{16}= [];
groups{2}= [];    groups{5}= [];        groups{8}=[];     groups{11}= [];       groups{14}= [];      groups{17}= [];
groups{3}= [];    groups{6}= [];        groups{9}=[];     groups{12}= [];       groups{15}= [];      groups{18}= [];


check=length([groups{:}]);

if check==nfiles
    cprintf('blue','all files are within a group');
else
    cprintf('red',' WARNING : NOT ALL FILES ARE IN GROUPS !!! ');
end

%% now re-arrange & group data + normalizations

clear norm;


for gg = 1:length(groups)
for ssgg=1:length(groups{gg})    
    
    if ssgg==1
gd.nuc.Alldata{gg}=nuc.Alldata{groups{gg}(ssgg)};             gd.mem.Alldata{gg}=mem.Alldata{groups{gg}(ssgg)};  
gd.nuc.meanAlldata{gg}= nuc.meanAlldata{groups{gg}(ssgg)};    gd.mem.meanAlldata{gg}= mem.meanAlldata{groups{gg}(ssgg)};
    else
gd.nuc.Alldata{gg}=cat(1,gd.nuc.Alldata{gg},nuc.Alldata{groups{gg}(ssgg)});                   gd.mem.Alldata{gg}=cat(1,gd.mem.Alldata{gg},mem.Alldata{groups{gg}(ssgg)});  
gd.nuc.meanAlldata{gg}=cat(1,gd.nuc.meanAlldata{gg},nuc.meanAlldata{groups{gg}(ssgg)});       gd.mem.meanAlldata{gg}=cat(1,gd.mem.meanAlldata{gg},mem.meanAlldata{groups{gg}(ssgg)});
   end
  
end
end

clear temptps; clear tempidxs;

for gg=1:length(groups)
temptps=unique(gd.nuc.meanAlldata{gg}(:,1)); 
for tt=1:length(temptps)
tempidxs{gg}{tt}=find(tt==gd.nuc.meanAlldata{gg}(:,1)); 
mem.tempidxs{gg}{tt}=find(tt==gd.mem.meanAlldata{gg}(:,1));
end
end

gd.nuc.tempmAd = [];  gd.nuc.tempstdAd=[]; gd.mem.tempmAd = [];  gd.mem.tempstdAd=[];

for gg=1:length(groups)
    temptps=unique(gd.nuc.Alldata{gg}(:,1));
for tt=1:length(temptps)

gd.nuc.tempmAd{gg}(tt,:)=[tt mean(gd.nuc.meanAlldata{gg}(tempidxs{gg}{tt},2))];
gd.mem.tempmAd{gg}(tt,:)=[tt mean(gd.mem.meanAlldata{gg}(mem.tempidxs{gg}{tt},2))];

end
end

clear gd.nuc.meanAlldata; clear gd.mem.meanAlldata;
clear gd.nuc.stdAlldata; clear gd.mem.stdAlldata;

gd.nuc.meanAlldata=gd.nuc.tempmAd;  gd.mem.meanAlldata=gd.mem.tempmAd;
norm.gd.nuc.meanAlldata=[]; norm.gd.mem.meanAlldata=[];


controls=[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]; %indicate controls for each day for each condition, this has to be manual

clear norm;
for gg=1:length(groups)

    
if length(gd.nuc.meanAlldata{gg}(:,end))<length(gd.nuc.meanAlldata{controls(gg)}(:,end));
   temp.diff =  length(gd.nuc.meanAlldata{controls(gg)}(:,end))-length(gd.nuc.meanAlldata{gg}(:,end));
   temp.length = length(gd.nuc.meanAlldata{gg}(:,end));
   temp.value=gd.nuc.meanAlldata{gg}(end,end);
   gd.nuc.meanAlldata{gg}(temp.length:temp.length+temp.diff,end)=temp.value;  
   gd.mem.meanAlldata{gg}(temp.length:temp.length+temp.diff,end)=temp.value;
end
    
   
norm.gd.nuc.meanAlldata{gg}=[gd.nuc.meanAlldata{gg}(:,1) gd.nuc.meanAlldata{gg}(:,end)./gd.nuc.meanAlldata{controls(gg)}(:,end)];
norm.gd.mem.meanAlldata{gg}=[gd.mem.meanAlldata{gg}(:,1) gd.mem.meanAlldata{gg}(:,end)./gd.mem.meanAlldata{controls(gg)}(:,end)];

end
% 

clear temptps; clear tempidxs; clear mem.tempidxs;

for gg=1:length(groups)
temptps=unique(gd.nuc.Alldata{gg}(:,1)); %times apply to both mem and nuc
for tt=1:length(temptps)
tempidxs{gg}{tt}=find(tt==gd.nuc.Alldata{gg}(:,1));
mem.tempidxs{gg}{tt}=find(tt==gd.mem.Alldata{gg}(:,1));
end
end

norm.gd.nuc.tempAd={}; norm.gd.mem.tempAd={}; 
gd.nuc.tempstdAd={}; gd.mem.tempstdAd={};
norm.gd.nuc.tempstdAd={}; norm.gd.mem.tempstdAd={};

for gg=1:length(groups)
    temptps=unique(gd.nuc.Alldata{gg}(:,1));
for tt=1:length(temptps)

gd.nuc.tempstdAd{gg}{tt}=std(gd.nuc.Alldata{gg}(tempidxs{gg}{tt},end))./2; % not norm
norm.gd.nuc.tempAd{gg}{tt}=gd.nuc.Alldata{gg}(tempidxs{gg}{tt},end)./gd.nuc.meanAlldata{gg}(tt,2); %norm
norm.gd.nuc.tempstdAd{gg}{tt}=std(norm.gd.nuc.tempAd{gg}{tt})./2; % norm

gd.mem.tempstdAd{gg}{tt}=std(gd.mem.Alldata{gg}(mem.tempidxs{gg}{tt},end))./2; %not norm
norm.gd.mem.tempAd{gg}{tt}=gd.mem.Alldata{gg}(mem.tempidxs{gg}{tt},end)./gd.mem.meanAlldata{gg}(tt,2); %norm
norm.gd.mem.tempstdAd{gg}{tt}=std(norm.gd.mem.tempAd{gg}{tt})./2;

end
end


gd.nuc.stdAlldata={}; norm.gd.nuc.Alldata={}; norm.gd.nuc.stdAlldata={};
gd.mem.stdAlldata={}; norm.gd.mem.Alldata={}; norm.gd.mem.stdAlldata={};

for gg=1:length(groups)
    temptps=unique(gd.nuc.Alldata{gg}(:,1));
for tt=1:length(temptps)
    
    if tt==1

    gd.nuc.stdAlldata{gg}=gd.nuc.tempstdAd{gg}{tt};
    norm.gd.nuc.Alldata{gg}=norm.gd.nuc.tempAd{gg}{tt};
    norm.gd.nuc.stdAlldata{gg}=norm.gd.nuc.tempstdAd{gg}{tt};
    
    gd.mem.stdAlldata{gg}=gd.mem.tempstdAd{gg}{tt};
    norm.gd.mem.Alldata{gg}=norm.gd.mem.tempAd{gg}{tt};
    norm.gd.mem.stdAlldata{gg}=norm.gd.mem.tempstdAd{gg}{tt};
    
    else
   
    gd.nuc.stdAlldata{gg}=cat(1,gd.nuc.stdAlldata{gg},gd.nuc.tempstdAd{gg}{tt});
    norm.gd.nuc.Alldata{gg}=cat(1,norm.gd.nuc.Alldata{gg},norm.gd.nuc.tempAd{gg}{tt});
    norm.gd.nuc.stdAlldata{gg}=cat(1,norm.gd.nuc.stdAlldata{gg},norm.gd.nuc.tempstdAd{gg}{tt});
     
    gd.mem.stdAlldata{gg}=cat(1,gd.mem.stdAlldata{gg},gd.mem.tempstdAd{gg}{tt});
    norm.gd.mem.Alldata{gg}=cat(1,norm.gd.mem.Alldata{gg},norm.gd.mem.tempAd{gg}{tt});
    norm.gd.mem.stdAlldata{gg}=cat(1,norm.gd.mem.stdAlldata{gg},norm.gd.mem.tempstdAd{gg}{tt});
    
    end
    
end
end


for gg=1:length(groups)
   
    maxcolors(gg,:)=[min(gd.nuc.Alldata{gg}(:,end)) max(nuc.maxes)]; %this is for the color, not data
  
   norm.gd.nuc.Alldata{gg}=horzcat(gd.nuc.Alldata{gg},norm.gd.nuc.Alldata{gg});
   norm.gd.mem.Alldata{gg}=horzcat(gd.mem.Alldata{gg},norm.gd.mem.Alldata{gg});
   
end



%% multiply timepoints for real time and Z for it's spacing value

%put the time in minutes of each day

 d1time= 31; % put this in min
 d2time = 42; % put this in min
 d3time = 48;% put this in min


for gg = gd.days(1,:);
%this is for day 1
norm.gd.nuc.Alldata{gg}(:,1)=(norm.gd.nuc.Alldata{gg}(:,1)).*(d1time/60);
gd.nuc.Alldata{gg}(:,1)=(gd.nuc.Alldata{gg}(:,1)).*(d1time/60);
gd.nuc.meanAlldata{gg}(:,1)=(gd.nuc.meanAlldata{gg}(:,1)).*(d1time/60);
norm.gd.nuc.meanAlldata{gg}(:,1)=(norm.gd.nuc.meanAlldata{gg}(:,1)).*(d1time/60);

gg


norm.gd.mem.Alldata{gg}(:,1)=(norm.gd.mem.Alldata{gg}(:,1)).*(d1time/60);
gd.mem.Alldata{gg}(:,1)=(gd.mem.Alldata{gg}(:,1)).*(d1time/60);
gd.mem.meanAlldata{gg}(:,1)=(gd.mem.meanAlldata{gg}(:,1)).*(d1time/60);
norm.gd.mem.meanAlldata{gg}(:,1)=(norm.gd.mem.meanAlldata{gg}(:,1)).*(d1time/60);

end

for gg=gd.days(2,:)
%for day 2


norm.gd.nuc.Alldata{gg}(:,1)=(norm.gd.nuc.Alldata{gg}(:,1)).*(d2time/60);
gd.nuc.Alldata{gg}(:,1)=(gd.nuc.Alldata{gg}(:,1)).*(d2time/60);
gd.nuc.meanAlldata{gg}(:,1)=(gd.nuc.meanAlldata{gg}(:,1)).*(d2time/60);
norm.gd.nuc.meanAlldata{gg}(:,1)=(norm.gd.nuc.meanAlldata{gg}(:,1)).*(d2time/60);

gg


norm.gd.mem.Alldata{gg}(:,1)=(norm.gd.mem.Alldata{gg}(:,1)).*(d2time/60);
gd.mem.Alldata{gg}(:,1)=(gd.mem.Alldata{gg}(:,1)).*(d2time/60);
gd.mem.meanAlldata{gg}(:,1)=(gd.mem.meanAlldata{gg}(:,1)).*(d2time/60);
norm.gd.mem.meanAlldata{gg}(:,1)=(norm.gd.mem.meanAlldata{gg}(:,1)).*(d2time/60);
end

%for day 3
for gg = gd.days(3,:)

norm.gd.nuc.Alldata{gg}(:,1)=(norm.gd.nuc.Alldata{gg}(:,1)).*(d3time/60);
gd.nuc.Alldata{gg}(:,1)=(gd.nuc.Alldata{gg}(:,1)).*(d3time/60);
gd.nuc.meanAlldata{gg}(:,1)=(gd.nuc.meanAlldata{gg}(:,1)).*(d3time/60);
norm.gd.nuc.meanAlldata{gg}(:,1)=(norm.gd.nuc.meanAlldata{gg}(:,1)).*(d3time/60);

gg


norm.gd.mem.Alldata{gg}(:,1)=(norm.gd.mem.Alldata{gg}(:,1)).*(d3time/60);
gd.mem.Alldata{gg}(:,1)=(gd.mem.Alldata{gg}(:,1)).*(d3time/60);
gd.mem.meanAlldata{gg}(:,1)=(gd.mem.meanAlldata{gg}(:,1)).*(d3time/60);
norm.gd.mem.meanAlldata{gg}(:,1)=(norm.gd.mem.meanAlldata{gg}(:,1)).*(d3time/60);
end


%here's the correction of Z in µm so it's equal for different acquistions
% %of other experiments

for gg=1:length(groups)


norm.gd.nuc.Alldata{gg}(:,end-2)=(norm.gd.nuc.Alldata{gg}(:,end-2).*spacing);
norm.gd.mem.Alldata{gg}(:,end-2)=(norm.gd.mem.Alldata{gg}(:,end-2).*spacing);

gd.nuc.Alldata{gg}(:,end-1)=(gd.nuc.Alldata{gg}(:,end-1).*spacing);
gd.mem.Alldata{gg}(:,end-1)=(gd.mem.Alldata{gg}(:,end-1).*spacing);


end

%filter area
 filterarea = 300;
 
 for gg=1:length(groups)
     temp.idx=gd.nuc.Alldata{gg}(:,4)<filterarea;
    gd.nuc.Alldata{gg}(temp.idx,:)=[]; 
 end

 %filter MI
 filterMI= 400; %this is in case there are black pixels in the movies 
 
 for gg=1:length(groups)
     temp.idx=gd.nuc.Alldata{gg}(:,end)<filterMI;
    gd.nuc.Alldata{gg}(temp.idx,:)=[]; 
 end


%% Calculations for color grids

clear tempidx; clear tempidxs; clear templist; clear temptps; clear temprz; clear temprzidx; clear gd.nuc.rzfk;
clear tempdfc; clear gd.nuc.rdfcfk; clear tempdfcidx;


repetitions= 3; %less than this will be eliminated


for gg=1:length(groups)

temptps=unique(gd.nuc.Alldata{gg}(:,1)); %going to change this for now 
    for tt=1:length(temptps)
tempidxs{gg}{tt}=gd.nuc.Alldata{gg}(:,1)==temptps(tt);
tempbytime{gg}{tt}=gd.nuc.Alldata{gg}(tempidxs{gg}{tt},:);
tempbytime{gg}{tt}(:,end-1)=round(tempbytime{gg}{tt}(:,end-1)); %round the z values to all
tempbytime{gg}{tt}(:,end-2)=round(tempbytime{gg}{tt}(:,end-2)); %round the dfc values

rzvalues{gg}{tt}= unique(tempbytime{gg}{tt}(:,end-1)); %rounded z values to first decimal, original: end-1
dfcvalues{gg}{tt} = unique(tempbytime{gg}{tt}(:,end-2)); %original : end-2

temp.counts=histc(tempbytime{gg}{tt}(:,end-1),rzvalues{gg}{tt});
rzvalues{gg}{tt}=rzvalues{gg}{tt}(temp.counts>repetitions);

        for rzz=1:length(rzvalues{gg}{tt})
temprzidx= tempbytime{gg}{tt}(:,end-1)==rzvalues{gg}{tt}(rzz);
temprz{gg}{tt}(rzz,:)= [tt temptps(tt) rzvalues{gg}{tt}(rzz) mean(tempbytime{gg}{tt}(temprzidx,end))+std(tempbytime{gg}{tt}(temprzidx,end))./2]; %timepoint - timepoint in hours - rounded z - mean nuc bcat
temprzstd{gg}{tt}(rzz,:)= [tt temptps(tt) rzvalues{gg}{tt}(rzz) std(tempbytime{gg}{tt}(temprzidx,end))./2]; %timepoint - timepoint in hours - rounded z - mean nuc bcat
        end

        for dfcc=1:length(dfcvalues{gg}{tt})
tempdfcidx=tempbytime{gg}{tt}(:,end-2)==dfcvalues{gg}{tt}(dfcc);
tempdfc{gg}{tt}(dfcc,:) = [tt temptps(tt) dfcvalues{gg}{tt}(dfcc) maxPer(tempbytime{gg}{tt}(tempdfcidx,end),97)]; % timepoint - timepoint in hours - rounded DFC - max nuc bcat
tempdfcstd{gg}{tt}(dfcc,:) = [tt temptps(tt) dfcvalues{gg}{tt}(dfcc) std(tempbytime{gg}{tt}(tempdfcidx,end))./2]; % timepoint - timepoint in hours - rounded DFC - mean nuc bcat

        end
    end
gd.nuc.rzfk{gg}=cat(1,temprz{gg}{:}); %rounded zs for kymograph         %timepoint - timepoint in hours - rounded z - mean nuc bcat
gd.nuc.rdfcfk{gg}=cat(1,tempdfc{gg}{:}); %rounded dfc for kymograph     % timepoint - timepoint in hours - rounded DFC - mean nuc bcat

gd.nuc.rzfk_std{gg}=cat(1,temprzstd{gg}{:}); %rounded zs for kymograph         %timepoint - timepoint in hours - rounded z - mean nuc bcat
gd.nuc.rdfcfk_std{gg}=cat(1,tempdfcstd{gg}{:}); %rounded dfc for kymograph     % timepoint - timepoint in hours - rounded DFC - mean nuc bcat


end


clear tempidx; clear tempidxs; clear templist; clear temptps; clear temprz; clear temprzidx;
clear tempdfc; clear tempdfcidx; clear tempZs; clear tempbyZ; 

for gg=1:length(groups)
    

   temp.list=unique(gd.nuc.rzfk{gg}(:,3));
   temp.counts=histc(gd.nuc.Alldata{gg}(:,8),temp.list); %este está mal, tiene que contar los de gd.nuc.Alldata
   tempZs{gg}= temp.list(temp.counts>repetitions);
   
   
   for zz=1:length(tempZs{gg})
   tempidxs{gg}{zz}=round(gd.nuc.Alldata{gg}(:,end-1))==tempZs{gg}(zz);
   tempbyZ{gg}{zz}=gd.nuc.Alldata{gg}(tempidxs{gg}{zz},:);

   tempbyZ{gg}{zz}(:,end-1)=round(tempbyZ{gg}{zz}(:,end-1)); %round the z values to all
   tempbyZ{gg}{zz}(:,end-2)=round(tempbyZ{gg}{zz}(:,end-2)); %round the dfc values
   
   dfcvaluesinZ{gg}{zz}=unique(tempbyZ{gg}{zz}(:,end-2));
   
        for dfcc=1:length(dfcvaluesinZ{gg}{zz})
            tempdfcidx=tempbyZ{gg}{zz}(:,end-2)==dfcvaluesinZ{gg}{zz}(dfcc);
            tempdfc{gg}{zz}(dfcc,:) = [zz tempZs{gg}(zz) dfcvaluesinZ{gg}{zz}(dfcc) mean(tempbyZ{gg}{zz}(tempdfcidx,end))]; % z - z ?? - rounded DFC - mean nuc bcat
        end 
   
   end
   gd.nuc.rdfc_rzfk{gg}=cat(1,tempdfc{gg}{:}); %rounded zs for kymograph
end


for gg=1:length(groups)
   topZ(gg)=max(rmoutliers(gd.nuc.rzfk{gg}(:,3))); %esto está con remove outliers, 
end
topZ=round(max(topZ)/10)*10;

 %% Kymograph !!! for Z vs T vs Nuc Bcat levels !!!!!

deltaforZ=2; %if there's a space between two z points, then it will be a false signal

X = cat(1,gd.nuc.rzfk{:});
Y = cat(1,gd.nuc.rzfk{:});

X = unique(X(:,1));
Y = unique(Y(:,3));

Ycodex = [(1:length(Y))' Y];
Xcodex = [(1:length(X))' X];

temp.grid=[];

for gg=1:length(groups)
    temp.grid{gg}=zeros(length(X),length(Y));
for kk=1:length(gd.nuc.rzfk{gg})
    [temp.x, temp.y]=find(Ycodex(:,2)==gd.nuc.rzfk{gg}(kk,3));
temp.grid{gg}(gd.nuc.rzfk{gg}(kk,1),Ycodex(temp.x,1))=gd.nuc.rzfk{gg}(kk,end);


end
end

%to clear outliers
for gg=1:length(groups)
    nonsense=size(temp.grid{gg});
for kk=1:nonsense(1)
temp.idxs=find(temp.grid{gg}(kk,:)~=0);
counter=0;

for mm = 2:length(temp.idxs)
 temp.diff=temp.idxs(mm)-temp.idxs(mm-1);
if temp.diff >= deltaforZ
  
   temp.grid{gg}(kk,temp.idxs(mm):end)=0;
  
end


end
end
end

%close all;

 %% actual graph now
%smooth version!!!!!

smoothlvl=1; % how much smooth you like to test? 1 is no smooth

for gg=1:length(groups)
temp.Sgrid{gg}=movmean(temp.grid{gg},smoothlvl);
temp.Sgrid{gg}=movmean(temp.Sgrid{gg},smoothlvl,2);
end

tiles=reshape(1:length(groups),length(days),length(conditions))'; % it has to be like this for it to be in the order of the graph tiles

posforbar = round(length(conditions)/2); %position of the tile for the color bar
posford23 = tiles(:,2:3); %positions for all in day 2 and 3
posford23 = reshape(posford23,1,[]);
posforconditions = tiles(:,1)';
posford1 = tiles(:,1)';
clear cp;

makeblack=0

if makeblack==0
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4)/1)]);
lettercolor='k'; oppositeLC='w'; oppositetriplet=[1 1 1]; colortriplet = [0 0 0];
cp=colormap(parula);
cp(1,:)=[1,1,1];
colormap(cp);

else
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4)/1.2)],'color','black');
lettercolor='w'; oppositeLC='k'; oppositetriplet=[0 0 0]; colortriplet = [1 1 1];
cp=colormap(turbo);
cp(1,:)=[0,0,0]; 
cp(end,:)=[1,0,0];
colormap(cp);

end
%for black
t=tiledlayout(length(conditions),length(days),'TileSpacing','none','Padding','loose');
ylabel(t,'Z centroid (µm)','FontWeight','bold','color',lettercolor);
xlabel(t,'Hours','FontWeight','bold','FontSize',15,'color',lettercolor);

for gg=1:length(groups)
    nexttile
pplot(gg)=imagesc(temp.Sgrid{gg}');
pplot(gg).YData = Ycodex(:,2);

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

if any(gg==posford23)
set(gca,'Ytick',[],'YColor','none','box','off','YTickLabel',[],'layer','top');
pplot(gg).Parent.YColor=oppositetriplet;
 %set(gca,'Visible','off');
end

if gg==tiles(posforbar,3)
cb(gg)=colorbar;
ylabel(cb(gg),'Nuclear Bcat Mean Intensity (A.U.)','FontWeight','bold','color',lettercolor);
cb(gg).Color = colortriplet;
cb(gg).Box='off'

end

if any(gg==posforconditions)
    tempidx=find(gg==posforconditions);
    ylabel(conditions{tempidx},'FontSize',15);
end

if any(gg==posford1)
pplot(gg).XData=pplot(gg).XData*(d1time/60);
xlim([0.5 19]);
axis.YAxis.FontSize=20;
axis.XTickLabel = [];
set(gca,'box','off','Ycolor',lettercolor,'layer','top');
end

if any(gg == gd.days(2,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d2time/60);
     xlim([.5 22]);
     axis.XTickLabel = [];
end


if any(gg == gd.days(3,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d3time/60);
     xlim([0.5 23]);
     axis.XTickLabel = [];
end

  ylim([2 60]);
 axis.XAxis.FontSize=12;
end


if makeblack==0
filename=strcat(experiment,'-','Kymograph of Z-Time-Nuc');
 saveas(fig, fullfile(graphsRD, filename), 'png');
  saveas(fig, fullfile(graphsRD, filename), 'svg'); 
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
else
    fig.InvertHardcopy = 'off';
   filename=strcat(experiment,'-','Kymograph of Z-Time-Nuc in black');
 saveas(fig,fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
end

%% check the low Zs & histograms (optional)


%histograms
fig=figure('Position',[1,1,round(screensize(3)),round(screensize(4)/2.5)]);min_value = Inf;
max_value = -Inf;
max_count = -Inf;

% Find min, max values in the data and max count for Y limits
for gg=1:height(gd.days)
    for ssgg=1:length(gd.days)
        temp_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,9);
        min_value = min(min_value, min(temp_data));
        max_value = max(max_value, max(temp_data));

        [counts, ~] = histcounts(temp_data); % Compute histogram without plotting
        max_count = max(max_count, max(counts));
    end
end

xlimits = [min_value, max_value];
ylimits = [0, max_count];

for gg=1:height(gd.days)
    sgtitle('mean intensity');
    nexttile
    for ssgg=1:length(gd.days)
        histogram(gd.nuc.Alldata{gd.days(gg,ssgg)}(:,9),'DisplayStyle','stairs','LineWidth',3); hold on;
        xlim(xlimits);
        ylim(ylimits);
        legend(conditions)
    end
end
filename=strcat(experiment,'-','Histogram for MI');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
   saveas(fig,fullfile(graphsRD, filename), 'svg');


%this is for Zs
fig=figure('Position',[1,1,round(screensize(3)),round(screensize(4)/2.5)]);min_value = Inf;
max_value = -Inf;
max_count = -Inf;

% Find min, max values in the data and max count for Y limits
for gg=1:height(gd.days)
    for ssgg=1:length(gd.days)
        temp_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,8);
        min_value = min(min_value, min(temp_data));
        max_value = max(max_value, max(temp_data));

        [counts, ~] = histcounts(temp_data); % Compute histogram without plotting
        max_count = max(max_count, max(counts));
    end
end

xlimits = [min_value, max_value];
ylimits = [0, max_count];

for gg=1:height(gd.days)
    sgtitle('Z centroid');
    nexttile
    for ssgg=1:length(gd.days)
        histogram(gd.nuc.Alldata{gd.days(gg,ssgg)}(:,8),'DisplayStyle','stairs','LineWidth',3); hold on;
        xlim(xlimits);
        ylim(ylimits);
        legend(conditions)
    end
end
% 
filename=strcat(experiment,'-','Histogram for Z');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
   saveas(fig,fullfile(graphsRD, filename), 'svg');

%for area

fig=figure('Position',[1,1,round(screensize(3)),round(screensize(4)/2.5)]);min_value = Inf;
max_value = -Inf;
max_count = -Inf;

% Find min, max values in the data and max count for Y limits
for gg=1:height(gd.days)
    for ssgg=1:length(gd.days)
        temp_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,4);
        min_value = min(min_value, min(temp_data));
        max_value = max(max_value, max(temp_data));

        [counts, ~] = histcounts(temp_data); % Compute histogram without plotting
        max_count = max(max_count, max(counts));
    end
end

xlimits = [min_value, max_value];
ylimits = [0, max_count];

for gg=1:height(gd.days)
    sgtitle('area');
    nexttile
    for ssgg=1:length(gd.days)
        histogram(gd.nuc.Alldata{gd.days(gg,ssgg)}(:,4),'DisplayStyle','stairs','LineWidth',3); hold on;
        xlim(xlimits);
        ylim(ylimits);
        legend(conditions)
    end
end
% 
filename=strcat(experiment,'-','Histogram for area');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
   saveas(fig,fullfile(graphsRD, filename), 'svg');

   
%% 3D Scatter for checking area vs mean intensity (optional)

%Histograms per file in condition and day to check wtf is making the wrong thing
group2check=5;
figure;
for nn=1:length(groups{group2check});
    the3rdaxis=ones(length(nuc.Alldata{groups{group2check}(nn)}(:,end-1)),1)*nn;
     scatter3(nuc.Alldata{groups{group2check}(nn)}(:,end-1),nuc.Alldata{groups{group2check}(nn)}(:,end),the3rdaxis,20,'filled'); hold on; %just MI and Z
     xlabel('Z');
     ylabel('MI'); %
     zlabel('nfile');
end
   legend(string([(groups{group2check})]));
   yline(filterMI,'-r');
   
%% Normal scatters

fig=figure('Position',[1,1,round(screensize(3)),round(screensize(4)/2.5)]);min_value = Inf;
min_x_value = Inf;
max_x_value = -Inf;
min_y_value = Inf;
max_y_value = -Inf;

% Find 1st and 98th percentile values in the data
for gg=1:height(gd.days)
    for ssgg=1:length(gd.days)
        temp_x_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,4);
        temp_y_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,end);
        min_x_value = min(min_x_value, prctile(temp_x_data, 1));
        max_x_value = max(max_x_value, prctile(temp_x_data, 99));
        min_y_value = min(min_y_value, prctile(temp_y_data, 1));
        max_y_value = max(max_y_value, prctile(temp_y_data, 99));
    end
end

xlimits = [0, max_x_value];
ylimits = [filterMI, max_y_value];

for gg=1:height(gd.days)
    sgtitle('area');
    nexttile
    for ssgg=1:length(gd.days)
       scatter(gd.nuc.Alldata{gd.days(gg,ssgg)}(:,4),gd.nuc.Alldata{gd.days(gg,ssgg)}(:,end),5,'filled'); hold on;
       xlim(xlimits);
       ylim(ylimits);
       xlabel('Area');
       ylabel('MI');
       legend(conditions)
    end
end
% 
filename=strcat(experiment,'-','scatter for area and MI');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
   saveas(fig,fullfile(graphsRD, filename), 'svg');
   
   %for Zs and MI
   
fig=figure('Position',[1,1,round(screensize(3)),round(screensize(4)/2.5)]);min_value = Inf;
min_x_value = Inf;
max_x_value = -Inf;
min_y_value = Inf;
max_y_value = -Inf;

% Find 1st and 98th percentile values in the data
for gg=1:height(gd.days)
    for ssgg=1:length(gd.days)
        temp_x_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,end-1);
        temp_y_data = gd.nuc.Alldata{gd.days(gg,ssgg)}(:,end);
        min_x_value = min(min_x_value, prctile(temp_x_data, 1));
        max_x_value = max(max_x_value, prctile(temp_x_data, 99));
        min_y_value = min(min_y_value, prctile(temp_y_data, 1));
        max_y_value = max(max_y_value, prctile(temp_y_data, 99));
    end
end

xlimits = [0, 40];
ylimits = [filterMI, max_y_value];

for gg=1:height(gd.days)
    sgtitle('Z');
    nexttile
    for ssgg=1:length(gd.days)
       scatter(gd.nuc.Alldata{gd.days(gg,ssgg)}(:,end-1),gd.nuc.Alldata{gd.days(gg,ssgg)}(:,end),5,'filled'); hold on;
       xlim(xlimits);
       ylim(ylimits);
       xlabel('Z');
       ylabel('MI');
       legend(conditions)
    end
end
% 
filename=strcat(experiment,'-','scatter for Z and MI');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
   saveas(fig,fullfile(graphsRD, filename), 'svg');
   
   
%% Kymograph with DFC vs Hours vs Bcat ( )

X = cat(1,gd.nuc.rdfcfk{:});
Y = cat(1,gd.nuc.rdfcfk{:});

X = unique(X(:,1));
Y = unique(Y(:,3));

Ycodex = [(1:length(Y))' Y];
Xcodex = [(1:length(X))' X];

temp.grid=[]; temp.ngrid=[];

%
bins = 25; %has to be divisible by 5
newY=[1:bins:351];


for gg=1:length(groups)
    temp.grid{gg}=zeros(length(X),length(Y));
for kk=1:length(gd.nuc.rdfcfk{gg})

    [temp.x, temp.y]=find(Ycodex(:,2)==gd.nuc.rdfcfk{gg}(kk,3));

temp.grid{gg}(gd.nuc.rdfcfk{gg}(kk,1),Ycodex(temp.x,1))=gd.nuc.rdfcfk{gg}(kk,end);
%
for tt=1:height(temp.grid{gg})
for nn=1:length(newY)-1

temp.ngrid{gg}(tt,nn)=mean(nonzeros(temp.grid{gg}(tt,newY(nn):newY(nn+1))));

end
end

%

end
end


temp.savedgrid=temp.ngrid;

%% take out nans

%this is for taking out the complete NaN columns
for gg=1:length(groups)
gridsize=size(temp.ngrid{gg});

for rr=fliplr(1:gridsize(1))
if all(isnan(temp.ngrid{gg}(rr,:)))
    temp.ngrid{gg}(rr,:)=[];
end
end
end

%take out the weird empty nans

for gg=1:length(groups)
    [m, n]= size(temp.ngrid{gg});
    isnan_A = isnan(temp.ngrid{gg});

% Find the indices of NaN elements
[row_indices, col_indices] = find(isnan_A);

% Replace NaN elements with the average of their surrounding values
for i = 1:length(row_indices)
    row = row_indices(i);
    col = col_indices(i);
    
    % Determine neighboring indices
    row_neighbors = max(1, row-1):min(m, row+1);
    col_neighbors = max(1, col-1):min(n, col+1);
    
    % Extract neighboring values
    neighbors = temp.ngrid{gg}(row_neighbors, col_neighbors);
    
    % Calculate the average of non-NaN neighboring values
    mean_value = mean(neighbors(~isnan(neighbors)));
    
    % Replace the NaN value with the calculated average
   temp.ngrid{gg}(row, col) = mean_value;
end
   
end


%% for bridging from Nodal

bridge = 2; %how many x points you want it to be bridged between the firsts and lasts of days

%now put this for all grids
for dd=1:length(days)-1
for ii=1:length(gd.conditions)
    
    temp_size_L=size(temp.ngrid{gd.conditions(ii,dd)});
    temp_size_R=size(temp.ngrid{gd.conditions(ii,dd+1)});
    
    leftchunk=zeros(bridge,temp_size_L(2));
    rightchunk=zeros(bridge,temp_size_R(2));
    
    leftchunk=temp.ngrid{gd.conditions(ii,dd)}(end-(bridge-1):end,:);
    rightchunk=temp.ngrid{gd.conditions(ii,dd)+1}(1:bridge,:);
    
    temp_chunk=(leftchunk+rightchunk)/2;
    [gd.conditions(ii,dd), gd.conditions(ii,dd)+1]
    
    temp.ngrid{gd.conditions(ii,dd)}(end-(bridge-1):end,:) = temp_chunk;
    temp.ngrid{gd.conditions(ii,dd)+1}(1:bridge,:) = temp_chunk;
end
end

%% smooth graph 

smoothlevel = 5;
temp.Sgrid={};
for gg=1:length(groups)
temp.Sgrid{gg}=movmean(temp.ngrid{gg},smoothlevel);
temp.Sgrid{gg}=movmean(temp.Sgrid{gg},smoothlevel,2);

end

%% now finally make the graph
makeblack=0

if makeblack==0
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4)/1)]);
lettercolor='k'; oppositeLC='w'; oppositetriplet=[1 1 1]; colortriplet = [0 0 0];
colormap(parula);


else
fig=figure('Position',[1,1,round(screensize(3)/2.3),round(screensize(4)/1)],'color','black');
lettercolor='w'; oppositeLC='k'; oppositetriplet=[0 0 0]; colortriplet = [1 1 1];

colormap(cp);
end
%for black
t=tiledlayout(length(conditions),length(days),'TileSpacing','none','Padding','loose');
sgtitle('Kymograph Distance from Center','FontSize',20,'FontWeight','bold','color',lettercolor);
ylabel(t,'Distance from center (µm)','FontWeight','bold','color',lettercolor);
xlabel(t,'Hours','FontWeight','bold','FontSize',15,'color',lettercolor);

for gg=1:length(groups)
    nexttile

pplot(gg)=imagesc(temp.Sgrid{gg}','Interpolation','bilinear');
pplot(gg).YData = [0:bins:350]';

caxis([round(min(gd.nuc.meanAlldata{1}(:,end))) round(max(nuc.maxes)/100)*100]);

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

if gg==tiles(posforbar,3)
cb(gg)=colorbar;
ylabel(cb(gg),'Nuclear Bcat Mean Intensity (A.U.)','FontWeight','bold','color',lettercolor);
cb(gg).Color = colortriplet;
cb(gg).Box='off'
end

if any(gg==posforconditions)
    tempidx=find(gg==posforconditions);
    ylabel(conditions{tempidx},'FontSize',15);
end

if any(gg==posford1)
pplot(gg).XData=pplot(gg).XData*(d1time/60);
xlim([1 22]);
axis.YAxis.FontSize=18;
set(gca,'box','off','Ycolor',lettercolor,'layer','top');
%
axis.XTickLabel=[];
end

if any(gg == gd.days(2,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d2time/60);
     xlim([1 22]);   
axis.XTickLabel=[];
end

if any(gg == gd.days(3,:),'all')>0
    pplot(gg).XData=pplot(gg).XData*(d3time/60);
     xlim([0.5 22]);   
axis.XTickLabel=[];
end

 ylim([0 350]);
axis.XAxis.FontSize=12;
axis.LineWidth=1;
end


if makeblack==0
filename=strcat(experiment,'-','Kymograph of DFC-Time-Nuc');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
   saveas(fig,fullfile(graphsRD, filename), 'svg');
else
    fig.InvertHardcopy = 'off';
   filename=strcat(experiment,'-','Kymograph of DFC-Time-Nuc in black');
 saveas(fig,fullfile(graphsRD, filename), 'png');
  saveas(fig,fullfile(graphsRD, filename), 'svg');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
end

%% Kymograph with DFC vs Z vs Bcat 

X = cat(1,gd.nuc.rdfc_rzfk{:});
Y = cat(1,gd.nuc.rdfc_rzfk{:});

X = unique(X(:,3));
Y = unique(Y(:,1));

Ycodex = [(1:length(Y))' Y];
Xcodex = [(1:length(X))' X];

temp.grid=[]; temp.ngrid=[]; temp.idxs=[];

for gg=1:length(groups)
    temp.grid{gg}=zeros(length(X),length(Y));
for kk=1:length(gd.nuc.rdfc_rzfk{gg})
    [temp.x, temp.y]=find(Xcodex(:,2)==gd.nuc.rdfc_rzfk{gg}(kk,3));
if gd.nuc.rdfc_rzfk{gg}(kk,3)<1
    gd.nuc.rdfc_rzfk{gg}(kk,3)=1;
end
temp.grid{gg}(Xcodex(temp.x,1),gd.nuc.rdfc_rzfk{gg}(kk,1))=gd.nuc.rdfc_rzfk{gg}(kk,end);

end

end

%to clear outliers
for gg=1:length(groups)

    
for kk=1:length(temp.grid{gg})
temp.idxs=find(temp.grid{gg}(kk,:)~=0);
for mm = 2:length(temp.idxs)
 temp.diff=temp.idxs(mm)-temp.idxs(mm-1);
 
if temp.diff >= deltaforZ
  
   temp.grid{gg}(kk,temp.idxs(mm):end)=0;
  
end


end
end
end

for gg=1:length(groups)

    nonsense=size(temp.grid{gg});

for zz=1:nonsense(2)
for nn=1:length(newY)-1

temp.ngrid{gg}(zz,nn)=mean(nonzeros(temp.grid{gg}(newY(nn):newY(nn+1),zz)));

end
end
    
end



%% the actual graph
close all;

makeblack=0

if makeblack==0
fig=figure('Position',[1,1,round(screensize(3)/2),round(screensize(4)/1.2)]);
lettercolor='k'; oppositeLC='w'; oppositetriplet=[1 1 1]; colortriplet = [0 0 0];
colormap(cp);

else
fig=figure('Position',[1,1,round(screensize(3)/2),round(screensize(4)/1.2)],'color','black');
lettercolor='w'; oppositeLC='k'; oppositetriplet=[0 0 0]; colortriplet = [1 1 1];
colormap(cp);
end
%for black
t=tiledlayout(length(conditions),length(days),'TileSpacing','none','Padding','loose');
sgtitle('Kymograph Z-DFC-NucBcat','FontSize',20,'FontWeight','bold','color',lettercolor);
ylabel(t,'Z centroid (µm)','FontWeight','bold','color',lettercolor);
xlabel(t,'Distance from the center (µm)','FontWeight','bold','FontSize',15,'color',lettercolor);

for gg=1:length(groups)
    nexttile

pplot(gg)=imagesc(temp.ngrid{gg});


pplot(gg).XData = [0:bins:350];

caxis([round(min(gd.nuc.meanAlldata{1}(:,end))) round(max(nuc.maxes)/100)*100]);

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
 %set(gca,'Visible','off');
end

if gg==tiles(posforbar,3)
cb(gg)=colorbar;
ylabel(cb(gg),'Nuclear Bcat Mean Intensity (A.U.)','FontWeight','bold','color',lettercolor);
cb(gg).Color = colortriplet;
cb(gg).Box='off'
end

if any(gg==posforconditions)
    tempidx=find(gg==posforconditions);
    ylabel(conditions{tempidx},'FontSize',15);
end

if any(gg==posford1)
pplot(gg).XData=newY;
xlim([1 351]);
axis.YAxis.FontSize=18;
set(gca,'box','off','Ycolor',lettercolor,'layer','top');

end

if any(gg == gd.days(2,:),'all')>0
    pplot(gg).XData=newY;
     xlim([1 351]);   

end

if any(gg == gd.days(3,:),'all')>0
    pplot(gg).XData=newY;
     xlim([1 351]);   
end

  ylim([2 40]);
  xlim([1 351]);
 
axis.XAxis.FontSize=20;
axis.LineWidth=1;
axis.TickDir='out';
 
end


if makeblack==0
filename=strcat(experiment,'-','Kymograph of DFC-Z-Nuc');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
  saveas(fig,fullfile(graphsRD, filename), 'svg');
else
    fig.InvertHardcopy = 'off';
   filename=strcat(experiment,'-','Kymograph of DFC-Z-Nuc in black');
 saveas(fig,fullfile(graphsRD, filename), 'png');
  saveas(fig,fullfile(graphsRD, filename), 'svg');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
end
%%
