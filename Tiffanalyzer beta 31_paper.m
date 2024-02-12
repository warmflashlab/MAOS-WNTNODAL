%% 1 - Read files

clc; close all; clear; clear memory;

dataDir = '/Volumes//'; %this is the path where the images are
tablesDir = '/Volumes/';
EdataDir = dataDir;
addpath(genpath('~/Documents/GitHub/stemcells')); 


experiment= 'E23-'; %this will be the name of your experiment 
pathtofile =[dataDir];

files= dir(fullfile(pathtofile,'*.tif')); %this will find all the tifs in the path you put before

chlabel={'DAPI', 'ECAD', 'NANOG', 'NCAD'}; % label the cannels
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';

nfiles=length(files);

dataChannels = [1 2 3 4]; %if you have 4 channels, leave it like this

fullDF={};
for ii=1:nfiles
fullDF{ii}=[files(ii).folder '/' files(ii).name];
end


tiffs=cell(1,nfiles);
channels=cell(1,length(dataChannels));

%with Imread
for ii=1:nfiles
    for jj=1:length(dataChannels)

channels{1,jj}=imread(fullDF{ii},jj); 
    end
    tiffs{ii}=channels;
end

clear channels;

shortdash=9; %for long names

for ii=1:length(files)
 finddash=strfind(fullDF(ii),'-');
  conditions{ii}=fullDF{ii}(finddash{1}(shortdash)+1:end-4);
  originalconditions{ii}=fullDF{ii}(finddash{1}(shortdash)+1:end-4);
 
end
conditions

dataChannels = [1 2 3 4];
Alldata={};
resultsname=['results ' experiment];

resultsDir = fullfile(EdataDir,resultsname);
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

imagesRD=[resultsDir '/images/']; 
graphsRD=[resultsDir '/graphs/']; 
masksRD=[resultsDir '/masks/']; 
scatterplotsRD=[resultsDir '/scatterplots/']; 

if ~exist(imagesRD,'dir')
    mkdir(imagesRD);
end

if ~exist(graphsRD,'dir')
    mkdir(graphsRD);
end

if ~exist(masksRD,'dir')
    mkdir(masksRD);
end

if ~exist(scatterplotsRD,'dir')
    mkdir(scatterplotsRD);
end



%% 2 - check the scale bar


spacing = 2; %how much spacing between each slice

%pixtoum=0.325; %this is for the epi 10x
%pixtoum=0.414; %this is for 30x on LSM2
pixtoum=0.6125; %this is for the 20x SD

objfactor=pixtoum;
humbar=100/pixtoum;


SD20x=1; %if the pictures to analyze were not taken with the spinning disk at 20x, change it to zero.
         % 1 means that they were taken at 20x in the spinning disk. 
         
%scale bars are 100 um

%% 3 - find limits auto


limits={};

for cc=1:length(dataChannels)
    for jj=1:length(conditions)
    
        if weirdtiff==1
        limits{cc,jj}=stretchlim(mat2gray(tiffs{jj}{cc}),[0.01 0.995]); %test for weird
        else
         limits{cc,jj}=stretchlim(tiffs{jj}{cc},[0.01 0.995]); %siempre tiene que ir row=channel, column=image
        end
        
    end
end

limits=cell2mat(limits);%est� haciendo la matrix sin orden.
limits=limits';

limits(limits==1)=NaN;

limmin=zeros(1,length(conditions));
limmax=zeros(1,length(conditions));

for jj=1:length(conditions)
for ii=(1:2:length(limits(1,:)))
    limmin(jj,ii)=limits(jj,ii);
   limforall(ii)=min(limmin(:,ii));
end
end

for jj=1:length(conditions)
for ii=(2:2:length(limits(1,:)))
    limmax(jj,ii)=limits(jj,ii);
    limforall(ii)=max(limmax(:,ii));
    
end
end

for ii=(2:2:length(limits(1,:)))
   limits(:,ii)=limforall(ii);
end

if height(limits)>1
limits=max(limits);
end



adjimgs=cell(1,length(conditions));
adjcil=cell(1,length(dataChannels));


for ii=1:length(conditions)
    
    %currentlims=limits(ii,:);
    currentlims=limits;
    
%     just for this
limDAPI=[currentlims(1) currentlims(2)]; %need to automate this
limGFP=[currentlims(3) currentlims(4)];  %need to automate this
limRFP=[currentlims(5) currentlims(6)];  %need to automate this
limCY5=[currentlims(7) currentlims(8)];  %need to automate this
finallimits={limDAPI, limGFP, limRFP, limCY5};    

    
    for jj=1:length(dataChannels)
    if weirdtiff==1
        
         adjcil{jj}=mat2gray(tiffs{ii}{jj},finallimits{jj});
    else
        adjcil{jj}=imadjust(tiffs{ii}{jj},finallimits{jj}); 
    end
    end                                                 
    adjimgs{ii}=adjcil;
end
clear adjcil;


%% 4 - adjust images manually (if you want to set up your own lookup table)
% 
% %from Fiji:
limits={};

DAPI=[650 2200]; %DAPI
GFP=[650 2000]; %GFP
RFP=[600 1500]; %RFP
CY5=[630 2500]; %CY5

for ii=1:length(conditions)
    limits{1,ii}=DAPI;
    limits{2,ii}=GFP;
    limits{3,ii}=RFP;
    limits{4,ii}=CY5;
end

for ii=1:length(conditions)
    for jj=1:length(dataChannels)
    adjcil{1,jj}=uint8((2^8-1)*mat2gray(tiffs{ii}{jj},limits{jj,ii})); 
    end                                                
    adjimgs{ii}=adjcil;
end
clear adjcil;

 
%% 6 - Preview of each channel
 
% 　　　　 SELECT IF YOU WANT COLORBLIND OPTION !!!!!!!!!!! %
% 　　　　 SELECT IF YOU WANT COLORBLIND OPTION !!!!!!!!!!! %
% 　　　　 SELECT IF YOU WANT COLORBLIND OPTION !!!!!!!!!!! %
% 1 = yes , 0 = no;

colorblind=1;

% 　　　　  COLORBLIND OPTION ABOVE   !!!!!!!!!!! %
% 　　　　  COLORBLIND OPTION ABOVE   !!!!!!!!!!! %
% 　　　　  COLORBLIND OPTION ABOVE   !!!!!!!!!!! %


N = numel(conditions);
n = ceil(N/4);
m = ceil(N/n);

screensize = get( 0, 'Screensize' );
margin = 50;
fs = 30-15;
w = screensize(3)/2;
h = n*(screensize(3)/m + margin/2);
yt=100;
p=5;
divisor=length(conditions);
totalmat=length(conditions).*p;
allnum=[1:totalmat];
seriesmat=zeros(divisor,p);
series=zeros(divisor,p);
firstseries=[0:p-1];


xyforbar=[];
 titlestr = [experiment];
  
%DAPI
% DAPI (405) will always be in gray

fig=figure('Position', [1, 1, w, h],'color','black');
  sgtitle(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
for ii=1:length(conditions)
   
subplot_tight(n,m,ii)

 
 
  imshow(cat(3,adjimgs{ii}{1},adjimgs{ii}{1},adjimgs{ii}{1}));
    
        labelstr = ['\color{gray}'   chlabel{1}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
                
                 labelstr={' '}; %for now
                
       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
    
        if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
       
end

filename=strcat(experiment,chlabel{1});
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(imagesRD, filename), 'png');

%GFP

fig=figure('Position', [1, 1, w, h],'color','black');
for ii=1:length(conditions)
   
subplot_tight(n,m,ii)

 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},zeros(size(adjimgs{ii}{2}))));
  labelstr = ['\color{green}'   chlabel{2}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
                
                labelstr={' '}; %for now
                
 else
     imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},adjimgs{ii}{2})); %CYAN
     labelstr = ['\color{cyan}'   chlabel{2}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
 end
     
        
       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
    
        if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
    
    
end

if colorblind==0
filename=strcat(experiment,chlabel{2}, ' green ');
else
    filename=strcat(experiment,chlabel{2},' Cyan ');
end

fig.InvertHardcopy = 'off';
saveas(fig, fullfile(imagesRD, filename), 'png');


%RFP

fig=figure('Position', [1, 1, w, h],'color','black');
for ii=1:length(conditions)
   
subplot_tight(n,m,ii)

 if colorblind==0
  imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),zeros(size(adjimgs{ii}{3}))));
  labelstr = ['\color{red}'   chlabel{3}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
  labelstr={' '}; %for now
 else
     imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),adjimgs{ii}{3}));
     
      labelstr = ['\color{magenta}'   chlabel{3}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]                                   
 end
    

       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
    
        if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
       
end

if colorblind==0
filename=strcat(experiment,chlabel{3}), ' red ';
else
    filename=strcat(experiment,chlabel{3},' Magenta ');
end
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(imagesRD, filename), 'png');

%CY5

fig=figure('Position', [1, 1, w, h],'color','black');
for ii=1:length(conditions)
   
subplot_tight(n,m,ii)
 %asumiendo que tengas 4 canales, otherwise you'll be fucked. 
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{4})),zeros(size(adjimgs{ii}{4})),adjimgs{ii}{4}));
   labelstr = ['\color{blue}'   chlabel{4}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127
  
   labelstr={' '}; %for now
 else
      imshow(cat(3,adjimgs{ii}{4},adjimgs{ii}{4},zeros(size(adjimgs{ii}{4}))));
       labelstr = ['\color{yellow}'   chlabel{4}...
                    '\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127
    
 end
        
       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
                 
                     if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
                    
end
if colorblind==0
filename=strcat(experiment,chlabel{4}, ' blue ');
else
    filename=strcat(experiment,chlabel{4},' yellow ');
end

fig.InvertHardcopy = 'off';
saveas(fig, fullfile(imagesRD, filename), 'png');
pause; close all; %if you want to look at each figure, then comment the close all command. 


%% 7 - Preview each figure with white dapi and no DAPI on merged


% 　　　　 SELECT IF YOU WANT COLORBLIND OPTION !!!!!!!!!!! %
% 　　　　 SELECT IF YOU WANT COLORBLIND OPTION !!!!!!!!!!! %
% 　　　　 SELECT IF YOU WANT COLORBLIND OPTION !!!!!!!!!!! %
% 1 = yes , 0 = no;

colorblind=1;

% 　　　　  COLORBLIND OPTION ABOVE   !!!!!!!!!!! %
% 　　　　  COLORBLIND OPTION ABOVE   !!!!!!!!!!! %
% 　　　　  COLORBLIND OPTION ABOVE   !!!!!!!!!!! %



 N = numel(conditions);
n = ceil(N/4);
m = ceil(N/n);
%m = N/2;
screensize = get( 0, 'Screensize' );
margin = 50;
w = screensize(3);
h = n*(screensize(3)/m + margin/2);
yt=100;
fs = 30;


noDAPImerged={};
if colorblind==0
    
for ii=1:divisor
noDAPImerged{ii}=cat(3,adjimgs{ii}{3}, adjimgs{ii}{2}, adjimgs{ii}{4}); %
end

else
 %CMY
for ii=1:divisor
noDAPImerged{ii}=cat(3,adjimgs{ii}{3} + adjimgs{ii}{4}, adjimgs{ii}{2} + adjimgs{ii}{4} , adjimgs{ii}{2} + adjimgs{ii}{3}); %
end
    
end


seriesmat(1,:)=firstseries;
for ii=2:divisor
seriesmat(ii,:)=seriesmat(ii-1,:)+4;
end
allnum=reshape(allnum,[divisor,p])';

for ii=1:divisor
 
merged{ii}=cat(3,adjimgs{ii}{3}+adjimgs{ii}{4}.*.01, adjimgs{ii}{2}+adjimgs{ii}{4}.*.01, adjimgs{ii}{1}+adjimgs{ii}{4}.*.01); %this will later be used for masks.     
yforfig=size(adjimgs{ii}{1});    
fig=figure('Position', [1, 1, w, h/1.1 ],'Color','Black');
sgtitle(experiment,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
 
%merged
subplot_tight(1,5,1)

    imshow(noDAPImerged{ii});
        titlestr = strcat('Merged - ', conditions{ii});
       title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
   if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
      
       
subplot_tight(1,5,2)
  imshow(cat(3,adjimgs{ii}{1},adjimgs{ii}{1},adjimgs{ii}{1}));

  
    titlestr = chlabel{1};
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
       if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
    
subplot_tight(1,5,3)
 %488 CHANNEL
 
 titlestr = chlabel{2};
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},zeros(size(adjimgs{ii}{2}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','green');
 else
     imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},adjimgs{ii}{2}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','cyan');
 end
       
 if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
    
subplot_tight(1,5,4)
 %555 CHANNEL
 titlestr = chlabel{3};
 
 if colorblind==0
  imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),zeros(size(adjimgs{ii}{3}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','red');
 else
     imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),adjimgs{ii}{3}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','magenta');
 end
        
           if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
           end
       
subplot_tight(1,5,5)
 %647 CHANNEL
 titlestr = chlabel{4};
 
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{4})),zeros(size(adjimgs{ii}{4})),adjimgs{ii}{4}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','blue');
        filename=strcat(experiment,conditions{ii},'-',string(ii),' not merged WHITE DAPI and all channels');

 else
   imshow(cat(3,adjimgs{ii}{4},adjimgs{ii}{4},zeros(size(adjimgs{ii}{4}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','yellow'); 
        filename=strcat(experiment,conditions{ii},'-',string(ii),' not merged WHITE DAPI and all channels CMY');

 end
 
 if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
 end
       
                    set(gca,'color','black')


fig.InvertHardcopy = 'off';
filename=strcat(experiment,conditions{ii},'-',string(ii),' not merged WHITE DAPI and all channels');
saveas(fig, fullfile(imagesRD, filename), 'png');

end
pause
close all;



%% 8 - Configure how you want the preview

% FIRST indicate if you want colorblind:

colorblind=0;

%indicate which channel you want off or on in the merge;

DAPI=0;
GFP=1;
RFP=1;
CY5=0;

chtext=chlabel;
%all process starts from here:
chanvalue1=1; chanvalue2=1; chanvalue3=1; chanvalue4=1;

if DAPI<1
    chanvalue1=0.05;
    chtext{1}=' ';
end

if GFP<1
    chanvalue2=0.05;
    chtext{2}=' ';
end

if RFP<1
    chanvalue3=0.05;
    chtext{3}=' ';
end

if CY5<1
    chanvalue4=0.05;
    chtext{4}=' ';
end

if colorblind==0
    
for ii=1:divisor
 merged{ii}=cat(3,adjimgs{ii}{3}.*chanvalue3+adjimgs{ii}{1}.*chanvalue1, adjimgs{ii}{2}.*chanvalue2+adjimgs{ii}{1}.*chanvalue1, adjimgs{ii}{4}.*chanvalue4+adjimgs{ii}{1}.*chanvalue1); %asumiendo que tengas 4 canales, otherwise you'll be fucked.    
end

else
    
for ii=1:divisor
merged{ii}=cat(3,adjimgs{ii}{3}.*chanvalue3 + adjimgs{ii}{4}.*chanvalue4 +adjimgs{ii}{1}.*chanvalue1, adjimgs{ii}{2}.*chanvalue2 + adjimgs{ii}{4}.*chanvalue4 +adjimgs{ii}{1}.*chanvalue1 , adjimgs{ii}{2}.*chanvalue2 + adjimgs{ii}{3}.*chanvalue3 +adjimgs{ii}{1}.*chanvalue1); %
end

end


for ii=1:divisor
 
yforfig=size(adjimgs{ii}{1});    
fig=figure('Position', [1, 1, w, yforfig(1)./1.5 ],'Color','Black');
sgtitle([experiment '-' [chtext{:}]],'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
 
%merged
subplot_tight(1,5,1)

    imshow(merged{ii});
        titlestr = strcat('Merged - ', conditions{ii});
       title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
   if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
   end
      
   %405 CHANNEL here :
   
     
subplot_tight(1,5,2)

  imshow(cat(3,adjimgs{ii}{1},adjimgs{ii}{1},adjimgs{ii}{1}));

  
    titlestr = chlabel{1};
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
       if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
    
subplot_tight(1,5,3)

 %488 CHANNEL from here:
 
 titlestr = chlabel{2};
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},zeros(size(adjimgs{ii}{2}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','green');
 else
     imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},adjimgs{ii}{2}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','cyan');
 end
       
 if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
       
 % 555 CHANNEL from here:       
         
subplot_tight(1,5,4)
 titlestr = chlabel{3};
 
 if colorblind==0
  imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),zeros(size(adjimgs{ii}{3}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','red');
 else
     imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),adjimgs{ii}{3}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','magenta');
 end
        
           if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
           end
       
subplot_tight(1,5,5)
 %647 CHANNEL
 titlestr = chlabel{4};
 
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{4})),zeros(size(adjimgs{ii}{4})),adjimgs{ii}{4}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','blue');
        filename=strcat(experiment,conditions{ii},'-',string(ii),' not merged WHITE DAPI and all channels');

 else
   imshow(cat(3,adjimgs{ii}{4},adjimgs{ii}{4},zeros(size(adjimgs{ii}{4}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','yellow'); 
        filename=strcat(experiment,conditions{ii},'-',string(ii),' not merged WHITE DAPI and all channels CMY');

 end
 
 if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
 end
       
                    set(gca,'color','black')

%save


fig.InvertHardcopy = 'off';
saveas(fig, fullfile(imagesRD, filename), 'png');

end
pause
close all;

%% 8b - if you want DAPI with a color


% FIRST indicate if you want colorblind:

colorblind=1;

%indicate which channel you want off other than DAPI;
%one channel must be off for DAPI to switch.

DAPI=1; %this always has to be ON

%one of these has to be off:
GFP=1;
RFP=0;
CY5=1;


%all process starts from here:
chtext=chlabel;
chanvalue1=1; chanvalue2=1; chanvalue3=1; chanvalue4=1;
cDAPI=1; cGFP=2; cRFP=3; cCY5=4;

if DAPI<1
    chanvalue1=0.05;
    chtext{1}=' ';
end

if GFP<1
    cGFP=1;
    chtext{2}=' ';
end

if RFP<1
    cRFP=1;
    chtext{3}=' ';
end

if CY5<1
    cCY5=1;
    chtext{4}=' ';
end

if colorblind==0
    
for ii=1:divisor
 merged{ii}=cat(3,adjimgs{ii}{cRFP}.*chanvalue3, adjimgs{ii}{cGFP}.*chanvalue2, adjimgs{ii}{cCY5}.*chanvalue4); %asumiendo que tengas 4 canales, otherwise you'll be fucked.    
end

else
    
for ii=1:divisor
merged{ii}=cat(3,adjimgs{ii}{cRFP}.*chanvalue3 + adjimgs{ii}{cCY5}.*chanvalue4, adjimgs{ii}{cGFP}.*chanvalue2 + adjimgs{ii}{cCY5}.*chanvalue4, adjimgs{ii}{cGFP}.*chanvalue2 + adjimgs{ii}{cRFP}.*chanvalue3); %
end

end

for ii=1:divisor
 
yforfig=size(adjimgs{ii}{1});    
fig=figure('Position', [1, 1, w, yforfig(1)./1.5 ],'Color','Black');
sgtitle([experiment '-' [chtext{:}]],'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
 
%merged
subplot_tight(1,5,1)

    imshow(merged{ii});
        titlestr = strcat('Merged - ', conditions{ii});
       title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
   if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
   end
      
   %405 CHANNEL here :
   
     
subplot_tight(1,5,2)
  imshow(cat(3,adjimgs{ii}{1}.*.01,adjimgs{ii}{1}.*.01,adjimgs{ii}{1}.*.01));

  
    titlestr = chlabel{1};
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','black');
       if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
    
subplot_tight(1,5,3)

 %488 CHANNEL from here:
 
 titlestr = chlabel{cGFP};
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{cGFP})),adjimgs{ii}{cGFP},zeros(size(adjimgs{ii}{cGFP}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','green');
 else
     imshow(cat(3,zeros(size(adjimgs{ii}{cGFP})),adjimgs{ii}{cGFP},adjimgs{ii}{cGFP}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','cyan');
 end
       
 if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
       end
       
 % 555 CHANNEL from here:       
         
subplot_tight(1,5,4)
 titlestr = chlabel{cRFP};
 
 if colorblind==0
  imshow(cat(3,adjimgs{ii}{cRFP},zeros(size(adjimgs{ii}{cRFP})),zeros(size(adjimgs{ii}{cRFP}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','red');
 else
     imshow(cat(3,adjimgs{ii}{cRFP},zeros(size(adjimgs{ii}{cRFP})),adjimgs{ii}{cRFP}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','magenta');
 end
        
           if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
           end
       
subplot_tight(1,5,5)
 %647 CHANNEL
 titlestr = chlabel{cCY5};
 
 if colorblind==0
  imshow(cat(3,zeros(size(adjimgs{ii}{cCY5})),zeros(size(adjimgs{ii}{cCY5})),adjimgs{ii}{cCY5}));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','blue');
        filename=strcat(experiment,conditions{ii},'-',string(ii),' DAPI switched ');

 else
   imshow(cat(3,adjimgs{ii}{cCY5},adjimgs{ii}{cCY5},zeros(size(adjimgs{ii}{cCY5}))));
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','yellow'); 
        filename=strcat(experiment,conditions{ii},'-',string(ii),' DAPI switched CMY');

 end
 
 if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
 end
       
                    set(gca,'color','black')

%save


fig.InvertHardcopy = 'off';
saveas(fig, fullfile(imagesRD, filename), 'png');
end
pause
close all;

%% 9 - read Ilastik's masks (if you have tables from Ilastik skip until 9B)

%the mask should have the label 1 (yellow) should be background and label 2
%(blue) should be the cells !!!! this is important otherwise it won't work.

close all;
nffmasks={};
h5cp={}; h5files={};
    
nucchannel=1; %where's the nuc?
compressed=0; %is it compressed h5file?
identities=0;

%Ilastik exports masks as t,y,x,c . However, the only thing needed for this
%section to work is have the mask in this arrangement: (1,1024,1024,1). 
%if you have more than one channel mask or time, it won't let you,
%although it should, but there's some bug with reading and transposing.
    
if identities>0
   
for ii=1:length(conditions)
[nffmasks{ii,1},nffmasks{ii,2},nffmasks{ii,3}]=fileparts(fullDF{ii});
nffmasks{ii,4}='_Object Identities.h5';
nffmasks{ii,5}=strcat(nffmasks{ii,2},nffmasks{ii,4});
h5cp{ii}=fullfile(dataDir,nffmasks{ii,5});
h5IDs{ii}=h5read(h5cp{ii}, '/exported_data');
h5files{ii}=h5IDs{ii}>0;
%transform identities to segmask

h5files{ii}=squeeze(h5files{ii}==nucchannel);
h5files{ii}=h5files{ii}(:,:,1)'; 

   
end

else 
    
for ii=1:length(conditions)
[nffmasks{ii,1},nffmasks{ii,2},nffmasks{ii,3}]=fileparts(fullDF{ii});
nffmasks{ii,4}='_Simple Segmentation.h5';
nffmasks{ii,5}=strcat(nffmasks{ii,2},nffmasks{ii,4});
h5cp{ii}=fullfile(dataDir,nffmasks{ii,5});

%h5cp{ii}=[SSfiles(ii).folder '/' SSfiles(ii).name];
h5files{ii}=h5read(h5cp{ii}, '/exported_data');


if compressed>0
    h5files{ii}=squeeze(h5files{ii}==nucchannel);
   h5files{ii}=h5files{ii}(:,:,1)'; 
else
   
h5files{ii}=squeeze(h5files{ii}==nucchannel)';
end

end

end


%% 10 - Very Simple Cleanup

%use this cleanup if watershed it's deleting too much cells or if ilastik's segmentation is already fine

for ii=1:length(conditions)
cleanmask{ii}=imerode(h5files{ii},strel('disk',1));

end

for ii=1:length(conditions)
cleanmask{ii}=imdilate(cleanmask{ii},strel('disk',1));
end

 %less than this pixels will be eliminated
for ii=1:length(conditions)
cleanmask{ii}=bwareaopen(cleanmask{ii},50);
end

%% 11 - preview of the before/after cleanup mask     
   

N = numel(conditions);
n = ceil(N/4);
m = ceil(N/n);

%Antes
fig=figure('Position', [1, 1, w, h]); 
title('Before');
for ii=1:length(conditions)
   
subplot_tight(n,m,ii);   
 imshow(h5files{ii});

        labelstr = ['\color{Blue}' ' Before ' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold');
   
end


filename=strcat(experiment,' mask before process');
 saveas(fig, fullfile(masksRD, filename), 'png');

 
%Despues

fig=figure('Position', [1, screensize(4), w, h]); 
title('After');
for ii=1:length(conditions)
   
subplot_tight(n,m,ii);   
 imshow(cleanmask{ii});

        labelstr = ['\color{Blue}' ' After ' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold');
   
end

filename=strcat(experiment,' mask after process');
 saveas(fig, fullfile(masksRD, filename), 'png');
 
 %pause(1); close all;

 
 
 %% 12 - Verify mask overlays with image
 
     %revisa la mascara ! 

for ii=1:length(cleanmask)
bordermasks{ii}=cleanmask{ii}-imerode(cleanmask{ii}, strel('disk',2));    
end

fs=30;

for ii=1:length(conditions)
fig=figure('Position', [1, 1, w, h]); 
title('Mask overlay segmentation');

   
imshowpair(bordermasks{ii},merged{ii});

  
        labelstr = ['\color{white}' conditions{ii}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  yt , labelstr,'FontSize',fs,'FontWeight','bold');
       
       
filename=strcat(experiment,conditions{ii},' overlay mask',string(ii));
 saveas(fig, fullfile(masksRD, filename), 'png');
   
end
 
 pause; close all;
 
%% 9B Read tables from Ilastik

channels=length(dataChannels);
nuctable={};
nucfiles=dir(fullfile(tablesDir,'*.csv')); %make sure you looking for the same tables and that it matches the name/number


for ii=1:nfiles
  nuctable{ii}=readtable([nucfiles(ii).folder '/' nucfiles(ii).name]); 
   tiffsize{ii}=size(tiffs{ii}{1});
end

if size(nuctable)<1
    error('no tables in this folder');
else
    disp('number of tables found:');
    length(nuctable)
end

if length(nuctable)==nfiles
   cprintf('blue','same tables as tifs');
else
    error('not the same tables as tif, check directories');
end

%h5 files sizes is [x y c z t]

%% 10B arrange extracted data from tables
% 
 nuc.ilabel={}; nuc.xpos={}; nuc.ypos={}; nuc.zpos={}; nuc.mean={}; nuc.class={}; nuc.Alljuntits={}; nuc.Alldata={}; nuc.area={};
% 

for ii=1:nfiles

nuc.ilabel{ii}=nuctable{ii}.labelimage_oid;

nuc.xpos{ii}=nuctable{ii}.ObjectCenter_0;
nuc.ypos{ii}=nuctable{ii}.ObjectCenter_1;
nuc.zpos{ii}=nuctable{ii}.ObjectCenter_2;

nuc.class{ii}=nuctable{ii}.PredictedClass;
nuc.dfc{ii}= sqrt((nuc.xpos{ii}-tiffsize{ii}(1)/2).^2 + (nuc.ypos{ii}-tiffsize{ii}(2)/2).^2).*objfactor;
nuc.area{ii}=nuctable{ii}.ObjectArea;

%that +1 is so the time does not start with zero

nuc.mean{ii}{1}=nuctable{ii}.MeanIntensity_0; %0 is the first channel, 1 is the second

if channels>1
nuc.mean{ii}{2}=nuctable{ii}.MeanIntensity_1; 
end
if channels>2
nuc.mean{ii}{3}=nuctable{ii}.MeanIntensity_2; 
end
if channels>3
nuc.mean{ii}{4}=nuctable{ii}.MeanIntensity_3; 
end

ii

nuc.Alljuntits{ii}=[repmat(ii,length(nuc.ilabel{ii}),1) nuc.ilabel{ii} nuc.area{ii} nuc.xpos{ii} nuc.ypos{ii} nuc.dfc{ii} nuc.zpos{ii} [nuc.mean{ii}{:}] ];


end

description=' file number (1) - ilabel (2) - area (3) - xpos (4) - ypos (5) - distance from center (6) - zpos (7) - mean intensities (8-end) ';

unfiltered.nuc.Alljuntits=nuc.Alljuntits;

%% 11B Filters and re-grouping
 
% % AREA FILTER
% % DFC FILTER
% % Z FILTER
% 
 filterarea=100; %any area below this will be eliminated;
 filterdfc=500; %any distance from the center greater than this will be deleted
 filterzc=1; %any z-centroid below this will be eliminated
 
% 
% % AREA FILTER
% % DFC FILTER
% % Z FILTER
% 
% 

for ii=1:nfiles

       nuc.Alljuntits{ii}(nuc.Alljuntits{ii}(:,3)<filterarea,:)=[]; %area
       nuc.Alljuntits{ii}(nuc.Alljuntits{ii}(:,6)>filterdfc,:)=[]; %dfc
       nuc.Alljuntits{ii}(nuc.Alljuntits{ii}(:,7)<filterzc,:)=[]; %Z
      
end

% start arranging data
nuc.meanAlljuntits={}; nuc.stdAlljuntits={}; 
for ii=1:nfiles

nuc.meanAlljuntits{ii}=mean(nuc.Alljuntits{ii}(:,8:end)) ; %last row because last row is mean intensity
nuc.stdAlljuntits{ii}=std(nuc.Alljuntits{ii}(:,8:end))./2;

end

nuc.Alldata=nuc.Alljuntits;

divisor= nfiles/length(conditions);

%curate all the scatter data

curate= 99; % what percentage of all data want to extract the maximum from?
nuc.checkformax={}; %nuc.maxes={}; nuc.mins={};

%for mean intensity
for ii=1:nfiles
nuc.checkformax{ii}=sort(nuc.Alldata{ii}(:,8:end));

nuc.maxes{ii}=nuc.checkformax{ii}(round(height(nuc.checkformax{ii})*curate/100),:); 
nuc.mins{ii}=min(nuc.Alldata{ii}(:,8:end));

end
 
%% 13 - All data & normalization
tic
Alldata={};

nfiles=length(files);
%
% CHECK IF THE READTABLES OPTION IS ON OR OFF !!!!!!   
% CHECK IF THE READTABLES OPTION IS ON OR OFF !!!!!!   
% CHECK IF THE READTABLES OPTION IS ON OR OFF !!!!!!   
% CHECK IF THE READTABLES OPTION IS ON OR OFF !!!!!!   
%
%
readtables=1; %this is when you have an exported csv table form Ilastik


if readtables>0
  
     % tablesfiles=dir(fullfile(pathtofile,'*.csv'));   
   
    for ii=1:nfiles
      
      
    Alldata{ii}{1}.Area=nuc.Alldata{ii}(:,3);
    Alldata{ii}{1}.Centroid=[nuc.Alldata{ii}(:,4) nuc.Alldata{ii}(:,5)];
       
    Alldata{ii}{1}.MeanIntensity=nuc.Alldata{ii}(:,8);
    Alldata{ii}{2}.MeanIntensity=nuc.Alldata{ii}(:,9);
    Alldata{ii}{3}.MeanIntensity=nuc.Alldata{ii}(:,10);
    Alldata{ii}{4}.MeanIntensity=nuc.Alldata{ii}(:,11);
            
      
    end
else
        for ii=1:nfiles
        for jj=1:length(dataChannels)
        
        Alldata{ii}{jj}=regionprops(cleanmask{ii},tiffs{ii}{jj},'Area','Centroid','MeanIntensity');
        end
        end

end

ncells={}; Iave=cell(1,nfiles); stdIave={};

for ii=1:nfiles
for jj=1:length(dataChannels)
    

    
    if isempty([Alldata{ii}{jj}.MeanIntensity])==1
        ncells{ii}=1; Iave{ii}{jj}=1; stdIave{ii}{jj}=1;
    else
     ncells{ii}=length([Alldata{ii}{jj}.MeanIntensity]);
    Iave{ii}{jj}=mean([Alldata{ii}{jj}.MeanIntensity]);
     stdIave{ii}{jj}=std([Alldata{ii}{jj}.MeanIntensity])./2;
    end
    
    
    
end
end

centroids={};
for ii=1:nfiles
      
    if isempty([Alldata{ii}{1}.Centroid])==1
            centroids{ii}=[round(size(tiffs{1}{1})./2)]; %if there is no centroid, then make it at the center
    else
         centroids{ii}=cat(1,Alldata{ii}{1}.Centroid);
    end
   
end

normalizedIave={};  normalizedstdIave={};
normalizedIave=cell(1,nfiles); normalizedstdIave=cell(1,nfiles);
for ii=1:nfiles
for jj=1:length(chlabel)
    normalizedIave{ii}{jj}=Iave{ii}{jj}./Iave{ii}{1};
     normalizedstdIave{ii}{jj}=stdIave{ii}{jj}./Iave{ii}{1};
end
end
% 
AMIs={};
AMIs=cell(1,nfiles);
%just a matrix of the mean intensity
for ii=1:nfiles
    for cc=1:length(dataChannels)
   
        if isempty([Alldata{ii}{cc}.MeanIntensity])==1
            AMIs{ii}(:,cc)=[1];
        else
        AMIs{ii}(:,cc)=[Alldata{ii}{cc}.MeanIntensity];
        end
    end
end

nAMIs={};
nAMIs=cell(1,nfiles);
%just a matrix of the mean intensity normalized to DAPI
for ii=1:nfiles
    for cc=1:length(dataChannels)
nAMIs{ii}(:,cc)=AMIs{ii}(:,cc)./AMIs{ii}(:,1);
    end
end


chAMIs={};
for cc=1:length(dataChannels)
    for ii=1:nfiles
     
 chAMIs{cc}{ii}=AMIs{ii}(:,cc); %new one
    end
end

centerofimgs={};
for ii=1:nfiles
centerofimgs{ii}=ceil(size(tiffs{ii}{1})/2);
end

if readtables>0
%distance from the center
for ii=1:nfiles
dfcenter{ii} = nuc.Alljuntits{ii}(:,6);
end
else
    for ii=1:nfiles
   dfcenter{ii} = sqrt((centroids{ii}(:,1)-centerofimgs{ii}(2)).^2 + (centroids{ii}(:,2)-centerofimgs{ii}(1)).^2).*pixtoum; 
    end
end

Areas={}; TotalArea={};
for ii=1:nfiles
    Areas{ii}=[Alldata{ii}{1}.Area]';
    TotalArea{ii}=sum(Areas{ii});
end



%% 13b - group data obtained


groupdata=1; %you'll need to manually modify. if each condition is separate and there's no data to be grouped, put 0 and continue.


groupsnames={}; groupsindex={};

if groupdata==1
  %originalconditions=conditions;  

groupsnames = {'1','2'};


groupsindex{1}={}; %
groupsindex{2}={}; %

groups = length(groupsnames);


if length([groupsindex{:}]) ~= nfiles
    cprintf('red','The files manually indexed do not match the total number of files');
   
    disp('the files you indexed are:');
    length([groupsindex{:}])
    disp('the total files are:');
    nfiles
     disp('Is that what you intended?');
     %error;
end
 
%from this point, all should be automatic: 

gncells={}; gIave=cell(1,groups); 


%ncells in groups ! 
for gg= 1:groups
for ssgg=1:length(groupsindex{gg})
   

if ssgg==1
gncells{gg}=ncells{groupsindex{gg}{ssgg}};
else
gncells{gg}=gncells{gg}+ncells{groupsindex{gg}{ssgg}};
end

end
end

%Iave and stdIave in groups ! 

for gg= 1:groups
    clear tempv; 
for ssgg=1:length(groupsindex{gg})
 for cc=1:length(dataChannels)
    
if ssgg==1
gIave{gg}{cc}=Iave{groupsindex{gg}{ssgg}}{cc};
%gstdIave{gg}{cc}=stdIave{groupsindex{gg}{ssgg}}{cc};
tempv(ssgg)=stdIave{groupsindex{gg}{ssgg}}{cc};

else
gIave{gg}{cc}=gIave{gg}{cc}+Iave{groupsindex{gg}{ssgg}}{cc};
%gstdIave{gg}{cc}=gstdIave{gg}{cc}+stdIave{groupsindex{gg}{ssgg}}{cc};
tempv(ssgg)=stdIave{groupsindex{gg}{ssgg}}{cc};
end

if ssgg==length(groupsindex{gg})
gIave{gg}{cc}=gIave{gg}{cc}/ssgg;
%gstdIave{gg}{cc}=gstdIave{gg}{cc}/ssgg;

tempsizeot=size(tempv);

if tempsizeot(2)<2
    gstdIave{gg}{cc}=stdIave{groupsindex{gg}{ssgg}}{cc};
else
gstdIave{gg}{cc}=std(tempv)./2;
end

end

end
end
end


%centroids in groups

for gg=1:groups
for ssgg=1:length(groupsindex{gg})

if ssgg==1
gcentroids{gg}=centroids{groupsindex{gg}{ssgg}};
else
gcentroids{gg}=cat(1,gcentroids{gg},centroids{groupsindex{gg}{ssgg}});
end

end
end


gnormalizedIave={}; gnormalizedstdIave={};

for gg= 1:groups
 for cc=1:length(dataChannels)
    
gnormalizedIave{gg}{cc}=gIave{gg}{cc}/gIave{gg}{1};
gnormalizedstdIave{gg}{cc}=gstdIave{gg}{cc}/gstdIave{gg}{1};

end
end


%for GROUPS start here !!!!
gAMIs={}; gnAMIs={}; gdfcenter={};

for gg=1:groups
for ssgg=1:length(groupsindex{gg})

    
    if ssgg==1
gAMIs{gg}=AMIs{groupsindex{gg}{ssgg}};
gnAMIs{gg}=nAMIs{groupsindex{gg}{ssgg}};
gdfcenter{gg}=dfcenter{groupsindex{gg}{ssgg}};
    else
gAMIs{gg}=cat(1,gAMIs{gg},AMIs{groupsindex{gg}{ssgg}});
gnAMIs{gg}=cat(1,gnAMIs{gg},nAMIs{groupsindex{gg}{ssgg}});
gdfcenter{gg}=cat(1,gdfcenter{gg},dfcenter{groupsindex{gg}{ssgg}});
    end

end
end


%gropuing the grouped AMIs by channels ! 

gchAMIs={};
% 


for cc= 1:length(dataChannels)
for gg= 1:groups
% for nn= 1:gncells{gg}
gchAMIs{cc}{gg}=gAMIs{gg}(:,cc);
%end
end
end
%

gTotalArea={}; %total area of all objects in the mask
gAreas={};

for gg=1:groups
for ssgg=1:length(groupsindex{gg});
    
    if ssgg==1
           gTotalArea{gg}=TotalArea{groupsindex{gg}{ssgg}};
           gAreas{gg}=Areas{groupsindex{gg}{ssgg}};
    else
           gTotalArea{gg}=gTotalArea{gg}+TotalArea{groupsindex{gg}{ssgg}};
           if readtables==1
               gAreas{gg}=cat(2,gAreas{gg},Areas{groupsindex{gg}{ssgg}}); % this could be 1 or 2
           else
           gAreas{gg}=cat(1,gAreas{gg},Areas{groupsindex{gg}{ssgg}}); % this could be 1 or 2
           end
    end
end
end

ncells=gncells;
Iave=gIave;
stdIave=gstdIave;
oricentroids=centroids; centroids=gcentroids;
normalizedIave=gnormalizedIave;
normalizedstdIave=gnormalizedstdIave;
oriAMIs=AMIs; AMIs=gAMIs;
nAMIs=gnAMIs;
chAMIs=gchAMIs;
oricenterofimgs=centerofimgs; %centerofimgs=gcenterofimgs;
dfcenter=gdfcenter;
conditions=groupsnames;
nfiles=groups;
divisor=groups;
oriTotalArea=TotalArea; oriAreas=Areas;
TotalArea=gTotalArea; 
Areas=gAreas;

cprintf('blue','Data grouped');

%group tables 

if readtables>0
for gg = 1:length(groupsindex)
for ssgg=1:length(groupsindex{gg})    
    
    if ssgg==1
gd.nuc.Alldata{gg}=nuc.Alldata{groupsindex{gg}{ssgg}};            
gd.nuc.meanAlldata{gg}= nuc.meanAlljuntits{groupsindex{gg}{ssgg}};    
    else
gd.nuc.Alldata{gg}=cat(1,gd.nuc.Alldata{gg},nuc.Alldata{groupsindex{gg}{ssgg}});                  
gd.nuc.meanAlldata{gg}=cat(1,gd.nuc.meanAlldata{gg},nuc.meanAlljuntits{groupsindex{gg}{ssgg}});
gd.nuc.meanAlldata{gg}=mean(gd.nuc.meanAlldata{gg});
    end
  
end

gd.nuc.stdAlldata{gg}=std(gd.nuc.Alldata{gg}(:,8:end))./2;

end

maxZ={}; percentile=99;
% to change the Z to um
for gg=1:groups
gd.nuc.originalZ{gg}=gd.nuc.Alldata{gg}(:,7);
gd.nuc.Alldata{gg}(:,7)=(gd.nuc.Alldata{gg}(:,7).*spacing);
tempZ=gd.nuc.Alldata{gg}(:,7);

n2look=round(length(gd.nuc.Alldata{gg}(:,7))*(percentile/100));
maxZ{gg}=round(tempZ(n2look));
end
maxZ

maxZ=cell2mat(maxZ);
maxZ=max(maxZ);
end

end

 %% 15 - Mean intensities in double bars 
 

cmap=colormap('lines'); close;

if divisor>8
    cmap=colormap('turbo'); close;
    cidxs=round(linspace(10,length(colormap)-10,divisor));
    
    cidxs(2:2:end)=fliplr(cidxs(2:2:end));
    
    newcmap=cmap(cidxs,:); close;
    cmap=newcmap;
end

%cmap=distinguishable_colors(divisor,'k');

rot=0; %for vertical text = 90, for horizontal text on top of bar = 0.

if rot==0
    fs=15;
    %fs=30-15;
else 
    fs=30-(divisor.*2);
end

aIave=zeros(nfiles,length(dataChannels));
stdaIave=zeros(nfiles,length(dataChannels));
anIave=zeros(nfiles,length(dataChannels));
stdnaIave=zeros(nfiles,length(dataChannels));

for ii=1:nfiles
aIave(ii,:)=cell2mat([Iave{ii}]);
stdaIave(ii,:)=cell2mat([stdIave{ii}]);
anIave(ii,:)=cell2mat([normalizedIave{ii}]);
stdnaIave(ii,:)=cell2mat([normalizedstdIave{ii}]);
end

xmax=6; % this has always to be the same
ymax=max(ceil(max(aIave))+ceil(max(stdaIave)));

cxbar=categorical(chlabel(2:4));
%cxbar=reordercats(cxbar,{'SOX2' 'BRA' 'HOXB4'});

Iavelist=string(round(aIave,2));
Iavestdlist=string(round(stdaIave,2));
Iaveplusstd=strcat(Iavelist(:),' +/- ',Iavestdlist(:));
ytextp=aIave(:,2:4)+stdaIave(:,2:4);

% FIGUREEEEEEEEE in white
fig=figure('Position', [1, 1, w/2, h]); 

b=bar([aIave(:,2)'; aIave(:,3)'; aIave(:,4)']); 

Xb=[];
for bb=1:length(b)
    Xb=[Xb;b(bb).XData+b(bb).XOffset];
end

for bb=1:length(b)
    b(bb).FaceColor=cmap(bb,:);
end

hold on; 
errorbar(Xb(:)', [aIave(:,2)' aIave(:,3)' aIave(:,4)'], [stdaIave(:,2)' stdaIave(:,3)' stdaIave(:,4)'],'.','Color',[.25 .25 .25]);
ylim([0 ymax]);
xlim([0 length(chlabel)]);
legend(conditions,'Location','northeastoutside');

     titlestr = ['Mean Intensity Average'];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
       
set(gca,'FontSize',30, 'xticklabel', chlabel(2:4))

filename=strcat(experiment,'Bars grouped Mean Intensity average');
saveas(fig, fullfile(graphsRD, filename), 'png');
 saveas(fig, fullfile(graphsRD, filename), 'svg');
  savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');

%FIGUREEEEEE in BLACK !!!!!!!

fig=figure('Position', [1, 1, w/2, h],'color', 'black'); 
b=bar([aIave(:,2)'; aIave(:,3)'; aIave(:,4)']); 
Xb=[];
for bb=1:length(b)
    Xb=[Xb;b(bb).XData+b(bb).XOffset];
end

for bb=1:length(b)
    b(bb).FaceColor=cmap(bb,:);
end
 
hold on; 
errorbar(Xb(:)', [aIave(:,2)' aIave(:,3)' aIave(:,4)'], [stdaIave(:,2)' stdaIave(:,3)' stdaIave(:,4)'],'.','Color',[.75 .75 .75]);
ylim([0 ymax]);
xlim([0 length(chlabel)]);
legend(conditions,'Location','northeastoutside','TextColor','white');

     titlestr = ['Mean Intensity Average'];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');   

axis=gca
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
        
set(gca,'FontSize',30, 'xticklabel', chlabel(2:4),'color','black')
set(legend,'color','black')


filename=strcat(experiment,'Bars grouped Mean Intensity average in black');
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
 saveas(fig, fullfile(graphsRD, filename), 'svg');
 
 pause; close all;
 
%% new 14 - MI from center to edge of the colony.


colonysize=700; %total colony size

%cmap=distinguishable_colors(nfiles); % TEMPRAOALLLL

newX=[0:25:colonysize/2]; %radius is the top now

DaMIs={};
for ii=1:length(conditions)
DaMIs{ii}=[dfcenter{ii} AMIs{ii}]; %distance from center, channel 1, channel 2, channel 3, channel 4. 
end
% 
ViR={}; stdViR={}; %values in radius/range
for ii=1:length(conditions)
for nn=1:length(newX)-1
tempidx=find(DaMIs{ii}(:,1)>newX(nn) & DaMIs{ii}(:,1)<newX(nn+1));
templist=DaMIs{ii}(tempidx,:);

ViR{ii}(nn,:)=[newX(nn+1) mean(templist(:,2)) mean(templist(:,3)) mean(templist(:,4)) mean(templist(:,5))];
% 
 if isnan(ViR{ii}(nn,2:end))>0
     ViR{ii}(nn,:)=[newX(nn+1) min(DaMIs{ii}(:,2)) min(DaMIs{ii}(:,3)) min(DaMIs{ii}(:,4)) min(DaMIs{ii}(:,5)) ]; %this is in case there's no cells at the border which is just simply not possible, if this happens is due bad segmentation.
 end

stdViR{ii}(nn,:)=[newX(nn+1) std(templist(:,2))/2 std(templist(:,3))/2 std(templist(:,4))/2 std(templist(:,5))/2];
stderrormViR{ii}(nn,:)=[newX(nn+1) std(templist(:,2))/length(templist(:,2)) std(templist(:,3))/length(templist(:,3)) std(templist(:,4))/length(templist(:,4)) std(templist(:,5))/length(templist(:,5))];

% 

 if isnan(ViR{ii}(nn,2:end))>0
     ViR{ii}(nn,:)=[newX(nn+1) 0 0 0]; %this is in case there's no cells at the border which is just simply not possible, if this happens is due bad segmentation.
 end


end
end
% mean a cada columna y su std
for ii=1:length(conditions)
ViR{ii}=sortrows(cat(1,ViR{ii},[0 ViR{ii}(1,2:end)]));
stdViR{ii}=sortrows(cat(1,stdViR{ii},[0 stdViR{ii}(1,2:end)]));
end

%sacar el maximo de maximos Y
maxY={}; minY={};
for ii=1:length(conditions)
for cc=1:length(dataChannels)
maxY{ii}(cc+1)=max(ViR{ii}(:,cc+1)+stdViR{ii}(:,cc+1));
%minY{ii}(cc+1)=min(ViR{ii}(:,cc+1)-stdViR{ii}(:,cc+1));
minY{ii}(cc+1)=min(DaMIs{ii}(:,cc+1));
end
end
maxY=reshape(cell2mat(maxY),[5 length(conditions)])'; minY=reshape(cell2mat(minY),[5 length(conditions)])';
maxY=max(maxY); minY=min(minY);

for cc=1:length(dataChannels)
if maxY(cc+1)<1000
    maxY(cc+1)=1000;
end
end
% 

if colorblind==1
    titlecolors={'white','cyan','magenta','yellow'};
else
titlecolors={'white','green','red','blue'};
end

%let's try it make it auto in a loop for all channels
w=screensize(3);
h=screensize(4);
%channel=2;

for channel=1:length(dataChannels)

fig=figure('Position', [1, 1, w/2, h/1.5],'color','black'); 
for ii=1:length(conditions)
    set(gca,'FontName','Helvetica Neue','FontWeight','bold','color','black','box','off');
    
     titlestr = strcat(chlabel{channel});
        title(titlestr,'FontSize',fs+20,'FontName','Helvetica','FontWeight','bold','color',titlecolors{channel});

mpl(ii)=plot(ViR{ii}(:,1),ViR{ii}(:,channel+1),'LineWidth',5,'color',cmap(ii,:)); hold on; %this data has to have channel+1, the first value is the distance

fill([ViR{ii}(:,1); flipud(ViR{ii}(:,1))], [ViR{ii}(:,channel+1)-stdViR{ii}(:,channel+1); flipud(ViR{ii}(:,channel+1)+stdViR{ii}(:,channel+1))],mpl(ii).Color,'FaceAlpha',0.2,'EdgeColor','none'); hold on;
% 
xlim([0 colonysize/2]);
ylim([minY(channel+1) maxY(channel+1)]); 

   
axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
axis.XAxis.FontSize = fs+10;
axis.YAxis.FontSize = fs+10;
axis.YAxis.Exponent=0;
% axis.YTick=modYticks;
   xlabel('Distance from center to edge');
   ylabel('Mean Intensity (A.U.)');

set(legend,'color','black')
end

legend(mpl(:),conditions{:},'Location','Northeastoutside','TextColor','white');

filename=strcat(experiment,' Distance from center ',chlabel{channel}, ' in black');
fig.InvertHardcopy = 'off';
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
 saveas(fig, fullfile(graphsRD, filename), 'svg');
end

%in white, because I NEED IT tbh


for channel=1:length(dataChannels)
fig=figure('Position', [1, 1, w/2, h/1.5],'Color','white'); 
for ii=1:length(conditions)
grid minor
    set(gca,'FontName','Helvetica Neue','FontWeight','bold','box','off');
    
     titlestr = strcat(chlabel{channel});
        title(titlestr,'FontSize',fs+20,'FontName','Helvetica','FontWeight','bold','color',titlecolors{channel});

mpl(ii)=plot(ViR{ii}(:,1),ViR{ii}(:,channel+1),'LineWidth',10,'color',cmap(ii,:)); hold on; %this data has to have channel+1, the first value is the distance
fill([ViR{ii}(:,1); flipud(ViR{ii}(:,1))], [ViR{ii}(:,channel+1)-stdViR{ii}(:,channel+1); flipud(ViR{ii}(:,channel+1)+stdViR{ii}(:,channel+1))],mpl(ii).Color,'FaceAlpha',0.2,'EdgeColor','none'); hold on;

xlim([0 colonysize/2]);

ylim([0 maxY(channel+1)]); 

   
axis=gca;
axis.GridAlpha=0.8;
axis.XAxis.FontSize = fs+10;
axis.YAxis.FontSize = fs+10;
axis.YAxis.Exponent=0;
% axis.YTick=modYticks;
   xlabel('Distance from center to edge');
   ylabel('Mean Intensity (A.U.)');

%set(legend,'color','black')
end

legend(mpl(:),conditions{:},'Location','Northeastoutside');

filename=strcat(experiment,' Distance from center ',chlabel{channel}, ' in white');
fig.InvertHardcopy = 'off';
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
 saveas(fig, fullfile(graphsRD, filename), 'svg');
end

pause; close all;


%% 16B - same as before but not normalized


cmap=colormap('lines'); close; h=screensize(4);

if divisor>8
    cmap=colormap('turbo'); close;
    cidxs=round(linspace(10,length(colormap)-10,divisor));
    
    cidxs(2:2:end)=fliplr(cidxs(2:2:end));
    
    newcmap=cmap(cidxs,:); close;
    cmap=newcmap;
end

%if you want to group all the data/images: how many groups are?

grouping= 0 ; 

if grouping>0
if rem(divisor,grouping)>0
    disp('The number of groups you wish to generate is not divisible by the number of files, all will be treated as individual');
else

groups=divisor/grouping;

newcmap=zeros(divisor,3);
listofn=1:divisor;
listofdiv=rem(listofn,groups);
listofdiv=listofdiv<1;
listofn=listofn(listofdiv);

newcmap(1:groups,:)=cmap(1:groups,:);

for ii=listofn(1)+1:listofn(end)
   newcmap(ii,:)=newcmap(ii-groups,:);
end
cmap=newcmap;
end
end


rot=45; %for vertical text = 90, for horizontal text on top of bar = 0.
LW=1200/divisor;

if rot==0
    fs=20;
    %fs=30-15;
else 
    fs=30-10;
end

Iavelist=string(round(aIave,2));
Iavestdlist=string(round(stdaIave,2));
Iaveplusstd=strcat(Iavelist(:),' +/- ',Iavestdlist(:));
ytextp=aIave(:,2:4)+stdaIave(:,2:4);

% FIGUREEEEEEEEE for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP for GFP 
fig=figure('Position', [1, 1, w, h],'color','black'); 

for ii=1:divisor
Lb=line([ii ii], [0 aIave(ii,2)'],'LineWidth',LW,'Color',cmap(ii,:)); hold on;
Le=line([ii ii], [aIave(ii,2)' [aIave(ii,2)'+stdaIave(ii,2)']],'LineWidth',1,'Color','white'); hold on;
Leh=line([ii-0.1 ii+0.1], [[aIave(ii,2)'+stdaIave(ii,2)'] [aIave(ii,2)'+stdaIave(ii,2)']],'LineWidth',1,'Color','white'); hold on;

end
ymax=max(aIave(:,2)+stdaIave(:,2));

if ymax<1500
   ymax=1500;
   modYticks=linspace(0,ymax,5);
   ylim([0 ymax]);
else
    ymax=ceil(ymax/10)*10;
   % ymax=round(ymax,-2);
    modYticks=round(linspace(0,ymax,5),-2);
    ylim([0 ymax]);
end

    
xlim([0 divisor+1]);

     titlestr = strcat(chlabel{2});
        title(titlestr,'FontSize',fs+10,'FontName','Helvetica Neue','FontWeight','bold','color','green');
        
axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
axis.XAxis.FontSize = fs;
axis.YAxis.FontSize = fs+35;
%axis.YTick=modYticks; axis.YLim=[0 modYticks(end)];


set(gca,'FontName','Helvetica Neue','FontWeight','bold','xtick', 1:divisor,'xticklabel',conditions,'color','black','XTickLabelRotation',rot,...
'YGrid','off'); %used to have: 'YTickLabel',axis.YTick );
ylabeltext={'Mean Intensity (A.U.)'};

ylabelstuff=ylabel(ylabeltext,'FontSize',fs+70,'FontName','Helvetica Neue','FontWeight','bold'); 

filename=strcat(experiment,'L-Bars',chlabel{2}, 'Mean Intensity average own ymax in black not normalized');
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(graphsRD, filename), 'png');
savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');

%for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP for RFP 
fig=figure('Position', [1, 1, w, h],'color','black'); 

for ii=1:divisor
Lb=line([ii ii], [0 aIave(ii,3)'],'LineWidth',LW,'Color',cmap(ii,:)); hold on;
Le=line([ii ii], [aIave(ii,3)' [aIave(ii,3)'+stdaIave(ii,3)']],'LineWidth',1,'Color','w'); hold on;
Leh=line([ii-0.1 ii+0.1], [[aIave(ii,3)'+stdaIave(ii,3)'] [aIave(ii,3)'+stdaIave(ii,3)']],'LineWidth',1,'Color','w'); hold on;
end

ymin=min(aIave(:,3));
ymax=max(aIave(:,3)+stdaIave(:,3));

if ymax<10000 
    roundingnumber=-2;
end

if ymax<1000
   ymax=1000;
   roundingnumber=-2;
end


    ymax=ceil(ymax);
    ymax=round(ymax,roundingnumber);
    modYticks=round(linspace(0,ymax,5),-roundingnumber);
    ylim([0 ymax]);
    
xlim([0 divisor+1]);

     titlestr = strcat(chlabel{3});
        title(titlestr,'FontSize',fs+10,'FontName','Helvetica Neue','FontWeight','bold','color','red');
         
axis=gca
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
axis.XAxis.FontSize = fs;
axis.YAxis.FontSize = fs+35;
%axis.YTick=modYticks; axis.YLim=[0 modYticks(end)];
        
set(gca,'FontName','Helvetica Neue','FontWeight','bold','xtick', 1:divisor,'xticklabel',conditions,'color','black','XTickLabelRotation',rot,...
'YGrid','off'); %used to have: 'YTickLabel',axis.YTick );

ylabeltext={'Mean Intensity (A.U.)'};
ylabel(ylabeltext,'FontSize',fs+70,'FontName','Helvetica Neue','FontWeight','bold'); hold on;

filename=strcat(experiment,'L-Bars',chlabel{3}, 'Mean Intensity average own ymax in black not normalized');
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(graphsRD, filename), 'png');
savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');

%for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 %for CY5 

fig=figure('Position', [1, 1, w, h],'color','black'); 

for ii=1:divisor
Lb=line([ii ii], [0 aIave(ii,4)'],'LineWidth',LW,'Color',cmap(ii,:)); hold on;
Le=line([ii ii], [aIave(ii,4)' [aIave(ii,4)'+stdaIave(ii,4)']],'LineWidth',1,'Color','w'); hold on;
Leh=line([ii-0.1 ii+0.1], [[aIave(ii,4)'+stdaIave(ii,4)'] [aIave(ii,4)'+stdaIave(ii,4)']],'LineWidth',1,'Color','w'); hold on;
end

ymax=max(aIave(:,4)+stdaIave(:,4));
ymin=min(aIave(:,4));

if ymax<1500
   ymax=1500;
end


     ylim([0 ymax]);

xlim([0 divisor+1]);

    titlestr = strcat(chlabel{4});
        title(titlestr,'FontSize',fs+10,'FontName','Helvetica Neue','FontWeight','bold','color','blue');
       
axis=gca
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
axis.XAxis.FontSize = fs;
axis.YAxis.FontSize = fs+35;
%axis.YTick=modYticks; axis.YLim=[0 modYticks(end)];   

      
set(gca,'FontName','Helvetica Neue','FontWeight','bold','xtick', 1:divisor,'xticklabel',conditions,'color','black','XTickLabelRotation',rot,...
'YGrid','off'); %used to have: 'YTickLabel',axis.YTick );

ylabeltext={'Mean Intensity (A.U.)'};
ylabelstuff=ylabel(ylabeltext,'FontSize',fs+70,'FontName','Helvetica Neue','FontWeight','bold'); 


filename=strcat(experiment,'L-Bars',chlabel{4}, 'Mean Intensity average own ymax in black not normalized');
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(graphsRD, filename), 'png');
savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');

 pause;close all;

%% violin plots, ya la neta

newylims={}; eachGCA={}; vioPL=0; curatedAMIs={};
percentage=99;
minY=000;

rot=0;

fs=25; %change this depending on thenumber of conditions

if nfiles>6
    fs=15;
end

if length(conditions)>15
    howwidth=1.1;
    rot=45;
else
    howwidth=2;
end

LW=3; %line width
%colorsforthis=cmap(ii,:);
linecolor='k'; %this is the default, if you want to match body-color, use cmap 

for cc=1:length(dataChannels)
clear eachGCA; clear newylims; clear newxlims; clear maxY;

fig=figure('Position', [1, 1, w/howwidth, h/2]); 


for ii=1:length(conditions)

curatedAMIs{ii}(:,cc)=sort(AMIs{ii}(:,cc));
tempv=round((length(AMIs{ii}(:,cc))*(percentage/100)));
shrink=(ncells{ii}-tempv)/2;

if shrink<1
    shrink=3;
end

lengthfragment=length(shrink:ncells{ii}-shrink);
newcuratedAMIs{ii}(1:lengthfragment,cc)=curatedAMIs{ii}(shrink:ncells{ii}-shrink,cc);
    
if ii==nfiles
    vioPL=1;
end
    
if ii==1
%subtightplot(1,nfiles,1,[0.1 0]); %the first
subplot(1,nfiles,1);
%VPs=violin(AMIs{ii}(:,cc),'facecolor',cmap(ii,:),'plotlegend',0,'LineWidth',LW,'edgecolor',linecolor); %with all files
VPs=violin(newcuratedAMIs{ii}(:,cc),'facecolor',cmap(ii,:),'plotlegend',0,'LineWidth',LW,'edgecolor',linecolor); %with all files
%title(conditions{ii},'FontSize',30);
 xlabel(conditions{1},'FontSize',fs,'Fontweight','bold','Rotation',rot);
else 
 %   subtightplot(1,nfiles,ii,[0.1 0]); %the rest
    subplot(1,nfiles,ii);
    %VPs=violin(AMIs{ii}(:,cc),'facecolor',cmap(ii,:),'plotlegend',vioPL,'LineWidth',LW,'edgecolor',linecolor); %with all files
     VPs=violin(newcuratedAMIs{ii}(:,cc),'facecolor',cmap(ii,:),'plotlegend',vioPL,'LineWidth',LW,'edgecolor',linecolor); %with all files
   % title(conditions{ii},'FontSize',30);
     xlabel(conditions{ii},'FontSize',fs,'FontWeight','bold','Rotation',rot);
end

eachGCA{ii}=gca;
newylims{ii}=[eachGCA{ii}.YLim];
newxlims{ii}=[eachGCA{ii}.XLim];

maxY=max([newylims{:}]);
% 
end

vioPL=0;
for ii=1:length(conditions)
set(eachGCA{ii}, 'YLim',[minY maxY],'XLim',[min([newxlims{:}]) max([newxlims{:}])],'box','off');
 %xlabel(conditions{ii});
end

set(eachGCA{1},'box','off','TickLength',[0.01,0.025],'FontSize',fs,'XTick',[],'LineWidth',1.5); 
for ii=2:length(conditions)
    set(eachGCA{ii},'box','off','YTickLabel',[],'XTick',[],'LineWidth',1.5);
   
    eachGCA{ii}.YAxis.Color = [1 1 1];
end
set(gcf,'color','w');
sgtitle(chlabel{cc},'FontSize',45,'FontWeight','bold'); %super title cannot be use due to subtightplot

filename=strcat(experiment,' - ',chlabel{cc}, 'Violin plots');

saveas(fig, fullfile(graphsRD, filename), 'svg');
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');

end


%pause; close  all;

%VIOLIN PLOTS IN BLACK FROM HERE

linecolor='w'; facealpha=0.9;

for cc=1:length(dataChannels)
clear eachGCA; clear newylims; clear newxlims; clear maxY;

fig=figure('Position', [1, 1, w/howwidth, h/2],'color','black'); 


for ii=1:length(conditions)

if ii==nfiles
    vioPL=1;
end
    
if ii==1

subplot(1,nfiles,1);
VPs=violin(newcuratedAMIs{ii}(:,cc),'facecolor',cmap(ii,:),'plotlegend',0,'LineWidth',LW,'edgecolor',cmap(ii,:),'facealpha',facealpha,'Boxcolor','k'); %use one letter for boxcolor

 xlabel(conditions{1},'FontSize',fs,'Fontweight','bold','Rotation',rot);
else 

    subplot(1,nfiles,ii);
    VPs=violin(newcuratedAMIs{ii}(:,cc),'facecolor',cmap(ii,:),'plotlegend',vioPL,'LineWidth',LW,'edgecolor',cmap(ii,:),'facealpha',facealpha,'Boxcolor','k'); %use one letter for boxcolor
 
     xlabel(conditions{ii},'FontSize',fs,'FontWeight','bold','Rotation',rot);
end

eachGCA{ii}=gca;
newylims{ii}=[eachGCA{ii}.YLim];
newxlims{ii}=[eachGCA{ii}.XLim];

maxY=max([newylims{:}]);
% 
end

vioPL=0;
for ii=1:length(conditions)
set(eachGCA{ii},'YLim',[minY maxY],'XLim',[min([newxlims{:}]) max([newxlims{:}])],'box','off','color','k','Ycolor','w','Xcolor','w');
 %xlabel(conditions{ii});
end

set(eachGCA{1},'box','off','TickLength',[0.01,0.025],'FontSize',fs,'XTick',[],'LineWidth',1.5); 
for ii=2:length(conditions)
    set(eachGCA{ii},'box','off','YTickLabel',[],'XTick',[],'LineWidth',1.5);
    eachGCA{ii}.YAxis.Color = [0 0 0];
end
set(gcf,'color','k');

sgtitle(chlabel{cc},'FontSize',45,'FontWeight','bold','color','w'); %super title cannot be use due to subtightplot

filename=strcat(experiment,' - ',chlabel{cc}, 'Violin plots in black');
fig.InvertHardcopy = 'off';
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
 

end

pause; close all;

%% 18 - Histograms


fs=30;

nfsortedAD={}; minsortedAD={}; maxsortedAD={};
for ii=1:length(conditions)
for jj=1:length(chlabel)
   nfsortedAD{ii}{jj}=newcuratedAMIs{ii}(:,jj);
   minsortedAD{ii}{jj}= min(newcuratedAMIs{ii}(:,jj));
   maxsortedAD{ii}{jj}= max(newcuratedAMIs{ii}(:,jj));
end
end

clear maxyh; clear minyh;
for ii=1:length(conditions)

maxyh(ii,:)=[maxsortedAD{ii}{:}];
minyh(ii,:)=[minsortedAD{ii}{:}];

end
minyh=min(minyh);

if height(maxyh)>1
maxyh=max(maxyh);
    if any(maxyh(2:4)<maxyh(1))
maxyh(maxyh<maxyh(1))=maxyh(1);
    end
end
%this will only work for 4 channels ! 

fig=figure('Position', [1, 1, w, h/2]);

subplot(1,length(dataChannels),1);
for ii=1:length(conditions)


[NN,edges]=histcounts(nfsortedAD{ii}{1}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;

 edges=[edges(1) edges];
 NN=[0 NN];
 
 plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;

xlabel('Intensity','FontSize',fs);
ylabel('Probability Density Estimate','FontSize',fs);

legend(conditions);
hold on;


 titlestr = [ chlabel{1} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold');
       hold on;
end
ylimit=max(NN);

subplot(1,length(dataChannels),2);
for ii=1:length(conditions)

%subplot_tight(n,m,2) %GFP

[NN,edges]=histcounts(nfsortedAD{ii}{2}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;

 edges=[edges(1) edges];
 NN=[0 NN];

plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;
xlabel('Intensity','FontSize',fs);
ylabel('Probability Density Estimate','FontSize',fs);

legend(conditions);
hold on;

 titlestr = [ chlabel{2} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold');
       hold on;
end

subplot(1,length(dataChannels),3);
for ii=1:length(conditions)


[NN,edges]=histcounts(nfsortedAD{ii}{3}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;


 edges=[edges(1) edges];
 NN=[0 NN];



plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;
xlabel('Intensity','FontSize',fs); 
ylabel('Probability Density Estimate','FontSize',fs);

legend(conditions);
hold on;

 titlestr = [ chlabel{3} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold');
       hold on;
end

subplot(1,length(dataChannels),4);
for ii=1:length(conditions)

%subplot_tight(n,m,4) %CY5

[NN,edges]=histcounts(nfsortedAD{ii}{4}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;


 edges=[edges(1) edges];
 NN=[0 NN];


plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;

xlabel('Intensity','FontSize',fs);
ylabel('Probability Density Estimate','FontSize',fs);

legend(conditions);
hold on;

 titlestr = [ chlabel{4} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold');
       hold on;
end

filename=strcat(experiment,'histograms');
saveas(fig, fullfile(graphsRD, filename), 'svg');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
 saveas(fig, fullfile(graphsRD, filename), 'png');

fig=figure('Position', [1, 1, w, h/2],'color','black');


subplot(1,length(dataChannels),1);
for ii=1:length(conditions)
%subplot_tight(n,m,1) %DAPI

[NN,edges]=histcounts(nfsortedAD{ii}{1}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;

edges=[edges(1) edges];
 NN=[0 NN];

plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;
set(gca,'FontSize',fs-10,'color','black');
axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;


xlabel('Intensity','FontSize',fs,'color','white');
ylabel('PDF','FontSize',fs,'color','white');

legend(conditions,'TextColor','white');
hold on;


 titlestr = [ chlabel{1} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
       hold on;
end

subplot(1,length(dataChannels),2);
for ii=1:length(conditions)

%subplot_tight(n,m,2) %GFP

[NN,edges]=histcounts(nfsortedAD{ii}{2}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;

edges=[edges(1) edges];
 NN=[0 NN];

plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;
set(gca,'FontSize',fs-10,'color','black');
axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;

xlabel('Intensity','FontSize',fs);
ylabel('PDF','FontSize',fs);

legend(conditions,'TextColor','white');
hold on;

 titlestr = [ chlabel{2} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
       hold on;
end

subplot(1,length(dataChannels),3);
for ii=1:length(conditions)

%subplot_tight(n,m,3) %RFP 

[NN,edges]=histcounts(nfsortedAD{ii}{3}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;

edges=[edges(1) edges];
 NN=[0 NN];

plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;
set(gca,'FontSize',fs-10,'color','black');
axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;


xlabel('Intensity','FontSize',fs); 
ylabel('PDF','FontSize',fs);

legend(conditions,'TextColor','white');
hold on;

 titlestr = [ chlabel{3} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold','color','w');
       hold on;
end

subplot(1,length(dataChannels),4);
for ii=1:length(conditions)

%subplot_tight(n,m,4) %CY5

[NN,edges]=histcounts(nfsortedAD{ii}{4}, 'Normalization','pdf');
edges=edges(2:end)-(edges(2)-edges(1))/2;

edges=[edges(1) edges];
 NN=[0 NN];

plot(edges,NN,'LineWidth',3,'Color',cmap(ii,:)); hold on;
set(gca,'FontSize',fs-10,'color','black');
axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;


xlabel('Intensity','FontSize',fs);
ylabel('PDF','FontSize',fs);

legend(conditions,'TextColor','white');
hold on;

 titlestr = [ chlabel{4} ];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
        labelstr = ['\color{blue}'   chlabel{1}...
                    '\color{green}' chlabel{2}...
                     '\color{red}' chlabel{3}...
                    '\color{magenta}' chlabel{4}]; %Magenta RGB 229, 9, 127 = [1 0 1]
       text(margin,  size(tiffs{1},2)-1*margin , labelstr,'FontSize',fs,'FontWeight','bold','color','white');
       hold on;
end

filename=strcat(experiment,'histograms im black');
 fig.InvertHardcopy = 'off';
 saveas(fig, fullfile(graphsRD, filename), 'png');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');

 pause; close all;

  %% 19 - Scatter plots cell Density ! 


limitlines= 0;  
percentage=99.5;
channeltocurate=3;

%
treatcondition = 2;
controlcondition = 1;
channelformin = 1;
channelformax = [4,8];


if limitlines==1
lineslimits=[0,0,0,0];
end

%
cs=ceil(w*h/100000);
fs=30;
curatedAMIs={}; sccuratedAMIs={};

for ii=1:length(conditions)
curatedAMIs{ii}=sortrows(AMIs{ii},channeltocurate);
tempv=round((length(AMIs{ii})*(percentage/100)));
shrink=round((ncells{ii}-tempv)/2);

if shrink<1
    shrink=3;
end

lengthfragment=length(shrink:ncells{ii}-shrink);
sccuratedAMIs{ii}=curatedAMIs{ii}(shrink:ncells{ii}-shrink,:);
end

n=round(nfiles/2);
m=round(nfiles/n);

if n<3
    n=2; m=2;
end

 
clear maxyh; clear minyh;
for ii=1:length(conditions)

maxyh(ii,:)=[maxsortedAD{ii}{:}];
minyh(ii,:)=[minsortedAD{ii}{:}];

end

if height(maxyh)>1
maxyh=max(maxyh);
minyh=min(minyh);
end


R2o=cell(length(conditions));

fig=figure('Position', [1, 1, w, h]);

for ii=1:length(chlabel)-1
 for jj=ii+1:length(chlabel)

 

    for kk=1:length(conditions)
        subplot(m,n,kk);
      
scpl(kk)=scatter_kde(sccuratedAMIs{kk}(:,ii),sccuratedAMIs{kk}(:,jj),'filled','MarkerSize',cs);
xlabel(string(chlabel(ii))+' - Mean Intensity (A.U.)','FontSize',fs);
ylabel(string(chlabel(jj))+' - Mean Intensity (A.U.)','FontSize',fs);

if maxyh(jj)<1200 
    maxyh(jj)=1200; 
end

if maxyh(ii)<1200
    maxyh(ii)=1200;
end

% 
 xlim([minyh(ii) maxyh(ii)]);
 ylim([minyh(jj) maxyh(jj)]);


if limitlines==1
xline(lineslimits(ii),'--r','LineWidth',2);
yline(lineslimits(jj),'--r','LineWidth',2);
end

LRstuff=fitlm(AMIs{kk}(:,ii),AMIs{kk}(:,jj));
R2o{kk}=LRstuff.Rsquared.Ordinary;

titlestr = [conditions{kk} string(round(R2o{kk},3))];
title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
 set(gca,'FontSize',fs-10,'FontName','Arial');
 
if kk==nfiles
     cb=colorbar('EastOutside');   
     cb.Label.String = 'Probability Density Estimate';
end
 
%to check what's doing step by step
filename=strcat(experiment,chlabel{ii},chlabel{jj});
  
           saveas(fig, fullfile(scatterplotsRD, filename), 'png');
          saveas(fig, fullfile(scatterplotsRD, filename), 'svg');
 savefig(fig,fullfile(scatterplotsRD, [filename '.fig']),'compact');
           
           
    end
   
 end
 
end  


%in black
fig=figure('Position', [1, 1, w, h],'color','black');

for ii=1:length(chlabel)-1
 for jj=ii+1:length(chlabel)

 

    for kk=1:length(conditions)
        subplot(m,n,kk);
       

scatter_kde(AMIs{kk}(:,ii),AMIs{kk}(:,jj),'filled','MarkerSize',cs);

xlabel(string(chlabel(ii))+' - Mean Intensity (A.U.)','FontSize',fs);
ylabel(string(chlabel(jj))+' - Mean Intensity (A.U.)','FontSize',fs);

xlim([minyh(ii) maxyh(ii)]);
ylim([minyh(jj) maxyh(jj)]);

LRstuff=fitlm(AMIs{kk}(:,ii),AMIs{kk}(:,jj));
R2o{kk}=LRstuff.Rsquared.Ordinary;

titlestr = [conditions{kk} string(round(R2o{kk},3))];
title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
 set(gca,'FontSize',fs-10,'color','black');

axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
 
 
if kk==nfiles
     cb=colorbar('EastOutside');   
     cb.Label.String = 'Probability Density Estimate';
     cb.Color=[1 1 1];
     cb.Box='off';
end
 

%to check what's doing step by step
filename=strcat(experiment,chlabel{ii},chlabel{jj},' in black');
            fig.InvertHardcopy = 'off';
           saveas(fig, fullfile(graphsRD, filename), 'png');
            saveas(fig, fullfile(graphsRD, filename), 'svg');
 savefig(fig,fullfile(graphsRD, [filename '.fig']),'compact');
    end
 end
end  

close all;
 

%% 20-B Scatter with Z. Only with ilastik tables


fig=figure('Position', [1, 1, w, h]);
colormap(turbo);
for ii=1:channels-1
 for jj=ii+1:channels
     
if maxyh(ii)<1000
    maxyh(ii)=1000;
end
  
if maxyh(jj)<1000
    maxyh(jj)=1000;
end     

    for kk=1:groups
        subplot(m,n,kk);

scatter(gd.nuc.Alldata{kk}(:,ii+7),gd.nuc.Alldata{kk}(:,jj+7),cs,round(gd.nuc.Alldata{kk}(:,7)),'filled','MarkerEdgeAlpha',.1,'MarkerFaceAlpha',0.6);

xlabel(string(chlabel(ii))+' - Mean Intensity (A.U.)','FontSize',fs);
ylabel(string(chlabel(jj))+' - Mean Intensity (A.U.)','FontSize',fs);


xlim([500 maxyh(ii)]);
ylim([500 maxyh(jj)]);


caxis([1 maxZ]);

titlestr = [conditions{kk}];
title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
 set(gca,'FontSize',fs-10);
 
if kk==groups
     cb=colorbar('EastOutside');   
     cb.Label.String = 'Z-centroid in µm';
end
filename=strcat(experiment,chlabel{ii},chlabel{jj},'scatterplot-Z');  
           saveas(fig, fullfile(scatterplotsRD, filename), 'png');
           saveas(fig, fullfile(scatterplotsRD, filename), 'svg');
 savefig(fig,fullfile(scatterplotsRD, [filename '.fig']),'compact');
    end  
 end
end  


fig=figure('Position', [1, 1, w, h],'color','black');
colormap(turbo);
for ii=1:length(chlabel)-1
 for jj=ii+1:length(chlabel)

       
if maxyh(ii)<1000
    maxyh(ii)=1000;
end
  
if maxyh(jj)<1000
    maxyh(jj)=1000;
end     
    for kk=1:groups
        subplot(m,n,kk);

scatter(gd.nuc.Alldata{kk}(:,ii+7),gd.nuc.Alldata{kk}(:,jj+7),cs,round(gd.nuc.Alldata{kk}(:,7)),'filled','MarkerEdgeAlpha',.1,'MarkerFaceAlpha',0.7);

xlabel(string(chlabel(ii))+' - Mean Intensity (A.U.)','FontSize',fs);
ylabel(string(chlabel(jj))+' - Mean Intensity (A.U.)','FontSize',fs);

xlim([minyh(ii)-50 maxyh(ii)]);
ylim([minyh(jj)-50 maxyh(jj)]);


caxis([1 maxZ]);


titlestr = [conditions{kk}];
title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');
 set(gca,'FontSize',fs-10,'color','black');

axis=gca;
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
 
if kk==groups
     cb=colorbar('EastOutside');   
     cb.Label.String = 'Z-centroid in µm';
     cb.Color=[1 1 1];
     cb.Box='off';
end
 
filename=strcat(experiment,chlabel{ii},chlabel{jj},'-scatterplot Z in black');
            fig.InvertHardcopy = 'off';
           saveas(fig, fullfile(scatterplotsRD, filename), 'png');
           savefig(fig,fullfile(scatterplotsRD, [filename '.fig']),'compact');
    end
 end
end  


pause;
close all;


%% check single overlay                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq%% Check the single object MI overlay
 
 % need to correct to pull up the original data !!!!! is pulling the
 % grouped one !!!

%which channel do you wish to view?
channel=4;

%for which file or group(if groupeD)?
file=1;


%from here should be automatic:
fig=figure(1); 
titlestr = chlabel{channel};
        
%for ii=1:nfiles
for ii=file;    
%subplot_tight(n.*2,m,ii)


if channel==1;
tempimg=cat(3,zeros(size(adjimgs{ii}{1})),zeros(size(adjimgs{ii}{1})),adjimgs{ii}{1});
end

if channel==2;
tempimg=cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},zeros(size(adjimgs{ii}{2})));
end

if channel==3;
tempimg=cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),zeros(size(adjimgs{ii}{3})));
end

if channel==4;
tempimg=cat(3,zeros(size(adjimgs{ii}{4})),zeros(size(adjimgs{ii}{4})),adjimgs{ii}{4});
end


imshow(imfuse(bordermasks{ii},tempimg,'blend'),[]); hold on;
title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','red');
%when data has been grouped
if groupdata==1
text(oricentroids{ii}(:,1),oricentroids{ii}(:,2),string(round(oriAMIs{ii}(:,channel),2)),'color','white','FontSize',7);
else
text(centroids{ii}(:,1),centroids{ii}(:,2),string(round(AMIs{ii}(:,channel),2)),'color','white','FontSize',7);
end

hold on;


           if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
           end
           
           
if groupdata==1           
    figure(2); title(originalconditions{ii});
else
    figure(2); title(conditions{ii});
end
          
if channel==1;
imshow(cat(3,zeros(size(adjimgs{ii}{1})),zeros(size(adjimgs{ii}{1})),adjimgs{ii}{1})); hold on;
end

if channel==2;
imshow(cat(3,zeros(size(adjimgs{ii}{2})),adjimgs{ii}{2},zeros(size(adjimgs{ii}{2})))); hold on;
end

if channel==3;
imshow(cat(3,adjimgs{ii}{3},zeros(size(adjimgs{ii}{3})),zeros(size(adjimgs{ii}{3})))); hold on;
end

if channel==4;
imshow(cat(3,zeros(size(adjimgs{ii}{4})),zeros(size(adjimgs{ii}{4})),adjimgs{ii}{4})); hold on;
end
 hold on;
 
%text(centroids{ii}(:,1),centroids{ii}(:,2),string(round(AMIs{ii}(:,channel))),'color','white');
           if SD20x==1
       yx=size(adjimgs{ii}{1}); %remember that for images is yx
       xbar=[yx(2)*.05 yx(2)*.05+humbar];
       ybar=[yx(1)*.95 yx(1)*.95]; 
       hold on
       plot(xbar,ybar,'color','white','linewidth',4);
           end

     
end
%
pause;
close all;
 
%% Calculate Area % to Area of DAPI


%you'll need a minimum number for this one though.

%min1=100; %DAPI from histogram comparison . 
 min2=3000; %SOX2 from histogram comparison . 
 min3=500; %BRA from histogram comparison . 
 min4=500; %CDX2 from histogram comparison . 
% 

if groupdata == 1
   
ori.aboves=cell(1,length(originalconditions));
for ii=1:length(originalconditions)
ori.aboves{ii}(:,2)=oriAMIs{ii}(:,2)>min2; %GFP
ori.aboves{ii}(:,3)=oriAMIs{ii}(:,3)>min3; %RFP
ori.aboves{ii}(:,4)=oriAMIs{ii}(:,4)>min4; %CFP
end
ori.tracedIDs=cell(1,length(originalconditions));
for ii=1:length(originalconditions)
     ori.tracedIDs{ii}{2}=find(ori.aboves{ii}(:,2)==1);
     ori.tracedIDs{ii}{3}=find(ori.aboves{ii}(:,3)==1);
     ori.tracedIDs{ii}{4}=find(ori.aboves{ii}(:,4)==1);
end    

ori.sumsA=cell(1,length(originalconditions)); % a lo mejor tienes que uncomment this
for ii=1:length(originalconditions)
for jj=1:length(dataChannels)
ori.sumsA{ii}{jj}=0;
end
end    

for ii=1:length(originalconditions)
for jj=2:length(dataChannels)
tempidx=oriAreas{ii};
for kk=1:length(ori.tracedIDs{ii}{jj})
ori.sumsA{ii}{jj}(:,kk)=tempidx(ori.tracedIDs{ii}{jj}(kk));
end
end
end    
    
for ii=1:length(originalconditions)
    for jj=1:length(dataChannels)
    ori.sumsA{ii}{jj}=sum(ori.sumsA{ii}{jj});
end
end

ori.coveredA=cell(1,length(originalconditions));
for ii=1:length(originalconditions)
for    jj=1:length(dataChannels)
    ori.coveredA{ii}{jj}=ori.sumsA{ii}{jj}./oriTotalArea{ii}.*100;
end
end

g.coveredA={}; 
for gg = 1:groups
for ssgg=1:length(groupsindex{gg})

if ssgg==1
    g.coveredA{gg}=ori.coveredA{groupsindex{gg}{ssgg}};
else
    g.coveredA{gg}=cat(1,g.coveredA{gg},ori.coveredA{groupsindex{gg}{ssgg}});
end
end
end

% 

for ii=1:groups
g.mCA(ii,1)=mean([g.coveredA{ii}{:,2}]); g.stdmCA(ii,1)=std([g.coveredA{ii}{:,2}]);
g.mCA(ii,2)=mean([g.coveredA{ii}{:,3}]); g.stdmCA(ii,2)=std([g.coveredA{ii}{:,3}]);
g.mCA(ii,3)=mean([g.coveredA{ii}{:,4}]); g.stdmCA(ii,3)=std([g.coveredA{ii}{:,4}]);
end

mCA=g.mCA;
coveredA=g.coveredA;


else
%With M.I.

aboves=cell(1,length(conditions));
for ii=1:length(conditions)
aboves{ii}(:,2)=AMIs{ii}(:,2)>min2; %GFP
aboves{ii}(:,3)=AMIs{ii}(:,3)>min3; %RFP
aboves{ii}(:,4)=AMIs{ii}(:,4)>min4; %CFP
end
%group

%trace the ID of those objects that are above the minimum objects

tracedIDs=cell(1,length(conditions));
for ii=1:length(conditions)
     tracedIDs{ii}{2}=find(aboves{ii}(:,2)==1);
     tracedIDs{ii}{3}=find(aboves{ii}(:,3)==1);
     tracedIDs{ii}{4}=find(aboves{ii}(:,4)==1);
end

sumsA=cell(1,length(conditions));
for ii=1:divisor
for jj=1:length(dataChannels)
sumsA{ii}{jj}=0;
end
end

for ii=1:length(conditions)
for jj=2:length(dataChannels)
tempidx=Areas{ii};
for kk=1:length(tracedIDs{ii}{jj})
sumsA{ii}{jj}(:,kk)=tempidx(tracedIDs{ii}{jj}(kk));
end
end
end

     
for ii=1:length(conditions)
    for jj=1:length(dataChannels)
    sumsA{ii}{jj}=sum(sumsA{ii}{jj});
end
end


coveredA=cell(1,length(conditions));
for ii=1:length(conditions)
for    jj=1:length(dataChannels)
    coveredA{ii}{jj}=sumsA{ii}{jj}./TotalArea{ii}.*100;
end
end

for ii=1:length(conditions)
mCA(ii,1)=coveredA{ii}{2};
mCA(ii,2)=coveredA{ii}{3};
mCA(ii,3)=coveredA{ii}{4};
end


end 

writetable(array2table(mCA),fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'meancoveredA.xlsx')));

if groupdata==1
writetable(array2table(cat(1,mCA,g.stdmCA)),fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'mthenstdmCoveredA.xlsx')));
end


txtArea=string(mCA);
fid = fopen(fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'coveredArearesults.txt')), 'wt');
fprintf(fid, '%s\n', txtArea);
fclose(fid);

% PLOTs !!!!! 

ytextp=string(mCA+1);

cxbar=categorical(chlabel(2:4));
cxbar=reordercats(cxbar,chlabel(2:4)); %these names must match the labels from the first section

% FIGUREEEEEEEEE in white
fig=figure('Position', [1, 1, w/1.5, h/1.5]); 

b=bar(cxbar,mCA); hold on
  
Xb=[];
for bb=1:length(b)
    Xb=[Xb;[1:3]'+b(bb).XOffset];
end
sortedXb=reshape(sort(Xb),length(b),length(dataChannels)-1);


hold off


ylabel('%','FontSize',fs,'FontWeight','bold');
legend(conditions,'Location','northeastoutside');
 titlestr = ['Area covered to Area of DAPI  ' experiment];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
  
        
  
text(sort(Xb(:)),mCA(:)+fs.*.15,string([round(mCA(:))]),'Vert','top','Horiz','center','FontWeight','normal','FontSize',fs/2+divisor,'Rotation',rot);
        
        
set(gca,'FontSize',30, 'xticklabel', chlabel(2:4),'box','off')

filename=strcat(experiment,'Area Bars grouped');
saveas(fig, fullfile(graphsRD, filename), 'png');

%FIGUREEEEEE in BLACK !!!!!!!

fig=figure('Position', [1, 1, w/1.5, h/1.5],'color', 'black'); 
b=bar(cxbar,[mCA(:,1)'; mCA(:,2)'; mCA(:,3)']); 
ylabel('%','FontSize',fs,'FontWeight','bold');
legend(conditions,'Location','northeastoutside','TextColor','white');
 
     titlestr = ['Area covered to Area of DAPI  ' experiment];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none','color','white');   

              
Xb=[];
for bb=1:length(b)
    Xb=[Xb;[1:3]'+b(bb).XOffset];
    b(bb).FaceColor=cmap(bb,:);
end
text(sort(Xb(:)),mCA(:)+fs.*.15,string([round(mCA(:))]),'Vert','top','Horiz','center','FontWeight','normal','FontSize',fs/2+divisor,'Rotation',rot,'color','white');
    
axis=gca
axis.XColor='w';
axis.YColor='w';
axis.GridAlpha=0.8;
        
set(gca,'FontSize',30, 'xticklabel', chlabel(2:4),'color','black','box','off')
set(legend,'color','black')


filename=strcat(experiment,'Area Bars grouped in black');
fig.InvertHardcopy = 'off';
saveas(fig, fullfile(graphsRD, filename), 'png');


writematrix(mCA,fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'Areacovered.xlsx')));
fid = fopen(fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'Areacovered.txt')), 'wt');
fprintf(fid, '%s\n', mCA);
fclose(fid);

%% Co-expression of M.I.

%how many channels co-expression do you want? 2 o 3?


%min1=100; %DAPI from histogram comparison . (need to automate this)
 min2=2000; %SOX2 from histogram comparison . (need to automate this)
 min3=500; %BRA from histogram comparison . (need to automate this)
 min4=500; %CDX2 from histogram comparison . (need to automate this)

chtocheck=[2,3]; %indicate which channels you which to check their co-expression in ascendant order
nchtocheck=length(chtocheck); %length of channels to check

%THIS WILL BE BASED ON THE MINS YOU PREIVOUSLY INDICATED

mins=[min2,min3,min4];

tempdata={};
%cells above threshold in each condition
for ii=1:length(conditions)
tempdata{ii}(:,1)=AMIs{ii}(:,chtocheck(1));
tempdata{ii}(:,2)=AMIs{ii}(:,chtocheck(2));
end

above={};
for ii=1:length(conditions)
above{ii}=tempdata{ii}>[mins(chtocheck(1)-1),mins(chtocheck(2)-1)];
%above{ii}=tempdata{ii}>[min2,min3];
end

coexp={};
doublep={};
for ii=1:length(conditions)
    coexp{ii}=0;
    for jj=1:length(above{ii})
    doublep=sum(above{ii}(jj,:));
    if doublep>1
        coexp{ii}=1+coexp{ii};
    end
    end
end

yield={};
for ii=1:length(conditions)
    yield{2,ii}=conditions{ii};
    yield{3,ii}=round((coexp{ii}/ncells{ii})*100,1);
   
    yield{1,ii}=([chlabel{chtocheck(2)},'-', chlabel{chtocheck(1)}]);
end

disp('Co-expression of:')
disp(chlabel(chtocheck));
yield


writetable(cell2table(yield),fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'yield.xlsx')));

txtyield=string(yield);

fid = fopen(fullfile(resultsDir, strcat(experiment,'-',[chlabel{chtocheck}],'-','yield.txt')), 'wt');
fprintf(fid, '%s\n', txtyield);
fclose(fid);


%% 20 - Results

result=cell(length(dataChannels).*2+3,length(conditions));

for ii=1:length(conditions)
    for jj=1:length(dataChannels)
        
result{1,ii}=string(conditions{ii});
result{2,ii}="Mean Intensity";
result{7,ii}="Std";
        
result{jj+2,ii}=normalizedIave{ii}{jj};
result{jj+7,ii}=normalizedstdIave{ii}{jj};
    end
end
writetable(cell2table(result),fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'nIaveresults.xlsx')));
result=string(result);

fid = fopen(fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'nIaveresults.txt')), 'wt');
fprintf(fid, '%s\n', result);
fclose(fid);


result=cell(length(dataChannels).*2+3,length(conditions));

for ii=1:length(conditions)
    for jj=1:length(dataChannels)
        
result{1,ii}=string(conditions{ii});
result{2,ii}="Mean Intensity";
result{7,ii}="Std";
        
result{jj+2,ii}=Iave{ii}{jj};
result{jj+7,ii}=Iave{ii}{jj};
    end
end
writetable(cell2table(result),fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'Iaveresults.xlsx')));
result=string(result);

fid = fopen(fullfile(resultsDir, strcat(experiment,string([chlabel{:}]),'Iaveresults.txt')), 'wt');
fprintf(fid, '%s\n', result);
fclose(fid);



%% Save the data

filename=[experiment 'Alldata.mat'];
save([resultsDir '/' filename],'Alldata','-v7.3');


