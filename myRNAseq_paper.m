%% Find mat dataset and load it

clc; close all; clear; clear memory;

dataDir= '//RNAseq/DataSets/';
addpath(genpath('~/Documents/GitHub/stemcells')); 

files= dir(fullfile(dataDir,'*.mat'));

ds=files(2).name;
load(ds);

conditions={'1xn', '3xn', '10xn', '1xP', '3xP', '10xP', 'Control'};
screensize=get(0,'Screensize');


resultsDir = fullfile(dataDir,'/heatmaps');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

diffDir=fullfile(dataDir,'/differential');
if ~exist(diffDir,'dir')
    mkdir(diffDir,'dir');
end

ds='Sox2Micropattern.mat'; %this is just the name
bmp='BMP48h.mat';

%% for looking any gene

%example: getExpression(dataset, 'gene name', 'that starts with this letters') => getExpression(ds,'ZIC','StartsWith');

%% Differential expressions


%example: getExpression(dataset, 'gene name') => getExpression(ds,'ZIC');
clear idxs; clear dERd; clear dERn; clear dERids;

[dERd dERn dERids]=differentialExpression(ds,{1:3,7},5,10);  %data set, which data to compare, fold change, expression level

%sort alphabetically
  [sorteddERn idxs]=sort(dERn);
  dERd=dERd(idxs,:);
  dERn=sorteddERn;
 
  length(dERn)
 

%% Heatmap of differential

close all;

fs=30;
sizeXY= [screensize(3)/2.5, screensize(4)/2];

%nameforthis='Pluripotent markers';
if length(dERn)>50
    fs=16;
    sizeXY=[screensize(3)/2, screensize(4)];
end

fig=figure('Position',[1,1,sizeXY(1),sizeXY(2)],'color','black');
hh=heatmap(round(dERd),'Colormap',hot,'CellLabelColor','none','FontColor','white','FontSize',fs);
hhSt=struct(hh).Heatmap;
hhStGrid=struct(hhSt).Grid;
hhStGrid.ColorData=uint8([50; 50; 50; 256]);
title('Differential Expression');

hh.YDisplayLabels=dERn;
hh.XDisplayLabels=conditions;

axs=struct(gca);
cb=axs.Colorbar;

axs.XAxis.FontSize=30;
axs.YAxis.FontSize=fs;

if max(cb.Limits)>200
 cb.Ticks=[0,ceil((max(hh.ColorData(:))-min(hh.ColorData(:)))/2/100)*100, ceil(max(hh.ColorData(:))/100)*100];
 cb.Limits=[0 ceil(max(hh.ColorData(:))/100)*100];
hh.ColorLimits=cb.Limits;
end

if max(cb.Limits)<5
    cb.Ticks=[0, 5, 10];
    cb.Limits=[0 10];
    hh.ColorLimits=cb.Limits;
end

filename=strcat('all up to control',' heatmap differential',' in black');
fig.InvertHardcopy = 'off';
% 
  saveas(fig, fullfile(diffDir, filename), 'png');
  savefig(fig,fullfile(diffDir, [filename '.fig']),'compact');


%% import list from excel and do a heatmap

tpaperlist=['human gastrula cluster genes paper.xlsx'];
tpaperlist=readtable(tpaperlist);

paperlistheaders=tpaperlist.Properties.VariableNames;
paperlist=table2cell(tpaperlist);

%% list and heatmap
%build the list with values
clear genedata; clear genenames;

%pluripotent
%genelist={'SOX2','POU5F1','NANOG','MYCN', 'KLF4','E-CAD',}; %pluripotent
%%Epithelial at the NSB 8.5-9.5
%genelist={'CDH1','LAMB1','KRT10','CLDN3','CLDN4','CRB3','TJP1','TJP2'}; nameforthis= 'Epithelial markers no KRT18';

%Mesenchymal markers 
%genelist={'CDH2','SNAI1','VIM','LEF1','ZEB1','ZEB2'}; nameforthis='Mesenchymal markers no FN1';

%from Briscoe paper: 
%genelist={'TBXT','NKX1-2','TBX6','CDX1','CDX2','CDX4','FGF8','FGF17','WNT8A','WNT3A'}; nameforthis='e8.5 NMPs with no SOX2';

%neural progenitros
%genelist={'SOX2','SOX1','SOX3','NKX1-2','ZIC1','ZIC2','NESTIN','PAX6','GFAP','RDH11','RDH10','NEUROG1','NEUROD1','NKX2-2'}; nameforthis='Neural progenitors';

%mesoderm somitic progenitors
%genelist={'TBXT','TBX6','MSGN1','MYOD1','MEOX1','PAX7','LFNG','HES7','RIPPLY2','CYP26A1'}; nameforthis='Somite progenitors';

%all endoderm
%genelist={'FOXA2','SOX17','EOMES','CXCR4','CER1','DKK1','LEFTY1','SOX7','GATA4','GATA6','HHEX','PDGFRA','HNF4A'}; nameforthis='All Endoderm markers';

%general for sevearl
genelist={'SOX2','POU5F1','NANOG','MYCN', 'KLF4','E-CAD','SNAI1','VIM','LEF1','TWIST1','ZEB1','ZEB2','TBXT','NKX1-2','TBX6','CDX2','FGF8','FGF17','WNT8A','WNT3A'...
    'SOX2','SOX1','SOX3','NKX1-2','ZIC1','ZIC2','NESTIN','PAX6','FOXA2','SOX17','EOMES','CXCR4','CER1','DKK1','LEFTY1','SOX7','GATA4','GATA6','HHEX','PDGFRA','HNF4A'};
 
%brain markers
%genelist={'FOXG1','EN1','EN2','PAX5','LMX1A','TH'}; nameforthis='Brain markers';

%HOX genes (remember to put this on normal mode)
%genelist={'HOX','PHOX'}; nameforthis='HOX genes';

%AP axis genes
%genelist={'OTX2','SIX3','GBX2','CDX1','CDX2','CDX4','HOXB1','HOXB2','HOXB3','HOXB4'}; nameforthis='AP axis markers';

%WNT signaling (check this one in normal)
%genelist={'WNT','FZD','LRP1','LRP2','LRP3','LRP4','LRP5','LRP6','DKK','AXIN','BCAT'}; nameforthis='WNT signaling';

%BMP signaling (check with normal)
%genelist={'BMP','SMAD','NOGGIN'}; nameforthis='BMP signaling';

%TGF-B signaling (check with normal)
%genelist={'TFGB','ACV','NODAL','LEFTY','INHA','INHBA'}; nameforthis='TGF-ß signaling';

%Full TGF-ß signaling
%genelist={'TFGB','ACV','NODAL','LEFTY','INH','BMP','SMAD','NOGGIN','CHRD','ALK','CFC','FST'}; nameforthis='Full TGF-ß signaling';


%SHH-signaling (look with normal)
%genelist={'SHH','SMO','PTCH1','PTCH2','GLI1','GLI2','GLI3','GLI4'}; nameforthis='Shh signaling';

%FGF-signaling (look with startswith
%genelist={'FGF','MAPK'}; nameforthis='FGF signaling';

% %basic list from paper
%  genelist={'SOX2','OTX2','CDH1','TFAP2A','GATA3','DLX5','FST','SOX1','NKX1-2','TBXT','MSGN1','TBX6','MESP1','PDGFRA','LEFTY2','LHX1','GATA6','BMP4','HAND1','SNAI1','SNAI2','FOXF1',...
%      'POSTN','ANXA1','CHRD','NOTO','CDH2','FOXA2','SOX17','EOMES','MEF2C','PECAM1','RUNX1','HBE1','GATA1'};
% nameforthis='Endoderm';

%to preserve order:

clear templist; clear genedata; clear genenames;
for ii=1:length(genelist)
    if ii==1
    templist(1,:)=getExpression(ds,genelist{ii},'Precise');
    else
    templist(ii,:)=getExpression(ds,genelist{ii},'Precise');
    end
end
genedata=templist; genenames=genelist;



length(genenames)

%% now graph it with a heatmap in black

close all;

fs=30;
sizeXY= [screensize(3)/2.5, screensize(4)/2];

%nameforthis='Pluripotent markers';
if length(genenames)>20
    fs=16;
    sizeXY=[screensize(3)/5, screensize(4)];
end

fig=figure('Position',[1,1,sizeXY(1),sizeXY(2)],'color','black');
hh=heatmap(round(genedata),'Colormap',hot,'CellLabelColor','none','FontColor','white','FontSize',fs);
hhSt=struct(hh).Heatmap;
hhStGrid=struct(hhSt).Grid;
hhStGrid.ColorData=uint8([100; 100; 100; 256]);
title(nameforthis);

hh.YDisplayLabels=genenames;
hh.XDisplayLabels=conditions;

axs=struct(gca);
cb=axs.Colorbar;

if max(cb.Limits)>200
 cb.Ticks=[0,ceil((max(hh.ColorData(:))-min(hh.ColorData(:)))/2/100)*100, ceil(max(hh.ColorData(:))/100)*100];
 cb.Limits=[0 ceil(max(hh.ColorData(:))/100)*100];
hh.ColorLimits=cb.Limits;
end

if max(cb.Limits)<5
    cb.Ticks=[0, 5, 10];
    cb.Limits=[0 10];
    hh.ColorLimits=cb.Limits;
end

filename=strcat(nameforthis,' heatmap',' in black');
fig.InvertHardcopy = 'off';
% 
  saveas(fig, fullfile(resultsDir, filename), 'png');
  savefig(fig,fullfile(resultsDir, [filename '.fig']),'compact');

 
%% now graph it with a heatmap in white

close all;

fs=30;
sizeXY= [screensize(3)/2.5, screensize(4)/2];

%nameforthis='Pluripotent markers';
if length(genenames)>20
    fs=16;
    sizeXY=[screensize(3)/5, screensize(4)];
end
% 


fig=figure('Position',[1,1,sizeXY(1),sizeXY(2)]);
hh=heatmap(round(genedata),'CellLabelColor','none','FontColor','black','FontSize',fs,'FontName','Arial');
hh.Colormap(1,:)=[1 1 1];
hhSt=struct(hh).Heatmap;
hhStGrid=struct(hhSt).Grid;
hhStGrid.ColorData=uint8([100; 100; 100; 256]);
title(nameforthis);

hh.YDisplayLabels=genenames;
hh.XDisplayLabels=conditions;

axs=struct(gca);
cb=axs.Colorbar;

if max(cb.Limits)>200
 cb.Ticks=[0,ceil((max(hh.ColorData(:))-min(hh.ColorData(:)))/2/100)*100, ceil(max(hh.ColorData(:))/100)*100];
 cb.Limits=[0 ceil(max(hh.ColorData(:))/100)*100];

 
 hh.ColorLimits=cb.Limits;

end

if max(cb.Limits)<5
    cb.Ticks=[0, 5, 10];
    cb.Limits=[0 10];
    hh.ColorLimits=cb.Limits;
end

filename=strcat(nameforthis,' heatmap',' in white');
%fig.InvertHardcopy = 'off';
% 
  saveas(fig, fullfile(resultsDir, filename), 'png');
  savefig(fig,fullfile(resultsDir, [filename '.fig']),'compact');
  saveas(fig, fullfile(resultsDir, filename), 'svg');

 
