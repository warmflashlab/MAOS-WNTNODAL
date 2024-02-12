// Este es para leer el folder principal y de ahí lea subfolders, los subfolders tienen que tener 5 archivos: 4 tiffs y un txt, si no no va a jalar ! !!!!
//va a leer el nombre del subfolder y se lo va poner al archivo final.

//mainDir="/Volumes/MAOS 3/Spinning Disk/210921-E2130/Day1/Appended/";

mainDir="/Volumes/MAOS 2/Spinning Disk/live plate/hyperstacks/";

saveDir="/Volumes/MAOS 2/Spinning Disk/live plate/hyperstacks/";

// how many channels? put in cc
folders=getFileList(mainDir);


print(".start");
 for (i=0; i<folders.length; i++) {
        if (endsWith(folders[i], "/")) {
           subDir=folders[i];
//ahora mete todo aquí, pero primero ve cómo harás el savename

files=getFileList(mainDir+subDir);
filename1=files[1]; filename2=files[2]; filename3=files[3]; filename4=files[4];

filesname=File.getName(files[1]);
//makesure all files have that -2 in their name for partitioning the filename
//saveName="FIJI-hyperstack-"+substring(filesname, 0,lastIndexOf(filesname, "-2"))+".tif"; //this is the old savename
saveName="FIJI-video-"+substring(filesname, 0,lastIndexOf(filesname, "-N"))+"-"+substring(subDir, 0, lastIndexOf(subDir,"/"))+".tif";
print(i," of ", folders.length);

//saveName="FIJI-testMacro-Hyperstack-2113-CHIP-NMP10x-800-c5.tif";
// a partir de aquí empiezan las operaciones
"---Starting new operation"
cc=2; 
"---Opening files 1 and 2"
open(mainDir+subDir+filename1);
open(mainDir+subDir+filename2);
selectWindow(filename1);
getDimensions(width, height, channels, slices, frames);

if (slices>1) {
	newslices=slices;
}

if (frames>1) {
	newslices=frames;
}

//newslices=frames; // check that sometimes slices are frames or viceversa.
"---Transforming to Hyperstacks"
selectWindow(filename1);
run("Stack to Hyperstack...", "order=xyczt(default) channels=cc slices=newslices frames=1 display=Color");
selectWindow(filename2);
run("Stack to Hyperstack...", "order=xyczt(default) channels=cc slices=newslices frames=1 display=Color");

print("---First pair-stitching initiated");
print("---Please wait...");
run("Pairwise stitching", "first_image=[filename1] second_image=[filename2] fusion_method=[Linear Blending] fused_image=[filename1<->filename2] check_peaks=5 ignore compute_overlap x=0.0000 y=0.0000 z=0.0000 registration_channel_image_1=[Average all channels] registration_channel_image_2=[Average all channels]");

selectWindow(filename1);
close();
selectWindow(filename2);
close();

"---Opening files 3 and 4"
open(mainDir+subDir+filename3);
open(mainDir+subDir+filename4);
"---Transforming to Hyperstacks"
selectWindow(filename3);
run("Stack to Hyperstack...", "order=xyczt(default) channels=cc slices=newslices frames=1 display=Color");
selectWindow(filename4);
run("Stack to Hyperstack...", "order=xyczt(default) channels=cc slices=newslices frames=1 display=Color");

print("---Second pair-stitching initiated");
print("---Please wait...");

textp1="first_image=["+filename3+"] ";
textp2="second_image=["+filename4+"] ";
textp3="fusion_method=[Linear Blending] fused_image=[filename3<->filename4] check_peaks=5 ignore compute_overlap x=0.0000 y=0.0000 z=0.0000 registration_channel_image_1=[Average all channels] registration_channel_image_2=[Average all channels]";

run("Pairwise stitching", textp1+textp2+textp3);

selectWindow(filename3);
close();
selectWindow(filename4);
close();

print("---Last pair-stitching initiated");
print("---Please wait...");
run("Pairwise stitching", "first_image=[filename1<->filename2] second_image=[filename3<->filename4] fusion_method=[Linear Blending] fused_image=[allfused] check_peaks=5 ignore compute_overlap x=0.0000 y=0.0000 z=0.0000 registration_channel_image_1=[Average all channels] registration_channel_image_2=[Average all channels]");
print("---Saving...");
saveAs("Tiff",saveDir+saveName);
print("----Hyperstack saved");
run("Close All");
print(".");

        }
 }
 print("---All subfolders processed ");     