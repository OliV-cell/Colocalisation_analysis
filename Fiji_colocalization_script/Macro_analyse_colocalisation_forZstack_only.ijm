//==================================================================
//Semi-automatic colocalization analysis script
//
//With data export of indices and graphics
//
//This script will be adapted for Z-stack, another one was made for serial images
//
//==================================================================

getBoolean("Alright, brief introduction, create 3 folders, one to store the images \n" +
"another for coefficient data and a last one for graph data \n" +
"Also, depending on the number of images, this can take time, let's say a few minutes \n " + 
"Also, we calculate the threshold with Renyi entropy \n " + 
"One more thing, currently I'm still working on a version where Z-stacks are in the same lif file so separate them first! \n" + 
"It's pretty great!", "Alright", "Never mind");

//Directory of the folder with different files to analyze

folder_direction = getDirectory("Enter the link to the folder where the lif file is located");
folder_direction_image_well = getDirectory("Enter the link to the folder where images for analysis should be stored");
folder_direction_data = getDirectory("Enter the link to the folder where colocalization factor files will be stored");
folder_direction_data_graph = getDirectory("Enter the link to the folder where graph data files will be stored");

filename = getString("Enter the file name with type", "Ex: Dvirilis_hook_.lif");
name = getString("Enter the species + Genotype + probe name", "Ex: Dvirilis_hook_DINE");


run("Bio-Formats Importer", "open=[" + folder_direction + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT open_all_series");

message = getBoolean("Do you want to measure colocalization with ROIs", "Yes", "No");

setBatchMode(true);

if ( message == true ) {
	
selectWindow(filename);
	
new_name = "Zstack" + name;

rename(new_name);

//Selection of the region of interest for the image

waitForUser("Draw the region to analyze for this image");

roiManager("add");

run("Crop");

run("Split Channels");

ch1 = "C1-" + new_name;

close(ch1);

// Rename and select the window

ch3 = "C3-" + new_name;

selectWindow(ch3);

//Measure the number of images in the Zstack

size = nSlices; 

//The following loops select each image, rename it and save it for later

run("Stack to Images");

for (i = 1; i <= size; i++) {
	a = IJ.pad(i, 4);
	
    selectImage("C3-" + new_name + "-"+ a);
    
    image_name = "C3_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", folder_direction_image_well + image_name);
    
    close("C3-" + new_name + "-"+ a);
};

// Rename and select the window

ch2 = "C2-" + new_name;

selectWindow(ch2);

run("Stack to Images");

for (i = 1; i <= size; i++) {
	
	a = IJ.pad(i, 4);

    selectImage("C2-" + new_name + "-" + a);
    
    image_name = "C2_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", folder_direction_image_well + image_name);
    
    close("C2-" + new_name + "-"+ a);
};

run("Close All");

};

else {
	
selectWindow(filename);
	
new_name = "Zstack" + name;

rename(new_name);

run("Split Channels");

ch1 = "C1-" + new_name;

close(ch1);

// Rename and select the window

ch3 = "C3-" + new_name;

selectWindow(ch3);

size = nSlices; 

run("Stack to Images");

//The following loops select each image, rename it and save it for later

for (i = 1; i <= size; i++) {
	a = IJ.pad(i, 4);
	
    selectImage("C3-" + new_name + "-"+ a);
    
    image_name = "C3_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", folder_direction_image_well + image_name);
    
    close("C3-" + new_name + "-"+ a);
};

// Rename and select the window

ch2 = "C2-" + new_name;

selectWindow(ch2);

run("Stack to Images");

for (i = 1; i <= size; i++) {
	
	a = IJ.pad(i, 4);

    selectImage("C2-" + new_name + "-" + a);
    
    image_name = "C2_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", folder_direction_image_well + image_name);
    
    close("C2-" + new_name + "-"+ a);
};

run("Close All");

};


//Warning message

setBatchMode(false);

message = "Image stored, the next step may be slow, continue?";
yesLabel = "Let's go!";
noLabel = "Nope";
getBoolean(message, yesLabel, noLabel);

//Initialize median calculation

run("Set Measurements...", "limit redirect=None decimal=3");



for (i = 1; i <=size; i++) {
	
	//Opening the 2 images to analyze
	
	imagea = folder_direction_image_well + "C2_Zstack" + name + "Image_" + i + ".tif";
    imageb = folder_direction_image_well + "C3_Zstack" + name + "Image_" + i + ".tif";

	open(imagea);
    
    selectImage("C2_Zstack" + name + "Image_" + i + ".tif");
    
    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
    
    run("Measure");
    
    thresholda = getResult("MinThr", 0);
 
    open(imageb);
    
    selectImage("C3_Zstack" + name + "Image_" + i + ".tif");
    
    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
    
    run("Measure");
    
    thresholdb = getResult("MinThr", 1);
    
	// Opening the plugin and launching the calculation of Pearson overlap MM cytofluo ICA and CCF
	
	
	 run("JACoP ", "imga=[" + "C2_Zstack" + name + "Image_" + i + ".tif" + "] imgb=[" + "C3_Zstack" + name + "Image_" + i + ".tif" + "] thra=" + thresholda + "  thrb=" + thresholdb + " pearson overlap mm ccf=100 cytofluo ica");
	 //Closing unnecessary windows
	 
	 selectWindow("ICA A (C2_Zstack" + name + "Image_" + i + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_ICA_A" + name + "_Image_" + i + ".csv");
	 
	 close("ICA A (C2_Zstack" + name + "Image_" + i + ".tif)");
	 
	 run("Clear Results");
	 
	 selectWindow("ICA B (C3_Zstack" + name + "Image_" + i + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_ICA_B" + name + "_Image_" + i + ".csv");
	 
	 close("ICA B (C3_Zstack" + name + "Image_" + i + ".tif)");
	 
	 run("Clear Results");
	 
	 selectWindow("Cytofluorogram between C2_Zstack" + name + "Image_" + i + ".tif and C3_Zstack" +  name + "Image_" + i + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_Cytofluogram" + name + "_Image_" + i + ".csv");
	 
	 close("Cytofluorogram between C2_Zstack" + name + "Image_" + i + " and C3_Zstack" + name + "Image_" + i);
	 
	 run("Clear Results");
	 
	 selectWindow("Van Steensel's CCF between C2_Zstack" + name + "Image_" + i + ".tif and C3_Zstack" + name + "Image_" + i + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_VanSteensel" + name + "_Image_" + i + ".csv");
	 
	 close("Van Steensel's CCF between C2_Zstack" + name + "Image_" + i + " and C3_Zstack" + name + "Image_" + i);
	 
	 run("Clear Results");
	 
	 run("Close All");
	 
// --- JACoP WINDOW KILLER ---
// This uses JavaScript to close windows that ImageJ doesn't see
eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }");	 
	 
	 //Selection and storage of graph data for R
};


setBatchMode(false);


logContent = getInfo("log");

filePath = folder_direction_data + "/Coef_coloc_"  + name + ".csv"; 

File.saveString(logContent, filePath);

print("Log content saved to: " + filePath);

getBoolean("Finish", "Ok", "Oki!!");

//Following this, I created an R script to analyze the information from the file
//So use it, or at least the file reading function, I will refine it further
