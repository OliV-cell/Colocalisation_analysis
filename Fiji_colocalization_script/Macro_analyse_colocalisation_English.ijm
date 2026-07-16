//==================================================================
//Semi-automatic colocalization analysis script
//
//
//==================================================================

//Directory of the folder with different files to analyze
getBoolean("Alright, brief introduction, create 3 folders, one to store the images \n" +
"another for the data of coefficients and a last one for graph data \n"
"Also, depending on the number of images, this can take time, let's say a few minutes. " + 
"Also, we calculate the threshold by cutting off the weakest 50% (at the median). It's pretty great!", "Alright", "Never mind");

//Call of directories for the folders where images, coefficient data and graph data
//will be stored

folder_direction = getDirectory("Enter the link to the folder where the lif file is located");
folder_direction_image_well = getDirectory("Enter the link to the folder where images for analysis should be stored");
folder_direction_data = getDirectory("Enter the link to the folder where colocalization factor files will be stored");
folder_direction_data_graph = getDirectory("Enter the link to the folder where graph data files will be stored");

//Name of the lif file and name to use to store data under a name = name

filename = getString("Enter the file name with type", "Ex: Dvirilis_hook_.lif");
name = getString("Enter the species + Genotype + probe name", "Ex: Dvirilis_hook_DINE");

//Here, enter the image numbers to analyze
n_image = getString("Image numbers to analyze", "ex: 1,6,8 ...");

n_image = split(n_image, ",");

length = n_image.length;

//Selection and opening of our images of interest

run("Bio-Formats Importer", "open=[" + folder_direction + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT open_all_series");

// j = image number to analyze in the file
// i = length of the array and it allows us to find specifically [i] to get j 

message = getBoolean("Do you want to measure colocalization with ROIs", "Yes", "No");

if ( message == true) {

for (i = 0; i < length; i++) {

j = n_image[i];

a = IJ.pad(j, 3);

//Selection of the region of interest for each image

selectWindow(filename + " - Image" + a);

waitForUser("Draw the region to analyze for this image");

roiManager("add");

run("Crop");

rename(name + "Image_" + j);

run("Split Channels");

ch1 = "C1-" + name + "Image_" + j;

ch2 = "C2-" + name + "Image_" + j;

ch3 = "C3-" + name + "Image_" + j;

close(ch1);

// Rename and select the window

selectWindow(ch2);

saveAs("Tiff", folder_direction_image_well + ch2);
close(ch2);

// Rename and select the window

selectWindow(ch3);
saveAs("Tiff", folder_direction_image_well + ch3);
close(ch3);

roiManager("Delete");

	}
}

else {
	
	setBatchMode(true);
	
	for (i = 0; i < length; i++) {

j = n_image[i];

a = IJ.pad(j, 3);

selectWindow(filename + " - Image" + a);

rename(name + "Image_" + j);

run("Split Channels");

ch1 = "C1-" + name + "Image_" + j;

ch2 = "C2-" + name + "Image_" + j;

ch3 = "C3-" + name + "Image_" + j;

close(ch1);

// Rename and select the window

selectWindow(ch2);

saveAs("Tiff", folder_direction_image_well + ch2);
close(ch2);

// Rename and select the window

selectWindow(ch3);
saveAs("Tiff", folder_direction_image_well + ch3);
close(ch3);
	}
}

run("Close All");

message = "Image stored, the next step may be slow, continue?";
yesLabel = "Let's go!";
noLabel = "Nope";
getBoolean(message, yesLabel, noLabel);

//Initialize median calculation

run("Set Measurements...", "limit redirect=None decimal=3");

for (i = 0; i < length; i++) {
	
	j = n_image[i];
	
	//Opening the 2 images to analyze
	
	imagea = folder_direction_image_well + "C2-" + name + "Image_" + j + ".tif";
    imageb = folder_direction_image_well + "C3-" + name + "Image_" + j + ".tif";

    open(imagea);
    
    selectImage("C2-" + name + "Image_" + j + ".tif");
    
    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
    
    run("Measure");
    
    thresholda = getResult("MinThr", 0);
 
    open(imageb);
    
    selectImage("C3-" + name + "Image_" + j + ".tif");
    
    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
    
    run("Measure");
    
    thresholdb = getResult("MinThr", 1);
    
	// Opening the plugin and launching the calculation of Pearson overlap MM cytofluo ICA and CCF
	
	 run("JACoP ", "imga=[" + "C2-" + name + "Image_" + j + ".tif" + "] imgb=[" + "C3-" + name + "Image_" + j + ".tif" + "] thra=" + thresholda + "  thrb=" + thresholdb + " pearson overlap mm ccf=10[...]

	 run("Clear Results");
	 
	 //Saving of different factors
	 
	 //ICA A
	 
	 selectWindow("ICA A (C2-" + name + "Image_" + j + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_ICA_A" + name + "_Image_" + j + ".csv");
	 
	 close("ICA A (C2-" + name + "Image_" + j + ".tif)");
	 
	 run("Clear Results");
	 
	 //ICA B
	 
	 selectWindow("ICA B (C3-" + name + "Image_" + j + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_ICA_B" + name + "_Image_" + j + ".csv");
	 
	 close("ICA B (C3-" + name + "Image_" + j + ".tif)");
	 
	 run("Clear Results");
	 
	 //Cytofluorogram
	 
	 selectWindow("Cytofluorogram between C2-" + name + "Image_" + j + ".tif and C3-" + name + "Image_" + j + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_Cytofluogram" + name + "_Image_" + j + ".csv");
	 
	 close("Cytofluorogram between C2-" + name + "Image_" + j + " and C3-" + name + "Image_" + j);
	 
	 run("Clear Results");
	 
	 //Van Steensel's CCF
	 
	 selectWindow("Van Steensel's CCF between C2-" + name + "Image_" + j + ".tif and C3-" + name + "Image_" + j + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", folder_direction_data_graph + "/Result_VanSteensel" + name + "_Image_" + j + ".csv");
	 
	 close("Van Steensel's CCF between C2-" + name + "Image_" + j + " and C3-" + name + "Image_" + j);
	 
	 run("Clear Results");
	 
	 //Close all windows
	 
	 close("C2-" + name + "Image_" + j + ".tif");
	 
	 close("C3-" + name + "Image_" + j + ".tif");
	 
	 run("Clear Results");
	 
	 run("Close All");
	 
	 // --- JACoP WINDOW KILLER --- This small code was made with Gemini-3
// This uses JavaScript to close windows that ImageJ doesn't see
eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/0[...]

	 
};

setBatchMode(false);

//Import data into CSV readable file

logContent = getInfo("log");

filePath = folder_direction_data + "/Coef_coloc_"+name+".csv"; 

File.saveString(logContent, filePath);

print("Log content saved to: " + filePath);

print("\\Clear");

getBoolean("Done", "super", "super!!");

//Following this, I created an R script to analyze the information from the file
//So use it, or at least the file reading function, I will refine it further
