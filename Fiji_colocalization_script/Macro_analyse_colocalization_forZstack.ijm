//==================================================================
//Semi-automatic colocalization analysis script
//
//With data export of indices and graphics
//
//This script will be adapted for Z-stack, another one was made for serial images
//
//==================================================================

//Directory of the folder with different files to analyze

//Ext.CLIJx_imageJ2RichardsonLucyDeconvolution(Image input, Image kernel_input, Image destination, Number num_iterations);
folder_direction = getDirectory("Enter the link to the folder where the lif file is located");
folder_direction_image_well = getDirectory("Enter the link to the folder where images for analysis should be stored");
folder_direction_data_image_zstack = getDirectory("Enter the link to the folder where colocalization factor files will be stored for Voxel analysis");
folder_direction_data_graph_zstack = getDirectory("Enter the link to the folder where graph data files will be stored for Voxel analysis");

folder_direction_data_image_to_image = getDirectory("Enter the link to the folder where colocalization factor files will be stored for single image");
folder_direction_data_graph_image_to_image = getDirectory("Enter the link to the folder where graph data files will be stored for single image");

filename = getString("Enter the file name with type", "Ex: Dvirilis_hook_.lif");
name = getString("Enter the species + Genotype + probe name", "Ex: Dvirilis_hook_DINE");

// Series numbers to analyze
n_image = getString("Numbers of series to analyze", "ex: 1,2,3,4 ...");

n_image = replace(n_image, " ", "");
n_image = split(n_image, ",");
number_of_series = n_image.length;

//Loading CLIJ2
run("CLIJ2 Macro Extensions", "cl_device="); 

for ( n = 0; n < number_of_series; n++) {
	
	message = getBoolean("Do you want to measure colocalization with ROIs", "Yes", "No");
	
	//Analysis of voxels and volume
	
	j = n_image[n];
	
	series_number = IJ.pad(j, 3);
	
	run("Bio-Formats Importer", "open=[" + folder_direction + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + j);
	
	if ( message == true ) {
	
	selectWindow(filename + " - Series" + series_number);
	
	new_name = "Zstack" +"_Series_" + j + name;

	rename(new_name);

	//Selection of the region of interest for the image

	waitForUser("Draw the region to analyze for this image");

	roiManager("add");

	run("Crop");
	
	roiManager("reset");

	run("Split Channels");

	ch1 = "C1-" + new_name;

	close(ch1);
	
	// Rename and select the window
	
	ch3 = "C3-" + new_name;
	
	selectWindow(ch3);
	
	saveAs("Tiff", folder_direction_image_well + ch3);
	
	rename(ch3); 
	
	//Measure the number of images in the Zstack
	
	size = nSlices; 
	
	//The following loops select each image, rename it and save it for later
	
	run("Stack to Images");
	
	for (i = 1; i <= size; i++) {
		a = IJ.pad(i, 4);
		
	    selectImage("C3-" + new_name + "-"+ a);
	    
	    image_name = "C3_Zstack" + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", folder_direction_image_well + image_name);
	    
	    close("C3-" + new_name + "-"+ a);
	};
	
	// Rename and select the window
	
	ch2 = "C2-" + new_name;
	
	selectWindow(ch2);
	
	saveAs("Tiff", folder_direction_image_well + ch2);
	
	rename(ch2); 

	run("Stack to Images");
	
	
	for (i = 1; i <= size; i++) {
		
		a = IJ.pad(i, 4);
	
	    selectImage("C2-" + new_name + "-" + a);
	    
	    image_name = "C2_Zstack"  + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", folder_direction_image_well + image_name);
	    
	    close("C2-" + new_name + "-"+ a);
	};
	

	
	run("Close All");
	
	};
	
	else {
		

		
	serial_number = IJ.pad(j, 3);
	
	selectWindow(filename + " - Series" + serial_number);
		
	new_name = "Zstack" +"_Series_" + j + "_" + name;
	
	rename(new_name);
	
	run("Split Channels");
	
	ch1 = "C1-" + new_name;
	
	close(ch1);
	
	// Rename and select the window
	
	ch3 = "C3-" + new_name;
	
	selectWindow(ch3);
	
	saveAs("Tiff", folder_direction_image_well);
	
	rename(ch3); 
	
	run("Stack to Images");
	
	for (i = 1; i <= size; i++) {
		a = IJ.pad(i, 4);
		
	    selectImage("C3-" + new_name + "-"+ a);
	    
	    image_name = "C3_Zstack" + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", folder_direction_image_well + image_name);
	    
	    close("C3-" + new_name + "-"+ a);
	};
	
	// Rename and select the window
	
	ch2 = "C2-" + new_name;
	
	selectWindow(ch2);
	
	saveAs("Tiff", folder_direction_image_well);
	
	rename(ch2); 
	
	run("Stack to Images");
	
	for (i = 1; i <= size; i++) {
		
		a = IJ.pad(i, 4);
	
	    selectImage("C2-" + new_name + "-" + a);
	    
	    image_name = "C2_Zstack"  + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", folder_direction_image_well + image_name);
	    
	    close("C2-" + new_name + "-"+ a);
	};
	

	
	run("Close All");
	
	};	
	
	
	//Warning message
	
	message = "Image stored, the next step may be slow, continue?";
	yesLabel = "Let's go!";
	noLabel = "Nope";
	getBoolean(message, yesLabel, noLabel);
	
	//Initialize threshold calculation
	
	run("Set Measurements...", "limit redirect=None decimal=3");
	
	//3D analysis function JACoP and DIANA
	
	zstacka = folder_direction_image_well + "C2-Zstack" + "_Series_" + j + name + ".tif";
	zstackb = folder_direction_image_well + "C3-Zstack" + "_Series_" + j + name + ".tif";
	
	open(zstacka);
	
	selectImage("C2-Zstack" + "_Series_" + j + name + ".tif");
	    
	    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
	    
	    run("Measure");
	    
	    zthresholda = getResult("MinThr", 0);
	 
	    open(zstackb);
	    
	    selectImage("C3-Zstack" + "_Series_" + j + name + ".tif");
	    
	    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
	    
	    run("Measure");
	    
	    zthresholdb = getResult("MinThr", 1);
	    
		// Opening the plugin and launching the calculation of Pearson overlap MM cytofluo ICA and CCF
		
		 run("JACoP ", "imga=[" + "C2-Zstack" + "_Series_" + j + name + ".tif" + "] imgb=[" + "C3-Zstack" + "_Series_" + j + name + ".tif" + "] thra=" + zthresholda + "  thrb=" + zthresholdb + " pearson overlap mm ccf=100 cytofluo ica");
		 //Closing unnecessary windows
		 
		 selectWindow("ICA A (C2-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_zstack + "/Result_ICA_A" + "_Series_" + j + name + ".csv");
		 
		 close("ICA A (C2-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("ICA B (C3-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_zstack + "/Result_ICA_B" + "_Series_" + j + name + ".csv");
		 
		 close("ICA B (C3-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("Cytofluorogram between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_zstack + "/Result_Cytofluogram" + "_Series_" + j + name + ".csv");
		 
		 close("Cytofluorogram between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 run("Clear Results");
		 
		 selectWindow("Van Steensel's CCF between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_zstack + "/Result_VanSteensel" + "_Series_" + j + name + ".csv");
		 
		 close("Van Steensel's CCF between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 run("Clear Results");
		 
		 
		 logContent = getInfo("log");
	
		filePath = folder_direction_data_image_zstack + "/Coef_coloc_zstack"  + "_Series_" + j + name + ".csv"; 
		
		File.saveString(logContent, filePath);
		
		print("Log content saved to: " + filePath);
			 
		 
		 run("Close All");
		 
		 print("\\Clear");
		 
		 // --- JACoP WINDOW KILLER ---
	// This uses JavaScript to close windows that ImageJ doesn't see
	eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }"); 


	for (i = 1; i <=size; i++) {
		
		//Opening the 2 images to analyze
		
		imagea = folder_direction_image_well + "C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif";
	    imageb = folder_direction_image_well + "C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif";
	
		open(imagea);
	    
	    selectImage("C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif");
	    
	    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
	    
	    run("Measure");
	    
	    thresholda = getResult("MinThr", 0);
	 
	    open(imageb);
	    
	    selectImage("C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif");
	    
	    setAutoThreshold("RenyiEntropy dark 16-bit no-reset");
	    
	    run("Measure");
	    
	    thresholdb = getResult("MinThr", 1);
	    
		// Opening the plugin and launching the calculation of Pearson overlap MM cytofluo ICA and CCF
		
		 run("JACoP ", "imga=[" + "C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif" + "] imgb=[" + "C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif" + "] thra=" + thresholda + "  thrb=" + thresholdb + " pearson overlap mm ccf=100 cytofluo ica");
		 //Closing unnecessary windows
		 
		 selectWindow("ICA A (C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_image_to_image + "/Result_ICA_A" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("ICA A (C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("ICA B (C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_image_to_image + "/Result_ICA_B" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("ICA B (C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("Cytofluorogram between C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif and C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_image_to_image + "/Result_Cytofluogram" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("Cytofluorogram between C2_Zstack" + "_Series_" + j + name + "Image_" + i + " and C3_Zstack" + "_Series_" + j + name + "Image_" + i);
		 
		 run("Clear Results");
		 
		 selectWindow("Van Steensel's CCF between C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif and C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", folder_direction_data_graph_image_to_image + "/Result_VanSteensel" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("Van Steensel's CCF between C2_Zstack" + "_Series_" + j + name + "Image_" + i + " and C3_Zstack" + name + "Image_" + i);
		 
		 run("Clear Results");
		 
		 run("Close All");
		 
	// --- JACoP WINDOW KILLER ---
	// This uses JavaScript to close windows that ImageJ doesn't see
	eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }"); 
 
		 //Selection and storage of graph data for R
	};
	
	
	logContent = getInfo("log");
	
	filePath = folder_direction_data_image_to_image + "/Coef_coloc_"  + "_Series_" + j + name + ".csv"; 
	
	File.saveString(logContent, filePath);
	
	print("Log content saved to: " + filePath);
	
	
	print("\\Clear");

};

getBoolean("Done", "super", "super!!");

//Following this, I created an R script to analyze the information from the file
//So use it, or at least the file reading function, I will refine it further
