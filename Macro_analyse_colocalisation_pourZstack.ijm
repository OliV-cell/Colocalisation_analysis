//==================================================================
//Script d'analyse de colocalisation semi-automatique
//
//Avec exportation des données des indices et des graphiques 
//
//Se script vas être adapté au Z-stack, un autres à été fait pour les images à la chaine
//
//==================================================================

//Direction du dossier avec les différents fichier à analyser

//Ext.CLIJx_imageJ2RichardsonLucyDeconvolution(Image input, Image kernel_input, Image destination, Number num_iterations);
direction_dossier = getDirectory("Mettre le lien du dossier où le fichier lif se trouve");
direction_dossier_puit_image = getDirectory("Mettre le lien du dossier où les images pour l'analyse doit être entreposer");
direction_dossier_puit_data_image_zstack = getDirectory("Mettre le lien du dossier où le fichier des facteurs de colocalisation serons entreposé pour l'anaylse de Voxel");
direction_dossier_puit_data_graphique_zstack = getDirectory("Mettre le lien du dossier où le fichier des données des graphiques serons entreposé pour l'anaylse de Voxel");

direction_dossier_puit_data_image_a_image = getDirectory("Mettre le lien du dossier où le fichier des facteurs de colocalisation serons entreposé pour image seul");
direction_dossier_puit_data_graphique_image_a_image = getDirectory("Mettre le lien du dossier où le fichier des données des graphiques serons entreposé pour image seul");

nom_du_fichier = getString("Mettre le nom du fichier avec le type", "Ex : Dvirilis_croche_.lif");
name = getString("Mettre le nom de l'espèce + Génotype + sonde", "Ex : Dvirilis_croche_DINE");

// Nombre de série à analyser
n_image = getString("Numéros des séries à analyser", "ex : 1,2,3,4 ...");

n_image = replace(n_image, " ", "");
n_image = split(n_image, ",");
nombre_de_serie = n_image.length;

//chargement de CLIJ2
run("CLIJ2 Macro Extensions", "cl_device="); 

for ( n = 0; n < nombre_de_serie; n++) {
	
	message = getBoolean("Es-ce que vous voulez mesurer la colocalisation avec des ROI", "Oui", "Non");
	
	//Analyse des voxels et volume
	
	j = n_image[n];
	
	numero_liste = IJ.pad(j, 3);
	
	run("Bio-Formats Importer", "open=[" + direction_dossier + nom_du_fichier + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + j);
	
	if ( message == true ) {
	
	selectWindow(nom_du_fichier + " - Series" + numero_liste);
	
	nouveau_nom = "Zstack" +"_Series_" + j + name;

	rename(nouveau_nom);

	//Sélection de la zone d'intéret pour l'image

	waitForUser("Dessiner la région à analyser pour cette image");

	roiManager("add");

	run("Crop");
	
	roiManager("reset");

	run("Split Channels");

	ch1 = "C1-" + nouveau_nom;

	close(ch1);
	
	// Renommer et selectioner la fenetre
	
	ch3 = "C3-" + nouveau_nom;
	
	selectWindow(ch3);
	
	saveAs("Tiff", direction_dossier_puit_image + ch3);
	
	rename(ch3); 
	
	//on mesure le nombre d'image dans le Zstack
	
	taille = nSlices; 
	
	//Les boucle qui suit, sélection i image et le renome et le sauvegarde pour tantot
	
	run("Stack to Images");
	
	for (i = 1; i <= taille; i++) {
		a = IJ.pad(i, 4);
		
	    selectImage("C3-" + nouveau_nom + "-"+ a);
	    
	    nom = "C3_Zstack" + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", direction_dossier_puit_image + nom);
	    
	    close("C3-" + nouveau_nom + "-"+ a);
	};
	
	// Renommer et selectioner la fenetre
	
	ch2 = "C2-" + nouveau_nom;
	
	selectWindow(ch2);
	
	saveAs("Tiff", direction_dossier_puit_image + ch2);
	
	rename(ch2); 

	run("Stack to Images");
	
	
	for (i = 1; i <= taille; i++) {
		
		a = IJ.pad(i, 4);
	
	    selectImage("C2-" + nouveau_nom + "-" + a);
	    
	    nom = "C2_Zstack"  + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", direction_dossier_puit_image + nom);
	    
	    close("C2-" + nouveau_nom + "-"+ a);
	};
	

	
	run("Close All");
	
	};
	
	else {
		

		
	numero = IJ.pad(j, 3);
	
	selectWindow(nom_du_fichier + " - Series" + numero);
		
	nouveau_nom = "Zstack" +"_Series_" + j + "_" + name;
	
	rename(nouveau_nom);
	
	run("Split Channels");
	
	ch1 = "C1-" + nouveau_nom;
	
	close(ch1);
	
	// Renommer et selectioner la fenetre
	
	ch3 = "C3-" + nouveau_nom;
	
	selectWindow(ch3);
	
	saveAs("Tiff", direction_dossier_puit_image);
	
	rename(ch3); 
	
	run("Stack to Images");
	
	for (i = 1; i <= taille; i++) {
		a = IJ.pad(i, 4);
		
	    selectImage("C3-" + nouveau_nom + "-"+ a);
	    
	    nom = "C3_Zstack" + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", direction_dossier_puit_image + nom);
	    
	    close("C3-" + nouveau_nom + "-"+ a);
	};
	
	// Renommer et selectioner la fenetre
	
	ch2 = "C2-" + nouveau_nom;
	
	selectWindow(ch2);
	
	saveAs("Tiff", direction_dossier_puit_image);
	
	rename(ch2); 
	
	run("Stack to Images");
	
	for (i = 1; i <= taille; i++) {
		
		a = IJ.pad(i, 4);
	
	    selectImage("C2-" + nouveau_nom + "-" + a);
	    
	    nom = "C2_Zstack"  + "_Series_" + j + name + "Image_" + i;
	    
	    saveAs("Tiff", direction_dossier_puit_image + nom);
	    
	    close("C2-" + nouveau_nom + "-"+ a);
	};
	

	
	run("Close All");
	
	};	
	
	
	//message d'avertissement
	
	message = "Image entreposer, la prochaine étape peux être lente, continuer?";
	yesLabel = "C'est partie!";
	noLabel = "Nope";
	getBoolean(message, yesLabel, noLabel);
	
	//initialisation du calcule du treshold
	
	run("Set Measurements...", "limit redirect=None decimal=3");
	
	//Fonction analyse 3D Jacop et DIaNA
	
	zstacka = direction_dossier_puit_image + "C2-Zstack" + "_Series_" + j + name + ".tif";
	zstackb = direction_dossier_puit_image + "C3-Zstack" + "_Series_" + j + name + ".tif";
	
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
	    
		// ouverture du puling et lencement du calcule de pearson overlap mm cytofluo ica et ccf
		
		 run("JACoP ", "imga=[" + "C2-Zstack" + "_Series_" + j + name + ".tif" + "] imgb=[" + "C3-Zstack" + "_Series_" + j + name + ".tif" + "] thra=" + zthresholda + "  thrb=" + zthresholdb + " pearson overlap mm ccf=100 cytofluo ica");
		 //Fermeture des fenêtres non utile
		 
		 selectWindow("ICA A (C2-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_zstack + "/Resultat_ICA_A" + "_Series_" + j + name + ".csv");
		 
		 close("ICA A (C2-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("ICA B (C3-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_zstack + "/Resultat_ICA_B" + "_Series_" + j + name + ".csv");
		 
		 close("ICA B (C3-Zstack" + "_Series_" + j + name + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("Cytofluorogram between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_zstack + "/Resultat_Cytofluogram" + "_Series_" + j + name + ".csv");
		 
		 close("Cytofluorogram between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 run("Clear Results");
		 
		 selectWindow("Van Steensel's CCF between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_zstack + "/Resultat_VanSteensel" + "_Series_" + j + name + ".csv");
		 
		 close("Van Steensel's CCF between C2-Zstack" + "_Series_" + j + name + ".tif and C3-Zstack" + "_Series_" + j + name + ".tif");
		 
		 run("Clear Results");
		 
		 
		 logContent = getInfo("log");
	
		filePath = direction_dossier_puit_data_image_zstack + "/Coef_coloc_zstack"  + "_Series_" + j + name + ".csv"; 
		
		File.saveString(logContent, filePath);
		
		print("Log content saved to: " + filePath);
			 
		 
		 run("Close All");
		 
		 print("\\Clear");
		 
		 // --- TUEUR DE FENÊTRE JACOP ---
	// Ceci utilise du JavaScript pour fermer les fenêtres que ImageJ ne voit pas
	eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }");
		 
	// Analyse des volume avec DIANa
	
		open(zstacka);
	
		selectImage("C2-Zstack" + "_Series_" + j + name + ".tif");
		
Ext.CLIJ2_clear();
		
		Ext.CLIJ2_push("C2-Zstack" + "_Series_" + j + name + ".tif");

		Ext.CLIJ2_voronoiOtsuLabeling("C2-Zstack" + "_Series_" + j + name + ".tif", label_image_C2, 2, 2);
		
		Ext.CLIJ2_pull(label_image_C2);
		
		Ext.CLIJ2_push(label_image_C2);
		
		Ext.CLIJ2_mergeTouchingLabels(label_image_C2, label_merged_C2);
		
		Ext.CLIJ2_pull(label_merged_C2);
		
		selectWindow(label_merged_C2);
		
		rename("Mask_Objects_C2-Zstack" + "_Series_" + j + name + ".tif");
		
		open(zstackb);
	
		selectImage("C3-Zstack" + "_Series_" + j + name + ".tif");
		
		Ext.CLIJ2_push("C3-Zstack" + "_Series_" + j + name + ".tif");

		Ext.CLIJ2_voronoiOtsuLabeling("C3-Zstack" + "_Series_" + j + name + ".tif", label_image_C3, 2, 2);
		
		Ext.CLIJ2_pull(label_image_C3);
		
		Ext.CLIJ2_push(label_image_C3);
		
		Ext.CLIJ2_mergeTouchingLabels(label_image_C3, label_merged_C3);

		Ext.CLIJ2_pull(label_merged_C3);
		
		selectWindow(label_merged_C3);
		
		rename("Mask_Objects_C3-Zstack" + "_Series_" + j + name + ".tif");
	
		selectImage("C2-Zstack" + "_Series_" + j + name + ".tif");
		
		getVoxelSize(width, height, depth, unit);
		
		selectImage("Mask_Objects_C2-Zstack" + "_Series_" + j + name + ".tif");
		
		run("16-bit");
		
		setVoxelSize(width, height, depth, unit);
	
		selectImage("C3-Zstack" + "_Series_" + j + name + ".tif");
		
		getVoxelSize(width, height, depth, unit);
		
		selectImage("Mask_Objects_C3-Zstack" + "_Series_" + j + name + ".tif");
		
		run("16-bit");
	
		setVoxelSize(width, height, depth, unit);
		
		//run("DiAna_Analyse", "img1=C1-" + i +".tif img2=C2-" + i +".tif lab1=labelled-A lab2=labelled-B coloc");
		
		selectImage("C3-Zstack" + "_Series_" + j + name + ".tif");
		selectImage("C2-Zstack" + "_Series_" + j + name + ".tif");
	
		run("DiAna_Analyse", " img1=[" + "C2-Zstack" + "_Series_" + j + name + ".tif" + "] img2=[" + "C3-Zstack" + "_Series_" + j + name + ".tif" + "] lab1=[" + "Mask_Objects_C2-Zstack" + "_Series_" + j + name + ".tif" + "] lab2=[" + "Mask_Objects_C3-Zstack" + "_Series_" + j + name + ".tif" + "] coloc");
	
		selectWindow("ColocResults");
	
		saveAs("Results", direction_dossier_puit_data_graphique_zstack + "/Resultat_DIANA" + "_Series_" + j + name + ".csv");
	
		selectImage("coloc");
		
		rename("coloc" + "_Series_" + j);
		
		saveAs("Tiff", direction_dossier_puit_data_graphique_zstack);
		
		run("Close All");
		
		run("Collect Garbage");
		
		selectWindow("Resultat_DIANA" + "_Series_" + j + name + ".csv"); run("Close");
		
		
Ext.CLIJ2_clear();
		
		// --- TUEUR DE FENÊTRE JACOP ---
	// Ceci utilise du JavaScript pour fermer les fenêtres que ImageJ ne voit pas
	eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('DiAna (Distance Analysis') != -1) { frames[i].dispose(); } }");                         
	//boulcle de mesure pour image par image
	
	for (i = 1; i <=taille; i++) {
		
		//ouverture des 2 images à analyser
		
		imagea = direction_dossier_puit_image + "C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif";
	    imageb = direction_dossier_puit_image + "C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif";
	
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
	    
		// ouverture du puling et lencement du calcule de pearson overlap mm cytofluo ica et ccf
		
		
		 run("JACoP ", "imga=[" + "C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif" + "] imgb=[" + "C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif" + "] thra=" + thresholda + "  thrb=" + thresholdb + " pearson overlap mm ccf=100 cytofluo ica");
		 //Fermeture des fenêtres non utile
		 
		 selectWindow("ICA A (C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_image_a_image + "/Resultat_ICA_A" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("ICA A (C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("ICA B (C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_image_a_image + "/Resultat_ICA_B" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("ICA B (C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif)");
		 
		 run("Clear Results");
		 
		 selectWindow("Cytofluorogram between C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif and C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_image_a_image + "/Resultat_Cytofluogram" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("Cytofluorogram between C2_Zstack" + "_Series_" + j + name + "Image_" + i + " and C3_Zstack" + "_Series_" + j + name + "Image_" + i);
		 
		 run("Clear Results");
		 
		 selectWindow("Van Steensel's CCF between C2_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif and C3_Zstack" + "_Series_" + j + name + "Image_" + i + ".tif");
		 
		 Plot.showValues();
		 
		 selectWindow("Results");
		 
		 saveAs("Results", direction_dossier_puit_data_graphique_image_a_image + "/Resultat_VanSteensel" + "_Series_" + j + name + "_Image_" + i + ".csv");
		 
		 close("Van Steensel's CCF between C2_Zstack" + "_Series_" + j + name + "Image_" + i + " and C3_Zstack" + name + "Image_" + i);
		 
		 run("Clear Results");
		 
		 run("Close All");
		 
	// --- TUEUR DE FENÊTRE JACOP ---
	// Ceci utilise du JavaScript pour fermer les fenêtres que ImageJ ne voit pas
	eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }");
		 
		 
		 //Sélection et entreposage des données des graphs pour R
	};
	
	
	logContent = getInfo("log");
	
	filePath = direction_dossier_puit_data_image_a_image + "/Coef_coloc_"  + "_Series_" + j + name + ".csv"; 
	
	File.saveString(logContent, filePath);
	
	print("Log content saved to: " + filePath);
	
	
	print("\\Clear");

};

getBoolean("Terminer", "super", "super!!");

//A la suite de ca, j'ai fait un script R pour analyser les informations du fichier
//Donc utiliser le, ou du moins, la fonction de lecture du fichier, je vais encore la paufiné









