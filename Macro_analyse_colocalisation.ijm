//==================================================================
//Script d'analyse de colocalisation semi-automatique
//
//
//==================================================================

//Direction du dossier avec les différents fichier à analyser
getBoolean("Oki, petite intro, faite 3 dossier, un pour deposer les images \n" +
"un autre pour les donnees des facteur et un dernier pour les donnees des graphiques \n" +
"Aussi, tout dependament du nombre d'image, cela peux prendre du temps, disons quelques minutes. " + 
"Aussi, on calcule le threashold en coupant ls 50% les plus faible( à la medianne). Ces pas mal tout!", "Okidoki", "Bon ben laisse faire");

//appelle des directions pour les dossier oû les images, les données des coéficiens ainsi que des données des graphiques 
//serons entreposer

direction_dossier = getDirectory("Mettre le lien du dossier où le fichier lif se trouve");
direction_dossier_puit_image = getDirectory("Mettre le lien du dossier où les images pour l'analyse doit être entreposer");
direction_dossier_puit_data = getDirectory("Mettre le lien du dossier où le fichier des facteurs de colocalisation serons entreposé");
direction_dossier_puit_data_graphique = getDirectory("Mettre le lien du dossier où le fichier des données des graphiques serons entreposé");

//Nom du fichier lif et nom utiliser pour entreposer les données sous un nom = name

nom_du_fichier = getString("Mettre le nom du fichier avec le type", "Ex : Dvirilis_croche_.lif");
name = getString("Mettre le nom de l'espèce + Génotype + sonde", "Ex : Dvirilis_croche_DINE");

//ici, mettre les numéros d'images à analysé
n_image = getString("Numéros d'image à analyser", "ex : 1,6,8 ...");

n_image = split(n_image, ",");

longeur = n_image.length;

//sélection et ouverture de nos images d'intérêt

run("Bio-Formats Importer", "open=[" + direction_dossier + nom_du_fichier + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT open_all_series");

// j = numéro d'image à analyser dans le fichier
// i = longeur de l'array et ca permet de retourver spécifiquement [i] pour avoir j 

message = getBoolean("Es-ce que vous voulez mesurer la colocalisation avec des ROI", "Oui", "Non");

if ( message == true) {

for (i = 0; i < longeur; i++) {

j = n_image[i];

a = IJ.pad(j, 3);

//Sélection de la zone d'intéret pour chaque image

selectWindow(nom_du_fichier + " - Image" + a);

waitForUser("Dessiner la région à analyser pour cette image");

roiManager("add");

run("Crop");

rename(name + "Image_" + j);

run("Split Channels");

ch1 = "C1-" + name + "Image_" + j;

ch2 = "C2-" + name + "Image_" + j;

ch3 = "C3-" + name + "Image_" + j;

close(ch1);

// Renommer et selectioner la fenetre

selectWindow(ch2);

saveAs("Tiff", direction_dossier_puit_image + ch2);
close(ch2);

// Renommer et selectioner la fenetre

selectWindow(ch3);
saveAs("Tiff", direction_dossier_puit_image + ch3);
close(ch3);

roiManager("Delete");

	}
}

else {
	
	setBatchMode(true);
	
	for (i = 0; i < longeur; i++) {

j = n_image[i];

a = IJ.pad(j, 3);

selectWindow(nom_du_fichier + " - Image" + a);

rename(name + "Image_" + j);

run("Split Channels");

ch1 = "C1-" + name + "Image_" + j;

ch2 = "C2-" + name + "Image_" + j;

ch3 = "C3-" + name + "Image_" + j;

close(ch1);

// Renommer et selectioner la fenetre

selectWindow(ch2);

saveAs("Tiff", direction_dossier_puit_image + ch2);
close(ch2);

// Renommer et selectioner la fenetre

selectWindow(ch3);
saveAs("Tiff", direction_dossier_puit_image + ch3);
close(ch3);
	}
}

run("Close All");

message = "Image entreposer, la prochaine étape peux être lente, continuer?";
yesLabel = "C'est partie!";
noLabel = "Nope";
getBoolean(message, yesLabel, noLabel);

//initialisation du calcule de médianne

run("Set Measurements...", "limit redirect=None decimal=3");

for (i = 0; i < longeur; i++) {
	
	j = n_image[i];
	
	//ouverture des 2 images à analyser
	
	imagea = direction_dossier_puit_image + "C2-" + name + "Image_" + j + ".tif";
    imageb = direction_dossier_puit_image + "C3-" + name + "Image_" + j + ".tif";

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
    
	// ouverture du puling et lencement du calcule de pearson overlap mm cytofluo ica et ccf
	
	 run("JACoP ", "imga=[" + "C2-" + name + "Image_" + j + ".tif" + "] imgb=[" + "C3-" + name + "Image_" + j + ".tif" + "] thra=" + thresholda + "  thrb=" + thresholdb + " pearson overlap mm ccf=100 cytofluo ica");

	 run("Clear Results");
	 
	 //Sauvegarde des différent facteur
	 
	 //ICA A
	 
	 selectWindow("ICA A (C2-" + name + "Image_" + j + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_ICA_A" + name + "_Image_" + j + ".csv");
	 
	 close("ICA A (C2-" + name + "Image_" + j + ".tif)");
	 
	 run("Clear Results");
	 
	 //ICA A
	 
	 selectWindow("ICA B (C3-" + name + "Image_" + j + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_ICA_B" + name + "_Image_" + j + ".csv");
	 
	 close("ICA B (C3-" + name + "Image_" + j + ".tif)");
	 
	 run("Clear Results");
	 
	 //Cyrofluogram
	 
	 selectWindow("Cytofluorogram between C2-" + name + "Image_" + j + ".tif and C3-" + name + "Image_" + j + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_Cytofluogram" + name + "_Image_" + j + ".csv");
	 
	 close("Cytofluorogram between C2-" + name + "Image_" + j + " and C3-" + name + "Image_" + j);
	 
	 run("Clear Results");
	 
	 //Van steensels CCF
	 
	 selectWindow("Van Steensel's CCF between C2-" + name + "Image_" + j + ".tif and C3-" + name + "Image_" + j + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_VanSteensel" + name + "_Image_" + j + ".csv");
	 
	 close("Van Steensel's CCF between C2-" + name + "Image_" + j + " and C3-" + name + "Image_" + j);
	 
	 run("Clear Results");
	 
	 //On se débarasse de toute les fenetres
	 
	 close("C2-" + name + "Image_" + j + ".tif");
	 
	 close("C3-" + name + "Image_" + j + ".tif");
	 
	 run("Clear Results");
	 
	 run("Close All");
	 
	 // --- TUEUR DE FENÊTRE JACOP --- Ce petit code à été fait avec Gemini-3
// Ceci utilise du JavaScript pour fermer les fenêtres que ImageJ ne voit pas
eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }");

	 
};

setBatchMode(false);

//Importation des données en fichier lisable par fichier csv

logContent = getInfo("log");

filePath = direction_dossier_puit_data + "/Coef_coloc_"+name+".csv"; 

File.saveString(logContent, filePath);

print("Log content saved to: " + filePath);

print("\\Clear");

getBoolean("Terminer", "super", "super!!");

//A la suite de ca, j'ai fait un script R pour analyser les informations du fichier
//Donc utiliser le, ou du moins, la fonction de lecture du fichier, je vais encore la paufiné






