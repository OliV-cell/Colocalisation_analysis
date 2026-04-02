//==================================================================
//Script d'analyse de colocalisation semi-automatique
//
//Avec exportation des données des indices et des graphiques 
//
//Se script vas être adapté au Z-stack, un autres à été fait pour les images à la chaine
//
//==================================================================

getBoolean("Oki, petite intro, faite 3 dossier, un pour déposer les images \n" +
"un autre pour les données des facteur et un dernier pour les données des graphiques \n" +
"Aussi, tout dépendament du nombre d'image, cela peux prendre du temps, disons quelques minutes \n " + 
"Aussi, on calcule le threashold en coupant ls 50% les plus faible( à la médianne). \n " + 
"Encore une dernière chose, live je travaille encore sur une version ou les Zstacks se trouve dans le même fichier lif donc séparer les avant! \n" + 
"Ces pas mal tout!", "Okidoki", "Bon ben laisse faire");

//Direction du dossier avec les différents fichier à analyser

direction_dossier = getDirectory("Mettre le lien du dossier où le fichier lif se trouve");
direction_dossier_puit_image = getDirectory("Mettre le lien du dossier où les images pour l'analyse doit être entreposer");
direction_dossier_puit_data = getDirectory("Mettre le lien du dossier où le fichier des facteurs de colocalisation serons entreposé");
direction_dossier_puit_data_graphique = getDirectory("Mettre le lien du dossier où le fichier des données des graphiques serons entreposé");

nom_du_fichier = getString("Mettre le nom du fichier avec le type", "Ex : Dvirilis_croche_.lif");
name = getString("Mettre le nom de l'espèce + Génotype + sonde", "Ex : Dvirilis_croche_DINE");


run("Bio-Formats Importer", "open=[" + direction_dossier + nom_du_fichier + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT open_all_series");

message = getBoolean("Es-ce que vous voulez mesurer la colocalisation avec des ROI", "Oui", "Non");

setBatchMode(true);

if ( message == true ) {
	
selectWindow(nom_du_fichier);
	
nouveau_nom = "Zstack" + name;

rename(nouveau_nom);

//Sélection de la zone d'intéret pour l'image

waitForUser("Dessiner la région à analyser pour cette image");

roiManager("add");

run("Crop");

run("Split Channels");

ch1 = "C1-" + nouveau_nom;

close(ch1);

// Renommer et selectioner la fenetre

ch3 = "C3-" + nouveau_nom;

selectWindow(ch3);

//on mesure le nombre d'image dans le Zstack

taille = nSlices; 

//Les boucle qui suit, sélection i image et le renome et le sauvegarde pour tantot

run("Stack to Images");

for (i = 1; i <= taille; i++) {
	a = IJ.pad(i, 4);
	
    selectImage("C3-" + nouveau_nom + "-"+ a);
    
    nom = "C3_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", direction_dossier_puit_image + nom);
    
    close("C3-" + nouveau_nom + "-"+ a);
};

// Renommer et selectioner la fenetre

ch2 = "C2-" + nouveau_nom;

selectWindow(ch2);

run("Stack to Images");

for (i = 1; i <= taille; i++) {
	
	a = IJ.pad(i, 4);

    selectImage("C2-" + nouveau_nom + "-" + a);
    
    nom = "C2_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", direction_dossier_puit_image + nom);
    
    close("C2-" + nouveau_nom + "-"+ a);
};

run("Close All");

};

else {
	
selectWindow(nom_du_fichier);
	
nouveau_nom = "Zstack" + name;

rename(nouveau_nom);

run("Split Channels");

ch1 = "C1-" + nouveau_nom;

close(ch1);

// Renommer et selectioner la fenetre

ch3 = "C3-" + nouveau_nom;

selectWindow(ch3);

taille = nSlices; 

run("Stack to Images");

//Les boucle qui suit, sélection i image et le renome et le sauvegarde pour tantot

for (i = 1; i <= taille; i++) {
	a = IJ.pad(i, 4);
	
    selectImage("C3-" + nouveau_nom + "-"+ a);
    
    nom = "C3_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", direction_dossier_puit_image + nom);
    
    close("C3-" + nouveau_nom + "-"+ a);
};

// Renommer et selectioner la fenetre

ch2 = "C2-" + nouveau_nom;

selectWindow(ch2);

run("Stack to Images");

for (i = 1; i <= taille; i++) {
	
	a = IJ.pad(i, 4);

    selectImage("C2-" + nouveau_nom + "-" + a);
    
    nom = "C2_Zstack" + name + "Image_" + i;
    
    saveAs("Tiff", direction_dossier_puit_image + nom);
    
    close("C2-" + nouveau_nom + "-"+ a);
};

run("Close All");

};


//message d'avertissement

setBatchMode(false);

message = "Image entreposer, la prochaine étape peux être lente, continuer?";
yesLabel = "C'est partie!";
noLabel = "Nope";
getBoolean(message, yesLabel, noLabel);

//initialisation du calcule de médianne

run("Set Measurements...", "limit redirect=None decimal=3");



for (i = 1; i <=taille; i++) {
	
	//ouverture des 2 images à analyser
	
	imagea = direction_dossier_puit_image + "C2_Zstack" + name + "Image_" + i + ".tif";
    imageb = direction_dossier_puit_image + "C3_Zstack" + name + "Image_" + i + ".tif";

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
    
	// ouverture du puling et lencement du calcule de pearson overlap mm cytofluo ica et ccf
	
	
	 run("JACoP ", "imga=[" + "C2_Zstack" + name + "Image_" + i + ".tif" + "] imgb=[" + "C3_Zstack" + name + "Image_" + i + ".tif" + "] thra=" + thresholda + "  thrb=" + thresholdb + " pearson overlap mm ccf=100 cytofluo ica");
	 //Fermeture des fenêtres non utile
	 
	 selectWindow("ICA A (C2_Zstack" + name + "Image_" + i + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_ICA_A" + name + "_Image_" + i + ".csv");
	 
	 close("ICA A (C2_Zstack" + name + "Image_" + i + ".tif)");
	 
	 run("Clear Results");
	 
	 selectWindow("ICA B (C3_Zstack" + name + "Image_" + i + ".tif)");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_ICA_B" + name + "_Image_" + i + ".csv");
	 
	 close("ICA B (C3_Zstack" + name + "Image_" + i + ".tif)");
	 
	 run("Clear Results");
	 
	 selectWindow("Cytofluorogram between C2_Zstack" + name + "Image_" + i + ".tif and C3_Zstack" +  name + "Image_" + i + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_Cytofluogram" + name + "_Image_" + i + ".csv");
	 
	 close("Cytofluorogram between C2_Zstack" + name + "Image_" + i + " and C3_Zstack" + name + "Image_" + i);
	 
	 run("Clear Results");
	 
	 selectWindow("Van Steensel's CCF between C2_Zstack" + name + "Image_" + i + ".tif and C3_Zstack" + name + "Image_" + i + ".tif");
	 
	 Plot.showValues();
	 
	 selectWindow("Results");
	 
	 saveAs("Results", direction_dossier_puit_data_graphique + "/Resultat_VanSteensel" + name + "_Image_" + i + ".csv");
	 
	 close("Van Steensel's CCF between C2_Zstack" + name + "Image_" + i + " and C3_Zstack" + name + "Image_" + i);
	 
	 run("Clear Results");
	 
	 run("Close All");
	 
// --- TUEUR DE FENÊTRE JACOP ---
// Ceci utilise du JavaScript pour fermer les fenêtres que ImageJ ne voit pas
eval("script", "importClass(java.awt.Frame); var frames = Frame.getFrames(); for (var i=0; i<frames.length; i++) { if (frames[i].getTitle().indexOf('Just Another Colocalisation Plugin v2.1.4 21/02/03') != -1) { frames[i].dispose(); } }");
	 
	 
	 //Sélection et entreposage des données des graphs pour R
};


setBatchMode(false);


logContent = getInfo("log");

filePath = direction_dossier_puit_data + "/Coef_coloc_"  + name + ".csv"; 

File.saveString(logContent, filePath);

print("Log content saved to: " + filePath);


//A la suite de ca, j'ai fait un script R pour analyser les informations du fichier
//Donc utiliser le, ou du moins, la fonction de lecture du fichier, je vais encore la paufiné









