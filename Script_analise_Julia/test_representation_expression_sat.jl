#################################################
# Script Julia pour analyser les données
# des transcript ADN satellite  
# et générer des visualisations.               
#################################################

# importer les bibliothèques nécessaires
using DataFrames
using CSV
using Plots
using Makie

# Charger les données depuis un fichier CSV
data_TTAC_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_photo_ovaire_R_ttac.csv", DataFrame)

data_CTAC_aveccaac_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAC_avec_caac.csv", DataFrame)

data_CTAC_avecctat_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAC_avec_ctat.csv", DataFrame)

data_CTAC_avecctat_antisens_rep3 = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAC_avec_ctat_rep3.csv", DataFrame)

data_CTAC_avecttac_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAC_avec_ttac.csv", DataFrame)

data_CTAC_avecttac_antisens_rep2 = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAC_avec_ttac_rep2.csv", DataFrame)

data_CTAC_avecctat_sens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAC_sens_avec_ctat_sens.csv", DataFrame)

data_TTACavecCTAC_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/TTAC_avec_CTAC.csv", DataFrame)

data_TTACavecCTAC_antisens_rep2 = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/TTAC_avec_CTAC_rep2.csv", DataFrame)

data_CTATavecCTAC_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAT_avec_CTAC.csv", DataFrame)

data_CTATavecCTAC_antisens_rep3 = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAT_avec_CTAC_rep3.csv", DataFrame)

data_CTAT_avecCTAC_sens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CTAT_sens_avec_CTAC_sens.csv", DataFrame)

data_CAACavecCTAC_antisens = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/CAAC_avec_CTAC.csv", DataFrame)
  

# ininatialiser un DataFrame pour stocker les résultats

name = ["CTAC_aveccaac_antisens","CTAC_avecttac_antisens","CTAC_avecttac_antisens_rep2",
        "CTAC_avecctat_antisens","CTAC_avecctat_antisens_rep3","TTACavecCTAC_antisens",
        "TTACavecCTAC_antisens_rep2","CTATavecCTAC_antisens","CTATavecCTAC_antisens_rep3",
        "CAACavecCTAC_antisens","CTAC_avecctat_sens","CTAT_avecCTAC_sens"]

# Créer un dictionnaire pour stocker tout les DataFrames avec leurs nom futurs comme clé d'utilisation
# pour la representation des données et les calcules futurs

Dictionaire_des_dataframe = Dict{String,DataFrame}()

for nom in name 
    data_temporaire = DataFrame(
        Null   = zeros(Int64, 10),
        Faible = zeros(Int64, 10),
        Moyen  = zeros(Int64, 10),
        Elever = zeros(Int64, 10)
    )
    Dictionaire_des_dataframe[nom] = data_temporaire
end

print(Dictionaire_des_dataframe["CTAC_avecctat_antisens_rep3"])

# Fonction pour compter les occurrences des catégories

function decompte(data_source, data_puit)

  for i in 3:12 
    index_row = i-2
    colone_source = data_source[:,i]
    #null
    data_puit[index_row,1] = count(x -> !ismissing(x) && (x == "rien"), colone_source)
    #faible
    data_puit[index_row,2] = count(x -> !ismissing(x) && (x == "faible"), colone_source)
    #moyen
    data_puit[index_row,3] = count(x -> !ismissing(x) && (x == "moyen"), colone_source)
    #elever
    data_puit[index_row,4] = count(x -> !ismissing(x) && (x == "elever"), colone_source)
  end
  
  return(data_puit)
end


for nom in name 

    # uilisation de eval et Meta.parse pour convertir
    # le nom de la variable en variable réelle contante au nom des donner initionales

    data_decompte = decompte(eval(Meta.parse("data_"*nom)), Dictionaire_des_dataframe[nom])

end

print(Dictionaire_des_dataframe["CTAC_avecctat_antisens_rep3"])

# Calculer le score d'occurence des transcript
# Mettre les donner en Float64, car les calcule retournent des Float64, étant des chiffres décimaux
# Construction d'un DataFrame pour stocker les scores d'occurence
score_occurence_data_combinaison = DataFrame(score_CTAC_aveccaac_antisens = zeros(Float64, 10),
                                               score_CTAC_avecttac_antisens = zeros(Float64, 10),
                                               score_CTAC_avecctat_antisens = zeros(Float64, 10),
                                               score_TTACavecCTAC_antisens = zeros(Float64, 10),
                                               score_CTATavecCTAC_antisens = zeros(Float64, 10),
                                               score_CAACavecCTAC_antisens = zeros(Float64, 10),
                                               score_CTAC_avecctat_sens = zeros(Float64, 10), 
                                               score_CTAT_avecCTAC_sens = zeros(Float64, 10),
                                               score_CTATavecCTAC_antisens_rep3 = zeros(Float64, 10), 
                                               score_CTAC_avecctat_antisens_rep3 = zeros(Float64, 10),
                                               score_CTAC_avecttat_antisens_rep2 = zeros(Float64, 10),
                                               score_TTACavecCTAC_antisens_rep2 = zeros(Float64, 10),
                                               stade = Float64.(1:10),
                                               row_names = ["stade 1","stade 2","stade 3",
                                                 "stade 4","stade 5","stade 6",
                                                 "stade 7","stade 8","stade 9",
                                                 "stade 10"])

for i in 1:10

  score_occurence_data_combinaison[i,1] = (Dictionaire_des_dataframe["CTAC_aveccaac_antisens"][i,2] + 2*Dictionaire_des_dataframe["CTAC_aveccaac_antisens"][i,3] + 3*Dictionaire_des_dataframe["CTAC_aveccaac_antisens"][i,4])
  score_occurence_data_combinaison[i,2] = (Dictionaire_des_dataframe["CTAC_avecttac_antisens"][i,2] + 2*Dictionaire_des_dataframe["CTAC_avecttac_antisens"][i,3] + 3*Dictionaire_des_dataframe["CTAC_avecttac_antisens"][i,4])
  score_occurence_data_combinaison[i,3] = (Dictionaire_des_dataframe["CTAC_avecctat_antisens"][i,2] + 2*Dictionaire_des_dataframe["CTAC_avecctat_antisens"][i,3] + 3*Dictionaire_des_dataframe["CTAC_avecctat_antisens"][i,4])
  score_occurence_data_combinaison[i,4] = (Dictionaire_des_dataframe["TTACavecCTAC_antisens"][i,2] + 2*Dictionaire_des_dataframe["TTACavecCTAC_antisens"][i,3] + 3*Dictionaire_des_dataframe["TTACavecCTAC_antisens"][i,4])
  score_occurence_data_combinaison[i,5] = (Dictionaire_des_dataframe["CTATavecCTAC_antisens"][i,2] + 2*Dictionaire_des_dataframe["CTATavecCTAC_antisens"][i,3] + 3*Dictionaire_des_dataframe["CTATavecCTAC_antisens"][i,4])
  score_occurence_data_combinaison[i,6] = (Dictionaire_des_dataframe["CAACavecCTAC_antisens"][i,2] + 2*Dictionaire_des_dataframe["CAACavecCTAC_antisens"][i,3] + 3*Dictionaire_des_dataframe["CAACavecCTAC_antisens"][i,4])
  score_occurence_data_combinaison[i,7] = (Dictionaire_des_dataframe["CTAC_avecctat_sens"][i,2] + 2*Dictionaire_des_dataframe["CTAC_avecctat_sens"][i,3] + 3*Dictionaire_des_dataframe["CTAC_avecctat_sens"][i,4])
  score_occurence_data_combinaison[i,8] = (Dictionaire_des_dataframe["CTAT_avecCTAC_sens"][i,2] + 2*Dictionaire_des_dataframe["CTAT_avecCTAC_sens"][i,3] + 3*Dictionaire_des_dataframe["CTAT_avecCTAC_sens"][i,4])
  score_occurence_data_combinaison[i,9] = (Dictionaire_des_dataframe["CTATavecCTAC_antisens_rep3"][i,2] + 2*Dictionaire_des_dataframe["CTATavecCTAC_antisens_rep3"][i,3] + 3*Dictionaire_des_dataframe["CTATavecCTAC_antisens_rep3"][i,4])
  score_occurence_data_combinaison[i,10] = (Dictionaire_des_dataframe["CTAC_avecctat_antisens_rep3"][i,2] + 2*Dictionaire_des_dataframe["CTAC_avecctat_antisens_rep3"][i,3] + 3*Dictionaire_des_dataframe["CTAC_avecctat_antisens_rep3"][i,4])
  score_occurence_data_combinaison[i,11] = (Dictionaire_des_dataframe["CTAC_avecttac_antisens_rep2"][i,2] + 2*Dictionaire_des_dataframe["CTAC_avecttac_antisens_rep2"][i,3] + 3*Dictionaire_des_dataframe["CTAC_avecttac_antisens_rep2"][i,4])
  score_occurence_data_combinaison[i,12] = (Dictionaire_des_dataframe["TTACavecCTAC_antisens_rep2"][i,2] + 2*Dictionaire_des_dataframe["TTACavecCTAC_antisens_rep2"][i,3] + 3*Dictionaire_des_dataframe["TTACavecCTAC_antisens_rep2"][i,4])

end

# Normaliser le score d'occurence par le nombre total d'occurence

score_occurence_data_normaliser_combinaison = score_occurence_data_combinaison

for i in 1:10
  
  score_occurence_data_normaliser_combinaison[i,1] = score_occurence_data_normaliser_combinaison[i,1]/(sum(Dictionaire_des_dataframe["CTAC_aveccaac_antisens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,2] = score_occurence_data_normaliser_combinaison[i,2]/(sum(Dictionaire_des_dataframe["CTAC_avecttac_antisens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,3] = score_occurence_data_normaliser_combinaison[i,3]/(sum(Dictionaire_des_dataframe["CTAC_avecctat_antisens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,4] = score_occurence_data_normaliser_combinaison[i,4]/(sum(Dictionaire_des_dataframe["TTACavecCTAC_antisens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,5] = score_occurence_data_normaliser_combinaison[i,5]/(sum(Dictionaire_des_dataframe["CTATavecCTAC_antisens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,6] = score_occurence_data_normaliser_combinaison[i,6]/(sum(Dictionaire_des_dataframe["CAACavecCTAC_antisens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,7] = score_occurence_data_normaliser_combinaison[i,7]/(sum(Dictionaire_des_dataframe["CTAC_avecctat_sens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,8] = score_occurence_data_normaliser_combinaison[i,8]/(sum(Dictionaire_des_dataframe["CTAT_avecCTAC_sens"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,9] = score_occurence_data_normaliser_combinaison[i,9]/(sum(Dictionaire_des_dataframe["CTATavecCTAC_antisens_rep3"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,10] = score_occurence_data_normaliser_combinaison[i,10]/(sum(Dictionaire_des_dataframe["CTAC_avecctat_antisens_rep3"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,11] = score_occurence_data_normaliser_combinaison[i,11]/(sum(Dictionaire_des_dataframe["CTAC_avecttac_antisens_rep2"][i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[i,12] = score_occurence_data_normaliser_combinaison[i,12]/(sum(Dictionaire_des_dataframe["TTACavecCTAC_antisens_rep2"][i,1:4]) +1)
  
end


# Visualiser les résultats

plot(1:10, score_occurence_data_normaliser_combinaison[:,1],
 xlabel="Stade", ylabel="Score d'occurence CTAC normalisé", 
 title="Score d'occurence CTAC normalisé par Stade", 
 legend=false, marker=:circle)
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,2])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,3])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,4])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,5])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,6])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,7])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,8])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,9])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,10])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,11])
    plot!(1:10, score_occurence_data_normaliser_combinaison[:,12])


