##############################################
#Script exploratoire de visualisation, des normalité ainsi
#que de la variance des données des graphiques
##############################################

using DataFrames
using CSV
using Plots
using Makie
using CairoMakie
using GLMakie
using Statistics
using GLM
using HypothesisTests
using Distributions
using Printf

#######################################################
#Chargement des données cytofluogtamme sens et antisens 
#######################################################

#sonde antisens
#numero des image à analyser pour chaque stade

numero_image_antisens_stade3_rep1 = (6, 21, 22, 23, 29, 32, 34, 36, 39, 47, 63, 66, 70, 72, 73, 78, 87, 90)
numero_image_antisens_stade4_rep1 = (9, 12, 17, 19, 20, 28, 39, 53)
numero_image_antisens_stade5_rep1 = (1, 4, 7, 9, 15, 16, 22, 23, 25, 26, 31, 32, 34, 37, 40, 42, 44, 45, 47, 57)
numero_image_antisens_stade6_rep1 = (8, 18, 27, 46, 59, 75, 96, 97, 108, 109)
numero_image_antisens_stade7_rep1 = (2, 3, 14, 24, 27, 30, 33, 35, 41, 43, 44, 50, 58, 60, 62, 68, 110)
numero_image_antisens_stade8_rep1 = (10, 49, 52, 55, 84, 89, 111)
numero_image_antisens_stade9_rep1 = (51, 61, 105)
numero_image_antisens_stade10_rep1 = (56, 67, 82, 83)

##############################################
#numero des image à analyser pour chaque stade
#pour la réplication 1

numero_image_antisens_stade3_rep2 = (5, 11, 19, 36, 43, 103)
numero_image_antisens_stade4_rep2 = (13, 77)
numero_image_antisens_stade5_rep2 = (4, 7, 10, 13, 19, 20, 51, 56, 59, 65, 69, 75, 83, 89, 91, 99, 102, 104)
numero_image_antisens_stade6_rep2 = (1, 9, 33, 35, 41, 51, 56, 64, 88, 91, 92, 94, 96)
numero_image_antisens_stade7_rep2 = (15, 18, 26, 29, 30, 34, 42, 46, 57, 68, 70, 74, 87)
numero_image_antisens_stade8_rep2 = (2, 26, 27, 31, 40, 47, 50, 52, 53, 54, 60, 67, 73, 86, 87, 89, 100)
numero_image_antisens_stade9_rep2 = (3, 6, 12, 22, 25, 32, 37, 39, 44, 45, 61, 63, 66, 72, 85, 90, 95)
numero_image_antisens_stade10_rep2 = (14, 16, 17, 38, 48, 49, 62, 81, 84)
##############################################
#création des dictionnaires pour stocker les données des graphiques

dictonaire_des_graphique_cytofluogramme_antisens_stade3 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade4 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade5 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade6 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade7 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade8 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade9 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_antisens_stade10 = Dict{String,DataFrame}()

###############################################
#boucle de lecture des données des graphiques pour les stdae 
#et concaténation des répetion dans les dictinnaire assocers

#pour le stade 3

for i in numero_image_antisens_stade3_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade3["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_rep1_stade3_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade3_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade3["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade3_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 4
for i in numero_image_antisens_stade4_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade4["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antosensboth_stade4_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade4_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade4["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade4_rep2_Image_" * string(i) * ".csv", DataFrame)
end
##############################################
#pour le stade 5
for i in numero_image_antisens_stade5_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade5["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antosensboth_stade5_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade5_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade5["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade5_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 6
for i in numero_image_antisens_stade6_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade6["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade6_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade6_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade6["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade6_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 7

for i in numero_image_antisens_stade7_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade7["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade7_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade7_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade7["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade7_rep2_Image_" * string(i) * ".csv", DataFrame)
end
################################################
#pour le stade 8
for i in numero_image_antisens_stade8_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade8["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade8_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade8_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade8["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade8_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 9
for i in numero_image_antisens_stade9_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade9["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade9_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade9_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade9["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade9_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 10
for i in numero_image_antisens_stade10_rep1
    dictonaire_des_graphique_cytofluogramme_antisens_stade10["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade10_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_antisens_stade10_rep2
    dictonaire_des_graphique_cytofluogramme_antisens_stade10["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade10_rep2_Image_" * string(i) * ".csv", DataFrame)
end


################################################
#Fonction pour mesurer la longueur maximale des
#données des cytofluogrammes
#pour uniformiser la longueur des données
#avant de faire la moyenne

function mesure_max(dictionaire, liste_numero_image, liste_numero_image_rep1)

    maximum_length = Array{Union{Missing,Float64}}(missing, length(dictionaire))

    maximum_length_rep1 = Array{Union{Missing,Float64}}(missing, length(dictionaire))

    for i in liste_numero_image

        n_injection = findfirst(x -> x == i, liste_numero_image)

        maximum_length[n_injection] = length(dictionaire["cytofluogramme_image_"*string(i)][!, :X0])

    end

    for j in liste_numero_image_rep1

        n_injection_rep1 = findfirst(x -> x == j, liste_numero_image_rep1)

        maximum_length_rep1[n_injection_rep1] = length(dictionaire["cytofluogramme_rep1_image_"*string(j)][!, :X0])
    end

    maximum_length = maximum(skipmissing(vcat(maximum_length, maximum_length_rep1)))

    return maximum_length

end


max_length_stade3 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade3, numero_image_antisens_stade3_rep1, numero_image_antisens_stade3_rep2)
max_length_stade4 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade4, numero_image_antisens_stade4_rep1, numero_image_antisens_stade4_rep2)
max_length_stade5 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade5, numero_image_antisens_stade5_rep1, numero_image_antisens_stade5_rep2)
max_length_stade6 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade6, numero_image_antisens_stade6_rep1, numero_image_antisens_stade6_rep2)
max_length_stade7 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade7, numero_image_antisens_stade7_rep1, numero_image_antisens_stade7_rep2)
max_length_stade8 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade8, numero_image_antisens_stade8_rep1, numero_image_antisens_stade8_rep2)
max_length_stade9 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade9, numero_image_antisens_stade9_rep1, numero_image_antisens_stade9_rep2)
max_length_stade10 = mesure_max(dictonaire_des_graphique_cytofluogramme_antisens_stade10, numero_image_antisens_stade10_rep1, numero_image_antisens_stade10_rep2)

##############################################
#Construiction des dataframes pour le calcul des moyennes

X0_data_antisens_stade3 = DataFrame()
Y0_data_antisens_stade3 = DataFrame()

X0_data_antisens_stade4 = DataFrame()
Y0_data_antisens_stade4 = DataFrame()

X0_data_antisens_stade5 = DataFrame()
Y0_data_antisens_stade5 = DataFrame()

X0_data_antisens_stade6 = DataFrame()
Y0_data_antisens_stade6 = DataFrame()

X0_data_antisens_stade7 = DataFrame()
Y0_data_antisens_stade7 = DataFrame()

X0_data_antisens_stade8 = DataFrame()
Y0_data_antisens_stade8 = DataFrame()

X0_data_antisens_stade9 = DataFrame()
Y0_data_antisens_stade9 = DataFrame()

X0_data_antisens_stade10 = DataFrame()
Y0_data_antisens_stade10 = DataFrame()

##############################################
#Fonction pour ajuster la taille des données des cytofluogrammes
#avant de faire la moyenne

function ajustement_taille(dictionaire, liste_numero_image, liste_numero_image_rep1, max_length, X0_data, Y0_data)

    for i in liste_numero_image

        if length(dictionaire["cytofluogramme_image_"*string(i)][!, :X0]) < max_length

            allowmissing!(dictionaire["cytofluogramme_image_"*string(i)], [:X0, :Y0])
            difference_length = max_length - length(dictionaire["cytofluogramme_image_"*string(i)][!, :X0])

            for j in 1:difference_length

                push!(dictionaire["cytofluogramme_image_"*string(i)][!, :X0], missing)
                push!(dictionaire["cytofluogramme_image_"*string(i)][!, :Y0], missing)
            end

        end
        n_injection = Symbol(findfirst(x -> x == i, liste_numero_image))

        X0_data[!, n_injection] = dictionaire["cytofluogramme_image_"*string(i)][!, :X0]
        Y0_data[!, n_injection] = dictionaire["cytofluogramme_image_"*string(i)][!, :Y0]
    end

    taille = length(X0_data[1, :])

    for i in liste_numero_image_rep1

        if length(dictionaire["cytofluogramme_rep1_image_"*string(i)][!, :X0]) < max_length

            allowmissing!(dictionaire["cytofluogramme_rep1_image_"*string(i)], [:X0, :Y0])
            difference_length = max_length - length(dictionaire["cytofluogramme_rep1_image_"*string(i)][!, :X0])

            for j in 1:difference_length

                push!(dictionaire["cytofluogramme_rep1_image_"*string(i)][!, :X0], missing)
                push!(dictionaire["cytofluogramme_rep1_image_"*string(i)][!, :Y0], missing)
            end

        end

        n_injection_rep1 = Int64(findfirst(x -> x == i, liste_numero_image_rep1))

        a_partir_de = Symbol((Int64(taille) + n_injection_rep1))

        X0_data[!, a_partir_de] = dictionaire["cytofluogramme_rep1_image_"*string(i)][!, :X0]
        Y0_data[!, a_partir_de] = dictionaire["cytofluogramme_rep1_image_"*string(i)][!, :Y0]
    end

end

#ajustement des taille et remplissage des dataframes

ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade3, numero_image_antisens_stade3_rep1, numero_image_antisens_stade3_rep2, max_length_stade3, X0_data_antisens_stade3, Y0_data_antisens_stade3)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade4, numero_image_antisens_stade4_rep1, numero_image_antisens_stade4_rep2, max_length_stade4, X0_data_antisens_stade4, Y0_data_antisens_stade4)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade5, numero_image_antisens_stade5_rep1, numero_image_antisens_stade5_rep2, max_length_stade5, X0_data_antisens_stade5, Y0_data_antisens_stade5)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade6, numero_image_antisens_stade6_rep1, numero_image_antisens_stade6_rep2, max_length_stade6, X0_data_antisens_stade6, Y0_data_antisens_stade6)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade7, numero_image_antisens_stade7_rep1, numero_image_antisens_stade7_rep2, max_length_stade7, X0_data_antisens_stade7, Y0_data_antisens_stade7)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade8, numero_image_antisens_stade8_rep1, numero_image_antisens_stade8_rep2, max_length_stade8, X0_data_antisens_stade8, Y0_data_antisens_stade8)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade9, numero_image_antisens_stade9_rep1, numero_image_antisens_stade9_rep2, max_length_stade9, X0_data_antisens_stade9, Y0_data_antisens_stade9)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_antisens_stade10, numero_image_antisens_stade10_rep1, numero_image_antisens_stade10_rep2, max_length_stade10, X0_data_antisens_stade10, Y0_data_antisens_stade10)

##############################################
#Calcul des moyennes des cytofluogrammes
#initialisation des dataframes pour stocker les moyennes

resultats_moyennes_antisens = Dict{Int64,DataFrame}()

for i in (3:10)

    X0_data = eval(Symbol("X0_data_antisens_stade$i"))
    Y0_data = eval(Symbol("Y0_data_antisens_stade$i"))

    x_moyens = [
        begin
            m = skipmissing(row)
            isempty(m) ? NaN : mean(m)
        end for row in eachrow(X0_data)
    ]

    y_moyens = [
        begin
            m = skipmissing(row)
            isempty(m) ? NaN : mean(m)
        end for row in eachrow(Y0_data)
    ]

    # 3. On crée le DataFrame final pour ce stade et on le stocke dans le dictionnaire
    DataFrame_temporaire = DataFrame(X0_moyen=x_moyens, Y0_moyen=y_moyens)

    filter!(row -> !isnan(row.X0_moyen) && !isnan(row.Y0_moyen), DataFrame_temporaire)

    resultats_moyennes_antisens[i] = DataFrame_temporaire
end

################################################
################################################
#sonde sens

#numero des image à analyser pour chaque stade

numero_image_sens_stade3 = (21, 34, 35, 36, 37, 52, 53, 55, 61, 75, 98)
numero_image_sens_stade4 = (13, 39, 44, 46, 48, 55)
numero_image_sens_stade5 = (12, 17, 22, 23, 29, 30, 37, 45, 51, 54, 55, 68, 69, 71, 74, 78, 80, 82, 83)
numero_image_sens_stade6 = (1, 3, 11, 22, 34, 49, 50, 63, 67, 70, 72, 81, 87, 93, 100, 101, 105)
numero_image_sens_stade7 = (5, 6, 7, 9, 11, 14, 16, 18, 20, 24, 28, 31, 38, 40, 54, 57, 58, 73, 76, 77, 84, 90, 91, 94)
numero_image_sens_stade8 = (2, 4, 10, 14, 19, 32, 33, 64, 65, 85, 104, 107, 113, 116)
numero_image_sens_stade9 = (8, 41, 59, 88, 92, 99, 103, 109, 111)
numero_image_sens_stade10 = (25, 26, 27, 102, 110, 114, 118)

##############################################
#numero des image à analyser pour chaque stade
#pour la réplication 1

numero_image_sens_stade3_rep1 = (4, 16, 20, 21, 27, 35, 39, 59, 66, 85, 91, 93)
numero_image_sens_stade4_rep1 = (1, 90, 92)
numero_image_sens_stade5_rep1 = (4, 10, 13, 13, 14, 16, 20, 21, 26, 30, 32, 33, 42, 48, 58, 66, 68, 69, 70, 72, 78, 83, 87, 90, 92, 94, 95)
numero_image_sens_stade6_rep1 = (5, 8, 9, 22, 34, 37, 41, 46, 60, 67, 79, 85)
numero_image_sens_stade7_rep1 = (2, 18, 19, 24, 25, 31, 38, 40, 44, 45, 50, 55, 56, 74, 75, 76, 83, 88, 89)
numero_image_sens_stade8_rep1 = (3, 6, 7, 12, 17, 20, 28, 29, 36, 43, 45, 48, 49, 51, 52, 57, 64, 70, 77, 80, 81, 82, 84, 86)
numero_image_sens_stade9_rep1 = (23, 53, 54, 61, 73)
numero_image_sens_stade10_rep1 = (11, 15, 47, 62, 63, 65)

##############################################
#création des dictionnaires pour stocker les données des graphiques

dictonaire_des_graphique_cytofluogramme_sens_stade3 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade4 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade5 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade6 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade7 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade8 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade9 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_sens_stade10 = Dict{String,DataFrame}()

###############################################
#boucle de lecture des données des graphiques pour les stdae 
#et concaténation des répetion dans les dictinnaire assocers

#pour le stade 3

for i in numero_image_sens_stade3
    dictonaire_des_graphique_cytofluogramme_sens_stade3["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_sade3_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade3_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade3["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade3_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 4
for i in numero_image_sens_stade4
    dictonaire_des_graphique_cytofluogramme_sens_stade4["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAC_sensboth_stade4_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade4_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade4["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade4_rep1_Image_" * string(i) * ".csv", DataFrame)
end
##############################################
#pour le stade 5
for i in numero_image_sens_stade5
    dictonaire_des_graphique_cytofluogramme_sens_stade5["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade5_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade5_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade5["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade5_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 6
for i in numero_image_sens_stade6
    dictonaire_des_graphique_cytofluogramme_sens_stade6["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade6_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade6_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade6["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade6_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 7
for i in numero_image_sens_stade7
    dictonaire_des_graphique_cytofluogramme_sens_stade7["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade7_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade7_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade7["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade7_rep1_Image_" * string(i) * ".csv", DataFrame)
end
################################################
#pour le stade 8
for i in numero_image_sens_stade8
    dictonaire_des_graphique_cytofluogramme_sens_stade8["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade8_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade8_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade8["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade8_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 9
for i in numero_image_sens_stade9
    dictonaire_des_graphique_cytofluogramme_sens_stade9["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade9_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade9_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade9["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade9_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 10
for i in numero_image_sens_stade10
    dictonaire_des_graphique_cytofluogramme_sens_stade10["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade10_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_sens_stade10_rep1
    dictonaire_des_graphique_cytofluogramme_sens_stade10["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade10_rep1_Image_" * string(i) * ".csv", DataFrame)
end

################################################
#Fonction pour mesurer la longueur maximale des
#données des cytofluogrammes
#pour uniformiser la longueur des données
#avant de faire la moyenne


max_length_sens_stade3 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade3, numero_image_sens_stade3, numero_image_sens_stade3_rep1)
max_length_sens_stade4 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade4, numero_image_sens_stade4, numero_image_sens_stade4_rep1)
max_length_sens_stade5 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade5, numero_image_sens_stade5, numero_image_sens_stade5_rep1)
max_length_sens_stade6 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade6, numero_image_sens_stade6, numero_image_sens_stade6_rep1)
max_length_sens_stade7 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade7, numero_image_sens_stade7, numero_image_sens_stade7_rep1)
max_length_sens_stade8 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade8, numero_image_sens_stade8, numero_image_sens_stade8_rep1)
max_length_sens_stade9 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade9, numero_image_sens_stade9, numero_image_sens_stade9_rep1)
max_length_sens_stade10 = mesure_max(dictonaire_des_graphique_cytofluogramme_sens_stade10, numero_image_sens_stade10, numero_image_sens_stade10_rep1)

##############################################
#Construiction des dataframes pour le calcul des moyennes

X0_data_sens_stade3 = DataFrame()
Y0_data_sens_stade3 = DataFrame()

X0_data_sens_stade4 = DataFrame()
Y0_data_sens_stade4 = DataFrame()

X0_data_sens_stade5 = DataFrame()
Y0_data_sens_stade5 = DataFrame()

X0_data_sens_stade6 = DataFrame()
Y0_data_sens_stade6 = DataFrame()

X0_data_sens_stade7 = DataFrame()
Y0_data_sens_stade7 = DataFrame()

X0_data_sens_stade8 = DataFrame()
Y0_data_sens_stade8 = DataFrame()

X0_data_sens_stade9 = DataFrame()
Y0_data_sens_stade9 = DataFrame()

X0_data_sens_stade10 = DataFrame()
Y0_data_sens_stade10 = DataFrame()

##############################################
#Fonction pour ajuster la taille des données des cytofluogrammes
#avant de faire la moyenne


#ajustement des taille et remplissage des dataframes

ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade3, numero_image_sens_stade3, numero_image_sens_stade3_rep1, max_length_sens_stade3, X0_data_sens_stade3, Y0_data_sens_stade3)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade4, numero_image_sens_stade4, numero_image_sens_stade4_rep1, max_length_sens_stade4, X0_data_sens_stade4, Y0_data_sens_stade4)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade5, numero_image_sens_stade5, numero_image_sens_stade5_rep1, max_length_sens_stade5, X0_data_sens_stade5, Y0_data_sens_stade5)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade6, numero_image_sens_stade6, numero_image_sens_stade6_rep1, max_length_sens_stade6, X0_data_sens_stade6, Y0_data_sens_stade6)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade7, numero_image_sens_stade7, numero_image_sens_stade7_rep1, max_length_sens_stade7, X0_data_sens_stade7, Y0_data_sens_stade7)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade8, numero_image_sens_stade8, numero_image_sens_stade8_rep1, max_length_sens_stade8, X0_data_sens_stade8, Y0_data_sens_stade8)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade9, numero_image_sens_stade9, numero_image_sens_stade9_rep1, max_length_sens_stade9, X0_data_sens_stade9, Y0_data_sens_stade9)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_sens_stade10, numero_image_sens_stade10, numero_image_sens_stade10_rep1, max_length_sens_stade10, X0_data_sens_stade10, Y0_data_sens_stade10)

##############################################
#Calcul des moyennes des cytofluogrammes
#initialisation des dataframes pour stocker les moyennes

resultats_moyennes_sens = Dict{Int64,DataFrame}()

for i in (3:10)

    X0_data = eval(Symbol("X0_data_sens_stade$i"))
    Y0_data = eval(Symbol("Y0_data_sens_stade$i"))

    x_moyens = [
        begin
            m = skipmissing(row)
            isempty(m) ? NaN : mean(m)
        end for row in eachrow(X0_data)
    ]

    y_moyens = [
        begin
            m = skipmissing(row)
            isempty(m) ? NaN : mean(m)
        end for row in eachrow(Y0_data)
    ]

    # 3. On crée le DataFrame final pour ce stade et on le stocke dans le dictionnaire
    DataFrame_temporaire = DataFrame(X0_moyen=x_moyens, Y0_moyen=y_moyens)

    filter!(row -> !isnan(row.X0_moyen) && !isnan(row.Y0_moyen), DataFrame_temporaire)

    resultats_moyennes_sens[i] = DataFrame_temporaire
end

####################################################################################
####################################################################################

#Vérification des normalité des données des X0 et Y0

#donner aberante X0 et Y0

fig_aberante_x0 = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_antisens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_aberante_x0[1, n_injection],
            title="Plot for x0  = CTAT of stage " * string(i),
            xlabel="Sonde antisens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_aberante_x0[2, n_injection-2],
            title="Plot for x0  = CTAT of stage " * string(i),
            xlabel="Sonde antisens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_aberante_x0[3, n_injection-4],
            title="Plot for x0  = CTAT of stage " * string(i),
            xlabel="Sonde antisens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_aberante_x0[4, n_injection-6],
            title="Plot for x0  = CTAT of stage " * string(i),
            xlabel="Sonde antisens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    Makie.scatter!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen])

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)

end

supertitle = Label(fig_aberante_x0[0, :],
    "Scatter plot for X0  = CTAT in stages 3 to 10 antisens",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

#chemin = joinpath("/Users", "verme", "Desktop", "Scatter_for_X0_sonde_antisens.png")
#save(chemin, fig_aberante_x0)


fig_aberante_y0 = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_antisens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_aberante_y0[1, n_injection],
            title="Plot for y0 = CTAC of stage " * string(i),
            xlabel="Sonde antisens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_aberante_y0[2, n_injection-2],
            title="Plot for y0 = CTAC of stage " * string(i),
            xlabel="Sonde antisens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_aberante_y0[3, n_injection-4],
            title="Plot for y0 = CTAC of stage " * string(i),
            xlabel="Sonde antisens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_aberante_y0[4, n_injection-6],
            title="Plot for y0 = CTAC of stage " * string(i),
            xlabel="Sonde antisens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    Makie.scatter!(ax, moyenne_cytofluogramme_stade[!, :Y0_moyen])

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)

end

supertitle = Label(fig_aberante_y0[0, :],
    "Scatter plot for Y0 = CTAC in stages 3 to 10 antisens",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "Scatter_for_Y0_sonde_antisens.png")
save(chemin, fig_aberante_y0)

##############################################

fig_aberante_x0 = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_sens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_aberante_x0[1, n_injection],
            title="Plot for x0 = CTAC of stage " * string(i),
            xlabel="Sonde sens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_aberante_x0[2, n_injection-2],
            title="Plot for x0 = CTAC of stage " * string(i),
            xlabel="Sonde sens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_aberante_x0[3, n_injection-4],
            title="Plot for x0 = CTAC of stage " * string(i),
            xlabel="Sonde sens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_aberante_x0[4, n_injection-6],
            title="Plot for x0 = CTAC of stage " * string(i),
            xlabel="Sonde sens CTAC",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    Makie.scatter!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen])

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)

end

supertitle = Label(fig_aberante_x0[0, :],
    "Scatter plot for X0 = CTAC in stages 3 to 10 sens",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "Scatter_for_X0_sonde_sens.png")
save(chemin, fig_aberante_x0)


fig_aberante_y0 = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_sens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_aberante_y0[1, n_injection],
            title="Plot for y0 = CTAT of stage " * string(i),
            xlabel="Sonde sens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_aberante_y0[2, n_injection-2],
            title="Plot for y0 = CTAT of stage " * string(i),
            xlabel="Sonde sens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_aberante_y0[3, n_injection-4],
            title="Plot for y0 = CTAT of stage " * string(i),
            xlabel="Sonde sens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_aberante_y0[4, n_injection-6],
            title="Plot for y0 = CTAT of stage " * string(i),
            xlabel="Sonde sens CTAT",
            ylabel="Moyenne des intensités",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    Makie.scatter!(ax, moyenne_cytofluogramme_stade[!, :Y0_moyen])

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)

end

supertitle = Label(fig_aberante_y0[0, :],
    "Scatter plot for Y0 = CTAT in stages 3 to 10 sens",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "Scatter_for_Y0_sonde_sens.png")
save(chemin, fig_aberante_y0)

####################################################################################
####################################################################################

# qqnorm pour X0

fig_qqnorm_x0_antisens = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_antisens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_qqnorm_x0_antisens[1, n_injection],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_qqnorm_x0_antisens[2, n_injection-2],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_qqnorm_x0_antisens[3, n_injection-4],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_qqnorm_x0_antisens[4, n_injection-6],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    qqplot!(ax, Normal, moyenne_cytofluogramme_stade[!, :X0_moyen], qqline=:identity)
    test = ShapiroWilkTest(sort(moyenne_cytofluogramme_stade[!, :X0_moyen]))
    p_value = pvalue(test)
    p_value_string = @sprintf("%.3e", p_value)
    text!(ax, p_value_string,
        position=(0.5, 0.5),
        align=(:center, :center),
        color=:black)

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)

end


supertitle = Label(fig_qqnorm_x0_antisens[0, :],
    "QQ plot of X0 = CTAT sonde antisens in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "qqnorm_for_X0_sonde_antisens.png")
save(chemin, fig_qqnorm_x0_antisens)

####################################################
#Test de fit avec Gamma, gaussinne inverse ou encore Weibull

#function ajustement_des_zéro(data)

#data.X0_moyen .= ifelse.(data.X0_moyen .<= 1e-4, data.X0_moyen .+ 1e-3, data.X0_moyen)
#data.Y0_moyen .= ifelse.(data.Y0_moyen .<= 1e-4, data.Y0_moyen .+ 1e-3, data.Y0_moyen)

#end

#for i in (3:10)
# ajustement_des_zéro(resultats_moyennes_sens[i])

#println(count(skipmissing(resultats_moyennes_sens[i][!, :Y0_moyen] .<= 1e-3)))
#end

#dist = fit(Weibull, resultats_moyennes_sens[3][!, :Y0_moyen])

#qqplot(resultats_moyennes_sens[3][!, :Y0_moyen], dist, qqline=:identity)

#ExactOneSampleKSTest(resultats_moyennes_sens[3][!, :Y0_moyen], dist)


#model = glm(@formula(Y0_moyen ~ X0_moyen), resultats_moyennes_sens[3], Gamma(), LogLink())

#summary(model)

# 2. Utiliser la fonction save de Makie
#save(chemin, fig_qqnorm_x0)


# qqnorm pour Y0

fig_qqnorm_y0_antisens = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_antisens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_qqnorm_y0_antisens[1, n_injection],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_qqnorm_y0_antisens[2, n_injection-2],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_qqnorm_y0_antisens[3, n_injection-4],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_qqnorm_y0_antisens[4, n_injection-6],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    qqplot!(ax, Normal, moyenne_cytofluogramme_stade[!, :Y0_moyen], qqline=:identity)
    test = ShapiroWilkTest(sort(moyenne_cytofluogramme_stade[!, :Y0_moyen]))
    p_value = pvalue(test)
    p_value_string = @sprintf("%.3e", p_value)
    text!(ax, p_value_string,
        position=(0.5, 0.5),
        align=(:center, :center),
        color=:black)

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :Y0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)
end

supertitle = Label(fig_qqnorm_y0_antisens[0, :],
    "QQ plot of Y0 = CTAC sonde antisens in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "qqnorm_for_Y0_sonde_antisens.png")
save(chemin, fig_qqnorm_y0_antisens)

#qqnorm pour X0 sens

fig_qqnorm_x0_sens = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_sens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_qqnorm_x0_sens[1, n_injection],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_qqnorm_x0_sens[2, n_injection-2],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_qqnorm_x0_sens[3, n_injection-4],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_qqnorm_x0_sens[4, n_injection-6],
            title="QQ plot for x0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    qqplot!(ax, Normal, moyenne_cytofluogramme_stade[!, :X0_moyen], qqline=:identity)
    test = ShapiroWilkTest(sort(moyenne_cytofluogramme_stade[!, :X0_moyen]))
    p_value = pvalue(test)
    p_value_string = @sprintf("%.3e", p_value)
    text!(ax, p_value_string,
        position=(0.5, 0.5),
        align=(:center, :center),
        color=:black)

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :X0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)

end


supertitle = Label(fig_qqnorm_x0_sens[0, :],
    "QQ plot of X0 = CTAC sonde sens in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "qqnorm_for_X0_sonde_sens.png")
save(chemin, fig_qqnorm_x0_sens)


# qqnorm pour Y0 sens

fig_qqnorm_y0_sens = Figure(size=(1200, 1600))

for i in (3:10)

    moyenne_cytofluogramme_stade = resultats_moyennes_sens[i]

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    if n_injection <= 2
        ax = Axis(fig_qqnorm_y0_sens[1, n_injection],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_qqnorm_y0_sens[2, n_injection-2],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))


    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_qqnorm_y0_sens[3, n_injection-4],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_qqnorm_y0_sens[4, n_injection-6],
            title="QQ plot for y0 of stage " * string(i),
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end


    qqplot!(ax, Normal, moyenne_cytofluogramme_stade[!, :Y0_moyen], qqline=:identity)
    test = ShapiroWilkTest(sort(moyenne_cytofluogramme_stade[!, :Y0_moyen]))
    p_value = pvalue(test)
    p_value_string = @sprintf("%.3e", p_value)
    text!(ax, p_value_string,
        position=(0.5, 0.5),
        align=(:center, :center),
        color=:black)

    #Makie.hist!(ax, moyenne_cytofluogramme_stade[!, :Y0_moyen],
    #    bins=500,
    #    strokewidth=0.5,
    #    strokecolor=:black,
    #    color=:blue)
end

supertitle = Label(fig_qqnorm_y0_sens[0, :],
    "QQ plot of Y0 = CTAT sonde sens in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "qqnorm_for_Y0_sonde_sens.png")
save(chemin, fig_qqnorm_y0_sens)

#################################################################
#Vérification des variance des X0 et Y0

for i in (3:10)

    print(LeveneTest(resultats_moyennes_sens[i][!, :X0_moyen], resultats_moyennes_sens[i][!, :Y0_moyen]))

end

for i in (3:10)

    print(LeveneTest(resultats_moyennes_antisens[i][!, :X0_moyen], resultats_moyennes_antisens[i][!, :Y0_moyen]))

end


levene_test_antisens_X0 = []
levene_test_antisens_Y0 = []
levene_test_sens_X0 = []
levene_test_sens_Y0 = []

for i in (3:10)

    for j in (3:10)

        if i != j

            push!(levene_test_antisens_X0, pvalue(LeveneTest(resultats_moyennes_antisens[i][!, :X0_moyen], resultats_moyennes_antisens[j][!, :X0_moyen])))

            push!(levene_test_antisens_Y0, pvalue(LeveneTest(resultats_moyennes_antisens[i][!, :Y0_moyen], resultats_moyennes_antisens[j][!, :Y0_moyen])))

            push!(levene_test_sens_X0, pvalue(LeveneTest(resultats_moyennes_sens[i][!, :X0_moyen], resultats_moyennes_sens[j][!, :X0_moyen])))

            push!(levene_test_sens_Y0, pvalue(LeveneTest(resultats_moyennes_sens[i][!, :Y0_moyen], resultats_moyennes_sens[j][!, :Y0_moyen])))

        end
    end
end

count(p -> p > 0.05, levene_test_antisens_X0)

count(p -> p > 0.05, levene_test_antisens_Y0)

count(p -> p > 0.05, levene_test_sens_X0)

count(p -> p > 0.05, levene_test_sens_Y0)

#decompte des 0

pourcentage_X0_antisens = []
pourcentage_Y0_antisens = []
pourcentage_X0_sens = []
pourcentage_Y0_sens = []

for i in (3:10)

    push!(pourcentage_X0_antisens, count(skipmissing(resultats_moyennes_antisens[i][!, :X0_moyen] .<= 0)) / length(resultats_moyennes_antisens[i][!, :X0_moyen]))

    push!(pourcentage_Y0_antisens, count(skipmissing(resultats_moyennes_antisens[i][!, :Y0_moyen] .<= 0)) / length(resultats_moyennes_antisens[i][!, :Y0_moyen]))

    push!(pourcentage_X0_sens, count(skipmissing(resultats_moyennes_sens[i][!, :X0_moyen] .<= 0)) / length(resultats_moyennes_sens[i][!, :X0_moyen]))

    push!(pourcentage_Y0_sens, count(skipmissing(resultats_moyennes_sens[i][!, :Y0_moyen] .<= 0)) / length(resultats_moyennes_sens[i][!, :Y0_moyen]))

end

println("Moyenne du pourcentage de 0 dans X0 antisens : ", round(mean(pourcentage_X0_antisens), digits=3))
println("Moyenne du pourcentage de 0 dans Y0 antisens : ", round(mean(pourcentage_Y0_antisens), digits=3))
println("Moyenne du pourcentage de 0 dans X0 sens : ", round(mean(pourcentage_X0_sens), digits=3))
println("Moyenne du pourcentage de 0 dans Y0 sens : ", round(mean(pourcentage_Y0_sens), digits=3))


#Vérification des covariance des X0 et Y0

cor_values = DataFrame(
    Stade=Int64[3, 4, 5, 6, 7, 8, 9, 10],
    Correlation=Float64[0, 0, 0, 0, 0, 0, 0, 0]
)


fig_cytofluogramme = Figure(resolution=(1200, 1600))

#Boucle for pour tracer les graphiques  

for i in (3:10)

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    moyenne_cytofluogramme_stade = resultats_moyennes_antisens[i]

    cor_values[n_injection, 2] = cor(moyenne_cytofluogramme_stade[!, :X0_moyen], moyenne_cytofluogramme_stade[!, :Y0_moyen])

    println("La corrélation pour le stade " * string(i) * " est de : " * string(round(cor_values[n_injection, 2], digits=3)))

    q2_CTAT_quantile = quantile(moyenne_cytofluogramme_stade.Y0_moyen, 0.5)
    q3_CTAT_quantile = quantile(moyenne_cytofluogramme_stade.Y0_moyen, 0.95)
    q2_CTAC_quantile = quantile(moyenne_cytofluogramme_stade.X0_moyen, 0.5)
    q3_CTAC_quantile = quantile(moyenne_cytofluogramme_stade.X0_moyen, 0.95)

    if n_injection <= 2
        ax = Axis(fig_cytofluogramme[1, n_injection],
            title="Correlation between X0 = CTAT and \n Y0 = CTAC antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_cytofluogramme[2, n_injection-2],
            title="Correlation between X0 = CTAT and \n Y0 = CTAC antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_cytofluogramme[3, n_injection-4],
            title="Correlation between X0 = CTAT and \n Y0 = CTAC antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_cytofluogramme[4, n_injection-6],
            title="Correlation between X0 = CTAT and \n Y0 = CTAC antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    end

    Makie.scatter!(ax,
        moyenne_cytofluogramme_stade[!, :X0_moyen],
        moyenne_cytofluogramme_stade[!, :Y0_moyen],
        color=:blue,
        markersize=5,
        alpha=0.6
    )

    #Ajout des R² et des PCC ainsi que là significativité des régressions
    #*** = p < 10^-99, extrêmement significatif
    #** = p < 10^-5, très significatif
    #* = p < 0.05, significatif
    #ns = p > 0.05, non significatif


    text!(ax,
        "PCC =" * string(round(cor(moyenne_cytofluogramme_stade[!, :X0_moyen], moyenne_cytofluogramme_stade[!, :Y0_moyen]), digits=3)),
        position=(maximum(moyenne_cytofluogramme_stade.X0_moyen) - 0.15 * maximum(moyenne_cytofluogramme_stade.X0_moyen), maximum(moyenne_cytofluogramme_stade.Y0_moyen) - 0.15 * maximum(moyenne_cytofluogramme_stade.Y0_moyen)),
        color=:black,
        align=(:right, :top),
        fontsize=12
    )


end

supertitle = Label(fig_cytofluogramme[0, :],
    "Correlation between X0 = CTAT and Y0 = CTAC antisens in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "Correlation_X0_Y0_antisens.png")

save(chemin, fig_cytofluogramme)


cor_values = DataFrame(
    Stade=Int64[3, 4, 5, 6, 7, 8, 9, 10],
    Correlation=Float64[0, 0, 0, 0, 0, 0, 0, 0]
)


fig_cytofluogramme = Figure(resolution=(1200, 1600))

#Boucle for pour tracer les graphiques  

for i in (3:10)

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    moyenne_cytofluogramme_stade = resultats_moyennes_sens[i]

    cor_values[n_injection, 2] = cor(moyenne_cytofluogramme_stade[!, :X0_moyen], moyenne_cytofluogramme_stade[!, :Y0_moyen])

    println("La corrélation pour le stade " * string(i) * " est de : " * string(round(cor_values[n_injection, 2], digits=3)))

    q2_CTAT_quantile = quantile(moyenne_cytofluogramme_stade.Y0_moyen, 0.5)
    q3_CTAT_quantile = quantile(moyenne_cytofluogramme_stade.Y0_moyen, 0.95)
    q2_CTAC_quantile = quantile(moyenne_cytofluogramme_stade.X0_moyen, 0.5)
    q3_CTAC_quantile = quantile(moyenne_cytofluogramme_stade.X0_moyen, 0.95)

    if n_injection <= 2
        ax = Axis(fig_cytofluogramme[1, n_injection],
            title="Correlation between X0 = CTAC and \n Y0 = CTAT antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_cytofluogramme[2, n_injection-2],
            title="Correlation between X0 = CTAC and \n Y0 = CTAT antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_cytofluogramme[3, n_injection-4],
            title="Correlation between X0 = CTAC and \n Y0 = CTAT antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_cytofluogramme[4, n_injection-6],
            title="Correlation between X0 = CTAC and \n Y0 = CTAT antisens for stage " * string(i),
            xlabel="X0",
            ylabel="Y0",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    end

    Makie.scatter!(ax,
        moyenne_cytofluogramme_stade[!, :X0_moyen],
        moyenne_cytofluogramme_stade[!, :Y0_moyen],
        color=:blue,
        markersize=5,
        alpha=0.6
    )


    #Ajout des R² et des PCC ainsi que là significativité des régressions
    #*** = p < 10^-99, extrêmement significatif
    #** = p < 10^-5, très significatif
    #* = p < 0.05, significatif
    #ns = p > 0.05, non significatif


    text!(ax,
        "PCC = " * string(round(cor(moyenne_cytofluogramme_stade[!, :X0_moyen], moyenne_cytofluogramme_stade[!, :Y0_moyen]), digits=3)),
        position=(maximum(moyenne_cytofluogramme_stade.X0_moyen) - 0.15 * maximum(moyenne_cytofluogramme_stade.X0_moyen), maximum(moyenne_cytofluogramme_stade.Y0_moyen) - 0.15 * maximum(moyenne_cytofluogramme_stade.Y0_moyen)),
        color=:black,
        align=(:right, :top),
        fontsize=12
    )


end

supertitle = Label(fig_cytofluogramme[0, :],
    "Correlation between X0 = CTAC and Y0 = CTAT sens in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

chemin = joinpath("/Users", "verme", "Desktop", "Correlation_X0_Y0_sens.png")

save(chemin, fig_cytofluogramme)


###############################################
#vérification des covariance entre X0 antisens


function verification_covariance_des_X0(dictionaire_des_moyennes, orientation)

    fig_correlation = Figure(size=(1200, 1200))

    for i in (3:10)

        for j in (3:10)

            if i != j && i < j

                if length(dictionaire_des_moyennes[i][!, :X0_moyen]) != length(dictionaire_des_moyennes[j][!, :X0_moyen])

                    longeur_minimal = min(length(dictionaire_des_moyennes[i][!, :X0_moyen]), length(dictionaire_des_moyennes[j][!, :X0_moyen]))

                    vecteur_temporaire_X0_i = dictionaire_des_moyennes[i][1:longeur_minimal, :X0_moyen]
                    vecteur_temporaire_X0_j = dictionaire_des_moyennes[j][1:longeur_minimal, :X0_moyen]

                else

                    vecteur_temporaire_X0_i = dictionaire_des_moyennes[i][!, :X0_moyen]
                    vecteur_temporaire_X0_j = dictionaire_des_moyennes[j][!, :X0_moyen]

                end

                ax = Axis(fig_correlation[i-2, j-2],
                    title="X0 stage " * string(i) * " vs " * string(j),
                    aspect=AxisAspect(1),
                    xticks=LinearTicks(5),
                    yticks=LinearTicks(5))

                hidedecorations!(ax)

                Makie.scatter!(ax,
                    vecteur_temporaire_X0_i,
                    vecteur_temporaire_X0_j,
                    color=:blue,
                    markersize=5,
                    alpha=0.6
                )

                text!(ax, string(round(cor(vecteur_temporaire_X0_i, vecteur_temporaire_X0_j), digits=3)),
                    position=(maximum(vecteur_temporaire_X0_i), maximum(vecteur_temporaire_X0_j)),
                    align=(:right, :top),
                    color=:black)

            elseif i != j && i > j
                ax_plus_petit = Axis(fig_correlation[i-2, j-2],
                    aspect=AxisAspect(1),
                    xticks=LinearTicks(5),
                    yticks=LinearTicks(5))
                hidedecorations!(ax_plus_petit)

                if length(dictionaire_des_moyennes[i][!, :X0_moyen]) != length(dictionaire_des_moyennes[j][!, :X0_moyen])

                    longeur_minimal = min(length(dictionaire_des_moyennes[i][!, :X0_moyen]), length(dictionaire_des_moyennes[j][!, :X0_moyen]))

                    vecteur_temporaire_X0_i = dictionaire_des_moyennes[i][1:longeur_minimal, :X0_moyen]
                    vecteur_temporaire_X0_j = dictionaire_des_moyennes[j][1:longeur_minimal, :X0_moyen]

                else

                    vecteur_temporaire_X0_i = dictionaire_des_moyennes[i][!, :X0_moyen]
                    vecteur_temporaire_X0_j = dictionaire_des_moyennes[j][!, :X0_moyen]

                end

                text!(ax_plus_petit, string(round(cor(vecteur_temporaire_X0_i, vecteur_temporaire_X0_j), digits=3)),
                    position=(0.5, 0.5),
                    align=(:center, :center),
                    color=:black)

            elseif i == j

                ax_identique = Axis(fig_correlation[i-2, j-2],
                    aspect=AxisAspect(1),
                    xticks=LinearTicks(5),
                    yticks=LinearTicks(5))

                hist!(ax_identique, dictionaire_des_moyennes[i][!, :X0_moyen],
                    bins=200,
                    color=:blue,
                    strokecolor=:black)

                hidedecorations!(ax_identique)

            end
        end
    end

    nom = string(Symbol(orientation))

    supertitle = Label(fig_correlation[0, :],
        "Correlation of X0 for the mean cytofluogramme in stages 3 to 10 for " * nom,
        fontsize=15, font="Arial", padding=(10, 10, 10, 10))

    nom_fichier = "correlation_X0_$(nom).png"

    chemin = joinpath(homedir(), "Desktop", nom_fichier)

    save(chemin, fig_correlation)

    return fig_correlation

end


verification_covariance_des_X0(resultats_moyennes_antisens, "antisens")

verification_covariance_des_X0(resultats_moyennes_sens, "sens")

#vérification des covariance entre Y0

function verification_covariance_des_Y0(dictionaire_des_moyennes, orientation)

    fig_correlation = Figure(size=(1200, 1200))

    for i in (3:10)

        for j in (3:10)

            if i != j && i < j

                if length(dictionaire_des_moyennes[i][!, :Y0_moyen]) != length(dictionaire_des_moyennes[j][!, :Y0_moyen])

                    longeur_minimal = min(length(dictionaire_des_moyennes[i][!, :Y0_moyen]), length(dictionaire_des_moyennes[j][!, :Y0_moyen]))

                    vecteur_temporaire_Y0_i = dictionaire_des_moyennes[i][1:longeur_minimal, :Y0_moyen]
                    vecteur_temporaire_Y0_j = dictionaire_des_moyennes[j][1:longeur_minimal, :Y0_moyen]

                else

                    vecteur_temporaire_Y0_i = dictionaire_des_moyennes[i][!, :Y0_moyen]
                    vecteur_temporaire_Y0_j = dictionaire_des_moyennes[j][!, :Y0_moyen]

                end

                ax = Axis(fig_correlation[i-2, j-2],
                    title="Y0 stage " * string(i) * " vs " * string(j),
                    aspect=AxisAspect(1),
                    xticks=LinearTicks(5),
                    yticks=LinearTicks(5))

                hidedecorations!(ax)

                Makie.scatter!(ax,
                    vecteur_temporaire_Y0_i,
                    vecteur_temporaire_Y0_j,
                    color=:blue,
                    markersize=5,
                    alpha=0.6
                )

                text!(ax, string(round(cor(vecteur_temporaire_Y0_i, vecteur_temporaire_Y0_j), digits=3)),
                    position=(maximum(vecteur_temporaire_Y0_i), maximum(vecteur_temporaire_Y0_j)),
                    align=(:right, :top),
                    color=:black)

            elseif i != j && i > j
                ax_plus_petit = Axis(fig_correlation[i-2, j-2],
                    aspect=AxisAspect(1),
                    xticks=LinearTicks(5),
                    yticks=LinearTicks(5))
                hidedecorations!(ax_plus_petit)

                if length(dictionaire_des_moyennes[i][!, :Y0_moyen]) != length(dictionaire_des_moyennes[j][!, :Y0_moyen])

                    longeur_minimal = min(length(dictionaire_des_moyennes[i][!, :Y0_moyen]), length(dictionaire_des_moyennes[j][!, :Y0_moyen]))

                    vecteur_temporaire_Y0_i = dictionaire_des_moyennes[i][1:longeur_minimal, :Y0_moyen]
                    vecteur_temporaire_Y0_j = dictionaire_des_moyennes[j][1:longeur_minimal, :Y0_moyen]

                else

                    vecteur_temporaire_Y0_i = dictionaire_des_moyennes[i][!, :Y0_moyen]
                    vecteur_temporaire_Y0_j = dictionaire_des_moyennes[j][!, :Y0_moyen]

                end

                text!(ax_plus_petit, string(round(cor(vecteur_temporaire_Y0_i, vecteur_temporaire_Y0_j), digits=3)),
                    position=(0.5, 0.5),
                    align=(:center, :center),
                    color=:black)

            elseif i == j

                ax_identique = Axis(fig_correlation[i-2, j-2],
                    aspect=AxisAspect(1),
                    xticks=LinearTicks(5),
                    yticks=LinearTicks(5))

                hist!(ax_identique, dictionaire_des_moyennes[i][!, :Y0_moyen],
                    bins=200,
                    color=:blue,
                    strokecolor=:black)

                hidedecorations!(ax_identique)

            end
        end
    end

    nom = string(Symbol(orientation))

    supertitle = Label(fig_correlation[0, :],
        "Correlation of Y0 for the mean cytofluogramme in stages 3 to 10 for " * nom,
        fontsize=15, font="Arial", padding=(10, 10, 10, 10))

    nom_fichier = "correlation_Y0_$(nom).png"

    chemin = joinpath(homedir(), "Desktop", nom_fichier)

    save(chemin, fig_correlation)

    return fig_correlation

end

verification_covariance_des_Y0(resultats_moyennes_antisens, "antisens")
verification_covariance_des_Y0(resultats_moyennes_sens, "sens")

#####################################

