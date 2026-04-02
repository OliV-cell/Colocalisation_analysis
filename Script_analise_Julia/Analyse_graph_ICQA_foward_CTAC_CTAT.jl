##############################################
#Script exploratoire de visualisation des graphique des test de 
#colocalisation de ICQA entre les sondes AAACTAC et AAACTAT
#sur les images Dvir48 en sens forward strand
##############################################
#importation des packages

using DataFrames
using CSV
using Plots
using Makie
using CairoMakie
using GLMakie
using Statistics
using GLM
using HypothesisTests

##############################################
#numero des image à analyser pour chaque stade

numero_image_stade3 = (21, 34, 35, 36, 37, 52, 53, 55, 61, 75, 98)
numero_image_stade4 = (13, 39, 44, 46, 48, 55)
numero_image_stade5 = (12, 17, 22, 23, 29, 30, 37, 45, 51, 54, 55, 68, 69, 71, 74, 78, 80, 82, 83)
numero_image_stade6 = (1, 3, 11, 22, 34, 49, 50, 63, 67, 70, 72, 81, 87, 93, 100, 101, 105)
numero_image_stade7 = (5, 6, 7, 9, 11, 14, 16, 18, 20, 24, 28, 31, 38, 40, 54, 57, 58, 73, 76, 77, 84, 90, 91, 94)
numero_image_stade8 = (2, 4, 10, 14, 19, 32, 33, 64, 65, 85, 104, 107, 113, 116)
numero_image_stade9 = (8, 41, 59, 88, 92, 99, 103, 109, 111)
numero_image_stade10 = (25, 26, 27, 102, 110, 114, 118)

##############################################
#numero des image à analyser pour chaque stade
#pour la réplication 1

numero_image_stade3_rep1 = (4, 16, 20, 21, 27, 35, 39, 59, 66, 85, 91, 93)
numero_image_stade4_rep1 = (1, 90, 92)
numero_image_stade5_rep1 = (4, 10, 13, 13, 14, 16, 20, 21, 26, 30, 32, 33, 42, 48, 58, 66, 68, 69, 70, 72, 78, 83, 87, 90, 92, 94, 95)
numero_image_stade6_rep1 = (5, 8, 9, 22, 34, 37, 41, 46, 60, 67, 79, 85)
numero_image_stade7_rep1 = (2, 18, 19, 24, 25, 31, 38, 40, 44, 45, 50, 55, 56, 74, 75, 76, 83, 88, 89)
numero_image_stade8_rep1 = (3, 6, 7, 12, 17, 20, 28, 29, 36, 43, 45, 48, 49, 51, 52, 57, 64, 70, 77, 80, 81, 82, 84, 86)
numero_image_stade9_rep1 = (23, 53, 54, 61, 73)
numero_image_stade10_rep1 = (11, 15, 47, 62, 63, 65)

##############################################
#création des dictionnaires pour stocker les données des graphiques

dictonaire_des_graphique_ICQA_stade3 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade4 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade5 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade6 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade7 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade8 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade9 = Dict{String,DataFrame}()
dictonaire_des_graphique_ICQA_stade10 = Dict{String,DataFrame}()

###############################################
#boucle de lecture des données des graphiques pour les stdae 
#et concaténation des répetion dans les dictinnaire assocers

#pour le stade 3

for i in numero_image_stade3
    dictonaire_des_graphique_ICQA_stade3["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_sade3_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade3_rep1
    dictonaire_des_graphique_ICQA_stade3["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade3_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 4
for i in numero_image_stade4
    dictonaire_des_graphique_ICQA_stade4["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAC_sensboth_stade4_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade4_rep1
    dictonaire_des_graphique_ICQA_stade4["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade4_rep1_Image_" * string(i) * ".csv", DataFrame)
end
##############################################
#pour le stade 5
for i in numero_image_stade5
    dictonaire_des_graphique_ICQA_stade5["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade5_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade5_rep1
    dictonaire_des_graphique_ICQA_stade5["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade5_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 6
for i in numero_image_stade6
    dictonaire_des_graphique_ICQA_stade6["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade6_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade6_rep1
    dictonaire_des_graphique_ICQA_stade6["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade6_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 7
for i in numero_image_stade7
    dictonaire_des_graphique_ICQA_stade7["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade7_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade7_rep1
    dictonaire_des_graphique_ICQA_stade7["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade7_rep1_Image_" * string(i) * ".csv", DataFrame)
end
################################################
#pour le stade 8
for i in numero_image_stade8
    dictonaire_des_graphique_ICQA_stade8["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade8_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade8_rep1
    dictonaire_des_graphique_ICQA_stade8["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade8_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 9
for i in numero_image_stade9
    dictonaire_des_graphique_ICQA_stade9["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade9_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade9_rep1
    dictonaire_des_graphique_ICQA_stade9["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade9_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 10
for i in numero_image_stade10
    dictonaire_des_graphique_ICQA_stade10["ICQA_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade10_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade10_rep1
    dictonaire_des_graphique_ICQA_stade10["ICQA_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_ICA_ADvir48_CTAC_CTAT_sensboth_stade10_rep1_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#Premiers jets de cyrofluogramme
#changer les truc manuellement pour chaque stade

fig_ICQA_stade4 = Figure()

ax = Axis(fig_ICQA_stade4[1, 1], title="ICQA des images Dvir48_CTAC_CTAT_sensboth_rep1_stade4")

for i in numero_image_stade4
    Makie.scatter!(ax,
        dictonaire_des_graphique_ICQA_stade4["ICQA_image_"*string(i)][!, :X],
        dictonaire_des_graphique_ICQA_stade4["ICQA_image_"*string(i)][!, :Y],
        color=:black,
        markersize=5,
        alpha=0.2)
end

for i in numero_image_stade4_rep1
    Makie.scatter!(ax,
        dictonaire_des_graphique_ICQA_stade4["ICQA_rep1_image_"*string(i)][!, :X],
        dictonaire_des_graphique_ICQA_stade4["ICQA_rep1_image_"*string(i)][!, :Y],
        color=:black,
        markersize=5,
        alpha=0.2)
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

        maximum_length[n_injection] = length(dictionaire["ICQA_image_"*string(i)][!, :X])

    end

    for j in liste_numero_image_rep1

        n_injection_rep1 = findfirst(x -> x == j, liste_numero_image_rep1)

        maximum_length_rep1[n_injection_rep1] = length(dictionaire["ICQA_rep1_image_"*string(j)][!, :X])
    end

    maximum_length = maximum(skipmissing(vcat(maximum_length, maximum_length_rep1)))

    return maximum_length

end


max_length_stade3 = mesure_max(dictonaire_des_graphique_ICQA_stade3, numero_image_stade3, numero_image_stade3_rep1)
max_length_stade4 = mesure_max(dictonaire_des_graphique_ICQA_stade4, numero_image_stade4, numero_image_stade4_rep1)
max_length_stade5 = mesure_max(dictonaire_des_graphique_ICQA_stade5, numero_image_stade5, numero_image_stade5_rep1)
max_length_stade6 = mesure_max(dictonaire_des_graphique_ICQA_stade6, numero_image_stade6, numero_image_stade6_rep1)
max_length_stade7 = mesure_max(dictonaire_des_graphique_ICQA_stade7, numero_image_stade7, numero_image_stade7_rep1)
max_length_stade8 = mesure_max(dictonaire_des_graphique_ICQA_stade8, numero_image_stade8, numero_image_stade8_rep1)
max_length_stade9 = mesure_max(dictonaire_des_graphique_ICQA_stade9, numero_image_stade9, numero_image_stade9_rep1)
max_length_stade10 = mesure_max(dictonaire_des_graphique_ICQA_stade10, numero_image_stade10, numero_image_stade10_rep1)

##############################################
#Construiction des dataframes pour le calcul des moyennes

X_data_stade3 = DataFrame()
Y_data_stade3 = DataFrame()

X_data_stade4 = DataFrame()
Y_data_stade4 = DataFrame()

X_data_stade5 = DataFrame()
Y_data_stade5 = DataFrame()

X_data_stade6 = DataFrame()
Y_data_stade6 = DataFrame()

X_data_stade7 = DataFrame()
Y_data_stade7 = DataFrame()

X_data_stade8 = DataFrame()
Y_data_stade8 = DataFrame()

X_data_stade9 = DataFrame()
Y_data_stade9 = DataFrame()

X_data_stade10 = DataFrame()
Y_data_stade10 = DataFrame()
##############################################
#Fonction pour ajuster la taille des données des cytofluogrammes
#avant de faire la moyenne

function ajustement_taille(dictionaire, liste_numero_image, liste_numero_image_rep1, max_length, X_data, Y_data)

    for i in liste_numero_image

        if length(dictionaire["ICQA_image_"*string(i)][!, :X]) < max_length

            allowmissing!(dictionaire["ICQA_image_"*string(i)], [:X, :Y])
            difference_length = max_length - length(dictionaire["ICQA_image_"*string(i)][!, :X])

            for j in 1:difference_length

                push!(dictionaire["ICQA_image_"*string(i)][!, :X], missing)
                push!(dictionaire["ICQA_image_"*string(i)][!, :Y], missing)
            end

        end
        n_injection = Symbol(findfirst(x -> x == i, liste_numero_image))

        X_data[!, n_injection] = dictionaire["ICQA_image_"*string(i)][!, :X]
        Y_data[!, n_injection] = dictionaire["ICQA_image_"*string(i)][!, :Y]
    end

    taille = length(X_data[1, :])

    for i in liste_numero_image_rep1

        if length(dictionaire["ICQA_rep1_image_"*string(i)][!, :X]) < max_length

            allowmissing!(dictionaire["ICQA_rep1_image_"*string(i)], [:X, :Y])
            difference_length = max_length - length(dictionaire["ICQA_rep1_image_"*string(i)][!, :X])

            for j in 1:difference_length

                push!(dictionaire["ICQA_rep1_image_"*string(i)][!, :X], missing)
                push!(dictionaire["ICQA_rep1_image_"*string(i)][!, :Y], missing)
            end

        end

        n_injection_rep1 = Int64(findfirst(x -> x == i, liste_numero_image_rep1))

        a_partir_de = Symbol((Int64(taille) + n_injection_rep1))

        X_data[!, a_partir_de] = dictionaire["ICQA_rep1_image_"*string(i)][!, :X]
        Y_data[!, a_partir_de] = dictionaire["ICQA_rep1_image_"*string(i)][!, :Y]
    end

end

#ajustement des taille et remplissage des dataframes

ajustement_taille(dictonaire_des_graphique_ICQA_stade3, numero_image_stade3, numero_image_stade3_rep1, max_length_stade3, X_data_stade3, Y_data_stade3)
ajustement_taille(dictonaire_des_graphique_ICQA_stade4, numero_image_stade4, numero_image_stade4_rep1, max_length_stade4, X_data_stade4, Y_data_stade4)
ajustement_taille(dictonaire_des_graphique_ICQA_stade5, numero_image_stade5, numero_image_stade5_rep1, max_length_stade5, X_data_stade5, Y_data_stade5)
ajustement_taille(dictonaire_des_graphique_ICQA_stade6, numero_image_stade6, numero_image_stade6_rep1, max_length_stade6, X_data_stade6, Y_data_stade6)
ajustement_taille(dictonaire_des_graphique_ICQA_stade7, numero_image_stade7, numero_image_stade7_rep1, max_length_stade7, X_data_stade7, Y_data_stade7)
ajustement_taille(dictonaire_des_graphique_ICQA_stade8, numero_image_stade8, numero_image_stade8_rep1, max_length_stade8, X_data_stade8, Y_data_stade8)
ajustement_taille(dictonaire_des_graphique_ICQA_stade9, numero_image_stade9, numero_image_stade9_rep1, max_length_stade9, X_data_stade9, Y_data_stade9)
ajustement_taille(dictonaire_des_graphique_ICQA_stade10, numero_image_stade10, numero_image_stade10_rep1, max_length_stade10, X_data_stade10, Y_data_stade10)

##############################################
#Calcul des moyennes des ICQA pour chaque stade
#initialisation des dataframes pour stocker les moyennes

resultats_moyennes = Dict{Int64,DataFrame}()

for i in (3:10)

    X_data = eval(Symbol("X_data_stade$i"))
    Y_data = eval(Symbol("Y_data_stade$i"))

    x_moyens = [
        begin
            m = skipmissing(row)
            isempty(m) ? NaN : mean(m)
        end for row in eachrow(X_data)
    ]

    y_moyens = [
        begin
            m = skipmissing(row)
            isempty(m) ? NaN : mean(m)
        end for row in eachrow(Y_data)
    ]

    # 3. On crée le DataFrame final pour ce stade et on le stocke dans le dictionnaire
    DataFrame_temporaire = DataFrame(X_moyen=x_moyens, Y_moyen=y_moyens)

    filter!(row -> !isnan(row.X_moyen) && !isnan(row.Y_moyen), DataFrame_temporaire)

    resultats_moyennes[i] = DataFrame_temporaire
end


##############################################
#Visualisation des cytofluogrammes des moyennes

cor_values = DataFrame(
    Stade=Int64[3, 4, 5, 6, 7, 8, 9, 10],
    Correlation=Float64[0, 0, 0, 0, 0, 0, 0, 0]
)

fig_ICQA = Figure(resolution=(1200, 1600))

for i in (3:10)

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    moyenne_ICQA_stade = resultats_moyennes[i]

    cor_values[n_injection, 2] = cor(moyenne_ICQA_stade[!, :X_moyen], moyenne_ICQA_stade[!, :Y_moyen])

    println("La corrélation pour le stade " * string(i) * " est de : " * string(round(cor_values[n_injection, 2], digits=3)))

    x_min, x_max = extrema(moyenne_ICQA_stade.X_moyen)

    if n_injection <= 2
        ax = Axis(fig_ICQA[1, n_injection],
            title="Mean ICA A for stage " * string(i),
            xlabel="Mean (Ai − a)(Bi − b)",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_ICQA[2, n_injection-2],
            title="Mean ICA A for stage " * string(i),
            xlabel="Mean (Ai − a)(Bi − b)",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_ICQA[3, n_injection-4],
            title="Mean ICA A for stage " * string(i),
            xlabel="Mean (Ai − a)(Bi − b)",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_ICQA[4, n_injection-6],
            title="Mean ICA A for stage " * string(i),
            xlabel="Mean (Ai − a)(Bi − b)",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5))

    end

    Makie.scatter!(ax,
        moyenne_ICQA_stade[!, :X_moyen],
        moyenne_ICQA_stade[!, :Y_moyen],
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
        "PCC =" * string(round(cor(moyenne_ICQA_stade[!, :X_moyen], moyenne_ICQA_stade[!, :Y_moyen]), digits=3)),
        position=(x_max, maximum(moyenne_ICQA_stade.Y_moyen) - 0.15 * maximum(moyenne_ICQA_stade.Y_moyen)),
        color=:black,
        align=(:right, :top),
        fontsize=12
    )
    hlines!(ax, 0,
        color=:black,
        linewidth=1,
        linestyle=:dot)


    vlines!(ax, 0,
        color=:black,
        linewidth=1,
        linestyle=:dot)

end

supertitle = Label(fig_ICQA[0, :],
    "Mean ICA A \n for the expression of AAACTAC and AAACTAT reverse LncRNA strand in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

##############################################

# 1. Définir le chemin proprement (joinpath assemble uniquement des bouts de texte)
chemin = joinpath("/Users", "verme", "Desktop", "ICQA_moyenne_reverse_lncRNAstrand_CTAC_CTAT.png")

# 2. Utiliser la fonction save de Makie
save(chemin, fig_ICQA)

##############################################
