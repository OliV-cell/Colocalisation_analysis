##############################################
#Script exploratoire de visualisation des graphique des test de 
#colocalisation de cytofluogramme entre les sondes AAACTAC et AAACTAT
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
using Printf
using Distributions

##############################################
#numero des image à analyser pour chaque stade

numero_image_stade3_rep1 = (6, 21, 22, 23, 29, 32, 34, 36, 39, 47, 63, 66, 70, 72, 73, 78, 87, 90)
numero_image_stade4_rep1 = (9, 12, 17, 19, 20, 28, 39, 53)
numero_image_stade5_rep1 = (1, 4, 7, 9, 15, 16, 22, 23, 25, 26, 31, 32, 34, 37, 40, 42, 44, 45, 47, 57)
numero_image_stade6_rep1 = (8, 18, 27, 46, 59, 75, 96, 97, 108, 109)
numero_image_stade7_rep1 = (2, 3, 14, 24, 27, 30, 33, 35, 41, 43, 44, 50, 58, 60, 62, 68, 110)
numero_image_stade8_rep1 = (10, 49, 52, 55, 84, 89, 111)
numero_image_stade9_rep1 = (51, 61, 105)
numero_image_stade10_rep1 = (56, 67, 82, 83)

##############################################
#numero des image à analyser pour chaque stade
#pour la réplication 1

numero_image_stade3_rep2 = (5, 11, 19, 36, 43, 103)
numero_image_stade4_rep2 = (13, 77)
numero_image_stade5_rep2 = (4, 7, 10, 13, 19, 20, 51, 56, 59, 65, 69, 75, 83, 89, 91, 99, 102, 104)
numero_image_stade6_rep2 = (1, 9, 33, 35, 41, 51, 56, 64, 88, 91, 92, 94, 96)
numero_image_stade7_rep2 = (15, 18, 26, 29, 30, 34, 42, 46, 57, 68, 70, 74, 87)
numero_image_stade8_rep2 = (2, 26, 27, 31, 40, 47, 50, 52, 53, 54, 60, 67, 73, 86, 87, 89, 100)
numero_image_stade9_rep2 = (3, 6, 12, 22, 25, 32, 37, 39, 44, 45, 61, 63, 66, 72, 85, 90, 95)
numero_image_stade10_rep2 = (14, 16, 17, 38, 48, 49, 62, 81, 84)
##############################################
#création des dictionnaires pour stocker les données des graphiques

dictonaire_des_graphique_cytofluogramme_stade3 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade4 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade5 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade6 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade7 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade8 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade9 = Dict{String,DataFrame}()
dictonaire_des_graphique_cytofluogramme_stade10 = Dict{String,DataFrame}()

###############################################
#boucle de lecture des données des graphiques pour les stdae 
#et concaténation des répetion dans les dictinnaire assocers

#pour le stade 3

for i in numero_image_stade3_rep1
    dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_rep1_stade3_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade3_rep2
    dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade3_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 4
for i in numero_image_stade4_rep1
    dictonaire_des_graphique_cytofluogramme_stade4["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antosensboth_stade4_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade4_rep2
    dictonaire_des_graphique_cytofluogramme_stade4["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade4_rep2_Image_" * string(i) * ".csv", DataFrame)
end
##############################################
#pour le stade 5
for i in numero_image_stade5_rep1
    dictonaire_des_graphique_cytofluogramme_stade5["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antosensboth_stade5_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade5_rep2
    dictonaire_des_graphique_cytofluogramme_stade5["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade5_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 6
for i in numero_image_stade6_rep1
    dictonaire_des_graphique_cytofluogramme_stade6["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade6_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade6_rep2
    dictonaire_des_graphique_cytofluogramme_stade6["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade6_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 7

for i in numero_image_stade7_rep1
    dictonaire_des_graphique_cytofluogramme_stade7["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade7_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade7_rep2
    dictonaire_des_graphique_cytofluogramme_stade7["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade7_rep2_Image_" * string(i) * ".csv", DataFrame)
end
################################################
#pour le stade 8
for i in numero_image_stade8_rep1
    dictonaire_des_graphique_cytofluogramme_stade8["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade8_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade8_rep2
    dictonaire_des_graphique_cytofluogramme_stade8["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade8_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 9
for i in numero_image_stade9_rep1
    dictonaire_des_graphique_cytofluogramme_stade9["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade9_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade9_rep2
    dictonaire_des_graphique_cytofluogramme_stade9["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade9_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#pour le stade 10
for i in numero_image_stade10_rep1
    dictonaire_des_graphique_cytofluogramme_stade10["cytofluogramme_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade10_rep1_Image_" * string(i) * ".csv", DataFrame)
end

for i in numero_image_stade10_rep2
    dictonaire_des_graphique_cytofluogramme_stade10["cytofluogramme_rep1_image_"*string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_antisens_Dvir_CTAC_CTAT/Entrepo_data_graph/Resultat_CytofluogramDvir48_CTAC_CTAT_antisensboth_stade10_rep2_Image_" * string(i) * ".csv", DataFrame)
end
###############################################
#Premiers jets de cyrofluogramme
#changer les truc manuellement pour chaque stade

fig_cytofluogramme_stade4 = Figure()

ax = Axis(fig_cytofluogramme_stade4[1, 1], title="Cytofluogrammes des images Dvir48_CTAC_CTAT_sensboth_rep1_stade4")

for i in numero_image_stade3_rep1
    Makie.scatter!(ax,
        dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_image_"*string(i)][!, :X0],
        dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_image_"*string(i)][!, :Y0],
        color=:black,
        markersize=5,
        alpha=1)


end

for i in numero_image_stade3_rep2
    Makie.scatter!(ax,
        dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_rep1_image_"*string(i)][!, :X0],
        dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_rep1_image_"*string(i)][!, :Y0],
        color=:black,
        markersize=5,
        alpha=1)


end


println(count(skipmissing(dictonaire_des_graphique_cytofluogramme_stade3["cytofluogramme_image_21"][!, :Y0] .<= 0)))

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


max_length_stade3 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade3, numero_image_stade3_rep1, numero_image_stade3_rep2)
max_length_stade4 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade4, numero_image_stade4_rep1, numero_image_stade4_rep2)
max_length_stade5 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade5, numero_image_stade5_rep1, numero_image_stade5_rep2)
max_length_stade6 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade6, numero_image_stade6_rep1, numero_image_stade6_rep2)
max_length_stade7 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade7, numero_image_stade7_rep1, numero_image_stade7_rep2)
max_length_stade8 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade8, numero_image_stade8_rep1, numero_image_stade8_rep2)
max_length_stade9 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade9, numero_image_stade9_rep1, numero_image_stade9_rep2)
max_length_stade10 = mesure_max(dictonaire_des_graphique_cytofluogramme_stade10, numero_image_stade10_rep1, numero_image_stade10_rep2)

##############################################
#Construiction des dataframes pour le calcul des moyennes

X0_data_stade3 = DataFrame()
Y0_data_stade3 = DataFrame()

X0_data_stade4 = DataFrame()
Y0_data_stade4 = DataFrame()

X0_data_stade5 = DataFrame()
Y0_data_stade5 = DataFrame()

X0_data_stade6 = DataFrame()
Y0_data_stade6 = DataFrame()

X0_data_stade7 = DataFrame()
Y0_data_stade7 = DataFrame()

X0_data_stade8 = DataFrame()
Y0_data_stade8 = DataFrame()

X0_data_stade9 = DataFrame()
Y0_data_stade9 = DataFrame()

X0_data_stade10 = DataFrame()
Y0_data_stade10 = DataFrame()

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

ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade3, numero_image_stade3_rep1, numero_image_stade3_rep2, max_length_stade3, X0_data_stade3, Y0_data_stade3)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade4, numero_image_stade4_rep1, numero_image_stade4_rep2, max_length_stade4, X0_data_stade4, Y0_data_stade4)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade5, numero_image_stade5_rep1, numero_image_stade5_rep2, max_length_stade5, X0_data_stade5, Y0_data_stade5)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade6, numero_image_stade6_rep1, numero_image_stade6_rep2, max_length_stade6, X0_data_stade6, Y0_data_stade6)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade7, numero_image_stade7_rep1, numero_image_stade7_rep2, max_length_stade7, X0_data_stade7, Y0_data_stade7)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade8, numero_image_stade8_rep1, numero_image_stade8_rep2, max_length_stade8, X0_data_stade8, Y0_data_stade8)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade9, numero_image_stade9_rep1, numero_image_stade9_rep2, max_length_stade9, X0_data_stade9, Y0_data_stade9)
ajustement_taille(dictonaire_des_graphique_cytofluogramme_stade10, numero_image_stade10_rep1, numero_image_stade10_rep2, max_length_stade10, X0_data_stade10, Y0_data_stade10)

##############################################
#Calcul des moyennes des cytofluogrammes
#initialisation des dataframes pour stocker les moyennes

resultats_moyennes = Dict{Int64,DataFrame}()

for i in (3:10)

    X0_data = eval(Symbol("X0_data_stade$i"))
    Y0_data = eval(Symbol("Y0_data_stade$i"))

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

    resultats_moyennes[i] = DataFrame_temporaire
end

##############################################
#Visualisation des cytofluogrammes des moyennes

cor_values = DataFrame(
    Stade=Int64[3, 4, 5, 6, 7, 8, 9, 10],
    Correlation=Float64[0, 0, 0, 0, 0, 0, 0, 0]
)

fig_cytofluogramme = Figure(resolution=(1200, 1600))

for i in (3:10)

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    moyenne_cytofluogramme_stade = resultats_moyennes[i]

    cor_values[n_injection, 2] = cor(moyenne_cytofluogramme_stade[!, :X0_moyen], moyenne_cytofluogramme_stade[!, :Y0_moyen])

    println("La corrélation pour le stade " * string(i) * " est de : " * string(round(cor_values[n_injection, 2], digits=3)))

    x_min, x_max = extrema(moyenne_cytofluogramme_stade.X0_moyen)

    q2_CTAT_quantile = quantile(moyenne_cytofluogramme_stade.Y0_moyen, 0.5)
    q3_CTAT_quantile = quantile(moyenne_cytofluogramme_stade.Y0_moyen, 0.95)
    q2_CTAC_quantile = quantile(moyenne_cytofluogramme_stade.X0_moyen, 0.5)
    q3_CTAC_quantile = quantile(moyenne_cytofluogramme_stade.X0_moyen, 0.95)

    if n_injection <= 2
        ax = Axis(fig_cytofluogramme[1, n_injection],
            title="Mean cytofluogramme for stage " * string(i),
            xlabel="Mean CTAT intensity",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_cytofluogramme[2, n_injection-2],
            title="Mean cytofluogramme for stage " * string(i),
            xlabel="Mean CTAT intensity",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_cytofluogramme[3, n_injection-4],
            title="Mean cytofluogramme for stage " * string(i),
            xlabel="Mean CTAT intensity",
            ylabel="Mean CTAC intensity",
            aspect=AxisAspect(2),
            xticks=LinearTicks(5),
            yticks=LinearTicks(5),
            xtickformat="{:.2e}",
            ytickformat="{:.2e}",
            xticklabelrotation=pi / 8)

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_cytofluogramme[4, n_injection-6],
            title="Mean cytofluogramme for stage " * string(i),
            xlabel="Mean CTAT intensity",
            ylabel="Mean CTAC intensity",
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

    text!(ax,
        "PCC = " * string(round(cor(moyenne_cytofluogramme_stade[!, :X0_moyen], moyenne_cytofluogramme_stade[!, :Y0_moyen]), digits=3)),
        position=(x_max, maximum(moyenne_cytofluogramme_stade.Y0_moyen) - 0.15 * maximum(moyenne_cytofluogramme_stade.Y0_moyen)),
        color=:black,
        align=(:right, :top),
        fontsize=12
    )

    #Ajout des lignes pour les quantiles
    #q2 = 50% = median
    #q3 = 95%
    #permet de visualiser la distribution des données

    hlines!(ax, q2_CTAT_quantile,
        color=:black,
        linewidth=1,
        linestyle=:dot)

    hlines!(ax, q3_CTAT_quantile,
        color=:black,
        linewidth=2,
        linestyle=:dot)

    vlines!(ax, q2_CTAC_quantile,
        color=:black,
        linewidth=1,
        linestyle=:dot)

    vlines!(ax, q3_CTAC_quantile,
        color=:black,
        linewidth=2,
        linestyle=:dot)

end

supertitle = Label(fig_cytofluogramme[0, :],
    "Mean Cytofluogrammes \n for the expression of AAACTAC and AAACTAT foward LncRNA strand in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

##############################################

# 1. Définir le chemin proprement (joinpath assemble uniquement des bouts de texte)
chemin = joinpath("/Users", "verme", "Desktop", "Cytofluogramme_moyenne_CTAC_CTAT_foward_lncRNA_strand.png")

# 2. Utiliser la fonction save de Makie
save(chemin, fig_cytofluogramme)

##############################################
