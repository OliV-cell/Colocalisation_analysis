##################################
#
#Script d'analyse des la oclocalisation de
#CTAC et CTAT sonde antisens en 3D
#
#####################################
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

#liste des numéro d'image

stade_3 = (1, 4, 6, 10, 15, 22, 30, 33, 35, 48, 59)
stade_4 = (2, 3, 17, 28, 44, 52, 58)
stade_5 = (1, 3, 4, 7, 9, 11, 13, 15, 26, 29, 32, 35, 46, 47, 57)
stade_6 = (2, 14, 19, 25, 27, 38, 50, 53)
stade_7 = (5, 8, 12, 20, 31, 34, 45)
stade_8 = (21, 23, 24, 37, 39, 42, 49)
stade_9 = (36, 40, 41, 51)
stade_10 = (54, 55, 56)

#initialsation des dictionnaire 
data_VanSteensel_3D_stade = Dict{Int,Dict{String,DataFrame}}()
data_Cytofluogram_3D_stade = Dict{Int,Dict{String,DataFrame}}()
data_ICQ_A_3D_stade = Dict{Int,Dict{String,DataFrame}}()
data_ICQ_B_3D_stade = Dict{Int,Dict{String,DataFrame}}()
data_DIANA_3D_stade = Dict{Int,Dict{String,DataFrame}}()


for i in (3:10)
    data_VanSteensel_3D_stade[i] = Dict{String,DataFrame}()
    data_Cytofluogram_3D_stade[i] = Dict{String,DataFrame}()
    data_ICQ_A_3D_stade[i] = Dict{String,DataFrame}()
    data_ICQ_B_3D_stade[i] = Dict{String,DataFrame}()
    data_DIANA_3D_stade[i] = Dict{String,DataFrame}()
end

#importation des données pour les différents stades

direction_data_graphique_3D = "/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_3D"

function importation(nom_graphique)

    nom = eval(Symbol("data_" * nom_graphique * "_3D_stade"))

    for i in stade_3
        nom[3][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade3.csv", DataFrame)
    end

    println("Importation des données pour le stade 3 terminée")

    for i in stade_4
        nom[4][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade4.csv", DataFrame)
    end

    println("Importation des données pour le stade 4 terminée")

    for i in stade_5
        nom[5][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade5.csv", DataFrame)
    end

    println("Importation des données pour le stade 5 terminée")

    for i in stade_6
        nom[6][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade6.csv", DataFrame)
    end

    println("Importation des données pour le stade 6 terminée")

    for i in stade_7
        nom[7][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade7.csv", DataFrame)
    end

    println("Importation des données pour le stade 7 terminée")

    for i in stade_8
        nom[8][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvri48_CTAC_CTAT_antisens_stade8.csv", DataFrame)
    end

    println("Importation des données pour le stade 8 terminée")

    for i in stade_9
        nom[9][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade9.csv", DataFrame)
    end

    println("Importation des données pour le stade 9 terminée")

    for i in stade_10
        nom[10][string(i)] = CSV.read(direction_data_graphique_3D * "/Resultat_" * nom_graphique * "_Series_" * string(i) * "Dvir48_CTAC_CTAT_antisens_stade10.csv", DataFrame)
    end

    println("Importation des données pour le stade 10 terminée")

end

importation("Cytofluogram")

################################################
#Fonction pour mesurer la longueur maximale des
#données des cytofluogrammes
#pour uniformiser la longueur des données
#avant de faire la moyenne

function mesure_max(dictionaire, liste_numero_image)

    maximum_length = Array{Union{Missing,Float64}}(missing, length(dictionaire))

    for i in liste_numero_image

        n_injection = findfirst(x -> x == i, liste_numero_image)

        maximum_length[n_injection] = length(dictionaire[string(i)][!, :X0])

    end

    maximum_length = maximum(skipmissing(maximum_length))

    return maximum_length

end

max_length_stade3 = mesure_max(data_Cytofluogram_3D_stade[3], stade_3)
max_length_stade4 = mesure_max(data_Cytofluogram_3D_stade[4], stade_4)
max_length_stade5 = mesure_max(data_Cytofluogram_3D_stade[5], stade_5)
max_length_stade6 = mesure_max(data_Cytofluogram_3D_stade[6], stade_6)
max_length_stade7 = mesure_max(data_Cytofluogram_3D_stade[7], stade_7)
max_length_stade8 = mesure_max(data_Cytofluogram_3D_stade[8], stade_8)
max_length_stade9 = mesure_max(data_Cytofluogram_3D_stade[9], stade_9)
max_length_stade10 = mesure_max(data_Cytofluogram_3D_stade[10], stade_10)

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

function ajustement_taille(dictionaire, liste_numero_image, max_length, X0_data, Y0_data)

    for i in liste_numero_image

        if length(dictionaire[string(i)][!, :X0]) < max_length

            allowmissing!(dictionaire[string(i)], [:X0, :Y0])
            difference_length = max_length - length(dictionaire[string(i)][!, :X0])

            for j in 1:difference_length

                push!(dictionaire[string(i)][!, :X0], missing)
                push!(dictionaire[string(i)][!, :Y0], missing)
            end

        end
        n_injection = Symbol(findfirst(x -> x == i, liste_numero_image))

        X0_data[!, n_injection] = dictionaire[string(i)][!, :X0]
        Y0_data[!, n_injection] = dictionaire[string(i)][!, :Y0]
    end

end

#ajustement des taille et remplissage des dataframes

ajustement_taille(data_Cytofluogram_3D_stade[3], stade_3, max_length_stade3, X0_data_stade3, Y0_data_stade3)
ajustement_taille(data_Cytofluogram_3D_stade[4], stade_4, max_length_stade4, X0_data_stade4, Y0_data_stade4)
ajustement_taille(data_Cytofluogram_3D_stade[5], stade_5, max_length_stade5, X0_data_stade5, Y0_data_stade5)
ajustement_taille(data_Cytofluogram_3D_stade[6], stade_6, max_length_stade6, X0_data_stade6, Y0_data_stade6)
ajustement_taille(data_Cytofluogram_3D_stade[7], stade_7, max_length_stade7, X0_data_stade7, Y0_data_stade7)
ajustement_taille(data_Cytofluogram_3D_stade[8], stade_8, max_length_stade8, X0_data_stade8, Y0_data_stade8)
ajustement_taille(data_Cytofluogram_3D_stade[9], stade_9, max_length_stade9, X0_data_stade9, Y0_data_stade9)
ajustement_taille(data_Cytofluogram_3D_stade[10], stade_10, max_length_stade10, X0_data_stade10, Y0_data_stade10)


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

CairoMakie.activate!() # Mieux pour sauvegarder de gros graphiques sans lag graphique

#Boucle for pour tracer les graphiques  

for (n_injection, i) in enumerate(3:10)

    fig_cytofluogramme = Figure(resolution=(800, 400))

    moyenne_cytofluogramme_stade = resultats_moyennes[i]
    x_col = moyenne_cytofluogramme_stade[!, :X0_moyen]
    y_col = moyenne_cytofluogramme_stade[!, :Y0_moyen]

    pcc = cor(x_col, y_col)
    cor_values[n_injection, 2] = pcc

    println("La corrélation pour le stade $i est de : $(round(pcc, digits=3))")

    x_min, x_max = extrema(x_col)
    y_max = maximum(y_col)

    q2_CTAT_quantile = quantile(y_col, 0.5)
    q3_CTAT_quantile = quantile(y_col, 0.95)
    q2_CTAC_quantile = quantile(x_col, 0.5)
    q3_CTAC_quantile = quantile(x_col, 0.95)

    ax = Axis(fig_cytofluogramme[1, 1],
        title="Mean cytofluogramme for stage $i",
        xlabel="Mean CTAC intensity",
        ylabel="Mean CTAT intensity",
        xticks=LinearTicks(5),
        yticks=LinearTicks(5),
        xtickformat="{:.2e}",
        ytickformat="{:.2e}",
        xticklabelrotation=pi / 8)

    Makie.scatter!(ax,
        x_col,
        y_col,
        color=:blue,
        markersize=5,
        alpha=0.6,
        rasterize=true # Réduit le poids si converti en PDF, mais aide aussi au rendu
    )

    text!(ax,
        "PCC = $(round(pcc, digits=3))",
        position=(x_max, y_max - 0.15 * y_max),
        color=:black,
        align=(:right, :top),
        fontsize=12
    )

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

    supertitle = Label(fig_cytofluogramme[0, :],
        "Mean 3D Cytofluogrammes \n for the expression of AAACTAC and AAACTAT forward LncRNA strand in stage $i",
        fontsize=15, font="Arial", padding=(10, 10, 10, 10))

    trim!(fig_cytofluogramme.layout)

    # 1. Définir le chemin proprement (joinpath assemble uniquement des bouts de texte)
    chemin = joinpath("/Users", "verme", "Desktop", "Graph_coloc_3D", "Cytofluogramme_colocalisation_Dvir48_CTAC_CTAT_sensforward_stade_$i.png")

    # 2. Utiliser la fonction save de Makie
    save(chemin, fig_cytofluogramme)

    empty!(fig_cytofluogramme) # Libérer la mémoire
end

function nettoyer_memoire_cytofluogramme()
    empty!(data_Cytofluogram_3D_stade)
    empty!(resultats_moyennes)
    empty!(X0_data_stade3)
    empty!(Y0_data_stade3)
    empty!(X0_data_stade4)
    empty!(Y0_data_stade4)
    empty!(X0_data_stade5)
    empty!(Y0_data_stade5)
    empty!(X0_data_stade6)
    empty!(Y0_data_stade6)
    empty!(X0_data_stade7)
    empty!(Y0_data_stade7)
    empty!(X0_data_stade8)
    empty!(Y0_data_stade8)
    empty!(X0_data_stade9)
    empty!(Y0_data_stade9)
    empty!(X0_data_stade10)
    empty!(Y0_data_stade10)
    empty!(cor_values)
    GC.gc()
    println("=> Mémoire libérée pour les cytofluogrammes avec succès !")
end

nettoyer_memoire_cytofluogramme()

################################################
#importation des données pour les différents stades
################################################
importation("VanSteensel")

max_length_stade3 = mesure_max(data_VanSteensel_3D_stade[3], stade_3)
max_length_stade4 = mesure_max(data_VanSteensel_3D_stade[4], stade_4)
max_length_stade5 = mesure_max(data_VanSteensel_3D_stade[5], stade_5)
max_length_stade6 = mesure_max(data_VanSteensel_3D_stade[6], stade_6)
max_length_stade7 = mesure_max(data_VanSteensel_3D_stade[7], stade_7)
max_length_stade8 = mesure_max(data_VanSteensel_3D_stade[8], stade_8)
max_length_stade9 = mesure_max(data_VanSteensel_3D_stade[9], stade_9)
max_length_stade10 = mesure_max(data_VanSteensel_3D_stade[10], stade_10)

#tout le monde est de la même taille, pas besoin de faire l'ajustement

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

#ajustement des taille et remplissage des dataframes de VanSteensel
#on prend le X3 qui est la moyenne des calcule de dephasage entre les images.
function ajustement_taille_vansteensel(dictionaire, liste_numero_image, max_length, X0_data, Y0_data)

    for i in liste_numero_image

        if length(dictionaire[string(i)][!, :X3]) < max_length

            allowmissing!(dictionaire[string(i)], [:X3, :Y3])
            difference_length = max_length - length(dictionaire[string(i)][!, :X3])

            for j in 1:difference_length

                push!(dictionaire[string(i)][!, :X3], missing)
                push!(dictionaire[string(i)][!, :Y3], missing)
            end

        end
        n_injection = Symbol(findfirst(x -> x == i, liste_numero_image))

        X0_data[!, n_injection] = dictionaire[string(i)][!, :X3]
        Y0_data[!, n_injection] = dictionaire[string(i)][!, :Y3]
    end

end

ajustement_taille_vansteensel(data_VanSteensel_3D_stade[3], stade_3, max_length_stade3, X_data_stade3, Y_data_stade3)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[4], stade_4, max_length_stade4, X_data_stade4, Y_data_stade4)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[5], stade_5, max_length_stade5, X_data_stade5, Y_data_stade5)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[6], stade_6, max_length_stade6, X_data_stade6, Y_data_stade6)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[7], stade_7, max_length_stade7, X_data_stade7, Y_data_stade7)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[8], stade_8, max_length_stade8, X_data_stade8, Y_data_stade8)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[9], stade_9, max_length_stade9, X_data_stade9, Y_data_stade9)
ajustement_taille_vansteensel(data_VanSteensel_3D_stade[10], stade_10, max_length_stade10, X_data_stade10, Y_data_stade10)

##############################################
#Calcul des moyennes des ICQB pour chaque stade
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

CairoMakie.activate!()

fig_VanSteensel = Figure(resolution=(600, 800))

for i in (3:10)

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    moyenne_VanSteensel_stade = resultats_moyennes[i]

    if n_injection <= 2
        ax = Axis(fig_VanSteensel[1, n_injection],
            title="Mean VanSteensel for stage " * string(i),
            xlabel="δ Shift",
            ylabel="PCF")

    elseif n_injection >= 3 && n_injection <= 4
        ax = Axis(fig_VanSteensel[2, n_injection-2],
            title="Mean VanSteensel for stage " * string(i),
            xlabel="δ Shift",
            ylabel="PCF")

    elseif n_injection >= 5 && n_injection <= 6
        ax = Axis(fig_VanSteensel[3, n_injection-4],
            title="Mean VanSteensel for stage " * string(i),
            xlabel="δ Shift",
            ylabel="PCF")

    elseif n_injection >= 7 && n_injection <= 8
        ax = Axis(fig_VanSteensel[4, n_injection-6],
            title="Mean VanSteensel for stage " * string(i),
            xlabel="δ Shift",
            ylabel="PCF")

    end

    for j in eval(Symbol("stade_" * string(i)))

        Makie.lines!(ax,
            data_VanSteensel_3D_stade[i][string(j)][!, :X3],
            data_VanSteensel_3D_stade[i][string(j)][!, :Y3],
            color=(:blue, 0.2),
            linestyle=:solid
        )

    end

    Makie.lines!(ax,
        resultats_moyennes[i][!, :X_moyen],
        resultats_moyennes[i][!, :Y_moyen],
        color=:black,
        linestyle=:solid,
        linewidth=3
    )

end

fig_VanSteensel

supertitle = Label(fig_VanSteensel[0, :],
    "Mean VanSteensel graph for the expression \n of AAACTAC and AAACTAT forward LncRNA strand in stages 3 to 10",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))

# Ajustement de la disposition et sauvegarde
trim!(fig_VanSteensel.layout)
chemin = joinpath("/Users", "verme", "Desktop", "Graph_coloc_3D", "VanSteensel_moyenne_forward_lncARN_CTAC_CTAT_TousStades.png")
save(chemin, fig_VanSteensel)

function nettoyer_memoire_VanSteensel()
    empty!(data_VanSteensel_3D_stade)
    empty!(resultats_moyennes)
    empty!(X0_data_stade3)
    empty!(Y0_data_stade3)
    empty!(X0_data_stade4)
    empty!(Y0_data_stade4)
    empty!(X0_data_stade5)
    empty!(Y0_data_stade5)
    empty!(X0_data_stade6)
    empty!(Y0_data_stade6)
    empty!(X0_data_stade7)
    empty!(Y0_data_stade7)
    empty!(X0_data_stade8)
    empty!(Y0_data_stade8)
    empty!(X0_data_stade9)
    empty!(Y0_data_stade9)
    empty!(X0_data_stade10)
    empty!(Y0_data_stade10)
    empty!(cor_values)
    GC.gc()
    println("=> Mémoire libérée pour les VanSteensel avec succès !")
end

nettoyer_memoire_VanSteensel()


##############################
importation("DIANA")



#importation("ICQ_A")
#importation("ICQ_B")




