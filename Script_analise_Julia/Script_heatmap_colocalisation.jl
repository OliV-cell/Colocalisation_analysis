##############################
# Script pour tracer les heatmaps de colocalisation
#############################

#importation des packages
using FilePathsBase
using DataFrames
using CSV
using Plots
using Makie
using CairoMakie
using GLMakie
using Statistics

#importation des données

direction_data_graphique_image_seule_3D = "/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image"

dictonaire_des_donnees_predict_stade = Dict{Int,DataFrame}()

for i in (3:10)
    dictonaire_des_donnees_predict_stade[i] = CSV.read(direction_data_graphique_image_seule_3D * "/df_predict_stade_" * string(i) * ".csv", DataFrame)
end

##############################################
#Visualisation des cytofluogrammes des moyennes
#modifier l'appelle de graphique pouir les heat map
#mettre le meme couleur de scale Resultat_ICA_B_Series_3test_Image_
#bleu doncer vers orange (style temoérature) 

#taille oocyte
oocyte_stade_3 = (3.8, 4.2, 2.3, 3.7)
oocyte_stade_4 = (3.8, 4.2, 2.3, 3.7)
oocyte_stade_5 = (4.2, 3.6, 4.3, 4.4)
oocyte_stade_6 = (4.8, 3.4, 3.2)
oocyte_stade_7 = (4.1, 3.97, 4.1, 4.7)
oocyte_stade_8 = (4.8, 4.3, 5, 4.88)
oocyte_stade_9 = (8.1, 6.9, 7.1, 6.2)
oocyte_stade_10 = (8.2, 8.7, 10, 7.2)

stades_dict = Dict(
    3 => oocyte_stade_3,
    4 => oocyte_stade_4,
    5 => oocyte_stade_5,
    6 => oocyte_stade_6,
    7 => oocyte_stade_7,
    8 => oocyte_stade_8,
    9 => oocyte_stade_9,
    10 => oocyte_stade_10
)

fig_heatmap = Figure(resolution=(600, 800))

#Boucle for pour tracer les graphiques  

for i in (3:10)

    n_injection = Int64(findfirst(x -> x == i, (3:10)))

    heatmap_stade = dictonaire_des_donnees_predict_stade[i]

    println("heatmap stade " * string(i) * " en traitement")

    row = cld(n_injection, 2)
    col = 2 * ((n_injection - 1) % 2) + 1

    ax = Axis(fig_heatmap[row, col],
        title="Heatmap de la colocalisation du stade " * string(i),
        ylabel="Décalage δ",
        xlabel="Profondeur (µm)",
        aspect=AxisAspect(2),
        xticks=LinearTicks(5),
        yticks=LinearTicks(5))

    hm = Makie.heatmap!(ax,
        heatmap_stade.Z,
        heatmap_stade.X3,
        heatmap_stade.fit,
        colormap=:inferno,
        colorscale=asinh,
        fxaa=true
    )

    Colorbar(fig_heatmap[row, col+1], hm,
        ticks=LinearTicks(4), ticklabelsize=10, width=15)
    colsize!(fig_heatmap.layout, col + 1, Auto(0.15))

    hlines!(ax, 0, color=:white, linewidth=1,
        linestyle=:dash)

    vlines!(ax, 50 - (mean(stades_dict[i]) / 2) - (std(stades_dict[i])), color=:white, linewidth=1,
        linestyle=:dash)

    vlines!(ax, 50 + (mean(stades_dict[i]) / 2) + (std(stades_dict[i])), color=:white, linewidth=1,
        linestyle=:dash)

end

trim!(fig_heatmap.layout)

# 1. Définir le chemin proprement (joinpath assemble uniquement des bouts de texte)
chemin = joinpath("/Users", "verme", "Desktop", "Heatmap_colocalisation_Dvir48_CTAC_CTAT_sensforward.png")

# 2. Utiliser la fonction save de Makie
save(chemin, fig_heatmap)

supertitle = Label(fig_heatmap[0, :],
    "Heatmap de la colocalisation",
    fontsize=15, font="Arial", padding=(10, 10, 10, 10))














