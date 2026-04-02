#######################################################
# Test analye de colocalisation 3D
#
#
#######################################################

# importer les bibliothèques nécessaires
using DataFrames
using CSV
using Plots
using Makie
using CairoMakie
using GLMakie

Dictionaire_des_dataframe_ICQ_A = Dict{String,DataFrame}()
Dictionaire_des_dataframe_ICQ_B = Dict{String,DataFrame}()
Dictionaire_des_dataframe_cytofluogramme = Dict{String,DataFrame}()
Dictionaire_des_dataframe_Vansten = Dict{String,DataFrame}()

# Charger les données depuis un fichier CSV
for i in 1:250
    Dictionaire_des_dataframe_ICQ_A["ICQ_A_"*string(i)] = CSV.read("/Users/verme/Desktop/test_coloc2/puis_data_graphique/Resultat_ICA_A_Series_3test_Image_" * string(i) * ".csv", DataFrame)
    Dictionaire_des_dataframe_ICQ_A["ICQ_A_"*string(i)][!, :Z] .= Float64(i)
    Dictionaire_des_dataframe_ICQ_B["ICQ_B_"*string(i)] = CSV.read("/Users/verme/Desktop/test_coloc2/puis_data_graphique/Resultat_ICA_B_Series_3test_Image_" * string(i) * ".csv", DataFrame)
    Dictionaire_des_dataframe_ICQ_B["ICQ_B_"*string(i)][!, :Z] .= Float64(i)
    Dictionaire_des_dataframe_cytofluogramme["Cytofluogramme_"*string(i)] = CSV.read("/Users/verme/Desktop/test_coloc2/puis_data_graphique/Resultat_Cytofluogram_Series_3test_Image_" * string(i) * ".csv", DataFrame)
    Dictionaire_des_dataframe_cytofluogramme["Cytofluogramme_"*string(i)][!, :Z] .= Float64(i)
    Dictionaire_des_dataframe_Vansten["Vansten_"*string(i)] = CSV.read("/Users/verme/Desktop/test_coloc2/puis_data_graphique/Resultat_VanSteensel_Series_3test_Image_" * string(i) * ".csv", DataFrame)
    Dictionaire_des_dataframe_Vansten["Vansten_"*string(i)][!, :Z] .= Float64(i)
end

fig_ICA_A = Figure()

ax = Axis3(fig_ICA_A[1, 1],
    aspect=(1, 1, 1),
    title="Tracés de l'évolution spatiale de l'indice de colocalisation ICQ_A",
    xlabel="Nomeros d'image",
    ylabel="Variance des canaux",
    zlabel="Indice de colocalisation")

for i in 1:250

    Makie.scatter!(ax,
        Dictionaire_des_dataframe_ICQ_A["ICQ_A_"*string(i)][!, :Z],
        Dictionaire_des_dataframe_ICQ_A["ICQ_A_"*string(i)][!, :X],
        Dictionaire_des_dataframe_ICQ_A["ICQ_A_"*string(i)][!, :Y],
        color=:black,
        markersize=5,
        alpha=Float64((i) / 250))

end

fig_ICA_A

record(fig_ICA_A, "/Users/verme/Desktop/rotation_coloc_ICA_A.mp4", 1:120; framerate=30) do frame
    # On fait tourner la caméra de quelques degrés à chaque frame
    ax.azimuth[] = 1.5pi + 0.3 * sin(2pi * frame / 120) # Petite oscillation
    # Ou une rotation complète :
    # ax.azimuth[] = 2pi * frame / 120 
end

saving_path = "/Users/verme/Desktop/Data_ovaire/Script_analise_Julia/graphique_3D_coloc.png"
save(saving_path, fig_ICA_A)

fig_ICA_B = Figure()

ax = Axis3(fig_ICA_B[1, 1], title="Tracés ICQ_B en accumulation")

for i in 1:250

    Makie.scatter!(ax,
        Dictionaire_des_dataframe_ICQ_B["ICQ_B_"*string(i)][!, :Z],
        Dictionaire_des_dataframe_ICQ_B["ICQ_B_"*string(i)][!, :X],
        Dictionaire_des_dataframe_ICQ_B["ICQ_B_"*string(i)][!, :Y],
        color=:black,
        markersize=5,
        alpha=Float64((i) / 250))

end

fig_ICA_B

fig_Cytofluogramme = Figure()
ax = Axis3(fig_ICA_A[1, 1],
    aspect=(1, 1, 1),
    title="Tracés de l'évolution spatiale du cytofluogramme",
    xlabel="Nomeros d'image",
    ylabel="Intensité du canal 1",
    zlabel="Intensité du canal 2")

for i in 1:250
    Makie.scatter!(ax,
        Dictionaire_des_dataframe_cytofluogramme["Cytofluogramme_"*string(i)][!, :Z],
        Dictionaire_des_dataframe_cytofluogramme["Cytofluogramme_"*string(i)][!, :X1],
        Dictionaire_des_dataframe_cytofluogramme["Cytofluogramme_"*string(i)][!, :Y1],
        color=:black,
        markersize=5,
        alpha=Float64((i) / 250))
end
fig_Cytofluogramme

saving_path = "/Users/verme/Desktop/graphique_3D_coloc_cytofluogramme.png"
save(saving_path, fig_Cytofluogramme)


fig_Vansten = Figure()
ax = Axis3(fig_Vansten[1, 1],
    aspect=(1, 1, 1),
    title="Tracés de l'évolution spatiale de l'indice de colocalisation VanSteensel",
    xlabel="Nomeros d'image",
    xticks=0:25:250,
    ylabel="Glissement des images en pixels",
    yticks=-100:25:100,
    zlabel="Indice de colocalisation",
    zticks=-0.1:0.05:0.3)
for i in 1:250
    Makie.scatter!(ax,
        Dictionaire_des_dataframe_Vansten["Vansten_"*string(i)][!, :Z],
        Dictionaire_des_dataframe_Vansten["Vansten_"*string(i)][!, :X0],
        Dictionaire_des_dataframe_Vansten["Vansten_"*string(i)][!, :Y0],
        color=:black,
        markersize=5)
end


display(fig_Vansten)

saving_path = "/Users/verme/Desktop/graphique_3D_coloc_Vansten.png"
save(saving_path, fig_Vansten)

print(Dictionaire_des_dataframe_Vansten["Vansten_1"])

GLMakie.activate!()

# On enregistre une rotation de 360 degrés
record(fig_Vansten, "/Users/verme/Desktop/rotation_coloc.mp4", 1:120; framerate=30) do frame
    # On fait tourner la caméra de quelques degrés à chaque frame
    ax.azimuth[] = 1.5pi + 0.3 * sin(2pi * frame / 120) # Petite oscillation
    # Ou une rotation complète :
    # ax.azimuth[] = 2pi * frame / 120 
end