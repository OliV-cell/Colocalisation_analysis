##################################
#
#Script d'analyse des la oclocalisation de
#CTAC et CTAT sonde antisens en 3D
#
#####################################

#faire test avec tweedie pour donner au stade 3 sans valeurs négatives
#idem pour autre stade et regarder les assomptions
#présentement le scat offre de quoi d'intéréssant mais est problématique
#regarder les autres méthodes

#importation des packages
using FilePathsBase
using DataFrames
using CSV
using Plots
using Makie
using CairoMakie
using GLMakie
using Statistics
using GLM
using HypothesisTests
using Glob
using RCall

#liste des numéro d'image

stade_3 = (1, 4, 6, 10, 15, 22, 30, 33, 35, 48, 59)
stade_4 = (2, 3, 17, 28, 44, 52, 58)
stade_5 = (1, 3, 4, 7, 9, 11, 13, 15, 26, 29, 32, 35, 46, 47, 57)
stade_6 = (2, 14, 19, 25, 27, 38, 50, 53)
stade_7 = (5, 8, 12, 20, 31, 34, 45)
stade_8 = (21, 23, 24, 37, 39, 42, 49)
stade_9 = (36, 40, 41, 51)
stade_10 = (54, 55, 56)

stades_dict = Dict(
    3 => stade_3,
    4 => stade_4,
    5 => stade_5,
    6 => stade_6,
    7 => stade_7,
    8 => stade_8,
    9 => stade_9,
    10 => stade_10
)

#initialisation des dictonnaires 

data_VanSteensel_image_seule_stade = Dict{Int,Dict{Int,DataFrame}}()
data_cytofluogramme_image_seule_stade = Dict{Int,Dict{Int,DataFrame}}()
data_ICQ_A_image_seule_stade = Dict{Int,Dict{Int,DataFrame}}()
data_ICQ_B_image_seule_stade = Dict{Int,Dict{Int,DataFrame}}()
data_DIANA_image_seule_stade = Dict{Int,Dict{Int,DataFrame}}()

#importation des données pour les différents stades

direction_data_graphique_image_seule_3D = "/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image"

# ======================================================== #
# FONCTIONS POUR IMPORTER LES DONNEES
function importer_stade(stade_num)
    println("\n=== Importation du stade $stade_num ===")
    stade_list = stades_dict[stade_num]

    for i in stade_list
        # Cytofluogramme
        fichiers_cyto = glob("Resultat_Cytofluogram_Series_$(i)Dvir48_CTAC_CTAT_antisens_stade*_Image_*.csv", direction_data_graphique_image_seule_3D)
        for n_fichier in 1:length(fichiers_cyto)
            get!(data_cytofluogramme_image_seule_stade, i, Dict{Int,DataFrame}())[n_fichier] = CSV.read(fichiers_cyto[n_fichier], DataFrame)
            data_cytofluogramme_image_seule_stade[i][n_fichier][!, :Z] .= Float64(n_fichier)
        end

        # VanSteensel
        fichiers_van = glob("Resultat_VanSteensel_Series_$(i)Dvir48_CTAC_CTAT_antisens_stade*_Image_*.csv", direction_data_graphique_image_seule_3D)
        for n_fichier in 1:length(fichiers_van)
            get!(data_VanSteensel_image_seule_stade, i, Dict{Int,DataFrame}())[n_fichier] = CSV.read(fichiers_van[n_fichier], DataFrame)
            data_VanSteensel_image_seule_stade[i][n_fichier][!, :Z] .= Float64(n_fichier)
        end

        # ICQ A
        fichiers_icq_a = glob("Resultat_ICQ_A_Series_$(i)Dvir48_CTAC_CTAT_antisens_stade*_Image_*.csv", direction_data_graphique_image_seule_3D)
        for n_fichier in 1:length(fichiers_icq_a)
            get!(data_ICQ_A_image_seule_stade, i, Dict{Int,DataFrame}())[n_fichier] = CSV.read(fichiers_icq_a[n_fichier], DataFrame)
            data_ICQ_A_image_seule_stade[i][n_fichier][!, :Z] .= Float64(n_fichier)
        end

        # ICQ B
        fichiers_icq_b = glob("Resultat_ICQ_B_Series_$(i)Dvir48_CTAC_CTAT_antisens_stade*_Image_*.csv", direction_data_graphique_image_seule_3D)
        for n_fichier in 1:length(fichiers_icq_b)
            get!(data_ICQ_B_image_seule_stade, i, Dict{Int,DataFrame}())[n_fichier] = CSV.read(fichiers_icq_b[n_fichier], DataFrame)
            data_ICQ_B_image_seule_stade[i][n_fichier][!, :Z] .= Float64(n_fichier)
        end

        println("=> Fichiers pour l'image $i (stade $stade_num) importés.")
    end
end

function nettoyer_memoire()
    empty!(data_VanSteensel_image_seule_stade)
    empty!(data_cytofluogramme_image_seule_stade)
    empty!(data_ICQ_A_image_seule_stade)
    empty!(data_ICQ_B_image_seule_stade)
    GC.gc()
    println("=> Mémoire libérée avec succès !")
end

function mesure_max(dictionaire, liste_numero_image)
    max_val = 0.0

    for j in 3:10

        println("=> Longueur maximale pour le stade $j")

        for i in liste_numero_image[j]
            if haskey(dictionaire, j)
                # On cherche le nombre de lignes maximum parmi tous les fichiers de l'image 'i'
                for df in values(dictionaire[j])
                    if hasproperty(df, :X0)
                        len = Float64(length(df[!, :X0]))
                        if len > max_val
                            max_val = len
                        end
                    end
                end
            end

            println("=> Longueur maximale pour l'image $i : $max_val")
        end

    end

    return max_val

end

function aligner_et_completer!(dictionaire, stade_num)
    liste_numero_image = stades_dict[stade_num]

    for i in liste_numero_image
        if !haskey(dictionaire, i) || isempty(dictionaire[i])
            continue
        end

        series = dictionaire[i]

        # 1. Identifier z_max en trouvant le Y0 MAXIMUM à X0=0, sous condition que Y0 <= 0.4 (pour ignorer les anomalies de bruit)
        max_y0_global = -Inf
        z_max = -1

        for (z, df) in series
            if hasproperty(df, :X3) && hasproperty(df, :Y3)
                idx_zero = findfirst(x -> !ismissing(x) && x == 0, df[!, :X3])

                if idx_zero !== nothing
                    valeur_y0_a_zero = df[idx_zero, :Y3]
                    if !ismissing(valeur_y0_a_zero) && valeur_y0_a_zero <= 0.5 && valeur_y0_a_zero > max_y0_global
                        max_y0_global = valeur_y0_a_zero
                        z_max = z
                    end
                end
            end
        end

        if z_max == -1
            println("=> ERREUR: Impossible de trouver un z_max valide pour l'image $i.")
            continue
        end

        # 2. Décalage pour aligner z_max sur Z=100 et supprimer Z < 0
        shift = 100 - z_max
        new_series = Dict{Int,DataFrame}()
        original_keys = sort(collect(keys(series)))

        for z in original_keys
            new_z = z + shift
            if new_z >= 0 && new_z <= 200
                df_copy = copy(series[z])
                df_copy[!, :Z] .= Float64(new_z)
                new_series[new_z] = df_copy
            end
        end

        # 3. Compléter depuis 0 et jusqu'à 200 en piochant aléatoirement dans les bordures
        current_max_z = isempty(new_series) ? -1 : maximum(keys(new_series))
        current_min_z = isempty(new_series) ? 201 : minimum(keys(new_series))

        if !isempty(new_series) && (current_max_z < 200 || current_min_z > 0)
            # Recherche de la frame avec la moyenne Y0 la plus basse de toute la série pour le padding
            lowest_mean_y0 = Inf
            best_z = -1

            for z in original_keys
                if hasproperty(series[z], :Y3)
                    y0_elements = skipmissing(series[z][!, :Y3])
                    if !isempty(y0_elements)
                        local_mean = mean(y0_elements)
                        if local_mean < lowest_mean_y0
                            lowest_mean_y0 = local_mean
                            best_z = z
                        end
                    end
                end
            end

            if best_z == -1
                best_z = original_keys[1]
            end

            pool = [best_z]

            # Remplissage par le haut (jusqu'à 200)
            if current_max_z < 200
                for pad_z in (current_max_z+1):200
                    chosen_z = rand(pool)
                    padded_df = copy(series[chosen_z])
                    padded_df[!, :Z] .= Float64(pad_z)
                    new_series[pad_z] = padded_df
                end
            end

            # Remplissage par le bas (depuis 0 jusqu'à current_min_z - 1)
            if current_min_z > 0
                for pad_z in 0:(current_min_z-1)
                    chosen_z = rand(pool)
                    padded_df = copy(series[chosen_z])
                    padded_df[!, :Z] .= Float64(pad_z)
                    new_series[pad_z] = padded_df
                end
            end
        end

        # 4. Appliquer les changements à la série
        dictionaire[i] = new_series
        println("=> Alignement Z=100 et padding aléatoire jusqu'à Z=200 terminés pour l'image $i.")
    end

    println("=> Opération terminée pour toutes les images du stade $stade_num.")
end

# ======================================================== #
# IMPORTATION DES DONNEES POUR LES DIFFERENTS STADES

## --- STADE 3 ---
importer_stade(3)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou
aligner_et_completer!(data_VanSteensel_image_seule_stade, 3)
#Tout le monde a la même longeur, pas d'ajustement nessesaire


# Fusionner tous les DataFrames de `data_VanSteensel_image_seule_stade` en un seul DataFrame final
Data_fusionner_stade_3 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)

    # 1. Delete positive parasites
    all_y3_pur = filter(y -> y <= 0.4, all_y3)
    all_y3_pur = filter(y -> y >= -0.05, all_y3_pur)

    poids_image = 1.0 / var(all_y3_pur)

    #on coupe les ourlier positif = bruit parasite
    #et les valeurs négatives = bruit de fond du a l'absence de déconvolution
    for (z, df) in dict_z
        # 1. Delete positive parasites from the spatial map
        df_copy = filter(row -> row.Y3 <= 0.4, df)
        df_copy = filter(row -> row.Y3 >= -0.05, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_3, df_copy)
    end
end

Data_fusionner_stade_3 = vcat(Data_fusionner_stade_3...)
select!(Data_fusionner_stade_3, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_3))

# on enleve les valeurs négatives (COMMENTED OUT TO KEEP BIOLOGY)
#minimum(Data_fusionner_stade_3[!, :Y3])
#for i in 1:length(Data_fusionner_stade_3[!, :Y3])
#    if Data_fusionner_stade_3[i, :Y3] < 0
#        Data_fusionner_stade_3[i, :Y3] = 0
#    end
#end

histogram(Data_fusionner_stade_3[!, :Y3])


@rput Data_fusionner_stade_3

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_3$ImageNumber = as.factor(Data_fusionner_stade_3$ImageNumber)

model_stade_3 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(30,30)),
                    data = Data_fusionner_stade_3,
                    method = "fREML",
                    select = TRUE,
                    weights = Poids,
                    gamma = 1.4,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )

summary(model_stade_3)
gam.check(model_stade_3)
k.check(model_stade_3)
summary.gam(model_stade_3)
"""

R"""
vis.gam(model_stade_3, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 40,                
        color = "topo",          
        main = "3D Surface of Colocalization")
"""

R"""
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_3, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""

Z_dense = range(minimum(Data_fusionner_stade_3.Z), maximum(Data_fusionner_stade_3.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_3.X3), maximum(Data_fusionner_stade_3.X3), length=100)

df_predict_stade = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade

R"""
prediction_stade_3 <- predict(model_stade_3, newdata = df_predict_stade, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade$fit <- prediction_stade_3$fit
df_predict_stade$se_fit <- prediction_stade_3$se.fit
"""
@rget df_predict_stade

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade[!, :Z] .= df_predict_stade[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_3.csv", df_predict_stade)

heatmap_stade_3 = Figure()
ax = Axis(heatmap_stade_3[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade.Z, df_predict_stade.X3, df_predict_stade.fit)
Colorbar(heatmap_stade_3[1, 2], hm)



heatmap_stade_3

nettoyer_memoire()

#========================================#
## --- STADE 4 ---
importer_stade(4)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou
aligner_et_completer!(data_VanSteensel_image_seule_stade, 4)
#Tout le monde a la même longeur, pas d'ajustement nessesaire

Data_fusionner_stade_4 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.45, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.05, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.45, df)
        df_copy = filter(row -> row.Y3 >= -0.05, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_4, df_copy)
    end
end

Data_fusionner_stade_4 = vcat(Data_fusionner_stade_4...)
select!(Data_fusionner_stade_4, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_4))
histogram(Data_fusionner_stade_4[!, :Y3])

@rput Data_fusionner_stade_4

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_4$ImageNumber = as.factor(Data_fusionner_stade_4$ImageNumber)

model_stade_4 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(30,30)),
                    data = Data_fusionner_stade_4,
                    method = "fREML",
                    weights = Poids,
                    gamma = 1.4,
                    select = TRUE,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )

summary(model_stade_4)
gam.check(model_stade_4)
k.check(model_stade_4)
summary.gam(model_stade_4)

"""

R"""

vis.gam(model_stade_4, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")

"""

R"""
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_4, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""


Z_dense = range(minimum(Data_fusionner_stade_4.Z), maximum(Data_fusionner_stade_4.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_4.X3), maximum(Data_fusionner_stade_4.X3), length=100)

df_predict_stade_4 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_4

R"""
prediction_stade_4 <- predict(model_stade_4, newdata = df_predict_stade_4, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_4$fit <- prediction_stade_4$fit
df_predict_stade_4$se_fit <- prediction_stade_4$se.fit
"""
@rget df_predict_stade_4

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_4[!, :Z] .= df_predict_stade_4[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_4.csv", df_predict_stade_4)

heatmap_stade_4 = Figure()
ax = Axis(heatmap_stade_4[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_4.Z, df_predict_stade_4.X3, df_predict_stade_4.fit)
Colorbar(heatmap_stade_4[1, 2], hm)

heatmap_stade_4

#modele n'est pas vraiment bon, on va essayer de le modifier
#refaire avec le poid et arranger le modele, il est dégeulasse
#on conserve le scat et relativement ok, viduellement le modèle fit bien les données

nettoyer_memoire()

#==========================================#
## --- STADE 5 ---
importer_stade(5)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou

aligner_et_completer!(data_VanSteensel_image_seule_stade, 5)
#Tout le monde a la même longeur, pas d'ajustement nessesaire

#représentation graphique des données

Data_fusionner_stade_5 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.45, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.01, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.45, df)
        df_copy = filter(row -> row.Y3 >= -0.01, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_5, df_copy)
    end
end
Data_fusionner_stade_5 = vcat(Data_fusionner_stade_5...)
select!(Data_fusionner_stade_5, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_5))
histogram(Data_fusionner_stade_5[!, :Y3])

@rput Data_fusionner_stade_5

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_5$ImageNumber = as.factor(Data_fusionner_stade_5$ImageNumber)

model_stade_5 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(35,35)),
                    data = Data_fusionner_stade_5,
                    method = "fREML",
                    weights = Poids,
                    select = TRUE,
                    gamma = 1.4,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )

summary(model_stade_5)
gam.check(model_stade_5)
k.check(model_stade_5)
summary.gam(model_stade_5)

"""

R"""
vis.gam(model_stade_5, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")

"""

R"""
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_5, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""

Z_dense = range(minimum(Data_fusionner_stade_5.Z), maximum(Data_fusionner_stade_5.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_5.X3), maximum(Data_fusionner_stade_5.X3), length=100)

df_predict_stade_5 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_5

R"""
prediction_stade_5 <- predict(model_stade_5, newdata = df_predict_stade_5, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_5$fit <- prediction_stade_5$fit
df_predict_stade_5$se_fit <- prediction_stade_5$se.fit
"""
@rget df_predict_stade_5

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_5[!, :Z] .= df_predict_stade_5[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_5.csv", df_predict_stade_5)

heatmap_stade_5 = Figure()
ax = Axis(heatmap_stade_5[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_5.Z, df_predict_stade_5.X3, df_predict_stade_5.fit)
Colorbar(heatmap_stade_5[1, 2], hm)

heatmap_stade_5

nettoyer_memoire()

#==========================================#
## --- STADE 6 ---
importer_stade(6)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou

aligner_et_completer!(data_VanSteensel_image_seule_stade, 6)
#Tout le monde a la même longeur, pas d'ajustement nessesaire

#représentation graphique des données

Data_fusionner_stade_6 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.5, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.05, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.5, df)
        df_copy = filter(row -> row.Y3 >= -0.05, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_6, df_copy)
    end
end
Data_fusionner_stade_6 = vcat(Data_fusionner_stade_6...)
select!(Data_fusionner_stade_6, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_6))
histogram(Data_fusionner_stade_6[!, :Y3])

@rput Data_fusionner_stade_6

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_6$ImageNumber = as.factor(Data_fusionner_stade_6$ImageNumber)

model_stade_6 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(30,30)),
                    data = Data_fusionner_stade_6,
                    method = "fREML",
                    weights = Poids,
                    select = TRUE,
                    gamma = 1.4,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )


summary(model_stade_6)
gam.check(model_stade_6)
k.check(model_stade_6)

summary.gam(model_stade_6)

"""

R"""
vis.gam(model_stade_6, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")
"""

R"""
library(performance)
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_6, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""


Z_dense = range(minimum(Data_fusionner_stade_6.Z), maximum(Data_fusionner_stade_6.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_6.X3), maximum(Data_fusionner_stade_6.X3), length=100)

df_predict_stade_6 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_6

R"""
prediction_stade_6 <- predict(model_stade_6, newdata = df_predict_stade_6, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_6$fit <- prediction_stade_6$fit
df_predict_stade_6$se_fit <- prediction_stade_6$se.fit
"""
@rget df_predict_stade_6

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_6[!, :Z] .= df_predict_stade_6[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_6.csv", df_predict_stade_6)

heatmap_stade_6 = Figure()
ax = Axis(heatmap_stade_6[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_6.Z, df_predict_stade_6.X3, df_predict_stade_6.fit)
Colorbar(heatmap_stade_6[1, 2], hm)

heatmap_stade_6

nettoyer_memoire()

#==========================================#
## --- STADE 7 ---
importer_stade(7)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou

aligner_et_completer!(data_VanSteensel_image_seule_stade, 7)
#Tout le monde a la même longeur, pas d'ajustement nessesaire

#représentation graphique des données

Data_fusionner_stade_7 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.5, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.1, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.5, df)
        df_copy = filter(row -> row.Y3 >= -0.1, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_7, df_copy)
    end
end
Data_fusionner_stade_7 = vcat(Data_fusionner_stade_7...)
select!(Data_fusionner_stade_7, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_7))
histogram(Data_fusionner_stade_7[!, :Y3])

@rput Data_fusionner_stade_7

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_7$ImageNumber = as.factor(Data_fusionner_stade_7$ImageNumber)

model_stade_7 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(30,30)),
                    data = Data_fusionner_stade_7,
                    method = "fREML",
                    weights = Poids,
                    select = TRUE,
                    gamma = 1.4,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )


summary(model_stade_7)
gam.check(model_stade_7)
k.check(model_stade_7)

summary.gam(model_stade_7)

"""

R"""
vis.gam(model_stade_7, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")
"""

R"""
library(performance)
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_7, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""


Z_dense = range(minimum(Data_fusionner_stade_7.Z), maximum(Data_fusionner_stade_7.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_7.X3), maximum(Data_fusionner_stade_7.X3), length=100)

df_predict_stade_7 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_7

R"""
prediction_stade_7 <- predict(model_stade_7, newdata = df_predict_stade_7, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_7$fit <- prediction_stade_7$fit
df_predict_stade_7$se_fit <- prediction_stade_7$se.fit
"""
@rget df_predict_stade_7

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_7[!, :Z] .= df_predict_stade_7[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_7.csv", df_predict_stade_7)

heatmap_stade_7 = Figure()
ax = Axis(heatmap_stade_7[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_7.Z, df_predict_stade_7.X3, df_predict_stade_7.fit)
Colorbar(heatmap_stade_7[1, 2], hm)

heatmap_stade_7

nettoyer_memoire()

#==========================================#
## --- STADE 8 ---
importer_stade(8)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou

aligner_et_completer!(data_VanSteensel_image_seule_stade, 8)
#Tout le monde a la même longeur, pas d'ajustement nessesaire

#représentation graphique des données

Data_fusionner_stade_8 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.5, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.1, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.5, df)
        df_copy = filter(row -> row.Y3 >= -0.1, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_8, df_copy)
    end
end
Data_fusionner_stade_8 = vcat(Data_fusionner_stade_8...)
select!(Data_fusionner_stade_8, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_8))
histogram(Data_fusionner_stade_8[!, :Y3])

@rput Data_fusionner_stade_8

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_8$ImageNumber = as.factor(Data_fusionner_stade_8$ImageNumber)

model_stade_8 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(30,30)),
                    data = Data_fusionner_stade_8,
                    method = "fREML",
                    weights = Poids,
                    select = TRUE,
                    gamma = 1.4,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )


summary(model_stade_8)
gam.check(model_stade_8)
k.check(model_stade_8)

summary.gam(model_stade_8)

"""

R"""
vis.gam(model_stade_8, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")
"""

R"""
library(performance)
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_8, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""

Z_dense = range(minimum(Data_fusionner_stade_8.Z), maximum(Data_fusionner_stade_8.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_8.X3), maximum(Data_fusionner_stade_8.X3), length=100)

df_predict_stade_8 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_8

R"""
prediction_stade_8 <- predict(model_stade_8, newdata = df_predict_stade_8, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_8$fit <- prediction_stade_8$fit
df_predict_stade_8$se_fit <- prediction_stade_8$se.fit
"""
@rget df_predict_stade_8

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_8[!, :Z] .= df_predict_stade_8[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_8.csv", df_predict_stade_8)

heatmap_stade_8 = Figure()
ax = Axis(heatmap_stade_8[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_8.Z, df_predict_stade_8.X3, df_predict_stade_8.fit)
Colorbar(heatmap_stade_8[1, 2], hm)

heatmap_stade_8

nettoyer_memoire()
#==========================================#
## --- STADE 9 ---
# importer_stade(9)
# ...

importer_stade(9)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou

aligner_et_completer!(data_VanSteensel_image_seule_stade, 9)
#Tout le monde a la même longeur, pas d'ajustement nessesaire

Data_fusionner_stade_9 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.40, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.05, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.5, df)
        df_copy = filter(row -> row.Y3 >= -0.05, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_9, df_copy)
    end
end
Data_fusionner_stade_9 = vcat(Data_fusionner_stade_9...)
select!(Data_fusionner_stade_9, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_9))
histogram(Data_fusionner_stade_9[!, :Y3])

@rput Data_fusionner_stade_9

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_9$ImageNumber = as.factor(Data_fusionner_stade_9$ImageNumber)

model_stade_9 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(30,30)),
                    data = Data_fusionner_stade_9,
                    method = "fREML",
                    weights = Poids,
                    select = TRUE,
                    gamma = 1.8,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )


summary(model_stade_9)
gam.check(model_stade_9)
k.check(model_stade_9)

summary.gam(model_stade_9)

"""

R"""
vis.gam(model_stade_9, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")
"""

R"""
library(performance)
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_9, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""
#représentation graphique des données
#le modele est à retoucher pour la variance

Z_dense = range(minimum(Data_fusionner_stade_9.Z), maximum(Data_fusionner_stade_9.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_9.X3), maximum(Data_fusionner_stade_9.X3), length=100)

df_predict_stade_9 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_9

R"""
prediction_stade_9 <- predict(model_stade_9, newdata = df_predict_stade_9, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_9$fit <- prediction_stade_9$fit
df_predict_stade_9$se_fit <- prediction_stade_9$se.fit
"""
@rget df_predict_stade_9

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_9[!, :Z] .= df_predict_stade_9[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_9.csv", df_predict_stade_9)

heatmap_stade_9 = Figure()
ax = Axis(heatmap_stade_9[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_9.Z, df_predict_stade_9.X3, df_predict_stade_9.fit)
Colorbar(heatmap_stade_9[1, 2], hm)

heatmap_stade_9
nettoyer_memoire()

#==========================================#
## --- STADE 10 ---
importer_stade(10)
#on regarde uniquement les données de VanSteensel
#les autres sont ici quelque peux redondante mais on les garde au cas ou

aligner_et_completer!(data_VanSteensel_image_seule_stade, 10)
#Tout le monde a la même longeur, pas d'ajustement nessesaire
Data_fusionner_stade_10 = DataFrame[]
for (img_num, dict_z) in data_VanSteensel_image_seule_stade

    all_y3 = vcat([df[!, :Y3] for df in values(dict_z)]...)
    all_y3_pur = filter(y -> y <= 0.5, all_y3)
    all_y3_min_pur = filter(y -> y >= -0.1, all_y3_pur)
    poids_image = 1.0 / var(all_y3_min_pur)

    for (z, df) in dict_z
        df_copy = filter(row -> row.Y3 <= 0.5, df)
        df_copy = filter(row -> row.Y3 >= -0.01, df_copy)
        df_copy[!, :ImageNumber] .= img_num
        df_copy[!, :Poids] .= poids_image
        push!(Data_fusionner_stade_10, df_copy)
    end
end

Data_fusionner_stade_10 = vcat(Data_fusionner_stade_10...)
select!(Data_fusionner_stade_10, Not([:X0, :Y0, :X1, :Y1, :X2, :Y2]))
println("Dimension du DataFrame fusionné: ", size(Data_fusionner_stade_10))
histogram(Data_fusionner_stade_10[!, :Y3])

@rput Data_fusionner_stade_10

R"""
library(mgcv)
library(nlme)

Data_fusionner_stade_10$ImageNumber = as.factor(Data_fusionner_stade_10$ImageNumber)

model_stade_10 = bam(Y3 ~ te(Z, X3, bs = c("fs", "fs"), k = c(25,25)),
                    data = Data_fusionner_stade_10,
                    method = "fREML",
                    select = TRUE,
                    weights = Poids,
                    gamma = 1.2,
                    discrete = TRUE,
                    optimizer = c("efs"),
                    family = scat(link="identity"),
                    control = gam.control(maxit = 1000)
                    )


summary(model_stade_10)
gam.check(model_stade_10)
k.check(model_stade_10)

summary.gam(model_stade_10)

"""

R"""
vis.gam(model_stade_10, 
        view = c("Z", "X3"),     
        theta = 45,              
        phi = 30,                
        color = "topo",          
        main = "3D Surface of Colocalization")
"""

R"""
library(performance)
library(DHARMa)
mod_sim_overlap <- simulateResiduals(fittedModel = model_stade_10, n = 250)

plot(mod_sim_overlap ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_overlap)
testOutliers(mod_sim_overlap)
testDispersion(mod_sim_overlap) 
testZeroInflation(mod_sim_overlap) 

"""

Z_dense = range(minimum(Data_fusionner_stade_10.Z), maximum(Data_fusionner_stade_10.Z), length=100)
X3_dense = range(minimum(Data_fusionner_stade_10.X3), maximum(Data_fusionner_stade_10.X3), length=100)

df_predict_stade_10 = DataFrame(
    Z=repeat(Z_dense, outer=100),
    X3=repeat(X3_dense, inner=100)
)

@rput df_predict_stade_10

R"""
prediction_stade_10 <- predict(model_stade_10, newdata = df_predict_stade_10, se.fit = TRUE, type = "link", discrete = TRUE)
df_predict_stade_10$fit <- prediction_stade_10$fit
df_predict_stade_10$se_fit <- prediction_stade_10$se.fit
"""
@rget df_predict_stade_10

#on rajoute lechelle en um, car 1 z-step = 0.5 um
df_predict_stade_10[!, :Z] .= df_predict_stade_10[!, :Z] * 0.5

CSV.write("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_image_par_image/df_predict_stade_10.csv", df_predict_stade_10)

heatmap_stade_10 = Figure()
ax = Axis(heatmap_stade_10[1, 1],
    xticks=(0:5:100, string.(0:5:100)),
    yticks=(-100:10:100, string.(-100:10:100)))
hm = Makie.heatmap!(ax, df_predict_stade_10.Z, df_predict_stade_10.X3, df_predict_stade_10.fit)
Colorbar(heatmap_stade_10[1, 2], hm)

heatmap_stade_10

#représentation graphique des données

