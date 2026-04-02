
##############################################
#Premier jet à conserver
##############################################
#essaie avec les données de Dvir48_CTAC_CTAT_sensboth_rep1
#avec le stade 5 pour commencer
numero_image = (4,10,13,14,16,20,21,26,30,32,33,42,48,58,66,68,69,70,72,78,83,87,90,92,94,95)

dictonaire_des_graphique_cytofluogramme = Dict{String,DataFrame}()

for i in numero_image

    dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_" * string(i)] = CSV.read("/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_sonde_sens_Dvir48_CTAC_CTAT/Donnée_des_graphiques/Resultat_CytofluogramDvir48_CTAC_CTAT_sensboth_stade5_rep1_Image_" * string(i) * ".csv", DataFrame)

end

fig_cytofluogramme = Figure()

ax = Axis(fig_cytofluogramme[1,1], title = "Cytofluogrammes des images Dvir48_CTAC_CTAT_sensboth_rep1_stade5")

for i in numero_image

    Makie.scatter!(ax, 
    dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:X0],
    dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:Y0],
    color = :black, 
    markersize = 5,
    alpha = 0.2)

end

function mesure_max(dictionaire,liste_numero_image)

    maximum_length = Array{Union{Missing,Float64}}(missing, length(dictionaire))

    for i  in liste_numero_image
        n_injection = findfirst(x -> x == i, numero_image)

        maximum_length[n_injection] = length(dictionaire["cytofluogramme_image_"*string(i)][!,:X0])

    end

    maximum_length = maximum(maximum_length)

    return maximum_length

end

max_length = mesure_max(dictonaire_des_graphique_cytofluogramme, numero_image)

X0_data = DataFrame()
Y0_data = DataFrame()

for i in numero_image

    if length(dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:X0]) < max_length

        allowmissing!(dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)], [:X0, :Y0])

        difference_length = max_length - length(dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:X0])

        for j in 1:difference_length

            push!(dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:X0], missing)
            push!(dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:Y0], missing)
        end


    end

    n_injection = Symbol(findfirst(x -> x == i, numero_image))

    X0_data[!, n_injection] = dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:X0]
    Y0_data[!, n_injection] = dictonaire_des_graphique_cytofluogramme["cytofluogramme_image_"*string(i)][!,:Y0]

end

print(X0_data[1,:])

moyenne_cytofluogramme = DataFrame(
    X0_moyen = Float64[0],
    Y0_moyen = Float64[0]
)

moyenne_cytofluogramme = repeat(moyenne_cytofluogramme, inner = Int64(max_length))

moyenne_cytofluogramme.X0_moyen = [
    let m = skipmissing(row)
        isempty(m) ? NaN : mean(m)
    end 
    for row in eachrow(X0_data)
]

moyenne_cytofluogramme.Y0_moyen = [
    let m = skipmissing(row)
        isempty(m) ? NaN : mean(m)
    end 
    for row in eachrow(Y0_data)
]


cor(moyenne_cytofluogramme[!,:X0_moyen], moyenne_cytofluogramme[!,:Y0_moyen])

model = lm(@formula(Y0_moyen ~ X0_moyen), moyenne_cytofluogramme)

fig_cytofluogramme = Figure()

ax = Axis(fig_cytofluogramme[1,1], title = "Cytofluogrammes des moyenne des intensité de Dvir48_CTAC_CTAT_sensboth_rep1_stade5")

Makie.scatter!(ax, 
    moyenne_cytofluogramme[!,:X0_moyen],
    moyenne_cytofluogramme[!,:Y0_moyen],
    color = :orange, 
    markersize = 5,
    alpha = 0.6
)
lines!(ax,
    moyenne_cytofluogramme[!,:X0_moyen],
    predict(model),
    color = :black,
    linestyle = (:solid,10),
    linewidth = 1
)
text!(ax,
    "R = " * string(round(cor(moyenne_cytofluogramme[!,:X0_moyen], moyenne_cytofluogramme[!,:Y0_moyen]), digits = 3)),
    position = (2000,1000),
    color = :black,
    align = (:left, :center),
    fontsize = 16
)

###############################################
#Graph seul
cor(moyenne_cytofluogramme_stade3[!,:X0_moyen], moyenne_cytofluogramme_stade3[!,:Y0_moyen])
model = lm(@formula(Y0_moyen ~ X0_moyen), moyenne_cytofluogramme_stade3)

x_min, x_max = extrema(moyenne_cytofluogramme_stade3.X0_moyen)
pred_df = DataFrame(X0_moyen = range(x_min, x_max, length=100))

pr = predict(model, pred_df, interval=:confidence, level=0.95)


fig_cytofluogramme_stade3 = Figure()

ax1 = Axis(fig_cytofluogramme_stade3[1,1], 
    title = "Cytofluogrammes des moyenne des intensité de Dvir48_CTAC_CTAT_sensboth_rep1_stade3",
    xlabel = "X0 Moyen",
    ylabel = "Y0 Moyen") 

ax2 = Axis(fig_cytofluogramme_stade3[2,1], 
    title = "Résidus du modèle linéaire",
    xlabel = "X0 Moyen",
    ylabel = "Résidus")

ax3 = Axis(fig_cytofluogramme_stade3[3,1], 
    title = "Résidus du modèle linéaire",
    xlabel = "X0 Moyen",
    ylabel = "Résidus")

band!(ax1,
    pred_df.X0_moyen,
    pr[:,2],
    pr[:,3],
    color = (:green, 0.3)
) +
Makie.scatter!(ax1, 
    moyenne_cytofluogramme_stade3[!,:X0_moyen],
    moyenne_cytofluogramme_stade3[!,:Y0_moyen],
    color = :blue, 
    markersize = 5,
    alpha = 0.6
) +
lines!(ax1,
    moyenne_cytofluogramme_stade3[!,:X0_moyen],
    predict(model),
    color = :black,
    linestyle = (:solid,10),
    linewidth = 1
) +
text!(ax1,
    "R = " * string(round(cor(moyenne_cytofluogramme_stade3[!,:X0_moyen], moyenne_cytofluogramme_stade3[!,:Y0_moyen]), digits = 3)),
    position = (x_max, maximum(moyenne_cytofluogramme_stade3.Y0_moyen)),
    color = :black,
    align = (:right, :top),
    fontsize = 16
)
