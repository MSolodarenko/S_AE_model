include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
#using PyPlot

using Interpolations
using Statistics
using JLD2
using ProgressMeter

#num_lambdas = 23

function create_plot(X::Vector{Float64},XLABEL::String,Y::Vector{Float64},YLABEL::String, IS_Y_PERCENTAGE::Bool=true, OCCUPATION::String="H", LOW_LIMIT=-Inf)
    num_lambdas = length(X)

    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    end

    YLIMS1 = max(minimum(Y), LOW_LIMIT)
    YLIMS2 = maximum(Y)
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
    DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))
    YPOS = mean(YTICKS[end-1:end])
    if maximum(Y[calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
        YPOS = mean(YTICKS[1:2])
    end
    if IS_Y_PERCENTAGE
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = scatter(X,Y,
                    #=regression=true,=#
                    color=COLOR,
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    yticks = YTICKS,
                    ylims = YLIMS )
    vline!([X[calibrated_lambda]], color="grey")
    if IS_Y_PERCENTAGE
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]*100; digits=DIGITS+1))%) "
    else
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]; digits=DIGITS+1))) "
    end

    annotate!([X[calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))
    return plt
end
function create_plots(X::Vector{Float64},XLABEL::String,Ys,YLABELs, IS_Y_PERCENTAGE::Bool=true, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)
    num_lambdas = length(X)

    plts = []
    half_range = maximum([ maximum(filter(!isnan,copy(Ys[i])))-minimum(filter(!isnan,copy(Ys[i]))) for i in 1:length(Ys)])/2.0
    for y_i = 1:length(Ys)
        Y = Ys[y_i]
        Ynanless = copy(Y)
        Ynanless = filter(!isnan,Ynanless)
        ns = num_lambdas - (length(Y) - length(Ynanless))
        YLABEL = YLABELs[y_i]

        YLIMS1 = max(middle(Ynanless[1:ns])-half_range, LOW_LIMIT)
        if IS_Y_PERCENTAGE
            YLIMS1 = max(YLIMS1, min(0.0, minimum(filter(!isnan,copy(Y)))))
        end
        YLIMS2 = YLIMS1+2*half_range
        YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
        YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
        YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))+1
        YPOS = mean(YTICKS[end-1:end])
        if maximum(filter(!isnan,copy(Y[calibrated_lambda:calibrated_lambda+4]))) > YTICKS[end-1]
            YPOS = mean(YTICKS[1:2])
        end
        if IS_Y_PERCENTAGE
            YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
        else
            YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
        end

        COLOR="blue"
        if OCCUPATION[y_i]=="W"
            COLOR="purple"
        elseif OCCUPATION[y_i]=="SP"
            COLOR="red"
        elseif OCCUPATION[y_i]=="EMP"
            COLOR="green"
        end
        plt = scatter(X,Y,
                        #=regression=true,=#
                        color = COLOR,
                        legend = false,
                        xlabel = XLABEL,
                        ylabel = YLABEL,
                        yticks = YTICKS,
                        ylims = YLIMS )
        vline!([X[calibrated_lambda]], color="grey")
        if IS_Y_PERCENTAGE
            TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]*100; digits=DIGITS+1))%) "
        else
            TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]; digits=DIGITS+1))) "
        end
        annotate!([X[calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))
        push!(plts, plt)
    end
    return plts
end
function create_combined_plot(X::Vector{Float64},XLABEL::String,Ys,YLABELs,YLABEL, IS_Y_PERCENTAGE::Bool=true, LEGENDPOS=:bottomright, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)
    num_lambdas = length(X)

    YLIMS1 = max(minimum(minimum.(Ys)), LOW_LIMIT)
    YLIMS2 = maximum(maximum.(Ys))
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
    DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))
    YPOS = mean(YTICKS[end-1:end])

    COLORS=[]
    for y_i = 1:length(Ys)
        if OCCUPATION[y_i]=="W"
            push!(COLORS,"purple")
        elseif OCCUPATION[y_i]=="SP"
            push!(COLORS,"red")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"green")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"yellow")
        else
            push!(COLORS,"blue")
        end

        if maximum(Ys[y_i][calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
            YPOS = mean(YTICKS[1:2])
        end
    end
    if IS_Y_PERCENTAGE
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = vline([X[calibrated_lambda]], color="grey",label="")
    if IS_Y_PERCENTAGE
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1))) "
    else
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1))) "
    end
    annotate!([X[calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))

    for y_i = 1:length(Ys)
        plt = scatter!(X,Ys[y_i],
                        #=regression=true,=#
                        color=COLORS[y_i],
                        legend=LEGENDPOS,
                        xlabel=XLABEL,
                        ylabel=YLABEL,
                        label=YLABELs[y_i],
                        yticks = YTICKS,
                        ylims = YLIMS )
    end

    return plt
end
function create_combined_plot(Xs,XLABEL::String,Ys,YLABELs,YLABEL, IS_Y_PERCENTAGE::Bool=true, LEGENDPOS=:bottomright, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)
    num_lambdas = length(Xs[1])

    YLIMS1 = max(minimum(minimum.(Ys)), LOW_LIMIT)
    YLIMS2 = maximum(maximum.(Ys))
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
    DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))
    YPOS = mean(YTICKS[end-1:end])

    COLORS=[]
    for y_i = 1:length(Ys)
        if OCCUPATION[y_i]=="W"
            push!(COLORS,"purple")
        elseif OCCUPATION[y_i]=="SP"
            push!(COLORS,"red")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"green")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"yellow")
        else
            push!(COLORS,"blue")
        end

        if maximum(Ys[y_i][calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
            YPOS = mean(YTICKS[1:2])
        end
    end
    if IS_Y_PERCENTAGE
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = vline([Xs[1][calibrated_lambda]], color="grey",label="")
    if IS_Y_PERCENTAGE
        TEXT = " Calibrated economy ($(round(Xs[1][calibrated_lambda]; digits=DIGITS+1))) "
    else
        TEXT = " Calibrated economy ($(round(Xs[1][calibrated_lambda]; digits=DIGITS+1))) "
    end
    annotate!([Xs[1][calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))

    for y_i = 1:length(Ys)
        plt = scatter!(Xs[y_i],Ys[y_i],
                        #=regression=true,=#
                        color=COLORS[y_i],
                        legend=LEGENDPOS,
                        xlabel=XLABEL,
                        ylabel=YLABEL,
                        label=YLABELs[y_i],
                        yticks = YTICKS,
                        ylims = YLIMS )
    end

    return plt
end

country = "italy"

LOCAL_DIR_SOURCE = "$(@__DIR__)/Results/Stationary/GE_lambda/General/"
if Sys.iswindows()
    LOCAL_DIR_SOURCE = "$(@__DIR__)\\Results\\Stationary\\GE_lambda\\General\\"
end
@load "$(LOCAL_DIR_SOURCE)SSS.jld2" SSS C_Ys Outputs Incomes Consumptions Rs Ws logcs loges giniWs giniEnts share_unbound means ginis avglogs varlogs avgs vars quantile_means TFPis TFPds mean_MPL var_MPL mean_MPK var_MPK Capital share_W_earnings_in_output share_SP_earnings_in_output share_EMP_earnings_in_output share_W_capital_income_in_output share_SP_capital_income_in_output share_EMP_capital_income_in_output

lambdas = zeros(length(SSS))
@showprogress for i = 1:length(SSS)
    lambdas[i] = SSS[i][5][1]
end
LAMBDA = 1.665907
calibrated_lambda = findfirst(x -> x==LAMBDA, lambdas)

LOCAL_DIR_SOURCE_fixed = "$(@__DIR__)/Results/Stationary/GE_lambda_fixed_occ/General/"
if Sys.iswindows()
    LOCAL_DIR_SOURCE_fixed = "\\Results\\Stationary\\GE_lambda_fixed_occ\\General\\"
end
@load "$(LOCAL_DIR_SOURCE_fixed)SSS_fixed.jld2" SSS_fixed_occ C_Ys_fixed_occ Outputs_fixed_occ Incomes_fixed_occ Consumptions_fixed_occ Rs_fixed_occ Ws_fixed_occ logcs_fixed_occ loges_fixed_occ giniWs_fixed_occ giniEnts_fixed_occ share_unbound_fixed_occ means_fixed_occ ginis_fixed_occ avglogs_fixed_occ varlogs_fixed_occ avgs_fixed_occ vars_fixed_occ quantile_means_fixed_occ TFPis_fixed_occ TFPds_fixed_occ mean_MPL_fixed_occ var_MPL_fixed_occ mean_MPK_fixed_occ var_MPK_fixed_occ Capital_fixed_occ share_W_earnings_in_output_fixed_occ share_SP_earnings_in_output_fixed_occ share_EMP_earnings_in_output_fixed_occ share_W_capital_income_in_output_fixed_occ share_SP_capital_income_in_output_fixed_occ share_EMP_capital_income_in_output_fixed_occ

LOCAL_DIR = "$(@__DIR__)/Results/Stationary/GE_lambda_comparison/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Stationary\\GE_lambda_comparison\\"
end
mkpath(LOCAL_DIR)

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)


function generate_plots(X,XLABEL,Y,Y_fixed_occ,YLABEL,PATHDIR,FILENAME,IS_PERCENTAGE::Bool=false)

    plt = create_combined_plot(X,XLABEL, [Y, Y_fixed_occ],["w occ mob", "wo occ mob"],YLABEL, IS_PERCENTAGE, :bottomright)
    display(plt)
    savefig(plt,"$(PATHDIR)$(country)_$(FILENAME).png")

    if !IS_PERCENTAGE
        plt = create_plot(X,XLABEL, (Y.-Y_fixed_occ)./Y,YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_diff.png")

        plt = create_plot(X,XLABEL, (Y.-(Y[calibrated_lambda]-Y_fixed_occ[calibrated_lambda]).-Y_fixed_occ)./Y,YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_did.png")
    else
        plt = create_plot(X,XLABEL, (Y.-Y_fixed_occ),YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_diff.png")

        plt = create_plot(X,XLABEL, (Y.-(Y[calibrated_lambda]-Y_fixed_occ[calibrated_lambda]).-Y_fixed_occ),YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_did.png")
    end
end

generate_plots(lambdas,"λ",C_Ys,C_Ys_fixed_occ,"Credit/Output",LOCAL_DIR_GENERAL,"lambda_credit_to_gdp",false)

generate_plots(lambdas,"λ",Outputs,Outputs_fixed_occ,"Output",LOCAL_DIR_GENERAL,"lambda_outputs",false)

generate_plots(lambdas,"λ",Incomes,Incomes_fixed_occ,"Income",LOCAL_DIR_GENERAL,"lambda_Incomes",false)

generate_plots(lambdas,"λ",Consumptions,Consumptions_fixed_occ,"Consumptions",LOCAL_DIR_GENERAL,"lambda_Consumptions",false)

generate_plots(lambdas,"λ",Capital[1,:],Capital_fixed_occ[1,:],"Capital",LOCAL_DIR_GENERAL,"lambda_Capital",false)
generate_plots(lambdas,"λ",Capital[2,:],Capital_fixed_occ[2,:],"Workers' Capital",LOCAL_DIR_GENERAL,"lambda_Capital_W",false)
generate_plots(lambdas,"λ",Capital[3,:],Capital_fixed_occ[3,:],"Sole Proprietors' Capital",LOCAL_DIR_GENERAL,"lambda_Capital_SP",false)
generate_plots(lambdas,"λ",Capital[4,:],Capital_fixed_occ[4,:],"Employers' Capital",LOCAL_DIR_GENERAL,"lambda_Capital_EMP",false)


#Interest rate and wage
generate_plots(lambdas,"λ",Rs,Rs_fixed_occ,"Interest Rate",LOCAL_DIR_GENERAL,"lambda_interest_rate",true)

generate_plots(lambdas,"λ",Ws,Ws_fixed_occ,"Wage",LOCAL_DIR_GENERAL,"lambda_wage",false)


LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#Credit/Output to: log-consumption, log-earnings
generate_plots(lambdas,"λ",logcs,logcs_fixed_occ,"Variance of log-consumption",LOCAL_DIR_INEQUALITY,"lambda_log_consumption",false)

generate_plots(lambdas,"λ",loges,loges_fixed_occ,"Variance of log-earnings",LOCAL_DIR_INEQUALITY,"lambda_log_earnings",false)

#Crdit/Output to Gini Workers and Entrepreneurs
generate_plots(lambdas,"λ",giniWs,giniWs_fixed_occ,"Gini for workers' income",LOCAL_DIR_INEQUALITY,"lambda_gini_workers",true)

generate_plots(lambdas,"λ",giniEnts,giniEnts_fixed_occ,"Gini for entrepreneurs' income",LOCAL_DIR_INEQUALITY,"lambda_gini_entrepreneurs",true)


for s = 1:4 # [1] = income, earnings, wealth, consumption
    # if s == 1
    stat_name = "Income"
    if s == 2
        stat_name = "Earnings"
    elseif s == 3
        stat_name = "Wealth"
    elseif s == 4
        stat_name = "Consumption"
    end

    CHOICE_NAMES = ["Households","Workers","Sole Proprietors","Employers","Entrepreneurs"]
    OCCS = ["H","W","SP","EMP"]

    LABELS=["1st","2nd","3rd","4th","5th"]
    for h=1:4
        try
            LABEL = "Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles) w occ mob"
            plt1 = create_combined_plot(lambdas,"λ", [quantile_means[1,h,s,:],quantile_means[2,h,s,:],quantile_means[3,h,s,:],quantile_means[4,h,s,:],quantile_means[5,h,s,:]],LABELS,LABEL, false, ["H","W","SP","EMP","ENT"])
            LABEL = "Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles) wo occ mob"
            plt2 = create_combined_plot(lambdas,"λ", [quantile_means_fixed_occ[1,h,s,:],quantile_means_fixed_occ[2,h,s,:],quantile_means_fixed_occ[3,h,s,:],quantile_means_fixed_occ[4,h,s,:],quantile_means_fixed_occ[5,h,s,:]],LABELS,LABEL, false, ["H","W","SP","EMP","ENT"])
            plt = plot(plt1,plt2,legend=false)
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)lambda_combined_mean_$(stat_name)_$(OCCS[h])_quantiles.png")

            LABEL = "Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles) diff"
            plot_mat = Array{Any}(undef,5)
            for iii=1:5
                plot_mat[iii] = (quantile_means[iii,h,s,:].-quantile_means_fixed_occ[iii,h,s,:])./quantile_means[iii,h,s,:]
            end
            plt = create_combined_plot(lambdas,"λ", plot_mat,LABELS,LABEL, true, ["H","W","SP","EMP","ENT"])
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)lambda_combined_mean_$(stat_name)_$(OCCS[h])_quantiles_diff.png")

            LABEL = "Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles) did"
            plot_mat = Array{Any}(undef,5)
            for iii=1:5
                plot_mat[iii] = (quantile_means[iii,h,s,:].-(quantile_means[iii,h,s,calibrated_lambda]-quantile_means_fixed_occ[iii,h,s,calibrated_lambda]).-quantile_means_fixed_occ[iii,h,s,:])./quantile_means[iii,h,s,:]
            end
            plt = create_combined_plot(lambdas,"λ", plot_mat,LABELS,LABEL, true, ["H","W","SP","EMP","ENT"])
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)lambda_combined_mean_$(stat_name)_$(OCCS[h])_quantiles_did.png")
        catch e
            println_sameline("Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles) have not been generated")
        end

    end
    # calculate mean
    #means[s,h,i]
    LABELS = ["Mean of $cn' $stat_name" for cn in CHOICE_NAMES]
    for h in 1:4
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        generate_plots(lambdas,"λ",means[s,h,:],means_fixed_occ[s,h,:],LABELS[h],LOCAL_DIR_INEQUALITY,"lambda_mean_$(choice_name)_$(stat_name)",false)
    end

    # calculate gini coefficent
    #ginis[s,h,i]
    LABELS = ["Gini of $cn' $stat_name" for cn in CHOICE_NAMES]
    for h in 1:4
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        generate_plots(lambdas,"λ",ginis[s,h,:],ginis_fixed_occ[s,h,:],LABELS[h],LOCAL_DIR_INEQUALITY,"lambda_gini_$(choice_name)_$(stat_name)",true)
    end
    # calculate variance of log-s
    #avgs[s,h,i]
    LABELS = ["Average of $cn' $stat_name" for cn in CHOICE_NAMES]
    for h in 1:4
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        generate_plots(lambdas,"λ",avgs[s,h,:],avgs_fixed_occ[s,h,:],LABELS[h],LOCAL_DIR_INEQUALITY,"lambda_avg_$(choice_name)_$(stat_name)",false)
    end

    #vars[s,h,i]
    LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    for h in 1:4
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        generate_plots(lambdas,"λ",vars[s,h,:],vars_fixed_occ[s,h,:],LABELS[h],LOCAL_DIR_INEQUALITY,"lambda_var_$(choice_name)_$(stat_name)",false)
    end

    # calculate variance of log-s
    #avglogs[s,h,i]
    LABELS = ["Average of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    for h in 1:4
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        generate_plots(lambdas,"λ",avglogs[s,h,:],avglogs_fixed_occ[s,h,:],LABELS[h],LOCAL_DIR_INEQUALITY,"lambda_avg_$(choice_name)_log$(stat_name)",false)
    end

    #varlogs[s,h,i]
    if s!=3
        LABELS = ["Variance of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    else
        LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    end
    for h in 1:4
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        generate_plots(lambdas,"λ",varlogs[s,h,:],varlogs_fixed_occ[s,h,:],LABELS[h],LOCAL_DIR_INEQUALITY,"lambda_var_$(choice_name)_log$(stat_name)",false)
    end
end

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)

# TFP_ideal for SP,EMP
generate_plots(lambdas,"λ",TFPis[1,:],TFPis_fixed_occ[1,:],"TFP for Entrepreneurs",LOCAL_DIR_PRODUCTIVITY,"lambda_tfp_ideal_ent",false)

generate_plots(lambdas,"λ",TFPis[2,:],TFPis_fixed_occ[2,:],"TFP for Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_tfp_ideal_sp",false)

generate_plots(lambdas,"λ",TFPis[3,:],TFPis_fixed_occ[3,:],"TFP for Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_tfp_ideal_emp",false)

# TFP_data for SP,EMP
generate_plots(lambdas,"λ",TFPds[2,:],TFPds_fixed_occ[2,:],"TFP for Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_tfp_data_sp",false)

generate_plots(lambdas,"λ",TFPds[3,:],TFPds_fixed_occ[3,:],"TFP for Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_tfp_data_emp",false)

# mean of MPL for ENT, SP,EMP
generate_plots(lambdas,"λ",mean_MPL[1,:],mean_MPL_fixed_occ[1,:],"Mean of MPL for Entrepreneurs",LOCAL_DIR_PRODUCTIVITY,"lambda_mean_mpl_ent",false)

generate_plots(lambdas,"λ",mean_MPL[2,:],mean_MPL_fixed_occ[2,:],"Mean of MPL for Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_mean_mpl_sp",false)

generate_plots(lambdas,"λ",mean_MPL[3,:],mean_MPL_fixed_occ[3,:],"Mean of MPL for Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_mean_mpl_emp",false)

# variance of MPL for ENT, SP,EMP
generate_plots(lambdas,"λ",var_MPL[1,:],var_MPL_fixed_occ[1,:],"Variance of MPL for Entrepreneurs",LOCAL_DIR_PRODUCTIVITY,"lambda_var_mpl_ent",false)

generate_plots(lambdas,"λ",var_MPL[2,:],var_MPL_fixed_occ[2,:],"Variance of MPL for Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_var_mpl_sp",false)

generate_plots(lambdas,"λ",var_MPL[3,:],var_MPL_fixed_occ[3,:],"Variance of MPL for Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_var_mpl_emp",false)

# mean of MPK for ENT, SP,EMP
generate_plots(lambdas,"λ",mean_MPK[1,:],mean_MPK_fixed_occ[1,:],"Mean of MPK for Entrepreneurs",LOCAL_DIR_PRODUCTIVITY,"lambda_mean_mpk_ent",false)

generate_plots(lambdas,"λ",mean_MPK[2,:],mean_MPK_fixed_occ[2,:],"Mean of MPK for Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_mean_mpk_sp",false)

generate_plots(lambdas,"λ",mean_MPK[3,:],mean_MPK_fixed_occ[3,:],"Mean of MPK for Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_mean_mpk_emp",false)

# variance of MPK for ENT, SP,EMP
generate_plots(lambdas,"λ",var_MPK[1,:],var_MPK_fixed_occ[1,:],"Variance of MPK for Entrepreneurs",LOCAL_DIR_PRODUCTIVITY,"lambda_var_mpk_ent",false)

generate_plots(lambdas,"λ",var_MPK[2,:],var_MPK_fixed_occ[2,:],"Variance of MPK for Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_var_mpk_sp",false)

generate_plots(lambdas,"λ",var_MPK[3,:],var_MPK_fixed_occ[3,:],"Variance of MPK for Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_var_mpk_emp",false)

# share of unbound ENT, SP,EMP
generate_plots(lambdas,"λ",share_unbound[1,:],share_unbound_fixed_occ[1,:],"Share of Unconstrained Entrepreneurs",LOCAL_DIR_PRODUCTIVITY,"lambda_share_ent_unbound",true)

generate_plots(lambdas,"λ",share_unbound[2,:],share_unbound_fixed_occ[2,:],"Share of Unconstrained Sole Proprietors",LOCAL_DIR_PRODUCTIVITY,"lambda_share_se_unbound",true)

generate_plots(lambdas,"λ",share_unbound[3,:],share_unbound_fixed_occ[3,:],"Share of Unconstrained Employers",LOCAL_DIR_PRODUCTIVITY,"lambda_share_emp_unbound",true)

#share of earnings in output
generate_plots(lambdas,"λ",share_W_earnings_in_output,share_W_earnings_in_output_fixed_occ,"Share of output as Workers' Earnings",LOCAL_DIR_PRODUCTIVITY,"lambda_share_of_output_W_earnings",true)

generate_plots(lambdas,"λ",share_SP_earnings_in_output,share_SP_earnings_in_output_fixed_occ,"Share of output as Sole Proprietors' Earnings",LOCAL_DIR_PRODUCTIVITY,"lambda_share_of_output_SP_earnings",true)

generate_plots(lambdas,"λ",share_EMP_earnings_in_output,share_EMP_earnings_in_output_fixed_occ,"Share of output as Employers' Earnings",LOCAL_DIR_PRODUCTIVITY,"lambda_share_of_output_EMP_earnings",true)

#share of capital income in output
generate_plots(lambdas,"λ",share_W_capital_income_in_output,share_W_capital_income_in_output_fixed_occ,"Share of output as Workers' Capital Income",LOCAL_DIR_PRODUCTIVITY,"lambda_share_of_output_W_capital_income",true)

generate_plots(lambdas,"λ",share_SP_capital_income_in_output,share_SP_capital_income_in_output_fixed_occ,"Share of output as Sole Proprietors' Capital Income",LOCAL_DIR_PRODUCTIVITY,"lambda_share_of_output_SP_capital_income",true)

generate_plots(lambdas,"λ",share_EMP_capital_income_in_output,share_EMP_capital_income_in_output_fixed_occ,"Share of output as Employers' Capital Income",LOCAL_DIR_PRODUCTIVITY,"lambda_share_of_output_EMP_capital_income",true)
