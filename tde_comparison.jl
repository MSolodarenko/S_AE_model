using JLD2
using Plots

using ProgressMeter
using SchumakerSpline
using Interpolations

function create_combined_plot(X,XLABEL::String,Ys,YLABELs,YLABEL,Y1s,Y2s, IS_Y_PERCENTAGE::Bool=false, LEGENDPOS=false, OCCUPATION=["ALL","ENT"])
    TICKSFONTSIZE = 12

    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=length(Ys[1]),length=NUM_XTICKS))))
    if length(X) <= 25
        if length(X) <= 10
            #XTICKS = Int64.(round.(collect(range(0;stop=T,step=1))))
            XTICKS = Int64.(round.(collect(range(1;stop=length(Ys[1]),step=1))))
        else
            #XTICKS = Int64.(round.(collect(range(0;stop=length(Y),step=5))))
            XTICKS = Int64.(round.(collect(range(0;stop=length(Ys[1]),step=5))))
        end
    end
    XTICKS = [1; XTICKS[2:end]]

    NUM_YTICKS = 7
    YLIMS1 = minimum([minimum.(Ys); minimum.(Y1s); minimum.(Y2s)])
    YLIMS2 = maximum([maximum.(Ys); maximum.(Y1s); maximum.(Y2s)])
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))

    COLORS=[]
    for y_i = 1:length(Ys)
        if OCCUPATION[y_i]=="W"
            push!(COLORS,"purple")
        elseif OCCUPATION[y_i]=="SP"
            push!(COLORS,"red")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"green")
        elseif OCCUPATION[y_i]=="ENT"
            # push!(COLORS,"yellow")
            push!(COLORS,RGB(255/255, 221/255, 0/255))
        else
            # push!(COLORS,"blue")
            push!(COLORS,RGB(0/255, 87/255, 183/255))
        end

    end

    if IS_Y_PERCENTAGE
        DIGITS = Int(max(2, round(log10(0.01/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = plot()
    for y_i = 1:length(Ys)
        if YLABELs[y_i] == "W"
            YLABELs[y_i] = "Workers"
        elseif YLABELs[y_i] == "SP"
            YLABELs[y_i] = "Sole Prop."
        elseif YLABELs[y_i] == "EMP"
            YLABELs[y_i] = "Employers"
        elseif YLABELs[y_i] == "ENT"
            YLABELs[y_i] = "Entrepreneurs"
        end
        #Y1line = ones(length(Ys[y_i])+1).*Y1s[y_i]
        Y1line = ones(length(Ys[y_i])).*Y1s[y_i]
        #Y2line = ones(length(Ys[y_i])+1).*Y2s[y_i]
        Y2line = ones(length(Ys[y_i])).*Y2s[y_i]
        #Y1line[10+2:end] .= NaN
        Y1line[10+1:end] .= NaN
        #Y2line[1:40] .= NaN
        Y2line[1:40-1] .= NaN
        #plot!(plt,collect([0; X]), collect.([[Y1s[y_i]; Ys[y_i]],Y1line,Y2line]),
        #plot!(plt,collect(X), collect.([Ys[y_i] ,Y1line,Y2line]),
        plot!(plt,collect(X), collect(Ys[y_i]),
                        #color=[COLORS[y_i] "green" "red"],
                        #color=[COLORS[y_i] COLORS[y_i] COLORS[y_i]],
                        color=COLORS[y_i],
                        #linestyle=[:solid :dot :dash],
                        linestyle=:solid,
                        legend=LEGENDPOS,
                        legendfontsize=ceil(Int64,TICKSFONTSIZE*0.72),
                        xlabel=XLABEL,
                        #label=[YLABELs[y_i] "" ""],
                        linewidth=2.5, thickness_scaling = 1,
                        xtickfontsize=TICKSFONTSIZE,
                        ytickfontsize=TICKSFONTSIZE,
                        xguidefontsize=TICKSFONTSIZE,
                        yguidefontsize=TICKSFONTSIZE,
                        label=YLABELs[y_i],
                        ylabel = YLABEL,
                        xticks = XTICKS,
                        yticks = YTICKS,
                        ylims = YLIMS )
        # text annotation
        #=if Y1s[y_i] != Y2s[y_i]
            if IS_Y_PERCENTAGE
                TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1s[y_i]*100;digits=DIGITS+1))%)"
            else
                TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1s[y_i];digits=DIGITS+1)))"
            end
            POS = Y1s[y_i]+YLIMMARGIN*1.25
            if Y1s[y_i] < Y2s[y_i]
                POS = Y1s[y_i]-YLIMMARGIN*1.25
                if POS < YLIMS1-YLIMMARGIN
                    POS = Y1s[y_i]+YLIMMARGIN*1.25
                end
            end
            annotate!([X[end]], POS, text(TEXT1, COLORS[y_i], :right, 7))

            if IS_Y_PERCENTAGE
                TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2s[y_i]*100;digits=DIGITS+1))%)"
            else
                TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2s[y_i];digits=DIGITS+1)))"
            end
            POS = Y2s[y_i]-YLIMMARGIN*1.25
            if POS < YLIMS1-YLIMMARGIN || Y1s[y_i] < Y2s[y_i]
                POS = Y2s[y_i]+YLIMMARGIN*1.25
            end
            annotate!([X[end]], POS, text(TEXT2, COLORS[y_i], :right, 7))
        end=#
    end
    return plt
end

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE\\General\\"
end
@load "$(LOCAL_DIR)trans_SSS.jld2" trans_SSS


LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE_fixed_occ/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE_fixed_occ\\General\\"
end
@load "$(LOCAL_DIR)trans_SSS_fixed.jld2" trans_SSS_fixed

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE_comparison/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE_comparison\\"
end
mkpath(LOCAL_DIR)

country = "Italy"
TIME_PERIODS = 50#100
Time = collect(1:TIME_PERIODS)

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)

function generate_plots(X,XLABEL,Y,Y_fixed_occ,YLABEL,PATHDIR,FILENAME,IS_PERCENTAGE::Bool=false, LEGENDPOS=false)

    plt = create_combined_plot(X,XLABEL, [Y[3], Y_fixed_occ[3]],["with mobility", "without mobility"],#=YLABEL=#"", [Y[1],Y_fixed_occ[1]],[Y[2],Y_fixed_occ[2]], IS_PERCENTAGE, LEGENDPOS)
    display(plt)
    savefig(plt,"$(PATHDIR)$(country)_$(FILENAME).png")

    # plt = create_combined_plot(X,XLABEL, [Y[3], Y_fixed_occ[3].+(Y[3][1]-Y_fixed_occ[3][1])],["with mobility", "without mobility"],YLABEL, [Y[1],Y_fixed_occ[1]+(Y[1]-Y_fixed_occ[1])],[Y[2],Y_fixed_occ[2]+(Y[1]-Y_fixed_occ[1])], IS_PERCENTAGE, :bottomright)
    # display(plt)
    # savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_adjusted.png")

    if !IS_PERCENTAGE
        # plt = create_plot(X,XLABEL, (Y[3].-Y_fixed_occ[3])./Y[3],YLABEL, true)
        # display(plt)
        # savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_diff.png")
        #
        # plt = create_plot(X,XLABEL, (Y[3].-(Y[3][1]-Y_fixed_occ[3][1]).-Y_fixed_occ[3])./Y[3],YLABEL, true)
        # display(plt)
        # savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_did.png")

        plt = create_combined_plot(X,XLABEL, [(Y[3].-Y[3][1])./Y[3][1], (Y_fixed_occ[3].-Y_fixed_occ[3][1])./Y_fixed_occ[3][1]],["with mobility", "without mobility"],#=YLABEL=#"", [0.0,0.0],[(Y[2].-Y[1])./Y[1],(Y_fixed_occ[2].-Y_fixed_occ[1])./Y_fixed_occ[1]], true, LEGENDPOS)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_growth_rate.png")
    else
        # plt = create_plot(X,XLABEL, (Y[3].-Y_fixed_occ[3]),YLABEL, true)
        # display(plt)
        # savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_diff.png")
        #
        # plt = create_plot(X,XLABEL, (Y[3].-(Y[3][1]-Y_fixed_occ[3][1]).-Y_fixed_occ[3]),YLABEL, true)
        # display(plt)
        # savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_did.png")

        plt = create_combined_plot(X,XLABEL, [Y[3].-Y[3][1], Y_fixed_occ[3].-Y_fixed_occ[3][1]],["with mobility", "without mobility"],#=YLABEL=#"", [0.0,0.0],[Y[2].-Y[1],Y_fixed_occ[2].-Y_fixed_occ[1]], true, LEGENDPOS)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_growth_rate.png")
    end
end

#Interest rate and wage
generate_plots(Time,"Time",trans_SSS[3],trans_SSS_fixed[3],"Interest Rate",LOCAL_DIR_GENERAL,"time_interest_rate",true)

generate_plots(Time,"Time",trans_SSS[4],trans_SSS_fixed[4],"Wage",LOCAL_DIR_GENERAL,"time_wage",false)

#Capital - All,W,SP,EMP
generate_plots(Time,"Time",[trans_SSS[5][1][1],trans_SSS[5][2][1],trans_SSS[5][3][1,:]],[trans_SSS_fixed[5][1][1],trans_SSS_fixed[5][2][1],trans_SSS_fixed[5][3][1,:]],"Capital",LOCAL_DIR_GENERAL,"time_capital",false)
generate_plots(Time,"Time",[trans_SSS[5][1][2],trans_SSS[5][2][2],trans_SSS[5][3][2,:]],[trans_SSS_fixed[5][1][2],trans_SSS_fixed[5][2][2],trans_SSS_fixed[5][3][2,:]],"Workers' Capital",LOCAL_DIR_GENERAL,"time_capital_W",false)
generate_plots(Time,"Time",[trans_SSS[5][1][3],trans_SSS[5][2][3],trans_SSS[5][3][3,:]],[trans_SSS_fixed[5][1][3],trans_SSS_fixed[5][2][3],trans_SSS_fixed[5][3][3,:]],"Sole Proprietors' Capital",LOCAL_DIR_GENERAL,"time_capital_SP",false)
generate_plots(Time,"Time",[trans_SSS[5][1][4],trans_SSS[5][2][4],trans_SSS[5][3][4,:]],[trans_SSS_fixed[5][1][4],trans_SSS_fixed[5][2][4],trans_SSS_fixed[5][3][4,:]],"Employers' Capital",LOCAL_DIR_GENERAL,"time_capital_EMP",false)

#Consumption
generate_plots(Time,"Time",trans_SSS[6],trans_SSS_fixed[6],"Consumption",LOCAL_DIR_GENERAL,"time_consumption",false)
#Credit     - All,SP,EMP,ENT
generate_plots(Time,"Time",[trans_SSS[7][1][1],trans_SSS[7][2][1],trans_SSS[7][3][1,:]],[trans_SSS_fixed[7][1][1],trans_SSS_fixed[7][2][1],trans_SSS_fixed[7][3][1,:]],"Credit",LOCAL_DIR_GENERAL,"time_credit",false)
#Credit-to-output

#Income
generate_plots(Time,"Time",trans_SSS[9],trans_SSS_fixed[9],"Income",LOCAL_DIR_GENERAL,"time_income",false)
#Output
generate_plots(Time,"Time",trans_SSS[10],trans_SSS_fixed[10],"Output",LOCAL_DIR_GENERAL,"time_output",false)

LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)

# income, earnings, wealth, consumption
for s = 1:4
    #if s==1
    stat_name = "Income"
    if s==2
        stat_name = "Earnings"
    elseif s==3
        stat_name = "Wealth"
    elseif s==4
        stat_name = "Consumption"
    end
    # All, W, SP, EMP, ENT
    CHOICE_NAMES = ["Households","Workers","Sole Proprietors","Employers","Entrepreneurs"]
    LOCAL_COLORS = ["H","W","SP","EMP","ENT"]
    for h = 1:4#:5
        choice_name = CHOICE_NAMES[h]
        choice_name_short = LOCAL_COLORS[h]

        #ginis
        generate_plots(Time,"Time",[trans_SSS[11][1][s,h],trans_SSS[11][2][s,h],trans_SSS[11][3][s,h,:]],[trans_SSS_fixed[11][1][s,h],trans_SSS_fixed[11][2][s,h],trans_SSS_fixed[11][3][s,h,:]],"Gini of $(choice_name)' $(stat_name)",LOCAL_DIR_INEQUALITY,"time_gini_$(stat_name)_$(choice_name_short)",true)

        #means
        generate_plots(Time,"Time",[trans_SSS[12][1][s,h],trans_SSS[12][2][s,h],trans_SSS[12][3][s,h,:]],[trans_SSS_fixed[12][1][s,h],trans_SSS_fixed[12][2][s,h],trans_SSS_fixed[12][3][s,h,:]],"Mean of $(choice_name)' $(stat_name)",LOCAL_DIR_INEQUALITY,"time_mean_$(stat_name)_$(choice_name_short)",false)

        #quantile mean
        try
            quantile_means_1 = trans_SSS[14][1][:,h,s]
            quantile_means_2 = trans_SSS[14][2][:,h,s]
            quantile_means = [trans_SSS[14][3][q,h,s,:] for q=1:5]

            quantile_means_fixed_1 = trans_SSS_fixed[14][1][:,h,s]
            quantile_means_fixed_2 = trans_SSS_fixed[14][2][:,h,s]
            quantile_means_fixed = [trans_SSS_fixed[14][3][q,h,s,:] for q=1:5]

            LABELS=["1st","2nd","3rd","4th","5th"]

            LABEL = "Mean of $(choice_name)' $(stat_name) (quantiles) w occ mob"
            plt1 = create_combined_plot(Time,"Time", quantile_means,LABELS,LABEL,quantile_means_1,quantile_means_2, false, :bottomright, ["H","W","SP","EMP","ENT"])
            LABEL = "Mean of $(choice_name)' $(stat_name) (quantiles) wo occ mob"
            plt2 = create_combined_plot(Time,"Time", quantile_means_fixed,LABELS,LABEL,quantile_means_fixed_1,quantile_means_fixed_2, false, :bottomright, ["H","W","SP","EMP","ENT"])
            plt = plot(plt1,plt2)
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_$(choice_name_short)_quantiles.png")

            LABEL = "Mean of $(choice_name)' $(stat_name) (quantiles) wo occ mob (adjusted)"
            quantile_means_fixed_adj = [quantile_means_fixed[q].+(quantile_means[q][1]-quantile_means_fixed[q][1]) for q=1:5]
            quantile_means_fixed_adj_1 = [quantile_means_fixed_1[q].+(quantile_means[q][1]-quantile_means_fixed[q][1]) for q=1:5]
            quantile_means_fixed_adj_2 = [quantile_means_fixed_2[q].+(quantile_means[q][1]-quantile_means_fixed[q][1]) for q=1:5]
            plt3 = create_combined_plot(Time,"Time", quantile_means_fixed_adj,LABELS,LABEL,quantile_means_fixed_adj_1,quantile_means_fixed_adj_2, false, :bottomright, ["H","W","SP","EMP","ENT"])
            plt = plot(plt1,plt3)
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_$(choice_name_short)_quantiles_adj.png")

            LABEL = "Mean of $(choice_name)' $(stat_name) (quantiles) diff"
            quantile_means_diff = [(quantile_means[q].-quantile_means_fixed[q])./quantile_means[q] for q=1:5]
            plts = [create_plot(Time,"Time", quantile_means_diff[q],LABELS[q], true, choice_name_short, false, 3, false ) for q = 1:5]
            plt = plot(plts[1],plts[2],plts[3],plts[4],plts[5],layout=(5,1))
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_$(choice_name_short)_quantiles_diff.png")

            LABEL = "Mean of $(choice_name)' $(stat_name) (quantiles) did"
            quantile_means_did = [(quantile_means[q].-(quantile_means[q][1]-quantile_means_fixed[q][1]).-quantile_means_fixed[q])./quantile_means[q] for q=1:5]
            plts = [create_plot(Time,"Time", quantile_means_did[q],LABELS[q], true, choice_name_short, false, 3, false ) for q = 1:5]
            plt = plot(plts[1],plts[2],plts[3],plts[4],plts[5],layout=(5,1))
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_$(choice_name_short)_quantiles_did.png")

        catch e

        end
    end
end

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)

#share of earnings in output
generate_plots(Time,"Time",trans_SSS[15],trans_SSS_fixed[15],"Share of output as Workers' Earnings",LOCAL_DIR_PRODUCTIVITY,"time_share_of_output_W_earnings",true)

generate_plots(Time,"Time",trans_SSS[16],trans_SSS_fixed[16],"Share of output as Sole Proprietors' Earnings",LOCAL_DIR_PRODUCTIVITY,"time_share_of_output_SP_earnings",true)

generate_plots(Time,"Time",trans_SSS[17],trans_SSS_fixed[17],"Share of output as Employers' Earnings",LOCAL_DIR_PRODUCTIVITY,"time_share_of_output_EMP_earnings",true)

#share of capital income in output
generate_plots(Time,"Time",trans_SSS[18],trans_SSS_fixed[18],"Share of output as Workers' Capital Income",LOCAL_DIR_PRODUCTIVITY,"time_share_of_output_W_capital_income",true)

generate_plots(Time,"Time",trans_SSS[19],trans_SSS_fixed[19],"Share of output as Sole Proprietors' Capital Income",LOCAL_DIR_PRODUCTIVITY,"time_share_of_output_SP_capital_income",true)

generate_plots(Time,"Time",trans_SSS[20],trans_SSS_fixed[20],"Share of output as Employers' Capital Income",LOCAL_DIR_PRODUCTIVITY,"time_share_of_output_EMP_capital_income",true)

# TFP_ideal, mean and var of MPK and MPL
occs = ["SP","EMP","ENT"]
occupations = ["Sole Proprietors","Employers","Entrepreneurs"]
for o in 1:length(occs)
    temp = [trans_SSS[21][1][o],trans_SSS[21][2][o],trans_SSS[21][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[21][1][o],trans_SSS_fixed[21][2][o],trans_SSS_fixed[21][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"TFP for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_tfp_$(occs[o])",false)

    temp = [trans_SSS[22][1][o],trans_SSS[22][2][o],trans_SSS[22][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[22][1][o],trans_SSS_fixed[22][2][o],trans_SSS_fixed[22][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Mean of MPL for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_mpl_mean_$(occs[o])",false)

    temp = [trans_SSS[23][1][o],trans_SSS[23][2][o],trans_SSS[23][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[23][1][o],trans_SSS_fixed[23][2][o],trans_SSS_fixed[23][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Variance of MPL for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_mpl_var_$(occs[o])",false)

    temp = [trans_SSS[24][1][o],trans_SSS[24][2][o],trans_SSS[24][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[24][1][o],trans_SSS_fixed[24][2][o],trans_SSS_fixed[24][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Mean of MPK for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_mpk_mean_$(occs[o])",false)

    temp = [trans_SSS[25][1][o],trans_SSS[25][2][o],trans_SSS[25][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[25][1][o],trans_SSS_fixed[25][2][o],trans_SSS_fixed[25][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Variance of MPK for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_mpk_var_$(occs[o])",false)
end

occs = ["W","SP","EMP"]
occupations = ["Workers","Sole Proprietors","Employers"]
for o in 1:length(occs)
    temp = [trans_SSS[26][1][o],trans_SSS[26][2][o],trans_SSS[26][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[26][1][o],trans_SSS_fixed[26][2][o],trans_SSS_fixed[26][3][o,:]]
    if o == 3
        generate_plots(Time,"Time",temp,temp_fixed_occ,"Avg m skill for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_m_skill_avg_$(occs[o])",false, :bottomright)
    else
        generate_plots(Time,"Time",temp,temp_fixed_occ,"Avg m skill for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_m_skill_avg_$(occs[o])",false)
    end

    temp = [trans_SSS[27][1][o],trans_SSS[27][2][o],trans_SSS[27][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[27][1][o],trans_SSS_fixed[27][2][o],trans_SSS_fixed[27][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Avg w skill for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_w_skill_avg_$(occs[o])",false)

    temp = [trans_SSS[28][1][o],trans_SSS[28][2][o],trans_SSS[28][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[28][1][o],trans_SSS_fixed[28][2][o],trans_SSS_fixed[28][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Variance m skill for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_m_skill_var_$(occs[o])",false)

    temp = [trans_SSS[29][1][o],trans_SSS[29][2][o],trans_SSS[29][3][o,:]]
    temp_fixed_occ = [trans_SSS_fixed[29][1][o],trans_SSS_fixed[29][2][o],trans_SSS_fixed[29][3][o,:]]
    generate_plots(Time,"Time",temp,temp_fixed_occ,"Variance w skill for $(occupations[o])",LOCAL_DIR_PRODUCTIVITY,"time_w_skill_var_$(occs[o])",false)
end
#
# #mean of m skill fixed
# trans_SSS_fixed[26][3][1:3,1]
# #variance of m skill fixed
# trans_SSS_fixed[28][3][1:3,1]
# #mean of w skill fixed
# trans_SSS_fixed[27][3][1:3,1]
# #variance of w skill fixed - WHY IT IS NOT EQUAL TO EACH OTHER???????
# trans_SSS_fixed[29][3][1:3,1]
#
#
# #m = u_i,alpha_m_i,zeta_i
# m_skill_distr = [permutedims(sum(trans_SSS_fixed[1][3][5][7][occ][1,:,:,:,:,:], dims=[1,5])[1,:,:,:,1], [1,3,2]) for occ in 1:3]
# sum.(m_skill_distr)
# m_skill_distr./sum.(m_skill_distr)
# #w = u_i,alpha_m_i,alpha_w_i
# w_skill_distr = [sum(trans_SSS_fixed[1][3][5][7][occ][1,:,:,:,:,:], dims=[1,3])[1,:,1,:,:] for occ in 1:3]
# sum.(w_skill_distr)
# w_skill_distr./sum.(w_skill_distr)
#
# denumerator = [sum(trans_SSS_fixed[1][3][5][7][occ][1,:,:,:,:,:]) for occ=1:3]
#
# avg_w_skill = [ sum(w_skill_distr[occ].*trans_SSS_fixed[1][1][4][2]) for occ in 1:3]./denumerator
# var_w_skill_1 = zeros(3)
# var_w_skill_2 = zeros(3)
# for occ = 1:3
#     var_w_skill_1[occ] = sum((w_skill_distr[occ]./denumerator[occ]).*(trans_SSS_fixed[1][1][4][2] .- avg_w_skill[occ]).^2)
#     var_w_skill_2[occ] = sum(w_skill_distr[occ].*(trans_SSS_fixed[1][1][4][2] .- avg_w_skill[occ].*denumerator[occ]).^2)/denumerator[occ]
# end
# (var_w_skill_1, var_w_skill_2)
#
#
#
# avg_m_skill_1 = [ sum(m_skill_distr[occ].*trans_SSS_fixed[1][1][4][1]) for occ in 1:3]./denumerator
# avg_m_skill_2 = [ sum((m_skill_distr[occ]./denumerator[occ]).*trans_SSS_fixed[1][1][4][1]) for occ in 1:3]
# (avg_m_skill_1,avg_m_skill_2)
# var_m_skill_1 = zeros(3)
# var_m_skill_2 = zeros(3)
# for occ = 1:3
#     var_m_skill_1[occ] = sum((m_skill_distr[occ]./denumerator[occ]).*(trans_SSS_fixed[1][1][4][1] .- avg_m_skill[occ]).^2)
#     var_m_skill_2[occ] = sum(m_skill_distr[occ].*(trans_SSS_fixed[1][1][4][1] .- avg_m_skill[occ].*denumerator[occ]).^2)/denumerator[occ]
# end
# (var_m_skill_1, var_m_skill_2)
