using JLD2
using Plots

using ProgressMeter
using SchumakerSpline
using Interpolations

function create_plot(X,XLABEL::String,Y,YLABEL::String,Y1,Y2, IS_Y_PERCENTAGE::Bool=false, OCCUPATION::String="H", TICKSFONTSIZE::Int64=9)
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    elseif OCCUPATION=="ENT"
        COLOR="brown"
    end

    NUM_YTICKS = 10
    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=length(Y),length=NUM_XTICKS))))
    if length(X) <= 25
        if length(X) <= 10
            #XTICKS = Int64.(round.(collect(range(0;stop=T,step=1))))
            XTICKS = Int64.(round.(collect(range(1;stop=length(Y),step=1))))
        else
            #XTICKS = Int64.(round.(collect(range(0;stop=length(Y),step=5))))
            XTICKS = Int64.(round.(collect(range(0;stop=length(Y),step=5))))
        end
    end
    XTICKS = [1; XTICKS[2:end]]

    YLIMS1 = minimum([Y; Y1; Y2])
    YLIMS2 = maximum([Y; Y1; Y2])
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    Y12low = minimum([Y1; Y2])
    Y12high = maximum([Y1; Y2])
    SHARElow = (Y12low-YLIMS1)/(YLIMS2-YLIMS1)
    SHAREmid = (Y12high-Y12low)/(YLIMS2-YLIMS1)
    SHAREup = (YLIMS2-Y12high)/(YLIMS2-YLIMS1)
    #display([SHAREup; SHAREmid; SHARElow])

    LOWLENGTH = Int64(round(NUM_YTICKS*SHARElow))
    YTICKSlow = [YLIMS1]
    if LOWLENGTH>1
        YTICKSlow = collect(range(YLIMS1; stop=Y12low, length=LOWLENGTH))
        STEP = abs(YTICKSlow[2]-YTICKSlow[1])
        YTICKSlow = YTICKSlow[1:end-1]
    end

    MIDLENGTH = Int64(round(NUM_YTICKS*SHAREmid))
    YTICKSmid = []
    if MIDLENGTH>1
        YTICKSmid = collect(range(Y12low; stop=Y12high, length=MIDLENGTH))
        STEP = abs(YTICKSmid[2]-YTICKSmid[1])
        YTICKSmid = YTICKSmid[2:end-1]
        YLIMS1 = Y12low-LOWLENGTH*STEP
        YTICKSlow = collect(range(Y12low; stop=YLIMS1, step=-STEP))[end:-1:1]
    end

    HIGHLENGTH = Int64(round(NUM_YTICKS*SHAREup))
    YTICKShigh = []
    if LOWLENGTH+MIDLENGTH<3
        YTICKShigh = collect(range(Y12high; stop=YLIMS2, length=HIGHLENGTH))
    else
        if HIGHLENGTH==0 && YLIMS2!=Y12high
            YTICKShigh = [YLIMS2]
        else
            YLIMS2 = Y12high+HIGHLENGTH*STEP
            YTICKShigh = collect(range(Y12high; stop=YLIMS2, step=STEP))
        end
    end
    #YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))
    YTICKS = [YTICKSlow; YTICKSmid; YTICKShigh]

    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    #=
    YPOS = mean(YTICKS[end-1:end])
    if maximum(Y[calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
        YPOS = mean(YTICKS[1:2])
    end
    =#
    if IS_Y_PERCENTAGE
        DIGITS = Int(max(2, round(log10(0.01/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    #Y1line = ones(length(Y)+1).*Y1
    Y1line = ones(length(Y)).*Y1
    #Y2line = ones(length(Y)+1).*Y2
    Y2line = ones(length(Y)).*Y2
    if length(Y) == 50
        #Y1line[10+2:end] .= NaN
        Y1line[10+1:end] .= NaN
        #Y2line[1:40] .= NaN
        Y2line[1:40-1] .= NaN
    elseif length(Y) == 10
        #Y1line[3+2:end] .= NaN
        Y1line[3+1:end] .= NaN
        #Y2line[1:7] .= NaN
        Y2line[1:8-1] .= NaN
    end
    #plt = plot(collect([0; X]), collect.([[Y1; Y],Y1line,Y2line]),
    plt = plot(collect(X), collect.([Y,Y1line,Y2line]),
                    #color=[COLOR "green" "red"],
                    color=[COLOR COLOR COLOR],
                    linestyle=[:solid :dot :dash],
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    xtickfontsize=TICKSFONTSIZE,
                    ytickfontsize=TICKSFONTSIZE,
                    xticks=XTICKS,
                    yticks = YTICKS,
                    ylims = YLIMS )
    # text annotation
    #=if Y1 != Y2
        if IS_Y_PERCENTAGE
            TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1*100;digits=DIGITS+1))%)"
        else
            TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1;digits=DIGITS+1)))"
        end
        POS = Y1+YLIMMARGIN*1.25
        if Y1 < Y2
            POS = Y1-YLIMMARGIN*1.25
            if POS < YLIMS1-YLIMMARGIN
                POS = Y1+YLIMMARGIN*1.25
            end
        end
        annotate!([X[end]], POS, text(TEXT1, COLOR, :right, 7))

        if IS_Y_PERCENTAGE
            TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2*100;digits=DIGITS+1))%)"
        else
            TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2;digits=DIGITS+1)))"
        end
        POS = Y2-YLIMMARGIN*1.25
        if POS < YLIMS1-YLIMMARGIN || Y1 < Y2
            POS = Y2+YLIMMARGIN*1.25
        end
        annotate!([X[end]], POS, text(TEXT2, COLOR, :right, 7))
    end=#
    return plt
end
function create_plot(X,XLABEL::String,Y,YLABEL::String, IS_Y_PERCENTAGE::Bool=false, OCCUPATION::String="H", IS_LABELX_SHOWN::Bool=true, NUM_YTICKS::Int64 = 10, IS_ZERO_INCLUDED::Bool=true)
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    elseif OCCUPATION=="ENT"
        COLOR="brown"
    end

    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=length(Y),length=NUM_XTICKS))))
    if length(X) <= 25
        if length(X) <= 10
            #XTICKS = Int64.(round.(collect(range(0;stop=T,step=1))))
            XTICKS = Int64.(round.(collect(range(1;stop=length(Y),step=1))))
        else
            #XTICKS = Int64.(round.(collect(range(0;stop=length(Y),step=5))))
            XTICKS = Int64.(round.(collect(range(0;stop=length(Y),step=5))))
        end
    end
    XTICKS = [1; XTICKS[2:end]]

    YLIMS1 = minimum([Y; 0.0])
    YLIMS2 = maximum([Y; 0.0])
    if !IS_ZERO_INCLUDED
        YLIMS1 = minimum(Y)
        YLIMS2 = maximum(Y)
    end
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))
    #=
    YPOS = mean(YTICKS[end-1:end])
    if maximum(Y[calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
        YPOS = mean(YTICKS[1:2])
    end
    =#
    if IS_Y_PERCENTAGE
        DIGITS = Int(max(2, round(log10(0.01/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    #plt = plot(collect([0; X]), collect.([[0.0; Y],ones(length(Y)+1).*0.0]),
    plot_mat = collect.([Y ,ones(length(Y)).*0.0])
    colors = [COLOR COLOR]
    linestyles = [:solid :dash]
    if !IS_ZERO_INCLUDED
        plot_mat = collect.([Y])
        colors = [COLOR]
        linestyles = [:solid]
    end
    plt = plot(collect(X), plot_mat,
                    #color=[COLOR "green" "red"],
                    color=colors,
                    linestyle=linestyles,
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    xticks = XTICKS,
                    yticks = YTICKS,
                    ylims = YLIMS )
    if !IS_LABELX_SHOWN
        plt = plot(collect(X), collect.([Y ,ones(length(Y)).*0.0]),
                        #color=[COLOR "green" "red"],
                        color=[COLOR COLOR],
                        linestyle=[:solid :dash],
                        legend=false,
                        ylabel=YLABEL,
                        yticks = YTICKS,
                        ylims = YLIMS )
    end
    return plt
end
function create_combined_plot(X,XLABEL::String,Ys,YLABELs,YLABEL,Y1s,Y2s, IS_Y_PERCENTAGE::Bool=false, LEGENDPOS=false, OCCUPATION=["SP","EMP"])
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

    NUM_YTICKS = 10
    YLIMS1 = minimum([minimum.(Ys); minimum.(Y1s); minimum.(Y2s)])
    YLIMS2 = maximum([maximum.(Ys); maximum.(Y1s); maximum.(Y2s)])
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=10))

    COLORS=[]
    for y_i = 1:length(Ys)
        if OCCUPATION[y_i]=="W"
            push!(COLORS,"purple")
        elseif OCCUPATION[y_i]=="SP"
            push!(COLORS,"red")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"green")
        elseif OCCUPATION[y_i]=="ENT"
            push!(COLORS,"brown")
        else
            push!(COLORS,"blue")
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
        plot!(plt,collect(X), collect.([Ys[y_i] ,Y1line,Y2line]),
                        #color=[COLORS[y_i] "green" "red"],
                        color=[COLORS[y_i] COLORS[y_i] COLORS[y_i]],
                        linestyle=[:solid :dot :dash],
                        legend=LEGENDPOS,
                        xlabel=XLABEL,
                        label=[YLABELs[y_i] "" ""],
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

TIME_PERIODS = 50#100
Time = collect(1:TIME_PERIODS)

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)

function generate_plots(X,XLABEL,Y,Y_fixed_occ,YLABEL,PATHDIR,FILENAME,IS_PERCENTAGE::Bool=false)

    plt = create_combined_plot(X,XLABEL, [Y[3], Y_fixed_occ[3]],["w occ mob", "wo occ mob"],YLABEL, [Y[1],Y_fixed_occ[1]],[Y[2],Y_fixed_occ[2]], IS_PERCENTAGE, :bottomright)
    display(plt)
    savefig(plt,"$(PATHDIR)$(country)_$(FILENAME).png")

    plt = create_combined_plot(X,XLABEL, [Y[3], Y_fixed_occ[3].+(Y[3][1]-Y_fixed_occ[3][1])],["w occ mob", "wo occ mob"],YLABEL, [Y[1],Y_fixed_occ[1]+(Y[1]-Y_fixed_occ[1])],[Y[2],Y_fixed_occ[2]+(Y[1]-Y_fixed_occ[1])], IS_PERCENTAGE, :bottomright)
    display(plt)
    savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_adjusted.png")

    if !IS_PERCENTAGE
        plt = create_plot(X,XLABEL, (Y[3].-Y_fixed_occ[3])./Y[3],YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_diff.png")

        plt = create_plot(X,XLABEL, (Y[3].-(Y[3][1]-Y_fixed_occ[3][1]).-Y_fixed_occ[3])./Y[3],YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_did.png")
    else
        plt = create_plot(X,XLABEL, (Y[3].-Y_fixed_occ[3]),YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_diff.png")

        plt = create_plot(X,XLABEL, (Y[3].-(Y[3][1]-Y_fixed_occ[3][1]).-Y_fixed_occ[3]),YLABEL, true)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_did.png")
    end
end

#Interest rate and wage
generate_plots(Time,"Time",trans_SSS[3],trans_SSS_fixed[3],"Interest Rate",LOCAL_DIR_GENERAL,"time_interest_rate",true)

generate_plots(Time,"Time",trans_SSS[4],trans_SSS_fixed[4],"Wage",LOCAL_DIR_GENERAL,"time_wage",false)

#Capital
generate_plots(Time,"Time",trans_SSS[5],trans_SSS_fixed[5],"Capital",LOCAL_DIR_GENERAL,"time_capital",false)
#Consumption
generate_plots(Time,"Time",trans_SSS[6],trans_SSS_fixed[6],"Consumption",LOCAL_DIR_GENERAL,"time_consumption",false)
#Credit     - All,SP,EMP,ENT
generate_plots(Time,"Time",[trans_SSS[7][1][1],trans_SSS[7][2][1],trans_SSS[7][3][1,:]],[trans_SSS_fixed[7][1][1],trans_SSS_fixed[7][2][1],trans_SSS_fixed[7][3][1,:]],"Credit",LOCAL_DIR_GENERAL,"time_credit",false)
#Credit-to-output
generate_plots(Time,"Time",trans_SSS[8],trans_SSS_fixed[8],"Credit-to-Output",LOCAL_DIR_GENERAL,"time_credit_to_output",false)

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
        generate_plots(Time,"Time",[trans_SSS[12][1][s,h],trans_SSS[12][2][s,h],trans_SSS[12][3][s,h,:]],[trans_SSS_fixed[12][1][s,h],trans_SSS_fixed[12][2][s,h],trans_SSS_fixed[12][3][s,h,:]],"Mean of $(choice_name)' $(stat_name)",LOCAL_DIR_INEQUALITY,"time_mean_$(stat_name)_$(choice_name_short)",true)

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
