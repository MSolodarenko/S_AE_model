using JLD2
using Plots

using ProgressMeter
using SchumakerSpline
using BasicInterpolators: LinearInterpolator

include("Functions/profit.jl")

country = "Italy"
TIME_PERIODS = 50#100
LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_fixed_occ/$(country)/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_fixed_occ\\$(country)\\"
end
mkpath(LOCAL_DIR)

@load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2

approx_object = ss_1[4]
global_approx_params = ss_1[1][47]
model_params = ss_1[5]

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)
open("$(LOCAL_DIR_GENERAL)table_stats.txt", "w") do f
    write(f, "hline \n")
    write(f, "Statistics                     & Calibrated Economy    & Developed Economy                           \n")
    write(f, "hline  \n")
    write(f, "Collateral constraint (lambda)& $(round(ss_1[5][1]; digits=2)) & $(round(ss_2[5][1]; digits=2))             \n")
    write(f, "Capital to GDP                 & $(round(ss_1[1][12]*100; digits=2))% & $(round(ss_2[1][12]*100; digits=2))%                         \n")
    write(f, "Credit to GDP                  & $(round(ss_1[1][13]*100; digits=2))% & $(round(ss_2[1][13]*100; digits=2))%                         \n")
    write(f, "multicolumn{3}{l}{textit{Occupational structure}}                        \n")
    write(f, "Workers                       & $(round(ss_1[1][14]*100; digits=2))% & $(round(ss_2[1][14]*100; digits=2))%                         \n")
    write(f, "Sole Proprietors               & $(round(ss_1[1][15]*100; digits=2))% & $(round(ss_2[1][15]*100; digits=2))%                         \n")
    write(f, "Employers                      & $(round(ss_1[1][16]*100; digits=2))% & $(round(ss_2[1][16]*100; digits=2))%                         \n")
    write(f, "multicolumn{3}{l}{textit{Income inequality: Between occupations}}        \n")
    write(f, "Variance of log-earnings       & $(round(ss_1[1][19]; digits=2)) & $(round(ss_2[1][19]; digits=2))                            \n")
    write(f, "multicolumn{3}{l}{textit{Income inequality: Within occupations}}         \n")
    write(f, "Gini for workers' income       & $(round(ss_1[1][20]; digits=2)) & $(round(ss_2[1][20]; digits=2))                            \n")
    write(f, "Gini for entrepreneurs' income & $(round(ss_1[1][21]; digits=2)) & $(round(ss_2[1][21]; digits=2))      					   \n")
    write(f, "hline    \n")

end

T = trans_res[1]
number_asset_grid = trans_res[9]
asset_grid = trans_res[10]
capital_s_distr_s = trans_res[5][7]
policy_s = trans_res[5][8]

lambda_1 = ss_1[5][1]
lambda_2 = ss_2[5][1]
lambda_s = trans_res[2]

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

function create_plot(X,XLABEL::String,Y,YLABEL::String, IS_Y_PERCENTAGE::Bool=false, OCCUPATION::String="H")
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
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=10))
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
    plt = plot(collect(X), collect.([Y ,ones(length(Y)).*0.0]),
                    #color=[COLOR "green" "red"],
                    color=[COLOR COLOR],
                    linestyle=[:solid :dash],
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    xticks = XTICKS,
                    yticks = YTICKS,
                    ylims = YLIMS )

    return plt
end

function create_combined_plot(X,XLABEL::String,Ys,YLABELs,YLABEL,Y1s,Y2s, IS_Y_PERCENTAGE::Bool=false, OCCUPATION=["SP","EMP"], LEGENDPOS=false)
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
#=
plt = create_plot(collect(1:T),"Time (Years)", lambda_s,"λ", lambda_1, lambda_2)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_lambdas.png")
=#
occ_W = ss_1[1][14]
occ_SP = ss_1[1][15]
occ_EMP = ss_1[1][16]
occ_ENT = occ_SP+occ_EMP

r_1 = ss_1[2]
r_2 = ss_2[2]
r_s = trans_res[3]
if length(r_s) != T
    r_s = r_s[1:T]
end
plt = create_plot(collect(1:T),"Time (Years)", r_s,""#="Interest rate"=#, r_1, r_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_interest_rates.png")

w_1 = ss_1[3]
w_2 = ss_2[3]
w_s = trans_res[4]
plt = create_plot(collect(1:T),"Time (Years)", w_s,""#="Wage"=#, w_1, w_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_wages.png")

# initialisation
Output_1 = sum([sum(ss_1[1][32][occ].*ss_1[1][5][occ]) for occ=1:3])
Output_2 = sum([sum(ss_2[1][32][occ].*ss_2[1][5][occ]) for occ=1:3])
Output_s = zeros(T)
Credit_to_Output_1 = ss_1[1][13]
Credit_to_Output_2 = ss_2[1][13]
Credit_to_Output_s = zeros(T)
Capital_1 = sum([sum(ss_1[1][3].*ss_1[1][5][occ]) for occ=1:3])
Capital_2 = sum([sum(ss_2[1][3].*ss_2[1][5][occ]) for occ=1:3])
Capital_s = zeros(T)
# All, SP, EMP, ENT
Credit = zeros(4,2)
Credit_s = zeros(4,T)
var_Credit = zeros(4,2)
var_Credit_s = zeros(4,T)
Income_1 = sum([sum((ss_1[1][23][occ] .- ones(size(ss_1[1][5][occ])).*ss_1[1][3]).*ss_1[1][5][occ]) for occ=1:3])
Income_2 = sum([sum((ss_2[1][23][occ] .- ones(size(ss_2[1][5][occ])).*ss_2[1][3]).*ss_2[1][5][occ]) for occ=1:3])
Income_s = zeros(T)
Consumption_1 = sum([sum((ss_1[1][23][occ] .- ss_1[1][4][occ]).*ss_1[1][5][occ]) for occ=1:3])
Consumption_2 = sum([sum((ss_2[1][23][occ] .- ss_2[1][4][occ]).*ss_2[1][5][occ]) for occ=1:3])
Consumption_s = zeros(T)

# income, earnings, wealth, consumption
# All, W, SP, EMP, ENT
# ss_1 ss_2
means = zeros(4,5,2)
ginis = zeros(4,5,2)
avglogs = zeros(4,5,2)
varlogs = zeros(4,5,2)
avgs = zeros(4,5,2)
vars = zeros(4,5,2)
coef_of_variation = zeros(4,5,2)

means_s = zeros(4,5,T)
ginis_s = zeros(4,5,T)
avglogs_s = zeros(4,5,T)
varlogs_s = zeros(4,5,T)
avgs_s = zeros(4,5,T)
vars_s = zeros(4,5,T)

##########################
calc_mean_quantile = true#false#
##########################
# 5 quantiles
# All, W, SP, EMP
# income, earnings, wealth, consumption
quantile_means = zeros(5,4,4,2)

quantile_means_s = zeros(5,4,4,T)

function quantile_mean(stat, distr)
    distr ./= sum(distr)

    stat_grid = vec(stat)
    stat_distr = vec(distr)
    ids = findall(x->x>0.0, stat_distr)
    ud_stat_grid = stat_grid[ids]
    ud_stat_distr = stat_distr[ids]
    try
        stat_grid = exp.(collect(range(log(1); stop=log(maximum(ud_stat_grid)+1-minimum(ud_stat_grid)), length=length(stat[:,1,1,1,1,1])))).+minimum(ud_stat_grid).-1
    catch e
        temp111 = findmin(ud_stat_grid)
        display([ud_stat_grid[temp111[2]], ud_stat_distr[temp111[2]]])
        id111 = findall(x->x<=0.0, stat_grid)
        display(sum(stat_distr[id111]))
        throw(e)
    end
    stat_distr = zeros(length(stat_grid))
    Threads.@threads for i in 1:length(stat_grid)
        iters1 = []
        iters2 = []
        if i == 1
            iters2 = findall(x->stat_grid[i]<=x<=stat_grid[i+1],ud_stat_grid)
        elseif i == length(stat_grid)
            iters1 = findall(x->stat_grid[i-1]<=x<=stat_grid[i],ud_stat_grid)
        else
            iters1 = findall(x->stat_grid[i-1]<x<=stat_grid[i],ud_stat_grid)
            iters2 = findall(x->stat_grid[i]<x<stat_grid[i+1],ud_stat_grid)
        end
        if !isempty(iters1)
            cand1 = ud_stat_grid[iters1]
            lottery_prob1 = (cand1 .- stat_grid[i-1])./(stat_grid[i]-stat_grid[i-1])
            cand_distr1 = ud_stat_distr[iters1]
            stat_distr[i] += sum(cand_distr1.*lottery_prob1)
        end
        if !isempty(iters2)
            cand2 = ud_stat_grid[iters2]
            lottery_prob2 = (stat_grid[i+1] .- cand2)./(stat_grid[i+1]-stat_grid[i])
            cand_distr2 = ud_stat_distr[iters2]
            stat_distr[i] += sum(cand_distr2.*lottery_prob2)
        end

    end
    cdf = cumsum(stat_distr)

    idx = unique(z -> cdf[z], 1:length(cdf))
    cdf = cdf[idx]
    stat_grid = stat_grid[idx]
    stat_distr = stat_distr[idx]

    p0 = plot(stat_distr)
    max_a_i = findlast(x->x<1.0-1e-4, cdf)
    try
        cdf = cdf[1:max_a_i]./cdf[max_a_i]
    catch e1
        try
            max_a_i = length(cdf)-1
            cdf = cdf[1:max_a_i]./cdf[max_a_i]
        catch e
            display(cdf)
            throw(e)
        end
    end
    stat_grid = stat_grid[1:max_a_i]
    stat_distr = stat_distr[1:max_a_i]

    idx = unique(z -> cdf[z], 1:length(cdf))
    cdf = cdf[idx]
    stat_grid = stat_grid[idx]
    stat_distr = stat_distr[idx]

    p0 = plot(stat_distr)
    max_a_i = findlast(x->x<1.0-1e-4, cdf)#length(cdf)#

    p1 = plot(cdf[1:max_a_i]./cdf[max_a_i])

    cdf_1 = []
    try
        #cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i],1:max_a_i; extrapolation=(Constant,Constant))
        cdf_1 = linear_interpolation(collect(cdf[1:max_a_i]./cdf[max_a_i]), collect(1:max_a_i), extrapolation_bc = Flat());
        p2 = plot(collect(range(0.0;stop=1.0,length=120)),cdf_1.(collect(range(0.0;stop=1.0,length=120))))
    catch e1
        try
            #cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i],1:max_a_i; extrapolation=(Constant,Constant))
            cdf_1 = linear_interpolation([0.0; collect(cdf[1:max_a_i]./cdf[max_a_i])], [1; collect(1:max_a_i)], extrapolation_bc = Flat());
            p2 = plot(collect(range(0.0;stop=1.0,length=120)),cdf_1.(collect(range(0.0;stop=1.0,length=120))))
        catch e2
            display(plot(p0,p1,legend=false))
            throw(e2)
        end
    end

    quantiles = collect(range(0.0;stop=1.0,length=6))
    #qwb = evaluate.(cdf_1,quantiles)    #quantile_wealth_bounds
    qwb = cdf_1.(quantiles)

    quantile_mean = zeros(5)
    for q = 1:5
        if qwb[q] > qwb[q+1]
            display(plot(p0,p1,p2,layout=(1,3),legend=false))
            display(qwb)
            throw(error)
        end
        _q = Int32(max(1,min(floor(qwb[q]),max_a_i)))
        q_ = Int32(max(1,min(ceil(qwb[q]),max_a_i)))
        _q1 = Int32(max(1,min(floor(qwb[q+1]),max_a_i)))
        q1_ = Int32(max(1,min(ceil(qwb[q+1]),max_a_i)))

        m_w_down_down = sum(stat_grid[_q:_q1] .* stat_distr[_q:_q1] / sum(stat_distr[_q:_q1]) )
        m_w_down_up = sum(stat_grid[_q:q1_] .* stat_distr[_q:q1_] / sum(stat_distr[_q:q1_]) )
        m_w_up_down = sum(stat_grid[q_:_q1] .* stat_distr[q_:_q1] / sum(stat_distr[q_:_q1]) )
        m_w_up_up = sum(stat_grid[q_:q1_] .* stat_distr[q_:q1_] / sum(stat_distr[q_:q1_]) )

        quantile_mean[q] = (m_w_down_down+m_w_down_up+m_w_up_down+m_w_up_up)/4
    end
    return quantile_mean
end

# Distribution of output to different channels of earnigs/income to occupations
share_W_earnings_in_output = zeros(2)
share_SP_earnings_in_output = zeros(2)
share_EMP_earnings_in_output = zeros(2)
share_W_capital_income_in_output = zeros(2)
share_SP_capital_income_in_output = zeros(2)
share_EMP_capital_income_in_output = zeros(2)

share_W_earnings_in_output_s = zeros(T)
share_SP_earnings_in_output_s = zeros(T)
share_EMP_earnings_in_output_s = zeros(T)
share_W_capital_income_in_output_s = zeros(T)
share_SP_capital_income_in_output_s = zeros(T)
share_EMP_capital_income_in_output_s = zeros(T)

# loop for _1 and _2
p = Progress(2*(1+5*(1+4*3)+4*4+1), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
if !calc_mean_quantile
    p = Progress(2*(1+5*(1+4*3)+1), dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=25)
end
Threads.@threads for i = 1:2
    ss = ss_1
    Output = Output_1
    if i == 2
        ss = ss_2
        Output = Output_2
    end

    lambda = ss[5][1]
    density_distr = ss[1][5]
    output = ss[1][32]
    capital_d = ss[1][26]
    labour_d = ss[1][29]
    credit = ss[1][27]

    eta = ss[5][5]
    theta = ss[5][6]

    z_m_nodes = ss[4][1]
    z_w_nodes = ss[4][2]

    unlimited_capital_choice = Array{Any}(undef,3)
    Threads.@threads for occ = 1:3
        unlimited_capital_choice[occ] = capital_d[occ].<(lambda.*(ss[1][3].*ones(size(capital_d[occ]))))
    end
    next!(p)

    Threads.@threads for h = 1:5 # [2] = All, W, SP, EMP, ENT
        #if h == 1 #All
        list_of_occs = [1,2,3]
        if h == 2 #W
            list_of_occs = [1]
        elseif h == 3 #SP
            list_of_occs = [2]
        elseif h == 4 #EMP
            list_of_occs = [3]
        elseif h == 5 #ENT
            list_of_occs = [2,3]
        end

        if h==1
            Credit[h,i] = sum([sum(density_distr[occ] .* credit[occ]) for occ=list_of_occs])
            var_Credit[h,i] = sum([sum(density_distr[occ] .* (credit[occ] .- Credit[h,i]).^2)  for occ=list_of_occs])/sum([sum(density_distr[occ]) for occ=list_of_occs])
            Credit[h,i] /= sum([sum(density_distr[occ]) for occ=list_of_occs])
        elseif h!=2
            Credit[h-1,i] = sum([sum(density_distr[occ] .* credit[occ]) for occ=list_of_occs])
            var_Credit[h-1,i] = sum([sum(density_distr[occ] .* (credit[occ] .- Credit[h-1,i]).^2)  for occ=list_of_occs])/sum([sum(density_distr[occ]) for occ=list_of_occs])
            Credit[h-1,i] /= sum([sum(density_distr[occ]) for occ=list_of_occs])
        end

        next!(p)

        Threads.@threads for s = 1:4 # [1] = income, earnings, wealth, consumption
            stat_distr = Array{Any}(undef,3)
            for occ = 1:3
                # if s == 1
                stat_distr[occ] = ss[1][23][occ] .- ones(size(ss[1][5][occ])).*ss[1][3]
                if s == 2
                    stat_distr[occ] = ss[1][24][occ]
                elseif s == 3
                    stat_distr[occ] = ones(size(ss[1][5][occ])).*ss[1][3]
                elseif s == 4
                    stat_distr[occ] = ss[1][23][occ] .- ss[1][4][occ]
                end
            end

            # calculate mean
            means[s,h,i] = sum([sum(stat_distr[occ] .* density_distr[occ]) for occ=list_of_occs]))/sum([sum(density_distr[occ]) for occ=list_of_occs])
            next!(p)

            # calculate gini coefficent
            density_distr_choice = []
            stat_choice = []
            if h==1
                density_distr_choice = ([density_distr[1]; density_distr[2]; density_distr[3]])./sum([sum(density_distr[occ]) for occ=1:3])
                stat_choice = [stat_distr[1]; stat_distr[2]; stat_distr[3]]
            elseif h==2
                density_distr_choice = (density_distr[1])./sum(density_distr[1])
                stat_choice = stat_distr[1]
            elseif h==3
                density_distr_choice = (density_distr[2])./sum(density_distr[2])
                stat_choice = stat_distr[2]
            elseif h==4
                density_distr_choice = (density_distr[3])./sum(density_distr[3])
                stat_choice = stat_distr[3]
            # elseif 2<=h && h>=4
            #     density_distr_choice = ([density_distr[occ] for occ=list_of_occs])./sum(sum([density_distr[occ] for occ=list_of_occs]))
            #     stat_choice = [stat_distr[occ] for occ=list_of_occs]
            elseif h==5
                density_distr_choice = ([density_distr[2]; density_distr[3]])./sum([sum(density_distr[occ]) for occ=2:3])
                stat_choice = [stat_distr[2]; stat_distr[3]]
            end
            density_distr_choice_vec = vec(density_distr_choice)
            stat_choice_vec = vec(stat_choice)
            index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
            density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
            stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]
            try
                ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))
            catch e
                ginis[s,h,i] = NaN
            end
            next!(p)

            avgs[s,h,i] = sum([sum(density_distr[occ].*max.(1e-12,stat_distr[occ])) for occ = list_of_occs])
            vars[s,h,i] = sum([sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).- avgs[s,h,i]).^2) for occ = list_of_occs])/sum([sum(density_distr[occ]) for occ = list_of_occs])
            avgs[s,h,i] /= sum([sum(density_distr[occ]) for occ = list_of_occs])
            coef_of_variation[s,h,i] = sqrt( vars[s,h,i] )/avgs[s,h,i]
            next!(p)

            if h!=5 && calc_mean_quantile
                try
                    quantile_means[:,h,s,i] .= quantile_mean(stat_choice, density_distr_choice)
                catch e
                    quantile_means[:,h,s,i] .= NaN
                end
                next!(p)
            end
        end
    end
    agg_c_e = ss[5][7]*sum(density_distr[3])
    share_W_earnings_in_output[i] = sum(ss[1][24][1] .* density_distr[1])/(Output-agg_c_e)
    share_SP_earnings_in_output[i] = sum(ss[1][24][2] .* density_distr[2])/(Output-agg_c_e)
    share_EMP_earnings_in_output[i] = sum(ss[1][24][3] .* density_distr[3])/(Output-agg_c_e)
    share_W_capital_income_in_output[i] = (ss[5][3]+ss[1][44])*sum(ones(size(ss[1][5][1])).*ss[1][3] .* density_distr[1])/(Output-agg_c_e)
    share_SP_capital_income_in_output[i] = (ss[5][3]+ss[1][44])*sum(ones(size(ss[1][5][2])).*ss[1][3] .* density_distr[2])/(Output-agg_c_e)
    share_EMP_capital_income_in_output[i] = (ss[5][3]+ss[1][44])*sum(ones(size(ss[1][5][3])).*ss[1][3] .* density_distr[3])/(Output-agg_c_e)
    next!(p)
end
println()
# print coeficient of variation
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
    LOCAL_COLORS = ["H","W","SP","EMP","ENT"]

    for h = 1:5
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        println("Coef.Variation of $(stat_name) for $(choice_name) - $(coef_of_variation[s,h,1]) - $(coef_of_variation[s,h,2])")

    end
end

#main computation loop
p = Progress(T*(1+5*(1+4*(2)+3)+1), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
if !calc_mean_quantile
    p = Progress(T*(1+5*(1+4*(2))+1), dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=25)
end
Threads.@threads for t=1:T
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

    beta = model_params[2]
    p_alpha = model_params[13]

    z_m_nodes = approx_object[1]
    z_w_nodes = approx_object[2]
    P_zeta = approx_object[3]
    P_u = approx_object[4]
    stat_P_u = approx_object[5]
    P_alpha = approx_object[6]
    number_u_nodes = approx_object[7]

    delta = model_params[3]
    gamma = model_params[4]
    eta = model_params[5]
    theta = model_params[6]
    c_e = model_params[7]

    z_m_nodes = approx_object[1]
    z_w_nodes = approx_object[2]
    number_u_nodes = approx_object[7]

    occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input = compute_income_profile_fixed_occ(asset_grid,number_asset_grid,r_s[t], w_s[t], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

    Output_s[t] = sum([sum(output[occ] .* capital_s_distr_s[occ][t,:,:,:,:,:]) for occ=1:3])

    Capital_s[t]= sum([sum(asset_grid .* capital_s_distr_s[occ][t,:,:,:,:,:]) for occ=1:3])

    agg_credit = sum([sum(credit[occ] .* capital_s_distr_s[occ][t,:,:,:,:,:]) for occ=1:3])

    Credit_to_Output_s[t] = agg_credit/Output_s[t]

    Income_s[t] = sum([sum((income[occ].- ones(size(capital_s_distr_s[occ][t,:,:,:,:,:])).*asset_grid) .* capital_s_distr_s[occ][t,:,:,:,:,:]) for occ=1:3])

    Consumption_s[t] = sum([sum((income[occ].-policy_s[occ][t,:,:,:,:,:]) .* capital_s_distr_s[occ][t,:,:,:,:,:]) for occ=1:3] )

    # income, earnings, wealth, consumption
    # All, W, SP, EMP, ENT
    # Time
    lambda = lambda_s[t]
    density_distr = Array{Any}(undef,3)
    Threads.@threads for occ = 1:3
        density_distr[occ] = capital_s_distr_s[occ][t,:,:,:,:,:]
    end
    output = output
    capital_d = capital_d
    labour_d = labour_d

    eta = ss_1[5][5]
    theta = ss_1[5][6]

    z_m_nodes = ss_1[4][1]
    z_w_nodes = ss_1[4][2]

    unlimited_capital_choice = Array{Any}(undef,3)
    Threads.@threads for occ = 1:3
        unlimited_capital_choice[occ] = capital_d[occ].<(lambda.*(asset_grid.*ones(size(capital_d[occ]))))
    end
    next!(p)

    Threads.@threads for h = 1:5 # [2] = All, W, SP, EMP, ENT
        #if h == 1 #All
        list_of_occs = [1,2,3]
        if h == 2 #W
            list_of_occs = [1]
        elseif h == 3 #SP
            list_of_occs = [2]
        elseif h == 4 #EMP
            list_of_occs = [3]
        elseif h == 5 #ENT
            list_of_occs = [2,3]
        end

        if h==1
            Credit[h,t] = sum([sum(density_distr[occ] .* credit[occ]) for occ=list_of_occs])
            var_Credit[h,t] = sum([sum(density_distr[occ] .* (credit[occ] .- Credit[h,t]).^2)  for occ=list_of_occs])/sum([sum(density_distr[occ]) for occ=list_of_occs])
            Credit[h,t] /= sum([sum(density_distr[occ]) for occ=list_of_occs])
        elseif h!=2
            Credit[h-1,t] = sum([sum(density_distr[occ] .* credit[occ]) for occ=list_of_occs])
            var_Credit[h-1,t] = sum([sum(density_distr[occ] .* (credit[occ] .- Credit[h-1,t]).^2)  for occ=list_of_occs])/sum([sum(density_distr[occ]) for occ=list_of_occs])
            Credit[h-1,t] /= sum([sum(density_distr[occ]) for occ=list_of_occs])
        end
        next!(p)

        Threads.@threads for s = 1:4 # [1] = income, earnings, wealth, consumption
            stat_distr = Array{Any}(undef,3)
            for occ = 1:3
                # if s == 1
                stat_distr[occ] = income[occ] .- ones(size(density_distr[occ])).*asset_grid
                if s == 2
                    stat_distr[occ] = earnings[occ]
                elseif s == 3
                    stat_distr[occ] = ones(size(density_distr[occ])).*asset_grid
                elseif s == 4
                    stat_distr[occ] = income[occ] .- policy_s[occ][t,:,:,:,:,:]
                end
            end

            # calculate mean
            means[s,h,t] = sum([sum(stat_distr[occ] .* density_distr[occ]) for occ=list_of_occs])/sum([sum(density_distr[occ]) for occ=list_of_occs])
            next!(p)

            # calculate gini coefficent
            density_distr_choice = []
            stat_choice = []
            if h==1
                density_distr_choice = ([density_distr[1]; density_distr[2]; density_distr[3]])./sum([sum(density_distr[occ]) for occ=1:3])
                stat_choice = [stat_distr[1]; stat_distr[2]; stat_distr[3]]
            elseif h==2
                density_distr_choice = (density_distr[1])./sum(density_distr[1])
                stat_choice = stat_distr[1]
            elseif h==3
                density_distr_choice = (density_distr[2])./sum(density_distr[2])
                stat_choice = stat_distr[2]
            elseif h==4
                density_distr_choice = (density_distr[3])./sum(density_distr[3])
                stat_choice = stat_distr[3]
            # elseif 2<=h && h>=4
            #     density_distr_choice = ([density_distr[occ] for occ=list_of_occs])./sum(sum([density_distr[occ] for occ=list_of_occs]))
            #     stat_choice = [stat_distr[occ] for occ=list_of_occs]
            elseif h==5
                density_distr_choice = ([density_distr[2]; density_distr[3]])./sum([sum(density_distr[occ]) for occ=2:3])
                stat_choice = [stat_distr[2]; stat_distr[3]]
            end
            density_distr_choice_vec = vec(density_distr_choice)
            stat_choice_vec = vec(stat_choice)
            index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
            density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
            stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]
            try
                ginis_s[s,h,t] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))
            catch e
                ginis_s[s,h,t] = NaN
            end
            next!(p)

            if h!=5 && calc_mean_quantile
                try
                    quantile_means[:,h,s,t] .= quantile_mean(stat_choice, density_distr_choice)
                catch e
                    quantile_means[:,h,s,t] .= NaN
                end
                next!(p)
            end
        end
    end

    agg_c_e = c_e*sum(density_distr[3])
    share_W_earnings_in_output_s[t] = sum(earnings[1] .* density_distr[1])/(Output_s[t]-agg_c_e)
    share_SP_earnings_in_output_s[t] = sum(earnings[2] .* density_distr[2])/(Output_s[t]-agg_c_e)
    share_EMP_earnings_in_output_s[t] = sum(earnings[3] .* density_distr[3])/(Output_s[t]-agg_c_e)
    share_W_capital_income_in_output_s[t] = (delta+r_s[t])*sum(ones(size(density_distr[1])).*asset_grid .* density_distr[1])/(Output_s[t]-agg_c_e)
    share_SP_capital_income_in_output_s[t] = (delta+r_s[t])*sum(ones(size(density_distr[2])).*asset_grid .* density_distr[2])/(Output_s[t]-agg_c_e)
    share_EMP_capital_income_in_output_s[t] = (delta+r_s[t])*sum(ones(size(density_distr[3])).*asset_grid .* density_distr[3])/(Output_s[t]-agg_c_e)

    next!(p)
end

Output_growth_rate_s = zeros(T)
Capital_growth_rate_s = zeros(T)
Income_growth_rate_s = zeros(T)
Consumption_growth_rate_s = zeros(T)

#create plots

plt = create_plot(collect(1:T),"Time (Years)", Output_s,""#="Output"=#, Output_1, Output_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_outputs.png")

plt = create_plot(collect(1:T),"Time (Years)", Capital_s./Output_s,""#="Capital/Output"=#, Capital_1/Output_1, Capital_2/Output_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_capital_to_outputs.png")
plt = create_plot(collect(1:T),"Time (Years)", Credit_s[1,:],""#="Credit"=#, Credit[1,1], Credit[1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_credits.png")
plt = create_plot(collect(1:T),"Time (Years)", Credit_to_Output_s,""#="Credit/Output"=#, Credit_to_Output_1, Credit_to_Output_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_credit_to_outputs.png")

#plt = create_plot(collect(1:T),"Time (Years)", Income_s,"Income", Income_1, Income_2, false)
plt = create_plot(collect(1:T),"Time (Years)", means_s[1,1,:],""#="Income"=#, means[1,1,1], means[1,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_incomes.png")

plt = create_plot(collect(1:T),"Time (Years)", means_s[4,1,:],""#="Consumption"=#, means[4,1,1], means[4,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_consumptions.png")
plt = create_plot(collect(1:T),"Time (Years)", means_s[2,1,:],""#="Earnings"=#, means[2,1,1], means[2,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", means_s[3,1,:],""#="Wealth"=#, means[3,1,1], means[3,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_wealths.png")

#Inequality
LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)

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
    LOCAL_COLORS = ["H","W","SP","EMP","ENT"]

    for h = 1#:5
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        # # calculate mean
        # #means[s,h,i]
        # plt = create_plot(collect(1:T),"Time (Years)", means_s[s,h,:],"Mean of $(CHOICE_NAMES[h])' $stat_name", means[s,h,1], means[s,h,2], false,LOCAL_COLORS[h])
        # display(plt)
        # savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_mean_$(choice_name)_$(stat_name).png")

        # # calculate gini coefficent
        # #ginis[s,h,i]
        # plt = create_plot(collect(1:T),"Time (Years)", ginis_s[s,h,:],""#="Gini of $(CHOICE_NAMES[h])' $stat_name"=#, ginis[s,h,1], ginis[s,h,2], false,LOCAL_COLORS[h])
        # display(plt)
        # savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_$(choice_name)_$(stat_name).png")

    end

    LEGENDPOS = false
    if s==4
        LEGENDPOS = :topright
    end
    # calculate mean
    #means[s,h,i]
    plt = create_combined_plot(collect(1:T),"Time (Years)",[means_s[s,2,:],means_s[s,3,:],means_s[s,4,:]],["W","SP","EMP"],""#="Mean $stat_name"=#,[means[s,2,1],means[s,3,1],means[s,4,1]],[means[s,2,2],means[s,3,2],means[s,4,2]], false, ["W","SP","EMP"], LEGENDPOS)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_mean_occupations_$(stat_name).png")

    LEGENDPOS = false
    if s==4
        LEGENDPOS = :right
    end
    # calculate gini coefficent
    #ginis[s,h,i]
    plt = create_combined_plot(collect(1:T),"Time (Years)",[ginis_s[s,2,:],ginis_s[s,3,:],ginis_s[s,4,:]],["W","SP","EMP"],""#="Gini $stat_name"=#,[ginis[s,2,1],ginis[s,3,1],ginis[s,4,1]],[ginis[s,2,2],ginis[s,3,2],ginis[s,4,2]], false, ["W","SP","EMP"], LEGENDPOS)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_occupations_$(stat_name).png")

end

if calc_mean_quantile
    LABELS=["1st","2nd","3rd","4th","5th"]
    stat_names = ["Income", "Earnings", "Wealth", "Consumption"]
    for stat = 1:4#[1,4]#
        stat_name = stat_names[stat]
        qm = quantile_means[:,:,stat,:]
        qm_s = quantile_means_s[:,:,stat,:]

        LEGENDPOS = false
        plt = create_combined_plot(collect(1:T),"Time (Years)", [qm_s[1,1,:],qm_s[2,1,:],qm_s[3,1,:],qm_s[4,1,:],qm_s[5,1,:]],LABELS,""#="Mean of $(stat_name) (quantiles)"=#, [qm[1,1,1],qm[2,1,1],qm[3,1,1],qm[4,1,1],qm[5,1,1]], [qm[1,1,2],qm[2,1,2],qm[3,1,2],qm[4,1,2],qm[5,1,2]], false,["H","W","SP","EMP","ENT"], LEGENDPOS)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_quantiles.png")

        LEGENDPOS = false
        plt = create_combined_plot(collect(1:T),"Time (Years)", [qm_s[1,2,:],qm_s[2,2,:],qm_s[3,2,:],qm_s[4,2,:],qm_s[5,2,:]],LABELS,""#="Mean of Workers' $(stat_name) (quantiles)"=#, [qm[1,2,1],qm[2,2,1],qm[3,2,1],qm[4,2,1],qm[5,2,1]], [qm[1,2,2],qm[2,2,2],qm[3,2,2],qm[4,2,2],qm[5,2,2]], false,["H","W","SP","EMP","ENT"], LEGENDPOS)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_w_quantiles.png")

        LEGENDPOS = false
        plt = create_combined_plot(collect(1:T),"Time (Years)", [qm_s[1,3,:],qm_s[2,3,:],qm_s[3,3,:],qm_s[4,3,:],qm_s[5,3,:]],LABELS,""#="Mean of Sole Proprietors' $(stat_name) (quantiles)"=#, [qm[1,3,1],qm[2,3,1],qm[3,3,1],qm[4,3,1],qm[5,3,1]], [qm[1,3,2],qm[2,3,2],qm[3,3,2],qm[4,3,2],qm[5,3,2]], false,["H","W","SP","EMP","ENT"], LEGENDPOS)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_sp_quantiles.png")

        LEGENDPOS = :topright
        plt = create_combined_plot(collect(1:T),"Time (Years)", [qm_s[1,4,:],qm_s[2,4,:],qm_s[3,4,:],qm_s[4,4,:],qm_s[5,4,:]],LABELS,""#="Mean of Employers' $(stat_name) (quantiles)"=#, [qm[1,4,1],qm[2,4,1],qm[3,4,1],qm[4,4,1],qm[5,4,1]], [qm[1,4,2],qm[2,4,2],qm[3,4,2],qm[4,4,2],qm[5,4,2]], false,["H","W","SP","EMP","ENT"], LEGENDPOS)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_emp_quantiles.png")
    end
end

# productivity
LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)

plt = create_plot(collect(1:T),"Time (Years)", share_W_earnings_in_output_s,""#="Share of output as Workers' Earnings"=#, share_W_earnings_in_output[1], share_W_earnings_in_output[2], true,"W")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_W_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", share_SP_earnings_in_output_s,""#="Share of output as Sole Proprietors' Earnings"=#, share_SP_earnings_in_output[1], share_SP_earnings_in_output[2], true,"SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_SP_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", share_EMP_earnings_in_output_s,""#="Share of output as Employers' Earnings"=#, share_EMP_earnings_in_output[1], share_EMP_earnings_in_output[2], true,"EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_EMP_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", share_W_capital_income_in_output_s,""#="Share of output as Workers' Capital Income"=#, share_W_capital_income_in_output[1], share_W_capital_income_in_output[2], true,"W")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_W_capital_income.png")
plt = create_plot(collect(1:T),"Time (Years)", share_SP_capital_income_in_output_s,""#="Share of output as Sole Proprietors' Capital Income"=#, share_SP_capital_income_in_output[1], share_SP_capital_income_in_output[2], true,"SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_SP_capital_income.png")
plt = create_plot(collect(1:T),"Time (Years)", share_EMP_capital_income_in_output_s,""#="Share of output as Employers' Capital Income"=#, share_EMP_capital_income_in_output[1], share_EMP_capital_income_in_output[2], true,"EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_EMP_capital_income.png")
