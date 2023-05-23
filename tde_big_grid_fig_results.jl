using JLD2
using Plots

using ProgressMeter
using SchumakerSpline
using BasicInterpolators: LinearInterpolator

include("Functions/profit.jl")

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [69,3,3,3,6,3]#[35,3,3,3,6,3]

CODENAME = #="SS_2642"#"SS_2065"=#"SS_2092"
CODENAME = "$(CODENAME)_$(GLOBAL_APPROX_PARAMS[1])"
# parameters of the model's economy (Italy)
                    #SS_2642 #SS_2065 #SS_2092 #prev.calibration
LAMBDA =            #=1.633951#1.405096=#1.665907 #1.513028
country = "Italy"
TIME_PERIODS = 50#100
LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_big_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_big_grid\\$(country)_$(CODENAME)\\"
end
mkpath(LOCAL_DIR)

@load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2

approx_object = ss_1[4]
global_approx_params = copy(GLOBAL_APPROX_PARAMS)
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
    write(f, "multicolumn{3}{l}{textit{Consumption inequality}}                        \n")
    write(f, "Var(log(c))                    & $(round(ss_1[1][17]; digits=2)) & $(round(ss_2[1][17]; digits=2))                            \n")
    write(f, "multicolumn{3}{l}{textit{Occupational one period transition rates}}      \n")
    write(f, "W to W                         & $(round(ss_1[1][18][1,1]*100; digits=2))% & $(round(ss_2[1][18][1,1]*100; digits=2))%                         \n")
    write(f, "W to SP                        & $(round(ss_1[1][18][1,2]*100; digits=2))% & $(round(ss_2[1][18][1,2]*100; digits=2))%                          \n")
    write(f, "SP to SP                       & $(round(ss_1[1][18][2,2]*100; digits=2))% & $(round(ss_2[1][18][2,2]*100; digits=2))%                         \n")
    write(f, "SP to EMP                      & $(round(ss_1[1][18][2,3]*100; digits=2))% & $(round(ss_2[1][18][2,3]*100; digits=2))%                         \n")
    write(f, "EMP to SP                      & $(round(ss_1[1][18][3,2]*100; digits=2))% & $(round(ss_2[1][18][3,2]*100; digits=2))%                          \n")
    write(f, "EMP to EMP                     & $(round(ss_1[1][18][3,3]*100; digits=2))% & $(round(ss_2[1][18][3,3]*100; digits=2))%                         \n")
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
Output_1 = sum(ss_1[1][32].*ss_1[1][5])
Output_2 = sum(ss_2[1][32].*ss_2[1][5])
Output_s = zeros(T)
Credit_to_Output_1 = ss_1[1][13]
Credit_to_Output_2 = ss_2[1][13]
Credit_to_Output_s = zeros(T)
Capital_1 = sum(ss_1[1][3].*ss_1[1][5])
Capital_2 = sum(ss_2[1][3].*ss_2[1][5])
Capital_s = zeros(T)
# All, SP, EMP, ENT
Credit = zeros(4,2)
Credit_s = zeros(4,T)
var_Credit = zeros(4,2)
var_Credit_s = zeros(4,T)
Income_1 = sum((ss_1[1][23] .- ones(size(ss_1[1][5])).*ss_1[1][3]).*ss_1[1][5])
Income_2 = sum((ss_2[1][23] .- ones(size(ss_2[1][5])).*ss_2[1][3]).*ss_2[1][5])
Income_s = zeros(T)
Consumption_1 = sum((ss_1[1][23] .- ss_1[1][4]).*ss_1[1][5])
Consumption_2 = sum((ss_2[1][23] .- ss_2[1][4]).*ss_2[1][5])
Consumption_s = zeros(T)
#Occupation
occ_Ws_1 = ss_1[1][14]
occ_Ws_2 = ss_2[1][14]
occ_Ws_s = zeros(T)
occ_SPs_1 = ss_1[1][15]
occ_SPs_2 = ss_2[1][15]
occ_SPs_s = zeros(T)
occ_EMPs_1 = ss_1[1][16]
occ_EMPs_2 = ss_2[1][16]
occ_EMPs_s = zeros(T)
occ_ENTs_1 = occ_SPs_1+occ_EMPs_1
occ_ENTs_2 = occ_SPs_2+occ_EMPs_2
occ_ENTs_s = zeros(T)
#Shares of unconstrained ENT, SP, EMP
unlimited_capital_choice_1 = ss_1[1][26].<(ss_1[5][1].*(ss_1[1][3].*ones(size(ss_1[1][26]))))
unlimited_capital_choice_2 = ss_2[1][26].<(ss_2[5][1].*(ss_2[1][3].*ones(size(ss_2[1][26]))))
unlimited_capital_choice_s = zeros(T,size(unlimited_capital_choice_1)...)
share_ENT_unbound_1 = sum( Float64.(ss_1[1][22].!=1.0) .* ss_1[1][5] .* unlimited_capital_choice_1 )/occ_ENTs_1
share_ENT_unbound_2 = sum( Float64.(ss_2[1][22].!=1.0) .* ss_2[1][5] .* unlimited_capital_choice_2 )/occ_ENTs_2
share_ENT_unbound_s = zeros(T)
share_SP_unbound_1 = sum( Float64.(ss_1[1][22].==2.0) .* ss_1[1][5] .* unlimited_capital_choice_1 )/occ_SPs_1
share_SP_unbound_2 = sum( Float64.(ss_2[1][22].==2.0) .* ss_2[1][5] .* unlimited_capital_choice_2 )/occ_SPs_2
share_SP_unbound_s = zeros(T)
share_EMP_unbound_1 = sum( Float64.(ss_1[1][22].==3.0) .* ss_1[1][5] .* unlimited_capital_choice_1 )/occ_EMPs_1
share_EMP_unbound_2 = sum( Float64.(ss_2[1][22].==3.0) .* ss_2[1][5] .* unlimited_capital_choice_2 )/occ_EMPs_2
share_EMP_unbound_s = zeros(T)

#Variances of log-consumption and log-earnings, Gini for W and ENT income
logcs_1 = ss_1[1][17]
logcs_2 = ss_2[1][17]
logcs_s = zeros(T)
loges_1 = ss_1[1][19]
loges_2 = ss_2[1][19]
loges_s = zeros(T)
giniWs_1 = ss_1[1][20]
giniWs_2 = ss_2[1][20]
giniWs_s = zeros(T)
giniEnts_1 = ss_1[1][21]
giniEnts_2 = ss_2[1][21]
giniEnts_s = zeros(T)

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

quantile_mean_wealth = zeros(5,4,2)
quantile_mean_consumption = zeros(5,4,2)
quantile_mean_earnings = zeros(5,4,2)
quantile_mean_income = zeros(5,4,2)

quantile_mean_wealth_s = zeros(5,4,T)
quantile_mean_consumption_s = zeros(5,4,T)
quantile_mean_earnings_s = zeros(5,4,T)
quantile_mean_income_s = zeros(5,4,T)

function quantile_mean(stat, distr)
    distr ./= sum(distr)

    ud_stat_grid = vec(stat)
    ud_stat_distr = vec(distr)
    stat_grid = exp.(collect(range(log(minimum(ud_stat_grid)+1); stop=log(maximum(ud_stat_grid)+1), length=length(stat[:,1,1,1,1,1])))).-1
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
        next!(p)
    end
    # grid = collect(1:length(ud_stat_grid))
    # stat_grid = []
    # stat_distr = []
    # while !isempty(grid)
    #     i = grid[1]
    #     #iters = findall(x->(ud_stat_grid[i]-1e-5<=x<=ud_stat_grid[i]+1e-5),ud_stat_grid)
    #     iters = findall(x->(x==ud_stat_grid[i]),ud_stat_grid)
    #     if isempty(stat_grid)
    #         push!(stat_grid,ud_stat_grid[i])
    #         push!(stat_distr,sum(ud_stat_distr[iters]))
    #     elseif stat_grid[end] < ud_stat_grid[i]
    #         push!(stat_grid,ud_stat_grid[i])
    #         push!(stat_distr,sum(ud_stat_distr[iters]))
    #     else
    #         ii = findfirst(x->ud_stat_grid[i]<x,stat_grid)
    #         insert!(stat_grid,ii,ud_stat_grid[i])
    #         insert!(stat_distr,ii,sum(ud_stat_distr[iters]))
    #     end
    #     splice!(ud_stat_grid,iters)
    #     splice!(ud_stat_distr,iters)
    #
    #     grid = collect(1:length(ud_stat_grid))
    #
    #     p1 = plot(ud_stat_grid,legend=false)
    #     p2 = plot(ud_stat_distr,legend=false)
    #     p3 = plot(stat_grid,legend=false)
    #     p4 = plot(stat_distr,legend=false)
    #     display(plot(p1,p2,p3,p4,layout=(2,2),legend=false))
    #
    # end

    # unsorted_stat_grid = vec(stat)
    # unsorted_stat_distr = vec(distr)
    #
    # perm = sortperm(unsorted_stat_grid)
    #
    # dirty_stat_grid = copy(unsorted_stat_grid[perm])
    # dirty_stat_distr = copy(unsorted_stat_distr[perm])
    #
    # stat_grid = []
    # stat_distr = []
    # i = 1
    # while i < length(dirty_stat_grid)
    #     iters = findall(x->x==dirty_stat_grid[i], dirty_stat_grid)
    #     push!(stat_grid,dirty_stat_grid[i])
    #     push!(stat_distr,sum(dirty_stat_distr[iters]))
    #     i = iters[end]+1
    # end
    cdf = cumsum(stat_distr)

    idx = unique(z -> cdf[z], 1:length(cdf))
    cdf = cdf[idx]
    stat_grid = stat_grid[idx]
    stat_distr = stat_distr[idx]

    p0 = plot(stat_distr)
    max_a_i = findlast(x->x<1.0-1e-4, cdf)
    cdf = cdf[1:max_a_i]./cdf[max_a_i]
    stat_grid = stat_grid[1:max_a_i]
    stat_distr = stat_distr[1:max_a_i]

    idx = unique(z -> cdf[z], 1:length(cdf))
    cdf = cdf[idx]
    stat_grid = stat_grid[idx]
    stat_distr = stat_distr[idx]

    p0 = plot(stat_distr)
    max_a_i = length(cdf)#findlast(x->x<1.0-1e-4, cdf)

    p1 = plot(cdf[1:max_a_i]./cdf[max_a_i])

    cdf_1 = []
    try
        #cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i],1:max_a_i; extrapolation=(Constant,Constant))
        cdf_1 = LinearInterpolator(collect(cdf[1:max_a_i]./cdf[max_a_i]), collect(1:max_a_i));
        p2 = plot(collect(range(0.0;stop=1.0,length=120)),cdf_1.(collect(range(0.0;stop=1.0,length=120))))
    catch e1
        try
            #cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i],1:max_a_i; extrapolation=(Constant,Constant))
            cdf_1 = LinearInterpolator([0.0; collect(cdf[1:max_a_i]./cdf[max_a_i])], [1; collect(1:max_a_i)]);
            p2 = plot(collect(range(0.0;stop=1.0,length=120)),cdf_1.(collect(range(0.0;stop=1.0,length=120))))
        catch e2
            display(plot(p0,p1,legend=false))
            throw(e)
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

#productivity
# SP, EMP, ENT
TFPis = zeros(3,2)
TFPds = zeros(3,2)
mean_MPL = zeros(3,2)
var_MPL = zeros(3,2)
mean_MPK = zeros(3,2)
var_MPK = zeros(3,2)
share_unbound = zeros(3,2)

TFPis_s = zeros(3,T)
TFPds_s = zeros(3,T)
mean_MPL_s = zeros(3,T)
var_MPL_s = zeros(3,T)
mean_MPK_s = zeros(3,T)
var_MPK_s = zeros(3,T)
share_unbound_s = zeros(3,T)

# W, SP,EMP, ENT
avg_w_skill = zeros(4,2)
avg_m_skill = zeros(4,2)
var_w_skill = zeros(4,2)
var_m_skill = zeros(4,2)

avg_w_skill_s = zeros(4,T)
avg_m_skill_s = zeros(4,T)
var_w_skill_s = zeros(4,T)
var_m_skill_s = zeros(4,T)

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
p = Progress(2*(5*(1+4)+4*4*number_asset_grid+3+4), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
if !calc_mean_quantile
    p = Progress(2*(5*(1+4)+3+4), dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=25)
end
for i = 1:2
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

    unlimited_capital_choice = capital_d.<(lambda.*(ss[1][3].*ones(size(capital_d))))

    for h = 1:5 # [2] = All, W, SP, EMP, ENT
        #if h == 1 #All
        choice = ones(size(ss[1][22]))
        if h == 2 #W
            choice = Float64.(ss[1][22].==1.0)
        elseif h == 3 #SP
            choice = Float64.(ss[1][22].==2.0)
        elseif h == 4 #EMP
            choice = Float64.(ss[1][22].==3.0)
        elseif h == 5 #ENT
            choice = Float64.(ss[1][22].!=1.0)
        end

        if h==1
            Credit[h,i] = sum(density_distr.*choice .* credit)
            var_Credit[h,i] = sum(density_distr.*choice .* (credit .- Credit[h,i]).^2)/sum(density_distr.*choice)
            Credit[h,i] /= sum(density_distr.*choice)
        elseif h!=2
            Credit[h-1,i] = sum(density_distr.*choice .* credit)
            var_Credit[h-1,i] = sum(density_distr.*choice .* (credit .- Credit[h-1,i]).^2)/sum(density_distr.*choice)
            Credit[h-1,i] /= sum(density_distr.*choice)
        end

        next!(p)

        for s = 1:4 # [1] = income, earnings, wealth, consumption
            # if s == 1
            stat_distr = ss[1][23] .- ones(size(ss[1][5])).*ss[1][3]
            if s == 2
                stat_distr = ss[1][24]
            elseif s == 3
                stat_distr = ones(size(ss[1][5])).*ss[1][3]
            elseif s == 4
                stat_distr = ss[1][23] .- ss[1][4]
            end

            # calculate mean
            means[s,h,i] = sum(stat_distr .* density_distr.*choice)/sum(density_distr.*choice)

            # calculate gini coefficent
            density_distr_choice = (density_distr.*choice)./sum(density_distr.*choice)
            stat_choice = stat_distr.*choice
            density_distr_choice_vec = vec(density_distr_choice)
            stat_choice_vec = vec(stat_choice)
            index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
            density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
            stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]
            ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))
            #=
            # calculate variance of log-s
            avglogs[s,h,i] = sum(density_distr.*choice.*log.(max.(1e-12,stat_distr))   )
            varlogs[s,h,i] = sum(density_distr.*choice.*(log.(max.(1e-12,stat_distr)).- avglogs[s,h,i]).^2)/sum(density_distr.*choice)
            avglogs[s,h,i] /= sum(density_distr.*choice)
            if s==3
                avglogs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
                varlogs[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avglogs[s,h,i]).^2)/sum(density_distr.*choice)
                avglogs[s,h,i] /= sum(density_distr.*choice)
            end
            =#
            avgs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
            vars[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avgs[s,h,i]).^2)/sum(density_distr.*choice)
            avgs[s,h,i] /= sum(density_distr.*choice)

            coef_of_variation[s,h,i] = sqrt( vars[s,h,i] )/avgs[s,h,i]

            next!(p)
        end
        #=
        if h>2
            denumerator = sum(choice.*density_distr)
            Y = sum(choice.*density_distr.*output)/denumerator
            K = sum(choice.*density_distr.*capital_d)/denumerator
            L = sum(choice.*density_distr.*labour_d)/denumerator
            TFPis[h-2,i] = Y/(K^eta*L^theta)
            TFPds[h-2,i] = Y/K^(eta/(theta+eta))

            c_1 = capital_d.^(-1.0)
            replace!(c_1,Inf=>0.0)
            mean_MPK[h-2,i] = sum(choice.*density_distr.*(eta.*output.*c_1))/denumerator
            var_MPK[h-2,i] = sum(choice.*density_distr.*(eta.*output.*c_1 .- mean_MPK[h-2,i]*denumerator).^2)/denumerator

            l_1 = labour_d.^(-1.0)
            replace!(l_1,Inf=>0.0)
            mean_MPL[h-2,i] = sum(choice.*density_distr.*(theta.*output.*l_1))/denumerator
            var_MPL[h-2,i] = sum(choice.*density_distr.*(theta.*output.*l_1 .- mean_MPL[h-2,i]*denumerator).^2)/denumerator

            share_unbound[h-2,i] = sum(choice.*density_distr.*unlimited_capital_choice)/sum(choice.*density_distr)

            next!(p)
        end
        =#

        #=
        if h>1
            #[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            #m = u_i,alpha_m_i,zeta_i
            m_skill_distr = permutedims(sum(density_distr.*choice, dims=[1,5])[1,:,:,:,1], [1,3,2])
            #w = u_i,alpha_m_i,alpha_w_i
            w_skill_distr = sum(density_distr.*choice, dims=[1,3])[1,:,1,:,:]

            avg_m_skill[h-1,i] = sum(m_skill_distr.*z_m_nodes)/sum(density_distr.*choice)
            avg_w_skill[h-1,i] = sum(w_skill_distr.*z_w_nodes)/sum(density_distr.*choice)

            var_m_skill[h-1,i] = sum(m_skill_distr.*(z_m_nodes .- avg_m_skill[h-1,i]*sum(density_distr.*choice)).^2)/sum(density_distr.*choice)
            var_w_skill[h-1,i] = sum(w_skill_distr.*(z_w_nodes .- avg_w_skill[h-1,i]*sum(density_distr.*choice)).^2)/sum(density_distr.*choice)

            next!(p)
        end
        =#

        if h!=5 && calc_mean_quantile
            quantile_mean_wealth[:,h,i] .= quantile_mean(ones(size(ss[1][5])).*ss[1][3], density_distr.*choice)
            next!(p)
            quantile_mean_income[:,h,i] .= quantile_mean(ss[1][23] .- ones(size(ss[1][5])).*ss[1][3], density_distr.*choice)
            next!(p)
            quantile_mean_earnings[:,h,i] .= quantile_mean(ss[1][24], density_distr.*choice)
            next!(p)
            quantile_mean_consumption[:,h,i] .= quantile_mean(ss[1][23] .- ss[1][4], density_distr.*choice)
            next!(p)
        end
    end
    agg_c_e = ss[5][7]*sum(density_distr.*Float64.(ss[1][22].==3.0))
    share_W_earnings_in_output[i] = sum(ss[1][24] .* density_distr.*Float64.(ss[1][22].==1.0))/(Output-agg_c_e)
    share_SP_earnings_in_output[i] = sum(ss[1][24] .* density_distr.*Float64.(ss[1][22].==2.0))/(Output-agg_c_e)
    share_EMP_earnings_in_output[i] = sum(ss[1][24] .* density_distr.*Float64.(ss[1][22].==3.0))/(Output-agg_c_e)
    share_W_capital_income_in_output[i] = (ss[5][3]+ss[1][44])*sum(ones(size(ss[1][5])).*ss[1][3] .* density_distr.*Float64.(ss[1][22].==1.0))/(Output-agg_c_e)
    share_SP_capital_income_in_output[i] = (ss[5][3]+ss[1][44])*sum(ones(size(ss[1][5])).*ss[1][3] .* density_distr.*Float64.(ss[1][22].==2.0))/(Output-agg_c_e)
    share_EMP_capital_income_in_output[i] = (ss[5][3]+ss[1][44])*sum(ones(size(ss[1][5])).*ss[1][3] .* density_distr.*Float64.(ss[1][22].==3.0))/(Output-agg_c_e)
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
p = Progress(T*(9+5*5+4*4*number_asset_grid+3+4), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
if !calc_mean_quantile
    p = Progress(T*(9+5*5+3+4), dt=0.5,
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

    occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input = compute_income_profile(asset_grid,number_asset_grid,r_s[t], w_s[t], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

    Output_s[t] = sum(output .* capital_s_distr_s[t,:,:,:,:,:])
    next!(p)

    Capital_s[t]= sum(asset_grid .* capital_s_distr_s[t,:,:,:,:,:])
    next!(p)

    agg_credit = sum(credit .* capital_s_distr_s[t,:,:,:,:,:])

    Credit_to_Output_s[t] = agg_credit/Output_s[t]
    next!(p)

    Income_s[t] = sum((income.- ones(size(capital_s_distr_s[t,:,:,:,:,:])).*asset_grid) .* capital_s_distr_s[t,:,:,:,:,:])
    next!(p)

    Consumption_s[t] = sum((income.-policy_s[t,:,:,:,:,:]) .* capital_s_distr_s[t,:,:,:,:,:] )
    next!(p)

    #occupation
    occ_Ws_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==1.0) )
    occ_SPs_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==2.0) )
    occ_EMPs_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==3.0) )
    occ_ENTs_s[t] = occ_SPs_s[t]+occ_EMPs_s[t]
    next!(p)

    #share of unconstrained ENT, SP,EMP
    #=
    unlimited_capital_choice_s[t,:,:,:,:,:] = capital_d.<(lambda_s[t].*(asset_grid.*ones(size(capital_d))))
    share_ENT_unbound_s[t] = sum( Float64.(occ_choice.!=1.0) .* capital_s_distr_s[t,:,:,:,:,:] .* unlimited_capital_choice_s[t,:,:,:,:,:] )/occ_ENTs_s[t]
    share_SP_unbound_s[t] = sum( Float64.(occ_choice.==2.0) .* capital_s_distr_s[t,:,:,:,:,:] .* unlimited_capital_choice_s[t,:,:,:,:,:] )/occ_SPs_s[t]
    share_EMP_unbound_s[t] = sum( Float64.(occ_choice.==3.0) .* capital_s_distr_s[t,:,:,:,:,:] .* unlimited_capital_choice_s[t,:,:,:,:,:] )/occ_EMPs_s[t]
    =#next!(p)

    #Variances of log-consumption and log-earnings, Gini for W and ENT income
    #=
    logcs_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,income.-policy_s[t,:,:,:,:,:])).^2 ) .- sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,income.-policy_s[t,:,:,:,:,:])) )^2
    loges_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,earnings)).^2 ) - sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,earnings)) )^2
    =#
    #=
    density_distr_workers = ( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==1.0) )./sum(( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==1.0) ))
    density_distr_workers_vec = vec(density_distr_workers)
    income_workers = income.*Float64.(occ_choice.==1.0)
    income_workers_vec = vec(income_workers)
    index_w_non_zero = findall(x-> x>1e-5,density_distr_workers_vec)
    density_distr_workers_vec_non_zero = density_distr_workers_vec[index_w_non_zero]
    income_workers_vec_non_zero = income_workers_vec[index_w_non_zero]
    index_w = sortperm(income_workers_vec_non_zero)
    ddwvs = density_distr_workers_vec_non_zero[index_w]
    iwvs = income_workers_vec_non_zero[index_w]
    S_w_income = cumsum(ddwvs.*iwvs)./sum(ddwvs.*iwvs)
    giniWs_s[t] = 1.0 - ( S_w_income[1]*ddwvs[1] + sum((S_w_income[2:end] .+ S_w_income[1:end-1]).*ddwvs[2:end]) )

    density_distr_entrepreneurs = (Float64.(occ_choice.!=1.0) .* capital_s_distr_s[t,:,:,:,:,:])./sum(Float64.(occ_choice.!=1.0) .* capital_s_distr_s[t,:,:,:,:,:])
    density_distr_entrepreneurs_vec = vec(density_distr_entrepreneurs)#reshape(density_distr_entrepreneurs, 1, :)
    income_entrepreneurs = income.*Float64.(occ_choice.!=1.0)
    income_entrepreneurs_vec = vec(income_entrepreneurs)#reshape(income_entrepreneurs,1,:)
    index_ent_non_zero = findall(x-> x>1e-5,density_distr_entrepreneurs_vec)
    density_distr_entrepreneurs_vec_non_zero = density_distr_entrepreneurs_vec[index_ent_non_zero]
    income_entrepreneurs_vec_non_zero = income_entrepreneurs_vec[index_ent_non_zero]
    index_ent = sortperm(income_entrepreneurs_vec_non_zero)
    ddevs = density_distr_entrepreneurs_vec_non_zero[index_ent]
    ievs = income_entrepreneurs_vec_non_zero[index_ent]
    S_ent_income = cumsum(ddevs.*ievs)./sum(ddevs.*ievs)
    giniEnts_s[t] = 1.0 - ( S_ent_income[1]*ddevs[1] + sum((S_ent_income[2:end] .+ S_ent_income[1:end-1]).*ddevs[2:end]) )
    =#

    next!(p)
    # income, earnings, wealth, consumption
    # All, W, SP, EMP, ENT
    # Time
    lambda = lambda_s[t]
    density_distr = capital_s_distr_s[t,:,:,:,:,:]
    output = output
    capital_d = capital_d
    labour_d = labour_d

    eta = ss_1[5][5]
    theta = ss_1[5][6]

    z_m_nodes = ss_1[4][1]
    z_w_nodes = ss_1[4][2]

    unlimited_capital_choice = capital_d.<(lambda.*(asset_grid.*ones(size(capital_d))))

    for h = 1:5 # [2] = All, W, SP, EMP, ENT
        #if h == 1 #All
        choice = ones(size(occ_choice))
        if h == 2 #W
            choice = Float64.(occ_choice.==1.0)
        elseif h == 3 #SP
            choice = Float64.(occ_choice.==2.0)
        elseif h == 4 #EMP
            choice = Float64.(occ_choice.==3.0)
        elseif h == 5 #ENT
            choice = Float64.(occ_choice.!=1.0)
        end

        if h==1
            Credit_s[h,t] = sum(capital_s_distr_s[t,:,:,:,:,:].*choice .* credit)
            var_Credit_s[h,t] = sum(capital_s_distr_s[t,:,:,:,:,:].*choice .* (credit .- Credit_s[h,t]).^2)/sum(capital_s_distr_s[t,:,:,:,:,:].*choice)
            Credit_s[h,t] /= sum(capital_s_distr_s[t,:,:,:,:,:].*choice)
        elseif h!=2
            Credit_s[h-1,t] = sum(capital_s_distr_s[t,:,:,:,:,:].*choice .* credit)
            var_Credit_s[h-1,t] = sum(capital_s_distr_s[t,:,:,:,:,:].*choice .* (credit .- Credit_s[h-1,t]).^2)/sum(capital_s_distr_s[t,:,:,:,:,:].*choice)
            Credit_s[h-1,t] /= sum(capital_s_distr_s[t,:,:,:,:,:].*choice)
        end
        next!(p)

        for s = 1:4 # [1] = income, earnings, wealth, consumption
            # if s == 1
            stat_distr = income .- ones(size(capital_s_distr_s[t,:,:,:,:,:])).*asset_grid
            if s == 2
                stat_distr = earnings
            elseif s == 3
                stat_distr = ones(size(capital_s_distr_s[t,:,:,:,:,:])).*asset_grid
            elseif s == 4
                stat_distr = income .- policy_s[t,:,:,:,:,:]
            end

            # calculate mean
            means_s[s,h,t] = sum(stat_distr .* density_distr.*choice)/sum(density_distr.*choice)

            # calculate gini coefficent
            density_distr_choice = (density_distr.*choice)./sum(density_distr.*choice)
            stat_choice = stat_distr.*choice
            density_distr_choice_vec = vec(density_distr_choice)
            stat_choice_vec = vec(stat_choice)
            index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
            density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
            stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]
            ginis_s[s,h,t] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))
            #=
            # calculate variance of log-s
            avglogs_s[s,h,t] = sum(density_distr.*choice.*log.(max.(1e-12,stat_distr))   )
            varlogs_s[s,h,t] = sum(density_distr.*choice.*(log.(max.(1e-12,stat_distr)).- avglogs_s[s,h,t]).^2)/sum(density_distr.*choice)
            avglogs_s[s,h,t] /= sum(density_distr.*choice)
            if s==3
                avglogs_s[s,h,t] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
                varlogs_s[s,h,t] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avglogs_s[s,h,t]).^2)/sum(density_distr.*choice)
                avglogs_s[s,h,t] /= sum(density_distr.*choice)
            end

            avgs_s[s,h,t] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
            vars_s[s,h,t] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avgs_s[s,h,t]).^2)/sum(density_distr.*choice)
            avgs_s[s,h,t] /= sum(density_distr.*choice)
            =#
            next!(p)
        end

        #=
        if h>2
            denumerator = sum(choice.*density_distr)
            Y = sum(choice.*density_distr.*output)/denumerator
            K = sum(choice.*density_distr.*capital_d)/denumerator
            L = sum(choice.*density_distr.*labour_d)/denumerator
            TFPis_s[h-2,t] = Y/(K^eta*L^theta)
            TFPds_s[h-2,t] = Y/K^(eta/(theta+eta))

            c_1 = capital_d.^(-1.0)
            replace!(c_1,Inf=>0.0)
            mean_MPK_s[h-2,t] = sum(choice.*density_distr.*(eta.*output.*c_1))/denumerator
            var_MPK_s[h-2,t] = sum(choice.*density_distr.*(eta.*output.*c_1 .- mean_MPK_s[h-2,t]*denumerator).^2)/denumerator

            l_1 = labour_d.^(-1.0)
            replace!(l_1,Inf=>0.0)
            mean_MPL_s[h-2,t] = sum(choice.*density_distr.*(theta.*output.*l_1))/denumerator
            var_MPL_s[h-2,t] = sum(choice.*density_distr.*(theta.*output.*l_1 .- mean_MPL_s[h-2,t]*denumerator).^2)/denumerator

            share_unbound_s[h-2,t] = sum(choice.*density_distr.*unlimited_capital_choice)/sum(choice.*density_distr)
            next!(p)
        end
        =#

        #=
        if h>1
            #[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            #m = u_i,alpha_m_i,zeta_i
            m_skill_distr = permutedims(sum(density_distr.*choice, dims=[1,5])[1,:,:,:,1], [1,3,2])
            #w = u_i,alpha_m_i,alpha_w_i
            w_skill_distr = sum(density_distr.*choice, dims=[1,3])[1,:,1,:,:]

            avg_m_skill_s[h-1,t] = sum(m_skill_distr.*z_m_nodes)/sum(density_distr.*choice)
            avg_w_skill_s[h-1,t] = sum(w_skill_distr.*z_w_nodes)/sum(density_distr.*choice)

            var_m_skill_s[h-1,t] = sum(m_skill_distr.*(z_m_nodes .- avg_m_skill_s[h-1,t]*sum(density_distr.*choice)).^2)/sum(density_distr.*choice)
            var_w_skill_s[h-1,t] = sum(w_skill_distr.*(z_w_nodes .- avg_w_skill_s[h-1,t]*sum(density_distr.*choice)).^2)/sum(density_distr.*choice)
            next!(p)
        end
        =#

        if h!=5 && calc_mean_quantile
            quantile_mean_wealth_s[:,h,t] .= quantile_mean(asset_grid.*ones(size(density_distr)), density_distr.*choice)
            next!(p)
            quantile_mean_income_s[:,h,t] .= quantile_mean(income .- ones(size(capital_s_distr_s[t,:,:,:,:,:])).*asset_grid, density_distr.*choice)
            next!(p)
            quantile_mean_earnings_s[:,h,t] .= quantile_mean(earnings, density_distr.*choice)
            next!(p)
            quantile_mean_consumption_s[:,h,t] .= quantile_mean(income .- policy_s[t,:,:,:,:,:], density_distr.*choice)
            next!(p)
        end
    end

    agg_c_e = c_e*sum(density_distr.*Float64.(occ_choice.==3.0))
    share_W_earnings_in_output_s[t] = sum(earnings .* density_distr.*Float64.(occ_choice.==1.0))/(Output_s[t]-agg_c_e)
    share_SP_earnings_in_output_s[t] = sum(earnings .* density_distr.*Float64.(occ_choice.==2.0))/(Output_s[t]-agg_c_e)
    share_EMP_earnings_in_output_s[t] = sum(earnings .* density_distr.*Float64.(occ_choice.==3.0))/(Output_s[t]-agg_c_e)
    share_W_capital_income_in_output_s[t] = (delta+r_s[t])*sum(ones(size(density_distr)).*asset_grid .* density_distr.*Float64.(occ_choice.==1.0))/(Output_s[t]-agg_c_e)
    share_SP_capital_income_in_output_s[t] = (delta+r_s[t])*sum(ones(size(density_distr)).*asset_grid .* density_distr.*Float64.(occ_choice.==2.0))/(Output_s[t]-agg_c_e)
    share_EMP_capital_income_in_output_s[t] = (delta+r_s[t])*sum(ones(size(density_distr)).*asset_grid .* density_distr.*Float64.(occ_choice.==3.0))/(Output_s[t]-agg_c_e)

    next!(p)
end

Output_growth_rate_s = zeros(T)
Capital_growth_rate_s = zeros(T)
Income_growth_rate_s = zeros(T)
Consumption_growth_rate_s = zeros(T)
# # SP, EMP, ENT
# mean_MPL_s = zeros(3,T)
# var_MPL_s = zeros(3,T)
# mean_MPK_s = zeros(3,T)
# var_MPK_s = zeros(3,T)
# # W, SP,EMP, ENT
# avg_w_skill_s = zeros(4,T)
# avg_m_skill_s = zeros(4,T)
# var_w_skill_s = zeros(4,T)
# var_m_skill_s = zeros(4,T)

#additional computation loop (for non-parallel computations)
#=
@showprogress for t=1:T
    if t > 1
        Output_growth_rate_s[t] = (Output_s[t]/Output_s[t-1])-1.0
        Capital_growth_rate_s[t] = (Capital_s[t]/Capital_s[t-1])-1.0
        Income_growth_rate_s[t] = (Income_s[t]/Income_s[t-1])-1.0
        Consumption_growth_rate_s[t] = (Consumption_s[t]/Consumption_s[t-1])-1.0
    else
        Output_growth_rate_s[1] = (Output_s[1]/Output_1)-1.0
        Capital_growth_rate_s[1] = (Capital_s[1]/Capital_1)-1.0
        Income_growth_rate_s[1] = (Income_s[1]/Income_1)-1.0
        Consumption_growth_rate_s[1] = (Consumption_s[1]/Consumption_1)-1.0
    end
end
=#

#create plots

plt = create_plot(collect(1:T),"Time (Years)", Output_s,""#="Output"=#, Output_1, Output_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_outputs.png")
#=
plt = create_plot(collect(1:T),"Time (Years)", Output_growth_rate_s,"Output growth rate", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_output_growth_rate.png")

plt = create_plot(collect(1:T),"Time (Years)", Capital_s,"Capital", Capital_1, Capital_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_capitals.png")

plt = create_plot(collect(1:T),"Time (Years)", Capital_growth_rate_s,"Capital growth rate", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_capital_growth_rate.png")
=#
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
#=
plt = create_plot(collect(1:T),"Time (Years)", Income_growth_rate_s,"Income growth rate", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_income_growth_rate.png")
=#
#plt = create_plot(collect(1:T),"Time (Years)", Consumption_s,"Consumption", Consumption_1, Consumption_2, false)
plt = create_plot(collect(1:T),"Time (Years)", means_s[4,1,:],""#="Consumption"=#, means[4,1,1], means[4,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_consumptions.png")
#=
plt = create_plot(collect(1:T),"Time (Years)", Consumption_growth_rate_s,"Consumption growth rate", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_consumption_growth_rate.png")
=#
plt = create_plot(collect(1:T),"Time (Years)", means_s[2,1,:],""#="Earnings"=#, means[2,1,1], means[2,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", means_s[3,1,:],""#="Wealth"=#, means[3,1,1], means[3,1,2], false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_wealths.png")


#Occupations
LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end
mkpath(LOCAL_DIR_OCCUPATION)
#plt = create_plot(collect(1:T),"Time (Years)", occ_Ws_s,"Share of Workers", occ_Ws_1, occ_Ws_2, true, "W")
plt = create_plot(collect(1:25),"Time (Years)", occ_Ws_s[1:25],""#="Share of Workers"=#, occ_Ws_1, occ_Ws_2, true, "W")
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_workers.png")
# plt = create_plot(collect(1:T),"Time (Years)", occ_ENTs_s,"Share of Entrepreneurs", occ_ENTs_1, occ_ENTs_2, true, "ENT")
# display(plt)
# savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_entrepreneurs.png")
#plt = create_plot(collect(1:T),"Time (Years)", occ_SPs_s,"Share of Sole Proprietors", occ_SPs_1, occ_SPs_2, true, "SP")
plt = create_plot(collect(1:25),"Time (Years)", occ_SPs_s[1:25],""#="Share of Sole Proprietors"=#, occ_SPs_1, occ_SPs_2, true, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_sole_proprietors.png")
#plt = create_plot(collect(1:T),"Time (Years)", occ_EMPs_s,"Share of Employers", occ_EMPs_1, occ_EMPs_2, true, "EMP")
plt = create_plot(collect(1:25),"Time (Years)", occ_EMPs_s[1:25],""#="Share of Employers"=#, occ_EMPs_1, occ_EMPs_2, true, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_employers.png")
#=
plt = create_combined_plot(collect(1:T),"Time (Years)",[occ_Ws_s,occ_SPs_s,occ_EMPs_s],["W","SP","EMP"],"Occupational shares",[occ_Ws_1,occ_SPs_1,occ_EMPs_1],[occ_Ws_2,occ_SPs_2,occ_EMPs_2], true, ["W","SP","EMP"], :right)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_shares_of_occupations.png")
=#
# plt1 = create_plot(collect(1:25),"Time (Years)", occ_Ws_s[1:25],"Share of Workers", occ_Ws_1, occ_Ws_2, true, "W", 5)
# plt2 = create_plot(collect(1:25),"Time (Years)", occ_SPs_s[1:25],"Share of Sole Proprietors", occ_SPs_1, occ_SPs_2, true, "SP", 5)
# plt3 = create_plot(collect(1:25),"Time (Years)", occ_EMPs_s[1:25],"Share of Employers", occ_EMPs_1, occ_EMPs_2, true, "EMP", 5)
# plt = plot(plt1,plt2,plt3, layout=(1,3))
# display(plt)
# savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_shares_of_occupations_zoomed_in.png")
# throw(error)
#Share of unconstrained ENT, SP,EMP
#=
plt = create_plot(collect(1:T),"Time (Years)", share_ENT_unbound_s,"Share of Unconstrained Entrepreneurs", share_ENT_unbound_1, share_ENT_unbound_2, true, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_unconstrained_entrepreneurs.png")
plt = create_plot(collect(1:T),"Time (Years)", share_SP_unbound_s,"Share of Unconstrained Sole Proprietors", share_SP_unbound_1, share_SP_unbound_2, true, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_unconstrained_sole_proprietors.png")
plt = create_plot(collect(1:T),"Time (Years)", share_EMP_unbound_s,"Share of Unconstrained Employers", share_EMP_unbound_1, share_EMP_unbound_2, true, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_unconstrained_employers.png")
=#

#Inequality
LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#Variance of log-consumption and log-earnings, Gini for W and ENT income
#=
plt = create_plot(collect(1:T),"Time (Years)", logcs_s,"Variance of log-consumption", logcs_1, logcs_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_log_consumption.png")
plt = create_plot(collect(1:T),"Time (Years)", loges_s,"Variance of log-earnings", loges_1, loges_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_log_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", giniWs_s,"Gini for Workers' income", giniWs_1, giniWs_2, false, "W")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_workers_income.png")
plt = create_plot(collect(1:T),"Time (Years)", giniEnts_s,"Gini for Entrepreneurs' income", giniEnts_1, giniEnts_2, false, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_entrepreneurs_income.png")
=#

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

        # calculate gini coefficent
        #ginis[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", ginis_s[s,h,:],""#="Gini of $(CHOICE_NAMES[h])' $stat_name"=#, ginis[s,h,1], ginis[s,h,2], false,LOCAL_COLORS[h])
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_$(choice_name)_$(stat_name).png")
        #=
        # calculate variance of log-s
        #avglogs[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", avglogs_s[s,h,:],"Average of $(CHOICE_NAMES[h])' Log-$stat_name", avglogs[s,h,1], avglogs[s,h,2], false,LOCAL_COLORS[h])
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_avg_$(choice_name)_log$(stat_name).png")
        #varlogs[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", varlogs_s[s,h,:],"Variance of $(CHOICE_NAMES[h])' Log-$stat_name", varlogs[s,h,1], varlogs[s,h,2], false,LOCAL_COLORS[h])
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_$(choice_name)_log$(stat_name).png")

        # calculate variance of s
        #avgs[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", avgs_s[s,h,:],"Average of $(CHOICE_NAMES[h])' $stat_name", avgs[s,h,1], avgs[s,h,2], false,LOCAL_COLORS[h])
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_avg_$(choice_name)_$(stat_name).png")
        #vars[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", vars_s[s,h,:],"Variance of $(CHOICE_NAMES[h])' $stat_name", vars[s,h,1], vars[s,h,2], false,LOCAL_COLORS[h])
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_$(choice_name)_$(stat_name).png")
        =#
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
    for stat = [1,4]#1:4 # wealth, income, earnings, consumption
        stat_name = "Wealth"
        qm = quantile_mean_wealth
        qm_s = quantile_mean_wealth_s
        if stat == 2
            stat_name = "Income"
            qm = quantile_mean_income
            qm_s = quantile_mean_income_s
        elseif stat == 3
            stat_name = "Earnings"
            qm = quantile_mean_earnings
            qm_s = quantile_mean_earnings_s
        elseif stat == 4
            stat_name = "Consumption"
            qm = quantile_mean_consumption
            qm_s = quantile_mean_consumption_s
        end

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
# TFP_ideal, TFP_data, mean and var of MPK and MPL, share_of_unconstrained SP,EMP,ENT
LABELS=["SP","EMP"]
#=
plt = create_combined_plot(collect(1:T),"Time (Years)", [TFPis_s[1,:],TFPis_s[2,:]],LABELS,"TFP", [TFPis[1,1],TFPis[2,1]], [TFPis[1,2],TFPis[2,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_tfp_ideal.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [mean_MPL_s[1,:],mean_MPL_s[2,:]],LABELS,"Mean of MPL", [mean_MPL[1,1],mean_MPL[2,1]], [mean_MPL[1,2],mean_MPL[2,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_mean_mpl.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [var_MPL_s[1,:],var_MPL_s[2,:]],LABELS,"Variance of MPL", [var_MPL[1,1],var_MPL[2,1]], [var_MPL[1,2],var_MPL[2,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_var_mpl.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [mean_MPK_s[1,:],mean_MPK_s[2,:]],LABELS,"Mean of MPK", [mean_MPK[1,1],mean_MPK[2,1]], [mean_MPK[1,2],mean_MPK[2,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_mean_mpk.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [var_MPK_s[1,:],var_MPK_s[2,:]],LABELS,"Variance of MPK", [var_MPK[1,1],var_MPK[2,1]], [var_MPK[1,2],var_MPK[2,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_var_mpk.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [share_unbound_s[1,:],share_unbound_s[2,:]],LABELS,"Share of Unconstrained Entrepreneurs", [share_unbound[1,1],share_unbound[2,1]], [share_unbound[1,2],share_unbound[2,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_share_unconstrained_sp_emp.png")
=#

# SP, EMP, ENT
#=
CHOICE_NAMES = ["Sole Proprietors","Employers","Entrepreneurs"]
LOCAL_COLORS = ["SP","EMP","ENT"]
for h=1:3
    choice_name = CHOICE_NAMES[h]
    if h==1
        choice_name = "SoleProprietors"
    end
    plt = create_plot(collect(1:T),"Time (Years)", TFPis_s[h,:],"TFP for $(CHOICE_NAMES[h])", TFPis[h,1], TFPis[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_tfp_ideal_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", TFPds_s[h,:],"TFP_data for $(CHOICE_NAMES[h])", TFPds[h,1], TFPds[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_tfp_data_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", mean_MPL_s[h,:],"Mean of MPL for $(CHOICE_NAMES[h])", mean_MPL[h,1], mean_MPL[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_mean_mpl_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_MPL_s[h,:],"Variance of MPL for $(CHOICE_NAMES[h])", var_MPL[h,1], var_MPL[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_mpl_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", mean_MPK_s[h,:],"Mean of MPK for $(CHOICE_NAMES[h])", mean_MPK[h,1], mean_MPK[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_mean_mpk_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_MPK_s[h,:],"Variance of MPK for $(CHOICE_NAMES[h])", var_MPK[h,1], var_MPK[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_mpk_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", share_unbound_s[h,:],"Share of Unconstrained $(CHOICE_NAMES[h])", share_unbound[h,1], share_unbound[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_unconstrained_$(choice_name).png")

end
=#
#skills
LABELS=["W", "SP","EMP", "ENT"]
#=
plt = create_combined_plot(collect(1:T),"Time (Years)", [avg_m_skill_s[1,:],avg_m_skill_s[2,:],avg_m_skill_s[3,:],avg_m_skill_s[4,:]],LABELS,"Average Managerial Skill", [avg_m_skill[1,1],avg_m_skill[2,1],avg_m_skill[3,1],avg_m_skill[4,1]], [avg_m_skill[1,2],avg_m_skill[2,2],avg_m_skill[3,2],avg_m_skill[4,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_avg_m_skill.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [avg_w_skill_s[1,:],avg_w_skill_s[2,:],avg_w_skill_s[3,:],avg_w_skill_s[4,:]],LABELS,"Average Working Skill", [avg_w_skill[1,1],avg_w_skill[2,1],avg_w_skill[3,1],avg_w_skill[4,1]], [avg_w_skill[1,2],avg_w_skill[2,2],avg_w_skill[3,2],avg_w_skill[4,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_avg_w_skill.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [var_m_skill_s[1,:],var_m_skill_s[2,:],var_m_skill_s[3,:],var_m_skill_s[4,:]],LABELS,"Variance of Managerial Skill", [var_m_skill[1,1],var_m_skill[2,1],var_m_skill[3,1],var_m_skill[4,1]], [var_m_skill[1,2],var_m_skill[2,2],var_m_skill[3,2],var_m_skill[4,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_var_m_skill.png")

plt = create_combined_plot(collect(1:T),"Time (Years)", [var_w_skill_s[1,:],var_w_skill_s[2,:],var_w_skill_s[3,:],var_w_skill_s[4,:]],LABELS,"Variance of Working Skill", [var_w_skill[1,1],var_w_skill[2,1],var_w_skill[3,1],var_w_skill[4,1]], [var_w_skill[1,2],var_w_skill[2,2],var_w_skill[3,2],var_w_skill[4,2]], false,LABELS)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_combined_var_w_skill.png")
CHOICE_NAMES = ["Workers","Sole Proprietors","Employers","Entrepreneurs"]
LOCAL_COLORS = ["W","SP","EMP","ENT"]
for h=1:4
    choice_name = CHOICE_NAMES[h]
    if h==2
        choice_name = "SoleProprietors"
    end
    plt = create_plot(collect(1:T),"Time (Years)", avg_m_skill_s[h,:],"Average Managerial Skill across $(CHOICE_NAMES[h])", avg_m_skill[h,1], avg_m_skill[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_avg_m_skill_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", avg_w_skill_s[h,:],"Average Working Skill across $(CHOICE_NAMES[h])", avg_w_skill[h,1], avg_w_skill[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_avg_w_skill_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_m_skill_s[h,:],"Variance of Managerial Skill across $(CHOICE_NAMES[h])", var_m_skill[h,1], var_m_skill[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_m_skill_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_w_skill_s[h,:],"Variance of Working Skill across $(CHOICE_NAMES[h])", var_w_skill[h,1], var_w_skill[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_w_skill_$(choice_name).png")

end
=#

CHOICE_NAMES = ["Households","Sole Proprietors","Employers","Entrepreneurs"]
LOCAL_COLORS = ["H","SP","EMP","ENT"]
# Credit for All, SP, EMP, ENT
#=
for h = 1:4
    choice_name = CHOICE_NAMES[h]
    if h==2
        choice_name = "SoleProprietors"
    end
    # calculate variance of s
    #avgs[s,h,i]
    plt = create_plot(collect(1:T),"Time (Years)", Credit_s[h,:],"Mean of $(CHOICE_NAMES[h])' credit", Credit[h,1], Credit[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_avg_$(choice_name)_credit.png")
    #vars[s,h,i]
    plt = create_plot(collect(1:T),"Time (Years)", var_Credit_s[h,:],"Variance of $(CHOICE_NAMES[h])' credit", var_Credit[h,1], var_Credit[h,2], false,LOCAL_COLORS[h])
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_$(choice_name)_credit.png")

end
=#

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

# using StatsPlots
# earnings_in_output = #=[Output_1.-Output_1;Output_s.-Output_1].*=#[[share_W_earnings_in_output[1];share_W_earnings_in_output_s] [share_SP_earnings_in_output[1];share_SP_earnings_in_output_s] [share_EMP_earnings_in_output[1];share_EMP_earnings_in_output_s]]
# plt = groupedbar(earnings_in_output,
#         bar_position = :stack,
#         xticks=(collect(1:5:T+1), collect(0:5:T)),
#         label=["Workers" "Sole Prop." "Employers"],
#         color=[:purple :red :green],
#         legend=:right,
#         xlabel="Time (Years)",
#         ylabel="Shares of Earnings in Output")
# display(plt)
# savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_earnings.png")
#
# capital_income_in_output = #=[Output_1.-Output_1;Output_s.-Output_1].*=#[[share_W_capital_income_in_output[1];share_W_capital_income_in_output_s] [share_SP_capital_income_in_output[1];share_SP_capital_income_in_output_s] [share_EMP_capital_income_in_output[1];share_EMP_capital_income_in_output_s]]
# plt = groupedbar(capital_income_in_output,
#         bar_position = :stack,
#         xticks=(collect(1:5:T+1), collect(0:5:T)),
#         label=["Workers" "Sole Prop." "Employers"],
#         color=[:purple :red :green],
#         legend=:right,
#         xlabel="Time (Years)",
#         ylabel="Shares of Capital Income in Output")
# display(plt)
# savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_of_output_capital_income.png")

#end

# check derivatives of factor prices
# deriv_r = diff(r_s)
# plt1 = plot(collect(3:T),deriv_r[2:end])
# deriv_w = diff(w_s)
# plt2 = plot(collect(3:T),deriv_w[2:end])
# plt3 = plot(collect(3:T),deriv_r[2:end].+deriv_w[2:end])
# plot(plt1,plt2,plt3, layout=(3,1))
#conclusion: in the first period both factor prices increase, after that interest rate falls and wage increases at the same speed

#                    1[1] 1[2] 1[3]
trans_SSS =[  [ss_1,ss_2,trans_res],
#                    2[1]     2[2]     2[3]
                    [lambda_1,lambda_2,lambda_s],
#                    3[1] 3[2] 3[3]
                    [r_1, r_2, r_s],
#                    4[1] 4[2] 4[3]
                    [w_1, w_2, w_s],
#                    5[1]      5[2]       5[3]
                    [Capital_1,Capital_2, Capital_s],
#                    6[1]          6[2]          6[3]
                    [Consumption_1,Consumption_2,Consumption_s],
#                    7[1]        7[2]        7[3]
                    [Credit[:,1],Credit[:,2],Credit_s],
#                    8[1]               8[2]               8[3]
                    [Credit_to_Output_1,Credit_to_Output_2,Credit_to_Output_s],
#                    9[1]     9[2]     9[3]
                    [Income_1,Income_2,Income_s],
#                    10[1]    10[2]    10[3]
                    [Output_1,Output_2,Output_s],
#                    11[1]        11[2]        11[3]
                    [ginis[:,:,1],ginis[:,:,2],ginis_s],
#                    12[1]        12[2]        12[3]
                    [means[:,:,1],means[:,:,2],means_s],
#                    13[1]           13[2]           13[3]
                    [var_Credit[:,1],var_Credit[:,2],var_Credit_s],
#                    14[1]                   14[2]                   14[3]
                    [quantile_means[:,:,:,1],quantile_means[:,:,:,2],quantile_means_s],
#                    15[1]                         15[2]                         15[3]
                    [share_W_earnings_in_output[1],share_W_earnings_in_output[2],share_W_earnings_in_output_s],
#                    16[1]                          16[2]                          16[3]
                    [share_SP_earnings_in_output[1],share_SP_earnings_in_output[2],share_SP_earnings_in_output_s],
#                    17[1]                           17[2]                           17[3]
                    [share_EMP_earnings_in_output[1],share_EMP_earnings_in_output[2],share_EMP_earnings_in_output_s],
#                    18[1]                               18[2]                               18[3]
                    [share_W_capital_income_in_output[1],share_W_capital_income_in_output[2],share_W_capital_income_in_output_s],
#                    19[1]                                19[2]                                19[3]
                    [share_SP_capital_income_in_output[1],share_SP_capital_income_in_output[2],share_SP_capital_income_in_output_s],
#                    20[1]                                 20[2]                                 20[3]
                    [share_EMP_capital_income_in_output[1],share_EMP_capital_income_in_output[2],share_EMP_capital_income_in_output_s]
                    ]

@save "$(LOCAL_DIR_GENERAL)trans_SSS.jld2" trans_SSS
