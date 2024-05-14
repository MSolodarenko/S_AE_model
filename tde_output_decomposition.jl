using JLD2
using Plots

using ProgressMeter
using SchumakerSpline
using BasicInterpolators: LinearInterpolator

include("Functions/profit.jl")

TIME_PERIODS = 50

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE\\General\\"
end
@load "$(LOCAL_DIR)trans_SSS.jld2" trans_SSS

# trans_res = trans_SSS[1][3]
density_distr_s = [trans_SSS[1][3][5][7][t,:,:,:,:,:] for t=1:TIME_PERIODS]
asset_grid = trans_SSS[1][3][10]
number_asset_grid = trans_SSS[1][3][9]
r_s = trans_SSS[1][3][3]
w_s = trans_SSS[1][3][4]
# ss_1 = trans_SSS[1][1]
# global_approx_params = trans_SSS[1][1][1][47]
number_zeta_nodes = trans_SSS[1][1][1][47][4]
number_alpha_m_nodes = trans_SSS[1][1][1][47][5]
number_alpha_w_nodes = trans_SSS[1][1][1][47][6]
lambda_s = trans_SSS[1][3][2]
# model_params = trans_SSS[1][1][5]
delta = trans_SSS[1][1][5][3]
gamma = trans_SSS[1][1][5][4]
eta = trans_SSS[1][1][5][5]
theta = trans_SSS[1][1][5][6]
c_e = trans_SSS[1][1][5][7]
# approx_object = trans_SSS[1][1][4]
z_m_nodes = trans_SSS[1][1][4][1]
z_w_nodes = trans_SSS[1][1][4][2]
number_u_nodes = trans_SSS[1][1][4][7]

# Conduct a Solow style decomposition exercise
# Let's assume that Y = A*M^gamma*K^eta*L^theta,
# where gamma+eta+theta=1
# and A will represent the production augmenting technology
# which will try to account for all inefficiencies
Output = zeros(TIME_PERIODS)
Management = zeros(TIME_PERIODS)
Capital = zeros(TIME_PERIODS)
Labour = zeros(TIME_PERIODS)
Technology = zeros(TIME_PERIODS)

MPM = zeros(TIME_PERIODS)
MPK = zeros(TIME_PERIODS)
MPL = zeros(TIME_PERIODS)

@showprogress for t = 1:TIME_PERIODS
    # t = 1

    occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input = compute_income_profile(asset_grid,number_asset_grid,r_s[t], w_s[t], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

    Output[t] = sum(output .* density_distr_s[t])
    Management[t] = sum(managerial_input .* density_distr_s[t])
    Capital[t] = sum(capital_d .* density_distr_s[t])
    Labour[t] = sum(labour_d .* density_distr_s[t])
    # A = Y/(M^gamma*K^eta*L^theta)
    Technology[t] = Output[t]/(Management[t]^gamma * Capital[t]^eta * Labour[t]^theta)

    # Y = A*M^gamma*K^eta*L^theta
    MPM[t] = gamma*Output[t]/Technology[t]
    MPK[t] = eta*Output[t]/Capital[t]
    MPL[t] = theta*Output[t]/Labour[t]
end

# log_change_Output = zeros(TIME_PERIODS)
# log_change_Management = zeros(TIME_PERIODS)
# log_change_Capital = zeros(TIME_PERIODS)
# log_change_Labour = zeros(TIME_PERIODS)
# log_change_Technology = zeros(TIME_PERIODS)
#
# @showprogress for t = 1:TIME_PERIODS
#     # t=2
#     log_change_Output[t] = log(Output[t]/Output[1])
#     log_change_Management[t] = gamma*log(Management[t]/Management[1])
#     log_change_Capital[t] = eta*log(Capital[t]/Capital[1])
#     log_change_Labour[t] = theta*log(Labour[t]/Labour[1])
#     log_change_Technology[t] = log(Technology[t]/Technology[1])
#
# end

# p11 = plot(1:TIME_PERIODS, log_change_Output, legend=false)
# p12 = plot(1:TIME_PERIODS, log_change_Management, legend=false)
# p13 = plot(1:TIME_PERIODS, log_change_Capital, legend=false)
# p14 = plot(1:TIME_PERIODS, log_change_Labour, legend=false)
# p15 = plot(1:TIME_PERIODS, log_change_Technology, legend=false)
# display(plot(p11,p12,p13,p14,p15, layout=(5,1)))


LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE_fixed_occ/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE_fixed_occ\\General\\"
end
@load "$(LOCAL_DIR)trans_SSS_fixed.jld2" trans_SSS_fixed

# trans_res = trans_SSS[1][3]
density_distr_s = [[trans_SSS_fixed[1][3][5][7][occ][t,:,:,:,:,:] for t=1:TIME_PERIODS] for occ=1:3]
asset_grid = trans_SSS_fixed[1][3][10]
number_asset_grid = trans_SSS_fixed[1][3][9]
r_s = trans_SSS_fixed[1][3][3]
w_s = trans_SSS_fixed[1][3][4]
# ss_1 = trans_SSS[1][1]
# global_approx_params = trans_SSS[1][1][1][47]
number_zeta_nodes = trans_SSS_fixed[1][1][1][47][4]
number_alpha_m_nodes = trans_SSS_fixed[1][1][1][47][5]
number_alpha_w_nodes = trans_SSS_fixed[1][1][1][47][6]
lambda_s = trans_SSS_fixed[1][3][2]
# model_params = trans_SSS[1][1][5]
delta = trans_SSS_fixed[1][1][5][3]
gamma = trans_SSS_fixed[1][1][5][4]
eta = trans_SSS_fixed[1][1][5][5]
theta = trans_SSS_fixed[1][1][5][6]
c_e = trans_SSS_fixed[1][1][5][7]
# approx_object = trans_SSS[1][1][4]
z_m_nodes = trans_SSS_fixed[1][1][4][1]
z_w_nodes = trans_SSS_fixed[1][1][4][2]
number_u_nodes = trans_SSS_fixed[1][1][4][7]

# Conduct a Solow style decomposition exercise
# Let's assume that Y = AM^gamma*K^eta*L^theta,
# where gamma+eta+theta=1
# and A will represent the production augmenting technology
# which will try to account for all inefficiencies
Output_fixed_occ = zeros(TIME_PERIODS)
Management_fixed_occ = zeros(TIME_PERIODS)
Capital_fixed_occ = zeros(TIME_PERIODS)
Labour_fixed_occ = zeros(TIME_PERIODS)
Technology_fixed_occ = zeros(TIME_PERIODS)

MPM_fixed_occ = zeros(TIME_PERIODS)
MPK_fixed_occ = zeros(TIME_PERIODS)
MPL_fixed_occ = zeros(TIME_PERIODS)

@showprogress for t = 1:TIME_PERIODS
    # t = 1

    occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input = compute_income_profile_fixed_occ(asset_grid,number_asset_grid,r_s[t], w_s[t], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

    Output_fixed_occ[t] = sum([sum(output[occ] .* density_distr_s[occ][t]) for occ=1:3])
    Management_fixed_occ[t] = sum([sum(managerial_input[occ] .* density_distr_s[occ][t]) for occ=1:3])
    Capital_fixed_occ[t] = sum([sum(capital_d[occ] .* density_distr_s[occ][t]) for occ=1:3])
    Labour_fixed_occ[t] = sum([sum(labour_d[occ] .* density_distr_s[occ][t]) for occ=1:3])
    # A = Y/(M^gamma*K^eta*L^theta)
    Technology_fixed_occ[t] = Output_fixed_occ[t]/(Management_fixed_occ[t]^gamma * Capital_fixed_occ[t]^eta * Labour_fixed_occ[t]^theta)

    # Y = A*M^gamma*K^eta*L^theta
    MPM_fixed_occ[t] = gamma*Output_fixed_occ[t]/Technology_fixed_occ[t]
    MPK_fixed_occ[t] = eta*Output_fixed_occ[t]/Capital_fixed_occ[t]
    MPL_fixed_occ[t] = theta*Output_fixed_occ[t]/Labour_fixed_occ[t]
end
throw(error)

res1 = ["Output","MFP","Management","Capital","Labour"]
res2 = round.(100 .*log.([Output[end]/Output[1],Technology[end]/Technology[1],(Management[end]/Management[1])^gamma,(Capital[end]/Capital[1])^eta,(Labour[end]/Labour[1])^theta]); digits=1)
res3 = round.(100 .*log.([Output_fixed_occ[end]/Output_fixed_occ[1],Technology_fixed_occ[end]/Technology_fixed_occ[1],(Management_fixed_occ[end]/Management_fixed_occ[1])^gamma,(Capital_fixed_occ[end]/Capital_fixed_occ[1])^eta,(Labour_fixed_occ[end]/Labour_fixed_occ[1])^theta]); digits=1)
res4 = round.(res2.-res3; digits=1)

println(res1)
println(res2)
println(res3)
println(res4)

# Contribution & Output  & MFP              & Management       & Capital          & Labour  \\ \hline
# Total        & +7.85\% & +2.25pp          & +0.52pp          & +4.17pp          & +0.91pp \\
# Direct       & +5.67pp & \textbf{+2.25pp} & \textbf{+0.19pp} & \textbf{+3.23pp} & +0.00pp \\
# Indirect     & +2.18pp & --0.00pp         & \textbf{+0.33pp} & \textbf{+0.94pp} & \textbf{+0.91pp}     \\ \hline

# management^gamma
generate_plots(Time,"Time",Management.^gamma,Management_fixed_occ.^gamma,"Management",LOCAL_DIR_PRODUCTIVITY,"time_management",false)

# capital^eta
generate_plots(Time,"Time",Capital.^eta,Capital_fixed_occ.^eta,"Capital",LOCAL_DIR_PRODUCTIVITY,"time_capital",false)

# labour^theta
generate_plots(Time,"Time",Labour.^theta,Labour_fixed_occ.^theta,"Labour",LOCAL_DIR_PRODUCTIVITY,"time_labour",false)

# economy efficiency
generate_plots(Time,"Time",Technology,Technology_fixed_occ,"Efficiency",LOCAL_DIR_PRODUCTIVITY,"time_efficiency",false,:bottomright)
# log_change_Output_fixed_occ = zeros(TIME_PERIODS)
# log_change_Management_fixed_occ = zeros(TIME_PERIODS)
# log_change_Capital_fixed_occ = zeros(TIME_PERIODS)
# log_change_Labour_fixed_occ = zeros(TIME_PERIODS)
# log_change_Technology_fixed_occ = zeros(TIME_PERIODS)
#
# @showprogress for t = 1:TIME_PERIODS
#     # t=2
#     log_change_Output_fixed_occ[t] = log(Output_fixed_occ[t]/Output_fixed_occ[1])
#     log_change_Management_fixed_occ[t] = gamma*log(Management_fixed_occ[t]/Management_fixed_occ[1])
#     log_change_Capital_fixed_occ[t] = eta*log(Capital_fixed_occ[t]/Capital_fixed_occ[1])
#     log_change_Labour_fixed_occ[t] = theta*log(Labour_fixed_occ[t]/Labour_fixed_occ[1])
#     log_change_Technology_fixed_occ[t] = log(Technology_fixed_occ[t]/Technology_fixed_occ[1])
#
# end

# p21 = plot(1:TIME_PERIODS, log_change_Output_fixed_occ, legend=false)
# p22 = plot(1:TIME_PERIODS, log_change_Management_fixed_occ, legend=false)
# p23 = plot(1:TIME_PERIODS, log_change_Capital_fixed_occ, legend=false)
# p24 = plot(1:TIME_PERIODS, log_change_Labour_fixed_occ, legend=false)
# p25 = plot(1:TIME_PERIODS, log_change_Technology_fixed_occ, legend=false)
# display(plot(p21,p22,p23,p24,p25, layout=(5,1)))

# display(plot(p11,p21,
#              p12,p22,
#              p13,p23,
#              p14,p24,
#              p15,p25, layout=(5,2)))

# diff_output = log_change_Output .- log_change_Output_fixed_occ
# diff_management = log_change_Management .- log_change_Management_fixed_occ
# diff_capital = log_change_Capital .- log_change_Capital_fixed_occ
# diff_labour = log_change_Labour .- log_change_Labour_fixed_occ
# diff_technology = log_change_Technology .- log_change_Technology_fixed_occ

# p31 = plot(1:TIME_PERIODS, diff_output, legend=false)
# p32 = plot(1:TIME_PERIODS, diff_management, legend=false)
# p33 = plot(1:TIME_PERIODS, diff_capital, legend=false)
# p34 = plot(1:TIME_PERIODS, diff_labour, legend=false)
# p35 = plot(1:TIME_PERIODS, diff_technology, legend=false)
#
# display(plot(p11,p21,p31,
#              p12,p22,p32,
#              p13,p23,p33,
#              p14,p24,p34,
#              p15,p25,p35, layout=(5,3)))

# plot the results
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

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE_comparison/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE_comparison\\"
end
mkpath(LOCAL_DIR)

country = "Italy"
TIME_PERIODS = 50#100
Time = collect(1:TIME_PERIODS)

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/Decomposition/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\Decomposition\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)

function generate_plots(X,XLABEL,Y,Y_fixed_occ,YLABEL,PATHDIR,FILENAME,IS_PERCENTAGE::Bool=false, LEGENDPOS=false)
    # LEGENDPOS = :bottomright
    plt = create_combined_plot(X,XLABEL, [Y, Y_fixed_occ],["with mobility", "without mobility"],#=YLABEL=#"", [Y[1],Y_fixed_occ[1]],[Y[end],Y_fixed_occ[end]], IS_PERCENTAGE, LEGENDPOS)
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

        plt = create_combined_plot(X,XLABEL, [log.(Y./Y[1]), log.(Y_fixed_occ./Y_fixed_occ[1])],["with mobility", "without mobility"],#=YLABEL=#"", [0.0,0.0],[log.(Y[end]./Y[1]),log.(Y_fixed_occ[end]./Y_fixed_occ[1])], true, LEGENDPOS)
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

        plt = create_combined_plot(X,XLABEL, [Y.-Y[1], Y_fixed_occ.-Y_fixed_occ[1]],["with mobility", "without mobility"],#=YLABEL=#"", [0.0,0.0],[Y[end].-Y[1],Y_fixed_occ[end].-Y_fixed_occ[1]], true, LEGENDPOS)
        display(plt)
        savefig(plt,"$(PATHDIR)$(country)_$(FILENAME)_growth_rate.png")
    end
end
throw(error)
# output
generate_plots(Time,"Time",Output,Output_fixed_occ,"Output",LOCAL_DIR_PRODUCTIVITY,"time_output",false)

# management^gamma
generate_plots(Time,"Time",Management.^gamma,Management_fixed_occ.^gamma,"Management",LOCAL_DIR_PRODUCTIVITY,"time_management",false)

# capital^eta
generate_plots(Time,"Time",Capital.^eta,Capital_fixed_occ.^eta,"Capital",LOCAL_DIR_PRODUCTIVITY,"time_capital",false)

# labour^theta
generate_plots(Time,"Time",Labour.^theta,Labour_fixed_occ.^theta,"Labour",LOCAL_DIR_PRODUCTIVITY,"time_labour",false)

# economy efficiency
generate_plots(Time,"Time",Technology,Technology_fixed_occ,"Efficiency",LOCAL_DIR_PRODUCTIVITY,"time_efficiency",false,:bottomright)

# marginal products
generate_plots(Time,"Time",MPM,MPM_fixed_occ,"Marginal product of management",LOCAL_DIR_PRODUCTIVITY,"time_mpm",false)
generate_plots(Time,"Time",MPK,MPK_fixed_occ,"Marginal product of capital",LOCAL_DIR_PRODUCTIVITY,"time_mpk",false)
generate_plots(Time,"Time",MPL,MPL_fixed_occ,"Marginal product of labour",LOCAL_DIR_PRODUCTIVITY,"time_mpl",false,:bottomright)
