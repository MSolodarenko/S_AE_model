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

include("Functions/profit.jl")


# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [69,             3,                3,                3,                 6,                    3]

country = "Italy"
# parameters of the model's economy (Italy)
LAMBDA =            1.665907
# parameters of the model's economy (Italy)
LOCAL_DIR = "$(@__DIR__)/Results/Stationary/GE_lambda_fixed_occ/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Stationary\\GE_lambda_fixed_occ\\"
end
global_approx_params = copy(GLOBAL_APPROX_PARAMS)

print_sameline("Loading data from Fixed_occ_shares/Lambda_grid")
#@load "$(LOCAL_DIR)SSS.jld2" SSS
SSS_names = readdir(LOCAL_DIR)
SSS_names = SSS_names[findall(x->occursin("SS_",x), SSS_names)]
SSS = Array{Any}(undef,length(SSS_names))
lambdas = zeros(length(SSS))
SS = []
@showprogress for i = 1:length(SSS)
    @load "$(LOCAL_DIR)$(SSS_names[i])" SS
    SSS[i] = copy(SS)
    lambdas[i] = SS[5][1]
end
temp_is = sortperm(lambdas)
lambdas = lambdas[temp_is]
SSS = SSS[temp_is]

num_lambdas = length(lambdas)#20# #instead of length(lambda)
calibrated_lambda = findfirst(x -> LAMBDA-1e-1<=x<=LAMBDA+1e-1, lambdas)

C_Ys = zeros(num_lambdas)
Outputs = zeros(num_lambdas)

Rs = zeros(num_lambdas)
Ws = zeros(num_lambdas)

# Distribution of output to different channels of earnigs/income to occupations
share_W_earnings_in_output = zeros(num_lambdas)
share_SP_earnings_in_output = zeros(num_lambdas)
share_EMP_earnings_in_output = zeros(num_lambdas)
share_W_capital_income_in_output = zeros(num_lambdas)
share_SP_capital_income_in_output = zeros(num_lambdas)
share_EMP_capital_income_in_output = zeros(num_lambdas)
# throw(error)

@showprogress for i = 1:num_lambdas
    # i=10
    #println("\n$(country) - $(i)/$(num_lambdas)")
    # ss_star = [res, r, w, approx_object, params]
    ss_star = copy(SSS[i])

    number_u_nodes = ss_star[4][7]
    number_asset_grid = ss_star[1][2]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

    delta = ss_star[5][3]
    gamma = ss_star[5][4]
    eta = ss_star[5][5]
    theta = ss_star[5][6]
    c_e = ss_star[5][7]

    z_m_nodes = ss_star[4][1]
    z_w_nodes = ss_star[4][2]

    asset_grid = ss_star[1][3]
    density_distr = ss_star[1][5]
    policy = ss_star[1][4]
    earnings = ss_star[1][24]
    capital_d = ss_star[1][26]
    labour_d = ss_star[1][29]
    output = ss_star[1][32]

    C_Ys[i] = ss_star[1][13]
    Rs[i] = ss_star[2]
    Ws[i] = ss_star[3]

    occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input = compute_income_profile_fixed_occ(asset_grid,number_asset_grid,Rs[i], Ws[i], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambdas[i], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

    wealth = Array{Any}(undef,3)
    for occ in 1:3
        wealth[occ] = ones(size(density_distr[occ])).*asset_grid
        income[occ] = income[occ] .- wealth[occ]
    end

    Outputs[i] = sum( [sum(output[occ] .* density_distr[occ]) for occ in 1:3] )

    agg_c_e = c_e*sum(density_distr[3])

    (Outputs[i]-agg_c_e)

    tA = sum([sum(earnings[1] .* density_distr[1]),
        sum(earnings[2] .* density_distr[2]),
        sum(earnings[3] .* density_distr[3]),
        (delta+Rs[i])*sum(ones(size(density_distr[1])).*asset_grid .* density_distr[1]),
        (delta+Rs[i])*sum(ones(size(density_distr[2])).*asset_grid .* density_distr[2]),
        (delta+Rs[i])*sum(ones(size(density_distr[3])).*asset_grid .* density_distr[3])])

    tA - (Outputs[i]-agg_c_e)

    agg_c_e = c_e*sum(density_distr[3])
    share_W_earnings_in_output[i] = sum(earnings[1] .* density_distr[1])/(Outputs[i]-agg_c_e)
    share_SP_earnings_in_output[i] = sum(earnings[2] .* density_distr[2])/(Outputs[i]-agg_c_e)
    share_EMP_earnings_in_output[i] = sum(earnings[3] .* density_distr[3])/(Outputs[i]-agg_c_e)
    share_W_capital_income_in_output[i] = (delta+Rs[i])*sum(ones(size(density_distr[1])).*asset_grid .* density_distr[1])/(Outputs[i]-agg_c_e)
    share_SP_capital_income_in_output[i] = (delta+Rs[i])*sum(ones(size(density_distr[2])).*asset_grid .* density_distr[2])/(Outputs[i]-agg_c_e)
    share_EMP_capital_income_in_output[i] = (delta+Rs[i])*sum(ones(size(density_distr[3])).*asset_grid .* density_distr[3])/(Outputs[i]-agg_c_e)

    #throw(error)
end


function create_plot(X::Vector{Float64},XLABEL::String,Y::Vector{Float64},YLABEL::String, IS_Y_PERCENTAGE::Bool=true, OCCUPATION::String="H", LOW_LIMIT=-Inf)
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
function create_combined_plot(X::Vector{Float64},XLABEL::String,Ys,YLABELs,YLABEL, IS_Y_PERCENTAGE::Bool=true, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)


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
                        legend=:outertopright,
                        xlabel=XLABEL,
                        ylabel=YLABEL,
                        label=YLABELs[y_i],
                        yticks = YTICKS,
                        ylims = YLIMS )
    end

    return plt
end

#throw(error)

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)

plt = create_plot(lambdas,"λ", Outputs,"Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Outputs.png")

#Interest rate and wage
plt = create_plot(C_Ys,"Credit/Output", Rs,"Interest Rate")
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_rate.png")
plt = create_plot(C_Ys,"Credit/Output", Ws,"Wage", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_wage.png")

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)

plts = create_plots(C_Ys,"Credit/Output", [share_W_earnings_in_output, share_SP_earnings_in_output, share_EMP_earnings_in_output, share_W_capital_income_in_output, share_SP_capital_income_in_output, share_EMP_capital_income_in_output],["Share of output as Workers' Earnings","Share of output as Sole Proprietors' Earnings","Share of output as Employers' Earnings","Share of output as Workers' Capital Income","Share of output as Sole Proprietors' Capital Income","Share of output as Employers' Capital Income"],true,["W","SP","EMP","W","SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_W_earnings.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_SP_earnings.png")
savefig(plts[3],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_EMP_earnings.png")
savefig(plts[4],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_W_capital_income.png")
savefig(plts[5],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_SP_capital_income.png")
savefig(plts[6],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_EMP_capital_income.png")

@load "$(LOCAL_DIR_GENERAL)SSS_fixed.jld2" share_W_earnings_in_output_fixed_occ share_SP_earnings_in_output_fixed_occ share_EMP_earnings_in_output_fixed_occ share_W_capital_income_in_output_fixed_occ share_SP_capital_income_in_output_fixed_occ share_EMP_capital_income_in_output_fixed_occ

share_EMP_earnings_in_output[1]
share_EMP_earnings_in_output_fixed_occ[1]
