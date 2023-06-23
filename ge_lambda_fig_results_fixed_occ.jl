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
Incomes = zeros(num_lambdas)
Consumptions = zeros(num_lambdas)

Rs = zeros(num_lambdas)
Ws = zeros(num_lambdas)

logcs = zeros(num_lambdas)
loges = zeros(num_lambdas)

giniWs = zeros(num_lambdas)
giniEnts = zeros(num_lambdas)

# ENT, SP, EMP
share_unbound = zeros(3,num_lambdas)

# income, earnings, wealth, consumption
# All, W, SP, EMP
# lambdas
means = zeros(4,4,num_lambdas)
ginis = zeros(4,4,num_lambdas)
avglogs = zeros(4,4,num_lambdas)
varlogs = zeros(4,4,num_lambdas)
avgs = zeros(4,4,num_lambdas)
vars = zeros(4,4,num_lambdas)

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
# 5 quantiles
# All, W, SP, EMP
# income, earnings, wealth, consumption
quantile_means = zeros(5,4,4,num_lambdas)

# All, W, SP, EMP
Capital = zeros(4,num_lambdas)

# ENT, SP, EMP
TFPis = zeros(3,num_lambdas)
TFPds = zeros(3,num_lambdas)

mean_MPL = zeros(3,num_lambdas)
var_MPL = zeros(3,num_lambdas)

mean_MPK = zeros(3,num_lambdas)
var_MPK = zeros(3,num_lambdas)

# Distribution of output to different channels of earnigs/income to occupations
share_W_earnings_in_output = zeros(num_lambdas)
share_SP_earnings_in_output = zeros(num_lambdas)
share_EMP_earnings_in_output = zeros(num_lambdas)
share_W_capital_income_in_output = zeros(num_lambdas)
share_SP_capital_income_in_output = zeros(num_lambdas)
share_EMP_capital_income_in_output = zeros(num_lambdas)
#throw(error)

@showprogress for i = 1:num_lambdas
    #println("\n$(country) - $(i)/$(num_lambdas)")
    # ss_star = [res, r, w, approx_object, params]
    ss_star = copy(SSS[i])

    number_u_nodes = SSS[i][4][7]
    number_asset_grid = ss_star[1][2]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

    eta = ss_star[5][5]
    theta = ss_star[5][6]

    z_m_nodes = ss_star[4][1]
    z_w_nodes = ss_star[4][2]

    asset_grid = ss_star[1][3]
    density_distr = ss_star[1][5]
    policy = ss_star[1][4]
    earnings = ss_star[1][24]
    capital_d = ss_star[1][26]
    labour_d = ss_star[1][29]
    output = ss_star[1][32]

    wealth = Array{Any}(undef,3)
    income = Array{Any}(undef,3)
    consumption = Array{Any}(undef,3)
    for occ in 1:3
        wealth[occ] = ones(size(density_distr[occ])).*asset_grid
        income[occ] = ss_star[1][23][occ] .- wealth[occ]
        consumption[occ] = income[occ] .+ wealth[occ] .- policy[occ]

        Capital[occ+1,i] = sum(wealth[occ] .* density_distr[occ])
    end

    C_Ys[i] = ss_star[1][13]
    Outputs[i] = sum( [sum(output[occ] .* density_distr[occ]) for occ in 1:3] )
    Incomes[i] = sum( [sum(income[occ] .* density_distr[occ]) for occ in 1:3] )
    Consumptions[i] = sum( [sum(consumption[occ] .* density_distr[occ]) for occ in 1:3] )
    Capital[1,i] = sum( Capital[2:4,i] )

    Rs[i] = ss_star[2]
    Ws[i] = ss_star[3]

    logcs[i] = ss_star[1][17]
    loges[i] = ss_star[1][19]

    giniWs[i] = ss_star[1][20]
    giniEnts[i] = ss_star[1][21]

    # [1] = income, earnings, wealth, consumption
    # [2] = All, W, SP, EMP
    # [3] = lambdas
    # = means, ginis, varlogs
    for s = 1:4 # [1] = income, earnings, wealth, consumption
        # if s == 1
        stat_distr = copy(income)
        if s == 2
            stat_distr = copy(earnings)
        elseif s == 3
            stat_distr = copy(wealth)
        elseif s == 4
            stat_distr = copy(consumption)
        end

        for h = 1:4 # [2] = All, W, SP, EMP
            temp_occ = 4
            if h == 1 #All
                # calculate mean
                means[s,h,i] = sum( [sum(stat_distr[occ] .* density_distr[occ]) for occ in 1:3])/sum( [sum(density_distr[occ]) for occ in 1:3])

                # calculate gini coefficent
                density_distr_choice = [density_distr[2]; density_distr[3]]./sum([density_distr[2]; density_distr[3]])
                stat_choice = [stat_distr[2]; stat_distr[3]]
                density_distr_choice_vec = vec(density_distr_choice)
                stat_choice_vec = vec(stat_choice)

                index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
                density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
                stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]

                ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))

                # calculate variance of log-s
                avglogs[s,h,i] = sum([ sum(density_distr[occ].*log.(max.(1e-12,stat_distr[occ]))   ) for occ in 1:3])
                varlogs[s,h,i] = sum([ sum(density_distr[occ].*(log.(max.(1e-12,stat_distr[occ])).-avglogs[s,h,i]).^2) for occ in 1:3])/sum([ sum(density_distr[occ]) for occ in 1:3])
                avglogs[s,h,i] /= sum([ sum(density_distr[occ]) for occ in 1:3])
                if s==3
                    avglogs[s,h,i] = sum([ sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   ) for occ in 1:3])
                    varlogs[s,h,i] = sum([ sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).-avglogs[s,h,i]).^2) for occ in 1:3])/sum([ sum(density_distr[occ]) for occ in 1:3])
                    avglogs[s,h,i] /= sum([ sum(density_distr[occ]) for occ in 1:3])
                end

                avgs[s,h,i] = sum([ sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   ) for occ in 1:3])
                vars[s,h,i] = sum([ sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).-avgs[s,h,i]).^2) for occ in 1:3])/sum([ sum(density_distr[occ]) for occ in 1:3])
                avgs[s,h,i] /= sum([ sum(density_distr[occ]) for occ in 1:3])

                try
                    quantile_means[:,h,s,i] .= quantile_mean([stat_distr[1];stat_distr[2];stat_distr[3]], [density_distr[1];density_distr[2];density_distr[3]])
                catch e
                    quantile_means[:,h,s,i] .= NaN
                end
            elseif h == 2 #W
                temp_occ = 1
            elseif h == 3 #SP
                temp_occ = 2
            elseif h == 4 #EMP
                temp_occ = 3
            end
            occ = temp_occ

            if h!=1
                # calculate mean
                means[s,h,i] = sum(stat_distr[occ] .* density_distr[occ])/sum(density_distr[occ])

                # calculate gini coefficent
                density_distr_choice = (density_distr[occ])./sum(density_distr[occ])
                stat_choice = stat_distr[occ]
                density_distr_choice_vec = vec(density_distr_choice)
                stat_choice_vec = vec(stat_choice)

                index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
                density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
                stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]

                ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))

                # calculate variance of log-s
                avglogs[s,h,i] = sum(density_distr[occ].*log.(max.(1e-12,stat_distr[occ]))   )
                varlogs[s,h,i] = sum(density_distr[occ].*(log.(max.(1e-12,stat_distr[occ])).-avglogs[s,h,i]).^2)/sum(density_distr[occ])
                avglogs[s,h,i] /= sum(density_distr[occ])
                if s==3
                    avglogs[s,h,i] = sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   )
                    varlogs[s,h,i] = sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).- avglogs[s,h,i]).^2)/sum(density_distr[occ])
                    avglogs[s,h,i] /= sum(density_distr[occ])
                end

                avgs[s,h,i] = sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   )
                vars[s,h,i] = sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).- avgs[s,h,i]).^2)/sum(density_distr[occ])
                avgs[s,h,i] /= sum(density_distr[occ])

                try
                    quantile_means[:,h,s,i] .= quantile_mean(stat_distr[occ], density_distr[occ])
                catch e
                    quantile_means[:,h,s,i] .= NaN
                end
            end

        end
    end

    unlimited_capital_choice = Array{Any}(undef,3)
    for occ in 1:3
        unlimited_capital_choice[occ] = capital_d[occ].<(lambdas[i].*(asset_grid.*ones(size(capital_d[occ]))))
    end

    for c = 1:4
        temp_occ = 4
        if c==1
            temp_occ = 1
        elseif c==2
            denumerator = sum([ sum(density_distr[occ]) for occ in 1:3])

            Y = sum([ sum(density_distr[occ].*output[occ]) for occ in 1:3])/denumerator
            K = sum([ sum(density_distr[occ].*capital_d[occ]) for occ in 1:3])/denumerator
            L = sum([ sum(density_distr[occ].*labour_d[occ]) for occ in 1:3])/denumerator
            TFPis[c-1,i] = Y/(K^eta*L^theta)
            TFPds[c-1,i] = Y/K^(eta/(theta+eta))

            mean_MPK[c-1,i] = sum([ sum(density_distr[occ].*(eta.*output[occ].*replace(capital_d[occ].^(-1.0),Inf=>0.0))) for occ in 1:3])/denumerator
            var_MPK[c-1,i] = sum([ sum(density_distr[occ].*(eta.*output[occ].*replace(capital_d[occ].^(-1.0),Inf=>0.0) .- mean_MPK[c-1,i]*denumerator).^2) for occ in 1:3])/denumerator

            mean_MPL[c-1,i] = sum([ sum(density_distr[occ].*(theta.*output[occ].*replace(labour_d[occ].^(-1.0),Inf=>0.0))) for occ in 1:3])/denumerator
            var_MPL[c-1,i] = sum([ sum(density_distr[occ].*(theta.*output[occ].*replace(labour_d[occ].^(-1.0),Inf=>0.0) .- mean_MPL[c-1,i]*denumerator).^2) for occ in 1:3])/denumerator

            #if c==2
            #    choice = ent_choice
            #end
            share_unbound[c-1,i] = sum([ sum(density_distr[occ].*unlimited_capital_choice[occ]) for occ in 1:3])/denumerator


        elseif c==3
            temp_occ = 2
        elseif c==4
            temp_occ = 3
        else
            throw("c is out of bounds")
        end
        occ = temp_occ

        if c!=2
            denumerator = sum(density_distr[occ])
            if c>1
                #if c==2
                #    choice = ones(size(choice))
                #end
                #if c==2
                #    denumerator = sum(choice.*density_distr)
                #end
                Y = sum(density_distr[occ].*output[occ])/denumerator
                K = sum(density_distr[occ].*capital_d[occ])/denumerator
                L = sum(density_distr[occ].*labour_d[occ])/denumerator
                TFPis[c-1,i] = Y/(K^eta*L^theta)
                TFPds[c-1,i] = Y/K^(eta/(theta+eta))

                c_1 = capital_d[occ].^(-1.0)
                replace!(c_1,Inf=>0.0)
                mean_MPK[c-1,i] = sum(density_distr[occ].*(eta.*output[occ].*c_1))/denumerator
                var_MPK[c-1,i] = sum(density_distr[occ].*(eta.*output[occ].*c_1 .- mean_MPK[c-1,i]*denumerator).^2)/denumerator

                l_1 = labour_d[occ].^(-1.0)
                replace!(l_1,Inf=>0.0)
                mean_MPL[c-1,i] = sum(density_distr[occ].*(theta.*output[occ].*l_1))/denumerator
                var_MPL[c-1,i] = sum(density_distr[occ].*(theta.*output[occ].*l_1 .- mean_MPL[c-1,i]*denumerator).^2)/denumerator

                #if c==2
                #    choice = ent_choice
                #end
                share_unbound[c-1,i] = sum(density_distr[occ].*unlimited_capital_choice[occ])/denumerator
            end
        end

    end

    agg_c_e = ss_star[5][7]*sum(density_distr[3])
    share_W_earnings_in_output[i] = sum(ss_star[1][24][1] .* density_distr[1])/(Outputs[i]-agg_c_e)
    share_SP_earnings_in_output[i] = sum(ss_star[1][24][2] .* density_distr[2])/(Outputs[i]-agg_c_e)
    share_EMP_earnings_in_output[i] = sum(ss_star[1][24][3] .* density_distr[3])/(Outputs[i]-agg_c_e)
    share_W_capital_income_in_output[i] = (ss_star[5][3]+ss_star[1][44])*sum(ones(size(ss_star[1][5][1])).*ss_star[1][3] .* density_distr[1])/(Outputs[i]-agg_c_e)
    share_SP_capital_income_in_output[i] = (ss_star[5][3]+ss_star[1][44])*sum(ones(size(ss_star[1][5][2])).*ss_star[1][3] .* density_distr[2])/(Outputs[i]-agg_c_e)
    share_EMP_capital_income_in_output[i] = (ss_star[5][3]+ss_star[1][44])*sum(ones(size(ss_star[1][5][3])).*ss_star[1][3] .* density_distr[3])/(Outputs[i]-agg_c_e)
    #throw(error)
end

#throw(error)

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

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)

plt = create_plot(lambdas,"λ", C_Ys,"Credit/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_credit_to_gdp.png")
plt = create_plot(lambdas,"λ", Outputs,"Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Outputs.png")
plt = create_plot(lambdas,"λ", Incomes,"Income", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Incomes.png")
plt = create_plot(lambdas,"λ", Consumptions,"Consumption", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Consumptions.png")

plt = create_plot(lambdas,"λ", Capital[1,:],"Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals.png")

plt = create_plot(lambdas,"λ", Capital[2,:],"Workers' Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals_W.png")
plt = create_plot(lambdas,"λ", Capital[3,:],"Sole Proprietors' Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals_SP.png")
plt = create_plot(lambdas,"λ", Capital[4,:],"Employers' Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals_EMP.png")

#Interest rate and wage
plt = create_plot(C_Ys,"Credit/Output", Rs,"Interest Rate")
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_rate.png")
plt = create_plot(C_Ys,"Credit/Output", Ws,"Wage", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_wage.png")


LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#Credit/Output to: log-consumption, log-earnings
plt = create_plot(C_Ys,"Credit/Output", logcs,"Variance of log-consumption", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_log_consumption.png")

plt = create_plot(C_Ys,"Credit/Output", loges,"Variance of log-earnings", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_log_earnings.png")

#Crdit/Output to Gini Workers and Entrepreneurs
plts = create_plots(C_Ys,"Credit/Output", [giniWs, giniEnts],["Gini for workers' income","Gini for entrepreneurs' income"],false,["W","ENT"])
for i in 1:length(plts)
    display(plts[i])
    if i==1
        savefig(plts[i],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_workers.png")
    elseif i==2
        savefig(plts[i],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_entrepreneurs.png")
    end
end

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
            LABEL = "Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles)"
            plt = create_combined_plot(C_Ys,"Credit/Output", [quantile_means[1,h,s,:],quantile_means[2,h,s,:],quantile_means[3,h,s,:],quantile_means[4,h,s,:],quantile_means[5,h,s,:]],LABELS,LABEL, false, ["H","W","SP","EMP","ENT"])
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_$(OCCS[h])_quantiles.png")
        catch e

        end
    end
    # calculate mean
    #means[s,h,i]
    LABELS = ["Mean of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [means[s,1,:], means[s,2,:], means[s,3,:], means[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", means[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name).png")
    end

    # calculate gini coefficent
    #ginis[s,h,i]
    LABELS = ["Gini of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], ginis[s,3,:], ginis[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", ginis[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name).png")
    end
    # calculate variance of log-s
    #avgs[s,h,i]
    LABELS = ["Average of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [avgs[s,1,:], avgs[s,2,:], avgs[s,3,:], avgs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avgs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name).png")
    end

    #vars[s,h,i]
    LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [vars[s,1,:], vars[s,2,:], vars[s,3,:], vars[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", vars[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name).png")
    end

    # calculate variance of log-s
    #avglogs[s,h,i]
    LABELS = ["Average of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], avglogs[s,3,:], avglogs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avglogs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end

    #varlogs[s,h,i]
    if s!=3
        LABELS = ["Variance of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    else
        LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    end
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], varlogs[s,3,:], varlogs[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", varlogs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
end

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)

# TFP_ideal for ENT
plt = create_plot(C_Ys,"Credit/Output",TFPis[1,:],"TFP for Entrepreneurs", false, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_ent.png")
# TFP_ideal for SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [TFPis[2,:], TFPis[3,:]],["TFP for Sole Proprietors","TFP for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_emp.png")
# TFP_data for SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [TFPds[2,:], TFPds[3,:]],["TFP_data for Sole Proprietors","TFP_data for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data_emp.png")


# mean of MPL for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",mean_MPL[1,:],"Mean of MPL for Entrepreneurs", false, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [mean_MPL[2,:], mean_MPL[3,:]],["Mean of MPL for Sole Proprietors","Mean of MPL for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [mean_MPL[2,:], mean_MPL[3,:]],LABELS,"Mean of MPL",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_mean_mpl.png")

# variance of MPL for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",var_MPL[1,:],"Variance of MPL for Entrepreneurs", false, "ENT", 0.0)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [var_MPL[2,:], var_MPL[3,:]],["Variance of MPL for Sole Proprietors","Variance of MPL for Employers"],false,["SP","EMP"],0.0)
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_MPL[2,:], var_MPL[3,:]],LABELS,"Variance of MPL",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_var_mpl.png")

# mean of MPK for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",mean_MPK[1,:],"Mean of MPK for Entrepreneurs", false, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [mean_MPK[2,:], mean_MPK[3,:]],["Mean of MPK for Sole Proprietors","Mean of MPK for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [mean_MPK[2,:], mean_MPK[3,:]],LABELS,"Mean of MPK",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_mean_mpk.png")

# variance of MPK for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",var_MPK[1,:],"Variance of MPK for Entrepreneurs", false, "ENT", 0.0)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [var_MPK[2,:], var_MPK[3,:]],["Variance of MPK for Sole Proprietors","Variance of MPK for Employers"],false,["SP","EMP"],0.0)
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_MPK[2,:], var_MPK[3,:]],LABELS,"Variance of MPK",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_var_mpk.png")

# share of unbound ENT, SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [share_unbound[1,:], share_unbound[2,:], share_unbound[3,:]],["Share of Unconstrained Entrepreneurs","Share of Unconstrained Sole Proprietors","Share of Unconstrained Employers"],true,["ENT","SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_ent_unbound.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_se_unbound.png")
savefig(plts[3],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_emp_unbound.png")

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

# Print statistics for calibrated economy

# income, earnings, wealth, consumption
# All, W, SP, EMP
# lambdas
#means = zeros(4,4,num_lambdas)
#ginis = zeros(4,4,num_lambdas)
#vars = zeros(4,4,num_lambdas)
for i=1:4

    if i==1
        stat_name="income"
    elseif i==2
        stat_name="earnings"
    elseif i==3
        stat_name="wealth"
    else  #i==4
        stat_name="consumption"
    end

    for j=[4]#1:4
        if j==1
            stat_type_name="Mean"
            stat_type = means
        elseif j==2
            stat_type_name="Variance"
            stat_type = vars
        elseif j==3
            stat_type_name="Varlog"
            stat_type = varlogs
        else  #i==4
            stat_type_name="Gini"
            stat_type = ginis
        end
        println("$(stat_type_name) of $(stat_name) for H,W,SP,EMP: $(round.(stat_type[i,:,calibrated_lambda];digits=2))")
    end
end

ss_star = SSS[10]
asset_grid = ss_star[1][3]
density_distr = ss_star[1][5]
policy = ss_star[1][4]
earnings = ss_star[1][24]

wealth = Array{Any}(undef,3)
income = Array{Any}(undef,3)
consumption = Array{Any}(undef,3)
for occ in 1:3
    wealth[occ] = ones(size(density_distr[occ])).*asset_grid
    income[occ] = ss_star[1][23][occ] .- wealth[occ]
    consumption[occ] = income[occ] .+ wealth[occ] .- policy[occ]
end

num_of_occs = 4

mean_consumption_to_income_wealth = Array{Any}(undef,num_of_occs)
var_consumption_to_income_wealth = Array{Any}(undef,num_of_occs)
std_consumption_to_income_wealth = Array{Any}(undef,num_of_occs)

mean_income_to_income_wealth = Array{Any}(undef,num_of_occs)
var_income_to_income_wealth = Array{Any}(undef,num_of_occs)
std_income_to_income_wealth = Array{Any}(undef,num_of_occs)

mean_consumption_to_income = Array{Any}(undef,num_of_occs)
var_consumption_to_income = Array{Any}(undef,num_of_occs)
std_consumption_to_income = Array{Any}(undef,num_of_occs)

mean_savings_to_income = Array{Any}(undef,num_of_occs)
var_savings_to_income = Array{Any}(undef,num_of_occs)
std_savings_to_income = Array{Any}(undef,num_of_occs)

mean_consumption_to_earnings = Array{Any}(undef,num_of_occs)
var_consumption_to_earnings = Array{Any}(undef,num_of_occs)
std_consumption_to_earnings = Array{Any}(undef,num_of_occs)

mean_savings_to_earnings = Array{Any}(undef,num_of_occs)
var_savings_to_earnings = Array{Any}(undef,num_of_occs)
std_savings_to_earnings = Array{Any}(undef,num_of_occs)

for i = 1:num_of_occs # H. W, SP, EMP
    if i == 1
        mean_consumption_to_income_wealth[i] = sum( [sum((consumption[occ]./ss_star[1][23][occ]).*density_distr[occ]) for occ=1:3])
        var_consumption_to_income_wealth[i] = sum( [sum((consumption[occ]./ss_star[1][23][occ] .- mean_consumption_to_income_wealth[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_consumption_to_income_wealth[i] = sqrt(var_consumption_to_income_wealth[i])


        mean_income_to_income_wealth[i] = sum( [sum((income[occ]./ss_star[1][23][occ]).*density_distr[occ]) for occ=1:3])
        var_income_to_income_wealth[i] = sum( [sum((income[occ]./ss_star[1][23][occ] .- mean_income_to_income_wealth[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_income_to_income_wealth[i] = sqrt(var_income_to_income_wealth[i])


        mean_consumption_to_income[i] = sum( [sum((consumption[occ]./income[occ]).*density_distr[occ]) for occ=1:3])
        var_consumption_to_income[i] = sum( [sum((consumption[occ]./income[occ] .- mean_consumption_to_income[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_consumption_to_income[i] = sqrt(var_consumption_to_income[i])

        mean_savings_to_income[i] = sum( [sum(((policy[occ].-wealth[occ])./income[occ]).*density_distr[occ]) for occ=1:3])
        var_savings_to_income[i] = sum( [sum(((policy[occ].-wealth[occ])./income[occ] .- mean_savings_to_income[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_savings_to_income[i] = sqrt(var_savings_to_income[i])

        mean_consumption_to_earnings[i] = sum( [sum((consumption[occ]./earnings[occ]).*density_distr[occ]) for occ=1:3])
        var_consumption_to_earnings[i] = sum( [sum((consumption[occ]./earnings[occ] .- mean_consumption_to_earnings[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_consumption_to_earnings[i] = sqrt(var_consumption_to_earnings[i])

        mean_savings_to_earnings[i] = sum( [sum(((policy[occ].-wealth[occ])./earnings[occ]).*density_distr[occ]) for occ=1:3])
        var_savings_to_earnings[i] = sum( [sum(((policy[occ].-wealth[occ])./earnings[occ] .- mean_savings_to_earnings[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_savings_to_earnings[i] = sqrt(var_savings_to_earnings[i])
    else

        if i==2
            choice = 1
        elseif i==3
            choice = 2
        elseif i==4
            choice = 3
        end

        mean_consumption_to_income_wealth[i] = sum((consumption[choice]./ss_star[1][23][choice]).*density_distr[choice])
        var_consumption_to_income_wealth[i] = sum((consumption[choice]./ss_star[1][23][choice] .- mean_consumption_to_income_wealth[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_consumption_to_income_wealth[i] /= sum(density_distr[choice])
        std_consumption_to_income_wealth[i] = sqrt(var_consumption_to_income_wealth[i])


        mean_income_to_income_wealth[i] = sum((income[choice]./ss_star[1][23][choice]).*density_distr[choice])
        var_income_to_income_wealth[i] = sum((income[choice]./ss_star[1][23][choice] .- mean_income_to_income_wealth[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_income_to_income_wealth[i] /= sum(density_distr[choice])
        std_income_to_income_wealth[i] = sqrt(var_income_to_income_wealth[i])


        mean_consumption_to_income[i] = sum((consumption[choice]./income[choice]).*density_distr[choice])
        var_consumption_to_income[i] = sum((consumption[choice]./income[choice] .- mean_consumption_to_income[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_consumption_to_income[i] /= sum(density_distr[choice])
        std_consumption_to_income[i] = sqrt(var_consumption_to_income[i])

        mean_savings_to_income[i] = sum(((policy[choice].-wealth[choice])./income[choice]).*density_distr[choice])
        var_savings_to_income[i] = sum(((policy[choice].-wealth[choice])./income[choice] .- mean_savings_to_income[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_savings_to_income[i] /= sum(density_distr[choice])
        std_savings_to_income[i] = sqrt(var_savings_to_income[i])

        mean_consumption_to_earnings[i] = sum((consumption[choice]./earnings[choice]).*density_distr[choice])
        var_consumption_to_earnings[i] = sum((consumption[choice]./earnings[choice] .- mean_consumption_to_earnings[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_consumption_to_earnings[i] /= sum(density_distr[choice])
        std_consumption_to_earnings[i] = sqrt(var_consumption_to_earnings[i])

        mean_savings_to_earnings[i] = sum(((policy[choice].-wealth[choice])./earnings[choice]).*density_distr[choice])
        var_savings_to_earnings[i] = sum(((policy[choice].-wealth[choice])./earnings[choice] .- mean_savings_to_earnings[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_savings_to_earnings[i] /= sum(density_distr[choice])
        std_savings_to_earnings[i] = sqrt(var_savings_to_earnings[i])
    end
end

println("Mean of the ratio of C/(Income+Wealth) for H,W,SP,EMP: $(round.(mean_consumption_to_income_wealth;digits=2))")
println("Std of the ratio of C/(Income+Wealth) for H,W,SP,EMP: $(round.(std_consumption_to_income_wealth;digits=2))")

println("Mean of the ratio of Income/(Income+Wealth) for H,W,SP,EMP: $(round.(mean_income_to_income_wealth;digits=2))")
println("Std of the ratio of Income/(Income+Wealth) for H,W,SP,EMP: $(round.(std_income_to_income_wealth;digits=2))")

println("Mean of the ratio of C/Income for H,W,SP,EMP: $(round.(mean_consumption_to_income;digits=2))")
println("Std of the ratio of C/Income for H,W,SP,EMP: $(round.(std_consumption_to_income;digits=2))")

println("Mean of the ratio of Savings/Income for H,W,SP,EMP: $(round.(mean_savings_to_income;digits=2))")
println("Std of the ratio of Savings/Income for H,W,SP,EMP: $(round.(std_savings_to_income;digits=2))")

println("Mean of the ratio of C/Earnings for H,W,SP,EMP: $(round.(mean_consumption_to_earnings;digits=2))")
println("Std of the ratio of C/Earnings for H,W,SP,EMP: $(round.(std_consumption_to_earnings;digits=2))")

println("Mean of the ratio of Savings/Earnings for H,W,SP,EMP: $(round.(mean_savings_to_earnings;digits=2))")
println("Std of the ratio of Savings/Earnings for H,W,SP,EMP: $(round.(std_savings_to_earnings;digits=2))")


SSS_fixed_occ = copy(SSS)
C_Ys_fixed_occ = copy(C_Ys)
Outputs_fixed_occ = copy(Outputs)
Incomes_fixed_occ = copy(Incomes)
Consumptions_fixed_occ = copy(Consumptions)
Rs_fixed_occ = copy(Rs)
Ws_fixed_occ = copy(Ws)
logcs_fixed_occ = copy(logcs)
loges_fixed_occ = copy(loges)
giniWs_fixed_occ = copy(giniWs)
giniEnts_fixed_occ = copy(giniEnts)
share_unbound_fixed_occ = copy(share_unbound)
means_fixed_occ = copy(means)
ginis_fixed_occ = copy(ginis)
avglogs_fixed_occ = copy(avglogs)
varlogs_fixed_occ = copy(varlogs)
avgs_fixed_occ = copy(avgs)
vars_fixed_occ = copy(vars)
quantile_means_fixed_occ = copy(quantile_means)
TFPis_fixed_occ = copy(TFPis)
TFPds_fixed_occ = copy(TFPds)
mean_MPL_fixed_occ = copy(mean_MPL)
var_MPL_fixed_occ = copy(var_MPL)
mean_MPK_fixed_occ = copy(mean_MPK)
var_MPK_fixed_occ = copy(var_MPK)
Capital_fixed_occ = copy(Capital)
share_W_earnings_in_output_fixed_occ = copy(share_W_earnings_in_output)
share_SP_earnings_in_output_fixed_occ = copy(share_SP_earnings_in_output)
share_EMP_earnings_in_output_fixed_occ = copy(share_EMP_earnings_in_output)
share_W_capital_income_in_output_fixed_occ = copy(share_W_capital_income_in_output)
share_SP_capital_income_in_output_fixed_occ = copy(share_SP_capital_income_in_output)
share_EMP_capital_income_in_output_fixed_occ = copy(share_EMP_capital_income_in_output)

@save "$(LOCAL_DIR_GENERAL)SSS_fixed.jld2" SSS_fixed_occ C_Ys_fixed_occ Outputs_fixed_occ Incomes_fixed_occ Consumptions_fixed_occ Rs_fixed_occ Ws_fixed_occ logcs_fixed_occ loges_fixed_occ giniWs_fixed_occ giniEnts_fixed_occ share_unbound_fixed_occ means_fixed_occ ginis_fixed_occ avglogs_fixed_occ varlogs_fixed_occ avgs_fixed_occ vars_fixed_occ quantile_means_fixed_occ TFPis_fixed_occ TFPds_fixed_occ mean_MPL_fixed_occ var_MPL_fixed_occ mean_MPK_fixed_occ var_MPK_fixed_occ Capital_fixed_occ share_W_earnings_in_output_fixed_occ share_SP_earnings_in_output_fixed_occ share_EMP_earnings_in_output_fixed_occ share_W_capital_income_in_output_fixed_occ share_SP_capital_income_in_output_fixed_occ share_EMP_capital_income_in_output_fixed_occ
