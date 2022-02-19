using JLD2
using Plots

using ProgressMeter

include("Functions/profit.jl")

country = "Italy"
TIME_PERIODS = 30#100
LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_big_grid_test/$(country)/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_big_grid_test\\$(country)\\"
end
mkpath(LOCAL_DIR)

@load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2

approx_object = ss_1[4]
global_approx_params = [69,3,3,3,6,3]
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

function create_plot(X,XLABEL::String,Y,YLABEL::String,Y1,Y2, IS_Y_PERCENTAGE::Bool=false, OCCUPATION::String="H")
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    end

    YLIMS1 = minimum([Y; Y1; Y2])
    YLIMS2 = maximum([Y; Y1; Y2])
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
    plt = plot(collect([0; X]), collect.([[Y1; Y],ones(length(Y)+1).*Y1,ones(length(Y)+1).*Y2]),
                    color=[COLOR "green" "red"],
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    yticks = YTICKS,
                    ylims = YLIMS )
    if IS_Y_PERCENTAGE
        TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1*100;digits=DIGITS+1))%)"
    else
        TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1;digits=DIGITS+1)))"
    end
    annotate!([X[end]], Y1+YLIMMARGIN*1.25, text(TEXT1, :green, :right, 7))

    if IS_Y_PERCENTAGE
        TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2*100;digits=DIGITS+1))%)"
    else
        TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2;digits=DIGITS+1)))"
    end
    annotate!([X[end]], Y2+YLIMMARGIN*1.25, text(TEXT2, :red, :right, 7))
    return plt
end

plt = create_plot(collect(1:T),"Time (Years)", lambda_s,"λ", lambda_1, lambda_2)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_lambdas.png")

r_1 = ss_1[2]
r_2 = ss_2[2]
r_s = trans_res[3]
if length(r_s) != T
    r_s = r_s[1:T]
end
plt = create_plot(collect(1:T),"Time (Years)", r_s,"Interest rate", r_1, r_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_interest_rates.png")

w_1 = ss_1[3]
w_2 = ss_2[3]
w_s = trans_res[4]
plt = create_plot(collect(1:T),"Time (Years)", w_s,"Wage", w_1, w_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_wages.png")

# initialisation
Output_1 = sum(ss_1[1][32].*ss_1[1][5])
Output_2 = sum(ss_2[1][32].*ss_2[1][5])
Output_s = zeros(T)
Credit_to_Output_1 = ss_1[1][13]
Credit_to_Output_2 = ss_2[1][13]
Credit_to_Output_s = zeros(T)
Credit_1 = Credit_to_Output_1*Output_1
Credit_2 = Credit_to_Output_2*Output_2
Credit_s = zeros(T)
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
unlimited_capital_choice_1 = ss_1[1][26].<(ss_1[5][1].*(asset_grid.*ones(size(ss_1[1][26]))))
unlimited_capital_choice_2 = ss_2[1][26].<(ss_2[5][1].*(asset_grid.*ones(size(ss_2[1][26]))))
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

means_s = zeros(4,5,T)
ginis_s = zeros(4,5,T)
avglogs_s = zeros(4,5,T)
varlogs_s = zeros(4,5,T)
avgs_s = zeros(4,5,T)
vars_s = zeros(4,5,T)

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

# loop for _1 and _2
for i = 1:2
    ss = ss_1
    if i == 2
        ss = ss_2
    end

    lambda = ss[5][1]
    density_distr = ss[1][5]
    output = ss[1][32]
    capital_d = ss[1][26]
    labour_d = ss[1][29]

    eta = ss[5][5]
    theta = ss[5][6]

    z_m_nodes = ss[4][1]
    z_w_nodes = ss[4][2]

    unlimited_capital_choice = capital_d.<(lambda.*(asset_grid.*ones(size(capital_d))))

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

        for s = 1:4 # [1] = income, earnings, wealth, consumption
            # if s == 1
            stat_distr = ss[1][23] .- ones(size(ss[1][5])).*asset_grid
            if s == 2
                stat_distr = ss[1][24]
            elseif s == 3
                stat_distr = ones(size(ss[1][5])).*asset_grid
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

            # calculate variance of log-s
            avglogs[s,h,i] = sum(density_distr.*choice.*log.(max.(1e-12,stat_distr))   )
            varlogs[s,h,i] = sum(density_distr.*choice.*(log.(max.(1e-12,stat_distr)).- avglogs[s,h,i]).^2)/sum(density_distr.*choice)
            avglogs[s,h,i] /= sum(density_distr.*choice)

            avgs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
            vars[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avgs[s,h,i]).^2)/sum(density_distr.*choice)
            avgs[s,h,i] /= sum(density_distr.*choice)
        end

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
        end

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

        end
    end
end

#main computation loop
p = Progress(T, dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=10)
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

    Credit_s[t] = sum(credit .* capital_s_distr_s[t,:,:,:,:,:])

    Credit_to_Output_s[t] = Credit_s[t]/Output_s[t]

    Income_s[t] = sum((income.- ones(size(capital_s_distr_s[t,:,:,:,:,:])).*asset_grid) .* capital_s_distr_s[t,:,:,:,:,:])

    Consumption_s[t] = sum((income.-policy_s[t,:,:,:,:,:]) .* capital_s_distr_s[t,:,:,:,:,:] )

    #occupation
    occ_Ws_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==1.0) )
    occ_SPs_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==2.0) )
    occ_EMPs_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice.==3.0) )
    occ_ENTs_s[t] = occ_SPs_s[t]+occ_EMPs_s[t]

    #share of unconstrained ENT, SP,EMP
    unlimited_capital_choice_s[t,:,:,:,:,:] = capital_d.<(lambda_s[t].*(asset_grid.*ones(size(capital_d))))
    share_ENT_unbound_s[t] = sum( Float64.(occ_choice.!=1.0) .* capital_s_distr_s[t,:,:,:,:,:] .* unlimited_capital_choice_s[t,:,:,:,:,:] )/occ_ENTs_s[t]
    share_SP_unbound_s[t] = sum( Float64.(occ_choice.==2.0) .* capital_s_distr_s[t,:,:,:,:,:] .* unlimited_capital_choice_s[t,:,:,:,:,:] )/occ_SPs_s[t]
    share_EMP_unbound_s[t] = sum( Float64.(occ_choice.==3.0) .* capital_s_distr_s[t,:,:,:,:,:] .* unlimited_capital_choice_s[t,:,:,:,:,:] )/occ_EMPs_s[t]

    #Variances of log-consumption and log-earnings, Gini for W and ENT income
    logcs_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,income.-policy_s[t,:,:,:,:,:])).^2 ) .- sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,income.-policy_s[t,:,:,:,:,:])) )^2
    loges_s[t] = sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,earnings)).^2 ) - sum( capital_s_distr_s[t,:,:,:,:,:].*log.(max.(1e-12,earnings)) )^2

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

            # calculate variance of log-s
            avglogs_s[s,h,t] = sum(density_distr.*choice.*log.(max.(1e-12,stat_distr))   )
            varlogs_s[s,h,t] = sum(density_distr.*choice.*(log.(max.(1e-12,stat_distr)).- avglogs_s[s,h,t]).^2)/sum(density_distr.*choice)
            avglogs_s[s,h,t] /= sum(density_distr.*choice)

            avgs_s[s,h,t] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
            vars_s[s,h,t] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avgs_s[s,h,t]).^2)/sum(density_distr.*choice)
            avgs_s[s,h,t] /= sum(density_distr.*choice)
        end

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
        end

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

        end

    end

    next!(p)
end

#create plots

plt = create_plot(collect(1:T),"Time (Years)", Output_s,"Output", Output_1, Output_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_outputs.png")
plt = create_plot(collect(1:T),"Time (Years)", Credit_s,"Credit", Credit_1, Credit_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_credits.png")
plt = create_plot(collect(1:T),"Time (Years)", Credit_to_Output_s,"Credit/Output", Credit_to_Output_1, Credit_to_Output_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_credit_to_outputs.png")
plt = create_plot(collect(1:T),"Time (Years)", Income_s,"Income", Income_1, Income_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_incomes.png")
plt = create_plot(collect(1:T),"Time (Years)", Consumption_s,"Consumption", Consumption_1, Consumption_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)time_consumptions.png")


#Occupations
LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end
mkpath(LOCAL_DIR_OCCUPATION)
plt = create_plot(collect(1:T),"Time (Years)", occ_Ws_s,"Share of Workers", occ_Ws_1, occ_Ws_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_workers.png")
plt = create_plot(collect(1:T),"Time (Years)", occ_ENTs_s,"Share of Entrepreneurs", occ_ENTs_1, occ_ENTs_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_entrepreneurs.png")
plt = create_plot(collect(1:T),"Time (Years)", occ_SPs_s,"Share of Sole Proprietors", occ_SPs_1, occ_SPs_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_sole_proprietors.png")
plt = create_plot(collect(1:T),"Time (Years)", occ_EMPs_s,"Share of Employers", occ_EMPs_1, occ_EMPs_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_employers.png")

#Share of unconstrained ENT, SP,EMP
plt = create_plot(collect(1:T),"Time (Years)", share_ENT_unbound_s,"Share of Unconstrained Entrepreneurs", share_ENT_unbound_1, share_ENT_unbound_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_unconstrained_entrepreneurs.png")
plt = create_plot(collect(1:T),"Time (Years)", share_SP_unbound_s,"Share of Unconstrained Sole Proprietors", share_SP_unbound_1, share_SP_unbound_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_unconstrained_sole_proprietors.png")
plt = create_plot(collect(1:T),"Time (Years)", share_EMP_unbound_s,"Share of Unconstrained Employers", share_EMP_unbound_1, share_EMP_unbound_2, true)
display(plt)
savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_share_of_unconstrained_employers.png")


#Inequality
LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#Variance of log-consumption and log-earnings, Gini for W and ENT income
plt = create_plot(collect(1:T),"Time (Years)", logcs_s,"Variance of log-consumption", logcs_1, logcs_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_log_consumption.png")
plt = create_plot(collect(1:T),"Time (Years)", loges_s,"Variance of log-earnings", loges_1, loges_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_log_earnings.png")
plt = create_plot(collect(1:T),"Time (Years)", giniWs_s,"Gini for Workers' income", giniWs_1, giniWs_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_workers_income.png")
plt = create_plot(collect(1:T),"Time (Years)", giniEnts_s,"Gini for Entrepreneurs' income", giniEnts_1, giniEnts_2, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_entrepreneurs_income.png")

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

    for h = 1:5
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end

        # calculate mean
        #means[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", means_s[s,h,:],"Mean of $(CHOICE_NAMES[h])' $stat_name", means[s,h,1], means[s,h,2], false)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_mean_$(choice_name)_$(stat_name).png")

        # calculate gini coefficent
        #ginis[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", ginis_s[s,h,:],"Gini of $(CHOICE_NAMES[h])' $stat_name", ginis[s,h,1], ginis[s,h,2], false)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_gini_$(choice_name)_$(stat_name).png")

        # calculate variance of log-s
        #avglogs[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", avglogs_s[s,h,:],"Average of $(CHOICE_NAMES[h])' Log-$stat_name", avglogs[s,h,1], avglogs[s,h,2], false)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_avg_$(choice_name)_log$(stat_name).png")
        #varlogs[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", varlogs_s[s,h,:],"Variance of $(CHOICE_NAMES[h])' Log-$stat_name", varlogs[s,h,1], varlogs[s,h,2], false)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_$(choice_name)_log$(stat_name).png")

        # calculate variance of s
        #avgs[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", avgs_s[s,h,:],"Average of $(CHOICE_NAMES[h])' $stat_name", avgs[s,h,1], avgs[s,h,2], false)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_avg_$(choice_name)_$(stat_name).png")
        #vars[s,h,i]
        plt = create_plot(collect(1:T),"Time (Years)", vars_s[s,h,:],"Variance of $(CHOICE_NAMES[h])' $stat_name", vars[s,h,1], vars[s,h,2], false)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_var_$(choice_name)_$(stat_name).png")

    end

end


# productivity
LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)
# TFP_ideal, TFP_data, mean and var of MPK and MPL, share_of_unconstrained SP,EMP,ENT

CHOICE_NAMES = ["Entrepreneurs","Sole Proprietors","Employers"]
for h=1:3
    choice_name = CHOICE_NAMES[h]
    if h==3
        choice_name = "SoleProprietors"
    end
    plt = create_plot(collect(1:T),"Time (Years)", TFPis_s[h,:],"TFP_ideal for $(CHOICE_NAMES[h])", TFPis[h,1], TFPis[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_tfp_ideal_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", TFPds_s[h,:],"TFP_data for $(CHOICE_NAMES[h])", TFPds[h,1], TFPds[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_tfp_data_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", mean_MPL_s[h,:],"Mean of MPL for $(CHOICE_NAMES[h])", mean_MPL[h,1], mean_MPL[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_mean_mpl_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_MPL_s[h,:],"Variance of MPL for $(CHOICE_NAMES[h])", var_MPL[h,1], var_MPL[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_mpl_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", mean_MPK_s[h,:],"Mean of MPK for $(CHOICE_NAMES[h])", mean_MPK[h,1], mean_MPK[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_mean_mpk_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_MPK_s[h,:],"Variance of MPK for $(CHOICE_NAMES[h])", var_MPK[h,1], var_MPK[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_mpk_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", share_unbound_s[h,:],"Share of Unconstrained $(CHOICE_NAMES[h])", share_unbound[h,1], share_unbound[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_share_unconstrained_$(choice_name).png")

end

#skills
CHOICE_NAMES = ["Workers","Sole Proprietors","Employers","Entrepreneurs"]
for h=1:4
    choice_name = CHOICE_NAMES[h]
    if h==2
        choice_name = "SoleProprietors"
    end
    plt = create_plot(collect(1:T),"Time (Years)", avg_m_skill_s[h,:],"Average Managerial Skill across $(CHOICE_NAMES[h])", avg_m_skill[h,1], avg_m_skill[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_avg_m_skill_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", avg_w_skill_s[h,:],"Average Working Skill across $(CHOICE_NAMES[h])", avg_w_skill[h,1], avg_w_skill[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_avg_w_skill_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_m_skill_s[h,:],"Variance of Managerial Skill across $(CHOICE_NAMES[h])", var_m_skill[h,1], var_m_skill[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_m_skill_$(choice_name).png")

    plt = create_plot(collect(1:T),"Time (Years)", var_w_skill_s[h,:],"Variance of Working Skill across $(CHOICE_NAMES[h])", var_w_skill[h,1], var_w_skill[h,2], false)
    display(plt)
    savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)time_var_w_skill_$(choice_name).png")

end



#end
