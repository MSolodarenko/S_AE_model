include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
#using PyPlot

using SchumakerSpline
using Statistics
using JLD2

country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Lambda_grid/$(country)/"
if Sys.iswindows()
    LOCAL_DIR = "\\Results\\Lambda_grid\\$(country)\\"
end

print_sameline("Loading data from Lambda_grid")
@load "$(LOCAL_DIR)SSS.jld2" SSS
lambdas = zeros(length(SSS))

calibrated_lambda = findfirst(x -> x==1.513028, lambdas)

Rs = zeros(length(lambdas))
Ws = zeros(length(lambdas))
K_Ys = zeros(length(lambdas))
C_Ys = zeros(length(lambdas))
occWs = zeros(length(lambdas))
occSEs = zeros(length(lambdas))
occEMPs = zeros(length(lambdas))
occENTs = zeros(length(lambdas))
logcs = zeros(length(lambdas))
loges = zeros(length(lambdas))
giniWs = zeros(length(lambdas))
giniEnts = zeros(length(lambdas))

TFPis = zeros(length(lambdas))
TFPds = zeros(length(lambdas))
Labour_productivity_s = zeros(length(lambdas))
Capital_productivity_s = zeros(length(lambdas))
Investment_productivity_s = zeros(length(lambdas))

Outputs = zeros(length(lambdas))
Incomes = zeros(length(lambdas))
Consumptions = zeros(length(lambdas))
Income_to_outputs = zeros(length(lambdas))

# income, earnings, wealth, consumption
# All, W, SP, EMP
# lambdas
means = zeros(4,4,length(lambdas))
ginis = zeros(4,4,length(lambdas))
avglogs = zeros(4,4,length(lambdas))
varlogs = zeros(4,4,length(lambdas))
avgs = zeros(4,4,length(lambdas))
vars = zeros(4,4,length(lambdas))

output_SPs  = zeros(length(lambdas))
output_EMPs = zeros(length(lambdas))

share_of_output_SPs  = zeros(length(lambdas))
share_of_output_EMPs = zeros(length(lambdas))

ratio_of_capitals = zeros(length(lambdas))
ratio_of_outputs  = zeros(length(lambdas))

interquartile_wealth = zeros(length(lambdas))
interquintile_wealth = zeros(length(lambdas))
interdecile_wealth = zeros(length(lambdas))
interquartile_wealth_w = zeros(length(lambdas))
interquintile_wealth_w = zeros(length(lambdas))
interdecile_wealth_w = zeros(length(lambdas))
interquartile_wealth_se = zeros(length(lambdas))
interquintile_wealth_se = zeros(length(lambdas))
interdecile_wealth_se = zeros(length(lambdas))
interquartile_wealth_emp = zeros(length(lambdas))
interquintile_wealth_emp = zeros(length(lambdas))
interdecile_wealth_emp = zeros(length(lambdas))

palma_wealth = zeros(length(lambdas))
palma_wealth_w = zeros(length(lambdas))
palma_wealth_se = zeros(length(lambdas))
palma_wealth_emp = zeros(length(lambdas))

ratio2020 = zeros(length(lambdas))
ratio2020_w = zeros(length(lambdas))
ratio2020_se = zeros(length(lambdas))
ratio2020_emp = zeros(length(lambdas))

gini_wealth = zeros(length(lambdas))
gini_wealth_w = zeros(length(lambdas))
gini_wealth_se = zeros(length(lambdas))
gini_wealth_emp = zeros(length(lambdas))

Capital = zeros(length(lambdas))

share_ent_unbound = zeros(length(lambdas))
share_se_unbound = zeros(length(lambdas))
share_emp_unbound = zeros(length(lambdas))

# H, W,SP,EMP
avg_w_skill = zeros(4,length(lambdas))
avg_m_skill = zeros(4,length(lambdas))

var_w_skill = zeros(4,length(lambdas))
var_m_skill = zeros(4,length(lambdas))

#throw(error)

for i = 1:length(lambdas)
    println("\n$(country) - $(i)/$(length(lambdas))")
    # ss_star = [res, r, w, approx_object, params]
    ss_star = copy(SSS[i])

    Rs[i] = ss_star[2]
    Ws[i] = ss_star[3]
    K_Ys[i] = ss_star[1][12]
    C_Ys[i] = ss_star[1][13]
    occWs[i] = ss_star[1][14]
    occSEs[i] = ss_star[1][15]
    occEMPs[i] = ss_star[1][16]
    occENTs[i] = occSEs[i] + occEMPs[i]
    logcs[i] = ss_star[1][17]
    loges[i] = ss_star[1][19]
    giniWs[i] = ss_star[1][20]
    giniEnts[i] = ss_star[1][21]

    TFPis[i] = ss_star[1][43][end-4]
    TFPds[i] = ss_star[1][43][end-3]
    Labour_productivity_s[i] = ss_star[1][43][end-2]
    Capital_productivity_s[i] = ss_star[1][43][end-1]
    Investment_productivity_s[i] = ss_star[1][43][end]

    z_m_nodes = ss_star[4][1]
    z_w_nodes = ss_star[4][2]

    asset_grid = ss_star[1][3]
    density_distr = ss_star[1][5]
    wealth = ones(size(density_distr)).*asset_grid
    policy = ss_star[1][4]
    income = ss_star[1][23] .- wealth
    earnings = ss_star[1][24]
    consumption = income .+ wealth .- policy
    occ_choice = ss_star[1][22]
    capital_d = ss_star[1][26]

    output = ss_star[1][32]

    Capital[i] = sum(asset_grid.*density_distr)

    Outputs[i] = sum(output .* density_distr)
    Incomes[i] = sum(income .* density_distr)
    Consumptions[i] = sum(consumption .* density_distr)#Incomes[i] - ss_star[1][7]
    Income_to_outputs[i] = Incomes[i]/Outputs[i]

    number_u_nodes = SSS[i][4][7]
    number_asset_grid = ss_star[1][2]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

    # [1] = income, earnings, wealth, consumption
    # [2] = All, W, SP, EMP
    # [3] = lambdas
    # = means, ginis, varlogs
    for s = 1:4 # [1] = income, earnings, wealth, consumption
        # if s == 1
        stat_distr = income
        if s == 2
            stat_distr = earnings
        elseif s == 3
            stat_distr = wealth
        elseif s == 4
            stat_distr = consumption
        end

        for h = 1:4 # [2] = All, W, SP, EMP
            choice = copy(occ_choice)
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                for a_i in 1:number_asset_grid
                    if h == 1 #All
                        choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
                    elseif h == 2 #W
                        if choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 1
                            choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
                        end
                    elseif h == 3 #SP
                        if choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 2
                            choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
                        else
                            choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
                        end
                    elseif h == 4 #EMP
                        if choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 3
                            choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
                        else
                            choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
                        end
                    end
                end
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
            varlogs[s,h,i] = sum(density_distr.*choice.*(log.(max.(1e-12,stat_distr)).-avglogs[s,h,i]).^2)/sum(density_distr.*choice)
            avglogs[s,h,i] /= sum(density_distr.*choice)
            if s==3
                avglogs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
                varlogs[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avglogs[s,h,i]).^2)/sum(density_distr.*choice)
                avglogs[s,h,i] /= sum(density_distr.*choice)
            end

            avgs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
            vars[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avgs[s,h,i]).^2)/sum(density_distr.*choice)
            avgs[s,h,i] /= sum(density_distr.*choice)
        end
    end

    w_choice   = copy(occ_choice)
    sp_choice  = copy(occ_choice)
    emp_choice = copy(occ_choice)
    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        for a_i in 1:number_asset_grid
            if w_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 1
                w_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
            end
            if sp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 2
                sp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
            else
                sp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
            end
            if emp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 3
                emp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
            else
                emp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
            end
        end
    end
    ent_choice = sp_choice.+emp_choice

    output_SPs[i]  = sum(ss_star[1][32].* density_distr.* sp_choice)
    output_EMPs[i] = sum(ss_star[1][32].* density_distr.*emp_choice)

    share_of_output_SPs[i]  = output_SPs[i]/Outputs[i]
    share_of_output_EMPs[i] = output_EMPs[i]/Outputs[i]

    ratio_of_capitals[i] = ss_star[1][6]

    asset_distr = sum(density_distr, dims=2:5)[:,1,1,1,1]
    asset_distr ./= sum(asset_distr)
    asset_distr_w = sum(density_distr.*w_choice, dims=2:5)[:,1,1,1,1]
    asset_distr_w ./= sum(asset_distr_w)
    asset_distr_se = sum(density_distr.*sp_choice, dims=2:5)[:,1,1,1,1]
    asset_distr_se ./= sum(asset_distr_se)
    asset_distr_emp = sum(density_distr.*emp_choice, dims=2:5)[:,1,1,1,1]
    asset_distr_emp ./= sum(asset_distr_emp)

    cdf = cumsum(asset_distr)
    cdf_w = cumsum(asset_distr_w)
    cdf_se = cumsum(asset_distr_se)
    cdf_emp = cumsum(asset_distr_emp)

    max_a_i = findlast(x->x<1.0-1e-5, cdf)

    cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Linear))
    cdf_w_1 = Schumaker(cdf_w[1:max_a_i]./cdf_w[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Linear))
    cdf_se_1 = Schumaker(cdf_se[1:max_a_i]./cdf_se[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Linear))
    cdf_emp_1 = Schumaker(cdf_emp[1:max_a_i]./cdf_emp[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Linear))

    #Q3-Q1
    interquartile_wealth[i] = evaluate(cdf_1,0.75) - evaluate(cdf_1,0.25)
    #Q8-Q2
    interquintile_wealth[i] = evaluate(cdf_1,0.8) - evaluate(cdf_1,0.2)
    #Q9-Q1
    interdecile_wealth[i] = evaluate(cdf_1,0.9) - evaluate(cdf_1,0.1)
    interquartile_wealth_w[i] = evaluate(cdf_w_1,0.75) - evaluate(cdf_w_1,0.25)
    interquintile_wealth_w[i] = evaluate(cdf_w_1,0.8) - evaluate(cdf_w_1,0.2)
    interdecile_wealth_w[i] = evaluate(cdf_w_1,0.9) - evaluate(cdf_w_1,0.1)
    interquartile_wealth_se[i] = evaluate(cdf_se_1,0.75) - evaluate(cdf_se_1,0.25)
    interquintile_wealth_se[i] = evaluate(cdf_se_1,0.8) - evaluate(cdf_se_1,0.2)
    interdecile_wealth_se[i] = evaluate(cdf_se_1,0.9) - evaluate(cdf_se_1,0.1)
    interquartile_wealth_emp[i] = evaluate(cdf_emp_1,0.75) - evaluate(cdf_emp_1,0.25)
    interquintile_wealth_emp[i] = evaluate(cdf_emp_1,0.8) - evaluate(cdf_emp_1,0.2)
    interdecile_wealth_emp[i] = evaluate(cdf_emp_1,0.9) - evaluate(cdf_emp_1,0.1)

    findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]

    #calc_share_of_wealth_between_shares_of_population
    function csowbsop(top_percent, bottom_percent, CDF_1, A_DISTR)
        a_top = evaluate(CDF_1, top_percent)
        a_bottom = evaluate(CDF_1, bottom_percent)
        #display([a_top, a_bottom])

        a_i_top = findnearest(asset_grid[1:max_a_i], a_top)
        a_i_bottom = findnearest(asset_grid[1:max_a_i], a_bottom)
        #display([a_i_top, a_i_bottom])

        #display([sum(asset_grid[a_i_bottom:a_i_top].*A_DISTR[a_i_bottom:a_i_top]),sum(asset_grid[1:max_a_i].*A_DISTR[1:max_a_i]), sum(asset_grid[a_i_bottom:a_i_top].*A_DISTR[a_i_bottom:a_i_top]) / sum(asset_grid[1:max_a_i].*A_DISTR[1:max_a_i])])
        return sum(asset_grid[a_i_bottom:a_i_top].*A_DISTR[a_i_bottom:a_i_top]) / sum(asset_grid[1:max_a_i].*A_DISTR[1:max_a_i])
    end

    palma_wealth[i] = csowbsop(1.0,0.9,cdf_1,asset_distr)/csowbsop(0.4,0.0,cdf_1,asset_distr)
    palma_wealth_w[i] = csowbsop(1.0,0.9,cdf_w_1,asset_distr_w)/csowbsop(0.4,0.0,cdf_w_1,asset_distr_w)
    palma_wealth_se[i] = csowbsop(1.0,0.9,cdf_se_1,asset_distr_se)/csowbsop(0.4,0.0,cdf_se_1,asset_distr_se)
    palma_wealth_emp[i] = csowbsop(1.0,0.9,cdf_emp_1,asset_distr_emp)/csowbsop(0.4,0.0,cdf_emp_1,asset_distr_emp)

    ratio2020[i] = csowbsop(1.0,0.8,cdf_1,asset_distr)/csowbsop(0.2,0.0,cdf_1,asset_distr)
    ratio2020_w[i] = csowbsop(1.0,0.8,cdf_w_1,asset_distr)/csowbsop(0.2,0.0,cdf_w_1,asset_distr)
    ratio2020_se[i] = csowbsop(1.0,0.8,cdf_se_1,asset_distr)/csowbsop(0.2,0.0,cdf_se_1,asset_distr)
    ratio2020_emp[i] = csowbsop(1.0,0.8,cdf_emp_1,asset_distr)/csowbsop(0.2,0.0,cdf_emp_1,asset_distr)

    #=
    gini_wealth[i] = sum([ cdf[y]*(1-cdf[y]) for y=1:max_a_i ]) / sum(asset_grid[1:max_a_i].*asset_distr[1:max_a_i])
    gini_wealth_w[i] = sum([ cdf_w[y]*(1-cdf_w[y]) for y=1:max_a_i ]) / sum(asset_grid[1:max_a_i].*asset_distr_w[1:max_a_i])
    gini_wealth_se[i] = sum([ cdf_se[y]*(1-cdf_se[y]) for y=1:max_a_i ]) / sum(asset_grid[1:max_a_i].*asset_distr_se[1:max_a_i])
    gini_wealth_emp[i] = sum([ cdf_emp[y]*(1-cdf_emp[y]) for y=1:max_a_i ]) / sum(asset_grid[1:max_a_i].*asset_distr_emp[1:max_a_i])
    =#
    gini_wealth[i] = sum([ asset_distr[y_i]*asset_distr[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr[1:max_a_i]))
    gini_wealth_w[i] = sum([ asset_distr_w[y_i]*asset_distr_w[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr_w[1:max_a_i]))
    gini_wealth_se[i] = sum([ asset_distr_se[y_i]*asset_distr_se[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr_se[1:max_a_i]))
    gini_wealth_emp[i] = sum([ asset_distr_emp[y_i]*asset_distr_emp[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr_emp[1:max_a_i]))

    unlimited_capital_choice = capital_d.<(lambdas[i].*(asset_grid.*ones(size(capital_d))))

    share_ent_unbound[i] = sum(ent_choice.*density_distr.*unlimited_capital_choice)/occENTs[i]
    share_se_unbound[i] = sum(sp_choice.*density_distr.*unlimited_capital_choice)/occSEs[i]
    share_emp_unbound[i] = sum(emp_choice.*density_distr.*unlimited_capital_choice)/occEMPs[i]

    for c = 1:4
        if c==1
            choice = ones(size(occ_choice))
        elseif c==2
            choice = w_choice
        elseif c==3
            choice = sp_choice
        elseif c==4
            choice = emp_choice
        else
            throw("c is out of bounds")
        end
        #[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
        #m = u_i,alpha_m_i,zeta_i
        m_skill_distr = permutedims(sum(density_distr.*choice, dims=[1,5])[1,:,:,:,1], [1,3,2])
        #w = u_i,alpha_m_i,alpha_w_i
        w_skill_distr = sum(density_distr.*choice, dims=[1,3])[1,:,1,:,:]

        avg_m_skill[c,i] = sum(m_skill_distr.*z_m_nodes)/sum(density_distr.*choice)
        avg_w_skill[c,i] = sum(w_skill_distr.*z_w_nodes)/sum(density_distr.*choice)

        var_m_skill[c,i] = sum(m_skill_distr.*(z_m_nodes .- avg_m_skill[c,i]*sum(density_distr.*choice)).^2)/sum(density_distr.*choice)
        var_w_skill[c,i] = sum(w_skill_distr.*(z_w_nodes .- avg_w_skill[c,i]*sum(density_distr.*choice)).^2)/sum(density_distr.*choice)

        #display([avg_w_skill[c,i],avg_m_skill[c,i], var_w_skill[c,i],var_m_skill[c,i]])
    end
    #throw(error)
end

asset_grid = SSS[1][1][3]
number_non_zero_asset_grid = 0
number_u_nodes = SSS[1][4][7]
number_zeta_nodes = global_approx_params[4]
number_alpha_m_nodes = global_approx_params[5]
number_alpha_w_nodes = global_approx_params[6]
for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
    global number_non_zero_asset_grid
    temp = findlast(round.(SSS[1][1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; digits=7) .!= 0.0)
    if isnothing(temp)
        temp = 0
    end
    if temp > number_non_zero_asset_grid
        number_non_zero_asset_grid = temp
    end
end
#policy = SSS[l][1][4]
num_lambdas = 20
num_lambdas_sp = findfirst(x -> x>0.75, C_Ys)-1

ratio_of_capitals = ratio_of_capitals./ratio_of_capitals[num_lambdas]
ratio_of_outputs = Outputs./Outputs[num_lambdas]

Capital = Capital./Capital[1]

#throw(error)

function create_plot(X::Vector{Float64},XLABEL::String,Y::Vector{Float64},YLABEL::String, IS_Y_PERCENTAGE::Bool=true, OCCUPATION::String="H")
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    end

    YLIMS1 = minimum(Y)
    YLIMS2 = maximum(Y[1:num_lambdas])
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
    plt = scatter(X[1:num_lambdas],Y[1:num_lambdas],
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
        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))
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
        plt = scatter(X[1:num_lambdas],Y[1:num_lambdas],
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

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)
#Interest rate and wage
plt = create_plot(C_Ys,"Credit/Output", Rs,"Interest Rate")
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_rate.png")
plt = create_plot(C_Ys,"Credit/Output", Ws,"Wage", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_wage.png")

#lambda to: Output, Consumption, Income/Output, Income
plt = create_plot(lambdas,"λ", Outputs,"Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Outputs.png")
plt = create_plot(lambdas,"λ", Consumptions,"Consumption", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Consumptions.png")
plt = create_plot(lambdas,"λ", Income_to_outputs,"Income/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Income_to_outputs.png")
plt = create_plot(lambdas,"λ", Incomes,"Income", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Incomes.png")
#Credit/Output to Income/Output
plt = create_plot(C_Ys,"Credit/Output", Income_to_outputs,"Income/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_Income_to_outputs.png")

#lambda to Credit/Output, Credit/Output to Capital/Output
plt = create_plot(lambdas,"λ", C_Ys,"Credit/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_credit_to_gdp.png")
plt = create_plot(C_Ys,"Credit/Output", K_Ys,"Capital/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_capital_to_gdp.png")

LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end
mkpath(LOCAL_DIR_OCCUPATION)
#Crdit/Output to share of W, SP, EMP
plts = create_plots(C_Ys,"Credit/Output", [occENTs, occSEs, occEMPs, occWs],["Share of Entrepreneurs","Share of Sole Proprietors","Share of Employers","Share of Workers"],true,["ENT","SP","EMP","W"])
display(plts[1])
savefig(plts[1],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_entrepreneurs.png")
display(plts[2])
savefig(plts[2],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_sole_proprietors.png")
display(plts[3])
savefig(plts[3],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_employers.png")
display(plts[4])
savefig(plts[4],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_workers.png")

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

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)
#TFP_ideal, TFP_data, Labour_productivity, Capital_productivity, Investment_productivity
plt = create_plot(C_Ys,"Credit/Output", TFPis,"TFP", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal.png")
plt = create_plot(C_Ys,"Credit/Output", TFPds,"TFP_data", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data.png")
plt = create_plot(C_Ys,"Credit/Output", Labour_productivity_s,"Output-to-labour", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_labour_productivity.png")
plt = create_plot(C_Ys,"Credit/Output", Capital_productivity_s,"Output-to-capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_capital_productivity.png")
plt = create_plot(C_Ys,"Credit/Output", Investment_productivity_s,"Expenses-to-revenue", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_investment_productivity.png")

# households' policy a_t+1
middle_point = Int64(round((1+num_lambdas)/2;digits=0))
middle_point1 = Int64(round((1+middle_point)/2;digits=0))
middle_point2 = Int64(round((middle_point+num_lambdas)/2;digits=0))
suitable_ls = [1,middle_point1,middle_point,middle_point2,num_lambdas]
#,u$(u),z$(zeta),aw$(alpha_w)
labels = permutedims(["λ=$(round(lambdas[l];digits=2)),am$(alpha_m)" for l in suitable_ls for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [#=Int64(round((1+number_alpha_m_nodes)/2;digits=0))=#1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]])
p2 = plot(asset_grid[1:number_non_zero_asset_grid],[SSS[l][1][4][1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for l in suitable_ls for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [#=Int64(round((1+number_alpha_m_nodes)/2;digits=0))=#1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]], label=labels, legend=:outertopleft)
plt = plot(p2, title="Policy functions (high vs low fixed effects)")
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_policy_hvsl_fixed_lambda.png")


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

    CHOICE_NAMES = ["Households","Workers","Sole Proprietors","Employers"]

    # calculate mean
    #means[s,h,i]
    LABELS = ["Mean of $cn' $stat_name" for cn in CHOICE_NAMES]
    ms = [means[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [means[s,1,:], means[s,2,:], ms, means[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [means[s,1,:], means[s,2,:], means[s,3,:], means[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name)_full.png")
    end


    # calculate gini coefficent
    #ginis[s,h,i]
    LABELS = ["Gini of $cn' $stat_name" for cn in CHOICE_NAMES]
    gs = [ginis[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], gs, ginis[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], ginis[s,3,:], ginis[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name)_full.png")
    end


    # calculate variance of log-s
    #avglogs[s,h,i]
    LABELS = ["Average of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    als = [avglogs[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], als, avglogs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], avglogs[s,3,:], avglogs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name)_full.png")
    end


    #varlogs[s,h,i]
    if s!=3
        LABELS = ["Variance of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    else
        LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    end
    vls = [varlogs[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], vls, varlogs[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], varlogs[s,3,:], varlogs[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name)_full.png")
    end
end

# output by SP
plt = create_plot(C_Ys,"Credit/Output", output_SPs,"Aggregate output by Sole Proprietors", false, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_output_sp.png")
# output by EMP
plt = create_plot(C_Ys,"Credit/Output", output_EMPs,"Aggregate output by Employers", false, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_output_emp.png")

# share of output by SP
plt = create_plot(C_Ys,"Credit/Output", share_of_output_SPs,"Share of output by Sole Proprietors", true, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_share_of_output_sp.png")
# share of output by EMP
plt = create_plot(C_Ys,"Credit/Output", share_of_output_EMPs,"Share of output by Employers", true, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_share_of_output_emp.png")

# ratio of capital restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_capitals,"Ratio of Restricted Capital to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_ratio_of_capital.png")
# ratio of output restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_outputs,"Ratio of Restricted Output to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_ratio_of_output.png")

# aggregate output ~~~ DOESN'T WORK
#=
using StatsPlots
xt = Plots.xticks(plt[1])
plt2 = groupedbar(C_Ys[1:num_lambdas],[output_SPs[1:num_lambdas] output_EMPs[1:num_lambdas]],
        bar_position = :stack,
        #bar_width=[abs(C_Ys[2]-C_Ys[1]); [abs(C_Ys[i+1]-C_Ys[i-1]) for i=2:num_lambdas-1]; abs(C_Ys[num_lambdas]-C_Ys[num_lambdas-1])]./2.0,
        bar_width=median([abs(C_Ys[i+1]-C_Ys[i-1]) for i=2:num_lambdas])./2.0,
        xticks=xt,
        label=["SP" "EMP"],
        color=["red" "green"])
scatter!(C_Ys[1:num_lambdas],[Outputs[1:num_lambdas] output_EMPs[1:num_lambdas]],#=regression=true,=#color=["blue" "green"],xlabel="Credit/Output",ylabel="Output",label=["Total" "EMP"],legend=:bottomright)
vline!([C_Ys[calibrated_lambda]], color="grey")
annotate!([C_Ys[calibrated_lambda]], Outputs[num_lambdas], text(" Calibrated economy ($(round(C_Ys[calibrated_lambda]; digits=2)),$(round(Outputs[num_lambdas]; digits=2))) ", :grey, :left, 7))
display(plt2)
savefig(plt2,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_Outputs_modified.png")
=#

# average and variance of working and managerial skills across occupations
CHOICE_NAMES = ["Households","Workers","Sole Proprietors","Employers"]
LABELS = ["Average Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[2,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,false,["H","W","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_avg_m_skill_$(occup).png")
end

LABELS = ["Average Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[2,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,false,["H","W","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_avg_w_skill_$(occup).png")
end

LABELS = ["Variance of Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[2,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_var_m_skill_$(occup).png")
end

LABELS = ["Variance of Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[2,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_var_w_skill_$(occup).png")
end

# cumulative growth & unconstrained
plt = create_plot(C_Ys,"Credit/Output", Capital.-1.0,"Cumulative growth of aggregate assets", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_cum_growth_agg_assets.png")

plt = create_plot(C_Ys,"Credit/Output", share_ent_unbound,"Share of unconstrained entrepreneurs", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_ent_unbound.png")

plt = create_plot(C_Ys,"Credit/Output", share_se_unbound,"Share of unconstrained Sole Proprietors", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_se_unbound.png")

plt = create_plot(C_Ys,"Credit/Output", share_emp_unbound,"Share of unconstrained employers", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_emp_unbound.png")

plts = create_plots(C_Ys,"Credit/Output", [share_ent_unbound, share_se_unbound, share_emp_unbound],["Share of Unconstrained Entrepreneurs","Share of Unconstrained Sole Proprietors","Share of Unconstrained Employers"],true,["ENT","SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_ent_unbound.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_se_unbound.png")
savefig(plts[3],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_emp_unbound.png")

# Print statistics for calibrated economy

# income, earnings, wealth, consumption
# All, W, SP, EMP
# lambdas
#means = zeros(4,4,length(lambdas))
#ginis = zeros(4,4,length(lambdas))
#vars = zeros(4,4,length(lambdas))
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

    println("Mean of $(stat_name) for H,W,SP,EMP: $(round.(means[i,:,calibrated_lambda];digits=2))")
    println("Variance of $(stat_name) for H,W,SP,EMP: $(round.(vars[i,:,calibrated_lambda];digits=2))")
    println("Varlog of $(stat_name) for H,W,SP,EMP: $(round.(varlogs[i,:,calibrated_lambda];digits=2))")
    println("Gini of $(stat_name) for H,W,SP,EMP: $(round.(ginis[i,:,calibrated_lambda];digits=2))")
end