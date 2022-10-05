##############
# Script to generate the probabilities of occupational mobility
##############

include("Functions/print_sameline.jl")
include("Functions/profit.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots

#using SchumakerSpline
#using Statistics
using JLD2
using ProgressMeter

CODENAME = "SS_2092_69"
country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_big_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_big_grid\\$(country)_$(CODENAME)\\"
end

# analytical part
function f()


    TIME_PERIODS = 50#100
    @load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2

    global_params = ss_1[1][46]
    global_approx_params = ss_1[1][47]
    model_params = ss_1[1][48]
    approx_object = ss_1[1][49]

    delta = model_params[3]
    gamma = model_params[4]
    eta = model_params[5]
    theta = model_params[6]
    c_e = model_params[7]

    number_u_nodes = approx_object[7]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]
    number_asset_grid = trans_res[9]
    asset_grid = trans_res[10]
    T = TIME_PERIODS

    z_m_nodes = approx_object[1]
    z_w_nodes = approx_object[2]
    number_u_nodes = approx_object[7]

    P_u = approx_object[4]
    p_alpha = model_params[13]
    P_zeta = approx_object[3]
    stat_P_u = approx_object[5]
    P_alpha = approx_object[6]

    lambda_s = trans_res[2]
    r_s = trans_res[3]
    w_s = trans_res[4]
    capital_s_distr_s = trans_res[5][7]
    policy_s = trans_res[5][8]

    function calculate_from_policies(policies)
        a1_indices = Array{Int64}(undef,T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_1 = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_2 = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

        p = Progress((number_asset_grid*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes), dt=0.5,
                     barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                     barlen=25)
        Threads.@threads for (t,(a_i,(u_i,(zeta_i,(alpha_m_i,alpha_w_i))))) in collect(Iterators.product(1:T,Iterators.product(1:number_asset_grid,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))))
            try
                a1 = policies[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                j_1 = sum(a1 .>= asset_grid)
                j = j_1 + 1

                a1_indices[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = j_1
                if j <= number_asset_grid && j_1 >= 1
                    lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
                    lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
                elseif j_1 == number_asset_grid
                    lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
                    lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
                end
            catch e
                println_sameline("Trouble is in loop for increasing grid for policies")
                throw(e)
            end
            next!(p)
        end

        return a1_indices, lottery_prob_1, lottery_prob_2
    end
    a1_indices, lottery_prob_1, lottery_prob_2 = calculate_from_policies(policy_s)
    println_sameline("calculate_from_policies - Done")

    init_density_distr = capital_s_distr_s[1,:,:,:,:,:]
    init_occ_choice,  = compute_income_profile(asset_grid,number_asset_grid,r_s[1], w_s[1], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[1], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
    init_w_choice = Float64.(init_occ_choice.==1.0)
    init_sp_choice = Float64.(init_occ_choice.==2.0)
    init_emp_choice = Float64.(init_occ_choice.==3.0)

    init_density_distr_w = init_density_distr.*init_w_choice
    init_density_distr_sp = init_density_distr.*init_sp_choice
    init_density_distr_emp = init_density_distr.*init_emp_choice

    density_distr_w = copy(init_density_distr_w)
    density_distr_sp = copy(init_density_distr_sp)
    density_distr_emp = copy(init_density_distr_emp)

    occ_trans_matrix = zeros(T,3,3)
    occ_trans_matrix[1,1,1] = 1.0
    occ_trans_matrix[1,2,2] = 1.0
    occ_trans_matrix[1,3,3] = 1.0

    absolute_occ_matrix = zeros(T,3,3)
    absolute_occ_matrix[1,1,1] = sum(density_distr_w)
    absolute_occ_matrix[1,2,2] = sum(density_distr_sp)
    absolute_occ_matrix[1,3,3] = sum(density_distr_emp)

    opposite_occ_trans_matrix = zeros(T,3,3)
    opposite_occ_trans_matrix[1,1,1] = 1.0
    opposite_occ_trans_matrix[1,2,2] = 1.0
    opposite_occ_trans_matrix[1,3,3] = 1.0

    mean_a_occ_trans_matrix = zeros(T,3,3)
    mean_z_m_occ_trans_matrix = zeros(T,3,3)
    mean_z_w_occ_trans_matrix = zeros(T,3,3)
    for i = 1:3
        dd_o = density_distr_w
        if i==2
            dd_o = density_distr_sp
        elseif i==3
            dd_o = density_distr_emp
        end
        mean_a_occ_trans_matrix[1,i,i] = sum(asset_grid .* dd_o)/sum(dd_o)
        m_skill_distr = permutedims(sum(dd_o, dims=[1,5])[1,:,:,:,1], [1,3,2])
        mean_z_m_occ_trans_matrix[1,i,i] = sum(m_skill_distr.*z_m_nodes)/sum(dd_o)
        w_skill_distr = sum(dd_o, dims=[1,3])[1,:,1,:,:]
        mean_z_w_occ_trans_matrix[1,i,i] = sum(w_skill_distr.*z_w_nodes)/sum(dd_o)
    end

    p = Progress((T-1)*((number_asset_grid*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes)+(number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes)+(number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes)), dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=100)
    for t=1:T-1

        density_distr_w_pr = copy(density_distr_w).*0.0
        density_distr_sp_pr = copy(density_distr_sp).*0.0
        density_distr_emp_pr = copy(density_distr_emp).*0.0

        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            for a_i in 1:number_asset_grid
                j_1 = a1_indices[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                density_distr_sp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                if j_1 != number_asset_grid
                    density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_sp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                end
                next!(p)
            end
        end

        density_distr_marginal_assets_w = sum(sum(sum(sum(density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        density_distr_marginal_assets_sp = sum(sum(sum(sum(density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        density_distr_marginal_assets_emp = sum(sum(sum(sum(density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
            temp_w = sum(density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_w = temp_w*P_u
            temp_sp = sum(density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_sp = temp_sp*P_u
            temp_emp = sum(density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_emp = temp_emp*P_u
            for zeta_prime_i in 1:number_zeta_nodes
                density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_w
                density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_sp
                density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_emp

                next!(p)
            end
        end

        Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_w
            density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_sp
            density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_emp

            next!(p)
        end

        occ_choice,  = compute_income_profile(asset_grid,number_asset_grid,r_s[t+1], w_s[t+1], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t+1], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
        w_choice = Float64.(occ_choice.==1.0)
        sp_choice = Float64.(occ_choice.==2.0)
        emp_choice = Float64.(occ_choice.==3.0)

        # conditional probability of households - P( ot = OCCj | o1 = OCCi )
        # W->W
        occ_trans_matrix[t+1, 1,1] = sum(w_choice.*density_distr_w_pr)/sum(init_density_distr_w)
        # W->SP
        occ_trans_matrix[t+1, 1,2] = sum(sp_choice.*density_distr_w_pr)/sum(init_density_distr_w)
        # W->EMP
        occ_trans_matrix[t+1, 1,3] = sum(emp_choice.*density_distr_w_pr)/sum(init_density_distr_w)

        # SP->W
        occ_trans_matrix[t+1, 2,1] = sum(w_choice.*density_distr_sp_pr)/sum(init_density_distr_sp)
        # W->SP
        occ_trans_matrix[t+1, 2,2] = sum(sp_choice.*density_distr_sp_pr)/sum(init_density_distr_sp)
        # W->EMP
        occ_trans_matrix[t+1, 2,3] = sum(emp_choice.*density_distr_sp_pr)/sum(init_density_distr_sp)

        # EMP->W
        occ_trans_matrix[t+1, 3,1] = sum(w_choice.*density_distr_emp_pr)/sum(init_density_distr_emp)
        # EMP->SP
        occ_trans_matrix[t+1, 3,2] = sum(sp_choice.*density_distr_emp_pr)/sum(init_density_distr_emp)
        # EMP->EMP
        occ_trans_matrix[t+1, 3,3] = sum(emp_choice.*density_distr_emp_pr)/sum(init_density_distr_emp)

        # conditional probability of households - P( o1 = OCCi & ot = OCCj )
        absolute_occ_matrix[t+1, 1,1] = sum(w_choice.*density_distr_w_pr)
        absolute_occ_matrix[t+1, 1,2] = sum(sp_choice.*density_distr_w_pr)
        absolute_occ_matrix[t+1, 1,3] = sum(emp_choice.*density_distr_w_pr)

        absolute_occ_matrix[t+1, 2,1] = sum(w_choice.*density_distr_sp_pr)
        absolute_occ_matrix[t+1, 2,2] = sum(sp_choice.*density_distr_sp_pr)
        absolute_occ_matrix[t+1, 2,3] = sum(emp_choice.*density_distr_sp_pr)

        absolute_occ_matrix[t+1, 3,1] = sum(w_choice.*density_distr_emp_pr)
        absolute_occ_matrix[t+1, 3,2] = sum(sp_choice.*density_distr_emp_pr)
        absolute_occ_matrix[t+1, 3,3] = sum(emp_choice.*density_distr_emp_pr)

        # conditional probability of households - P( o1 = OCCi | ot = OCCj )
        opposite_occ_trans_matrix[t+1, 1,1] = absolute_occ_matrix[t+1, 1,1]/sum(absolute_occ_matrix[t+1, :,1])
        opposite_occ_trans_matrix[t+1, 1,2] = absolute_occ_matrix[t+1, 1,2]/sum(absolute_occ_matrix[t+1, :,2])
        opposite_occ_trans_matrix[t+1, 1,3] = absolute_occ_matrix[t+1, 1,3]/sum(absolute_occ_matrix[t+1, :,3])

        opposite_occ_trans_matrix[t+1, 2,1] = absolute_occ_matrix[t+1, 2,1]/sum(absolute_occ_matrix[t+1, :,1])
        opposite_occ_trans_matrix[t+1, 2,2] = absolute_occ_matrix[t+1, 2,2]/sum(absolute_occ_matrix[t+1, :,2])
        opposite_occ_trans_matrix[t+1, 2,3] = absolute_occ_matrix[t+1, 2,3]/sum(absolute_occ_matrix[t+1, :,3])

        opposite_occ_trans_matrix[t+1, 3,1] = absolute_occ_matrix[t+1, 3,1]/sum(absolute_occ_matrix[t+1, :,1])
        opposite_occ_trans_matrix[t+1, 3,2] = absolute_occ_matrix[t+1, 3,2]/sum(absolute_occ_matrix[t+1, :,2])
        opposite_occ_trans_matrix[t+1, 3,3] = absolute_occ_matrix[t+1, 3,3]/sum(absolute_occ_matrix[t+1, :,3])

        # mean assets, managerial and working skills for ( o1 = OCCi & ot = OCCj )
        for (i,j) = collect(Iterators.product(1:3,1:3))
            dd_o = density_distr_w_pr
            if i==2
                dd_o = density_distr_sp_pr
            elseif i==3
                dd_o = density_distr_emp_pr
            end
            if j==1
                dd_o = w_choice.*dd_o
            elseif j==2
                dd_o = sp_choice.*dd_o
            elseif j==3
                dd_o = emp_choice.*dd_o
            end
            mean_a_occ_trans_matrix[t+1,i,j] = sum(asset_grid .* dd_o)/sum(dd_o)
            m_skill_distr = permutedims(sum(dd_o, dims=[1,5])[1,:,:,:,1], [1,3,2])
            mean_z_m_occ_trans_matrix[t+1,i,j] = sum(m_skill_distr.*z_m_nodes)/sum(dd_o)
            w_skill_distr = sum(dd_o, dims=[1,3])[1,:,1,:,:]
            mean_z_w_occ_trans_matrix[t+1,i,j] = sum(w_skill_distr.*z_w_nodes)/sum(dd_o)
        end

        density_distr_w = copy(density_distr_w_pr)
        density_distr_sp = copy(density_distr_sp_pr)
        density_distr_emp = copy(density_distr_emp_pr)

    end

    println_sameline("Final cycle - Done")

    return occ_trans_matrix, absolute_occ_matrix, opposite_occ_trans_matrix, mean_a_occ_trans_matrix,mean_z_m_occ_trans_matrix,mean_z_w_occ_trans_matrix
end
occ_trans_matrix, absolute_occ_matrix, opposite_occ_trans_matrix, mean_a_occ_trans_matrix,mean_z_m_occ_trans_matrix,mean_z_w_occ_trans_matrix = f()

# graphing part
LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end

occ_text = ["W","SP","EMP"]

plts_otm = Array{Any}(undef,3,3)
plts_aom = Array{Any}(undef,3,3)
plts_ootm = Array{Any}(undef,3,3)
plts_maotm = Array{Any}(undef,3,3)
plts_mzmotm = Array{Any}(undef,3,3)
plts_mzwotm = Array{Any}(undef,3,3)

for (i,j) = collect(Iterators.product(1:3,1:3))
    T = length(occ_trans_matrix[:,i,j])
    plts_otm[i,j] = plot(2:T, occ_trans_matrix[2:T,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_aom[i,j] = plot(1:T, absolute_occ_matrix[:,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_ootm[i,j] = plot(2:T, opposite_occ_trans_matrix[2:T,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

    plts_maotm[i,j] = plot(2:T, mean_a_occ_trans_matrix[2:T,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_mzmotm[i,j] = plot(2:T, mean_z_m_occ_trans_matrix[2:T,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_mzwotm[i,j] = plot(2:T, mean_z_w_occ_trans_matrix[2:T,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

end
# it shows the share of occupation1 at time 1 who become occupation2 at time t
plt_otm = plot(plts_otm[1,1],plts_otm[1,2],plts_otm[1,3], plts_otm[2,1],plts_otm[2,2],plts_otm[2,3], plts_otm[3,1],plts_otm[3,2],plts_otm[3,3], layout=(3,3))
display(plt_otm)
savefig(plt_otm,"$(LOCAL_DIR_OCCUPATION)time_occupational_mobility.png")

# it shows the share of households who were occupation1 at time 1 and become occupation2 at time t
plt_aom = plot(plts_aom[1,1],plts_aom[1,2],plts_aom[1,3], plts_aom[2,1],plts_aom[2,2],plts_aom[2,3], plts_aom[3,1],plts_aom[3,2],plts_aom[3,3], layout=(3,3))
display(plt_aom)
savefig(plt_aom,"$(LOCAL_DIR_OCCUPATION)time_occupational_conditional_distribution.png")

# it shows the share of occupation2 at time t who were occupation1 at time 1
plt_ootm = plot(plts_ootm[1,1],plts_ootm[1,2],plts_ootm[1,3], plts_ootm[2,1],plts_ootm[2,2],plts_ootm[2,3], plts_ootm[3,1],plts_ootm[3,2],plts_ootm[3,3], layout=(3,3))
display(plt_ootm)
savefig(plt_ootm,"$(LOCAL_DIR_OCCUPATION)time_opposite_occupational_mobility.png")

# it shows the mean assets of households who were occupation1 at time 1 and become occupation2 at time t
plt_maotm = plot(plts_maotm[1,1],plts_maotm[1,2],plts_maotm[1,3], plts_maotm[2,1],plts_maotm[2,2],plts_maotm[2,3], plts_maotm[3,1],plts_maotm[3,2],plts_maotm[3,3], layout=(3,3))
display(plt_maotm)
savefig(plt_maotm,"$(LOCAL_DIR_OCCUPATION)time_occupational_conditional_distribution_mean_a.png")

# it shows the mean z_m of households who were occupation1 at time 1 and become occupation2 at time t
plt_mzmotm = plot(plts_mzmotm[1,1],plts_mzmotm[1,2],plts_mzmotm[1,3], plts_mzmotm[2,1],plts_mzmotm[2,2],plts_mzmotm[2,3], plts_mzmotm[3,1],plts_mzmotm[3,2],plts_mzmotm[3,3], layout=(3,3))
display(plt_mzmotm)
savefig(plt_mzmotm,"$(LOCAL_DIR_OCCUPATION)time_occupational_conditional_distribution_mean_z_m.png")

# it shows the share of households who were occupation1 at time 1 and become occupation2 at time t
plt_mzwotm = plot(plts_mzwotm[1,1],plts_mzwotm[1,2],plts_mzwotm[1,3], plts_mzwotm[2,1],plts_mzwotm[2,2],plts_mzwotm[2,3], plts_mzwotm[3,1],plts_mzwotm[3,2],plts_mzwotm[3,3], layout=(3,3))
display(plt_mzwotm)
savefig(plt_mzwotm,"$(LOCAL_DIR_OCCUPATION)time_occupational_conditional_distribution_mean_z_w.png")
