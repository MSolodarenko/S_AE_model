##############
# Script to generate the probabilities of occupational mobility
##############

include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
#using Plots
#using PyPlot

#using SchumakerSpline
#using Statistics
using JLD2
using ProgressMeter

LOCAL_DIR = "$(@__DIR__)/Results/Stationary/SS/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Stationary\\SS\\"
end

LAMBDA =            1.665907#10.0#

@load "$(LOCAL_DIR)SS_lambda_$(round(LAMBDA;digits=2)).jld2" SS

global_params = SS[1][46]
global_approx_params = SS[1][47]
model_params = SS[1][48]
approx_object = SS[1][49]

number_u_nodes = approx_object[7]
number_zeta_nodes = global_approx_params[4]
number_alpha_m_nodes = global_approx_params[5]
number_alpha_w_nodes = global_approx_params[6]
number_asset_grid = SS[1][2]

P_u = approx_object[4]
p_alpha = model_params[13]
P_zeta = approx_object[3]
stat_P_u = approx_object[5]
P_alpha = approx_object[6]

density_distr = SS[1][5]
occ_choice = SS[1][22]

a1_indices = SS[1][39]
lottery_prob_1 = SS[1][40]
lottery_prob_2 = SS[1][41]

w_choice = Float64.(occ_choice.==1.0)
sp_choice = Float64.(occ_choice.==2.0)
emp_choice = Float64.(occ_choice.==3.0)

density_distr_w = density_distr.*w_choice
density_distr_sp = density_distr.*sp_choice
density_distr_emp = density_distr.*emp_choice

density_distr_w_pr = copy(density_distr_w).*0.0
density_distr_sp_pr = copy(density_distr_sp).*0.0
density_distr_emp_pr = copy(density_distr_emp).*0.0

p = Progress((number_asset_grid*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
    for a_i in 1:number_asset_grid
        j_1 = a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

        density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
        density_distr_sp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
        density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

        if j_1 != number_asset_grid
            density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            density_distr_sp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
        end
        next!(p)
    end
end

density_distr_marginal_assets_w = sum(sum(sum(sum(density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
density_distr_marginal_assets_sp = sum(sum(sum(sum(density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
density_distr_marginal_assets_emp = sum(sum(sum(sum(density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

p = Progress((number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
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

p = Progress((number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes), dt=0.5,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=25)
Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
    density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_w
    density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_sp
    density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_emp

    next!(p)
end

occ_trans_matrix = zeros(3,3)

# W->W
occ_trans_matrix[1,1] = sum(w_choice.*density_distr_w_pr)/sum(density_distr_w)
# W->SP
occ_trans_matrix[1,2] = sum(sp_choice.*density_distr_w_pr)/sum(density_distr_w)
# W->EMP
occ_trans_matrix[1,3] = sum(emp_choice.*density_distr_w_pr)/sum(density_distr_w)

# SP->W
occ_trans_matrix[2,1] = sum(w_choice.*density_distr_sp_pr)/sum(density_distr_sp)
# W->SP
occ_trans_matrix[2,2] = sum(sp_choice.*density_distr_sp_pr)/sum(density_distr_sp)
# W->EMP
occ_trans_matrix[2,3] = sum(emp_choice.*density_distr_sp_pr)/sum(density_distr_sp)

# EMP->W
occ_trans_matrix[3,1] = sum(w_choice.*density_distr_emp_pr)/sum(density_distr_emp)
# EMP->SP
occ_trans_matrix[3,2] = sum(sp_choice.*density_distr_emp_pr)/sum(density_distr_emp)
# EMP->EMP
occ_trans_matrix[3,3] = sum(emp_choice.*density_distr_emp_pr)/sum(density_distr_emp)

open("$(LOCAL_DIR)res_$(round(LAMBDA;digits=3)).txt", "a") do f
    write(f, "Occupational mobility probabilities:\n")
    write(f, "$(round.(occ_trans_matrix[1,:].*100;digits=2))\n")
    write(f, "$(round.(occ_trans_matrix[2,:].*100;digits=2))\n")
    write(f, "$(round.(occ_trans_matrix[3,:].*100;digits=2))\n")
end
println_sameline("Done - $(LAMBDA)")
