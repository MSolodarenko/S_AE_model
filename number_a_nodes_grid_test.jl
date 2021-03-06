include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2

print_sameline("Loading functions for steady_state procedure")
include("Functions/steady_state.jl")

# global parameters of the model's code
#                   1           2           3       4
#                gen_tol_x, gen_tol_f, distr_tol, val_tol
GLOBAL_PARAMS = [1e-6, 1e-4, 1e-9, 1e-10]#[1e-8, 1e-4, 1e-9, 1e-7]#[1e-8, 1e-4, 1e-12, 1e-9]#[1e-8, 1e-4, 1e-7, 1e-5]#

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [15,3,3,3,6,3]#[49,3,3,3,6,3]#[25,5,5,3,6,3]#

# parameters of the model's economy (Italy)
LAMBDA = 1.513028#1.548387
BETA = 0.917506#0.859906
DELTA = 0.1
GAMMA = 0.16
ETA = 0.3
THETA = 0.54

C_E = 0.007001#0.036321
RHO_M = 0.710266#0.849951
RHO_W = 0.96
SIGMA_EPS_M = 1.022857#0.942809
RHO_EPS_M_W = 0.361605#0.384691
SIGMA_ZETA = 0.091089#0.151709
P_ALPHA = 0.024445#0.043692
ETA_ALPHA = 5.896035#5.202754
PROB_NODE1_ALPHA = 0.39
MU_M_ALPHA = -4.817675#-3.909137
RHO_ALPHA_M_W = 0.181454#0.095903
SIGMA_ALPHA_W = 0.12
SIGMA_EPS_W = 0.021947#0.285519

CRRA = 1.0 # >0.0

R_ = -0.005338#-0.069128
W_ = 0.230294#0.393906

gen_tol_x = GLOBAL_PARAMS[1]
gen_tol_f = GLOBAL_PARAMS[2]

# interest rate bounds
r_min = -DELTA#-delta
r_max = 1/BETA-1#1/beta-1

# wage bounds #### NOT DONE!!! ########
w_min = 0.01#0.25
w_max = 0.99#0.5

optimal_r = R_#(r_min+r_max)/2#r#
optimal_w = W_#(w_min+w_max)/2#w#

text_output = false
fig_output = false
calc_add_results = false

country = "Italy"
#NUMBER_A_NODES = 49#15
guess_R = -0.008174259548335908#(r_min+r_max)/2#
guess_W = 0.23529682939411575#(w_min+w_max)/2#

MODEL_PARAMS = [LAMBDA, BETA, DELTA, GAMMA, ETA, THETA, C_E, RHO_M, RHO_W, SIGMA_EPS_M, RHO_EPS_M_W, SIGMA_ZETA, P_ALPHA, ETA_ALPHA, PROB_NODE1_ALPHA, MU_M_ALPHA, RHO_ALPHA_M_W, SIGMA_ALPHA_W, SIGMA_EPS_W, CRRA]
#=
open("Results/number_a_nodes_grid_results.txt","a") do io
    guess_R = -0.008174259548335908#(r_min+r_max)/2#
    guess_W = 0.23529682939411575#(w_min+w_max)/2#
    println_sameline(("NUMBER_A_NODES","time_length","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    println(io,("NUMBER_A_NODES","time_length","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    for NUMBER_A_NODES = 15:2:53
        try
            global_approx_params = [NUMBER_A_NODES,3,3,3,6,3]

            start_time = time()
            approx_object = build_skill_nodes(global_approx_params, MODEL_PARAMS)
            res = AllubErosa(guess_R,guess_W, GLOBAL_PARAMS, global_approx_params, MODEL_PARAMS, approx_object)

            time_length = round(time() - start_time; digits=0)
            len_r = round(res[10]; digits=7)
            len_w = round(res[11]; digits=7)
            cap_out=round(res[12]; digits=7)
            cre_out=round(res[13]; digits=7)
            occ_w = round(res[14]; digits=7)
            occ_se= round(res[15]; digits=7)
            occ_emp=round(res[16]; digits=7)

            println_sameline((NUMBER_A_NODES,time_length,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))
            println(io,(NUMBER_A_NODES,time_length,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))
        catch e
            println_sameline((NUMBER_A_NODES,e))
            println(io,(NUMBER_A_NODES,"failed to converge in time"))
        end
    end
end
=#
open("Results/number_a_nodes_ss_grid_results.txt","a") do io
    guess_R = -0.008174259548335908#(r_min+r_max)/2#
    guess_W = 0.23529682939411575#(w_min+w_max)/2#
    println_sameline(("NUMBER_A_NODES","time_length","r","w","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    println(io,("NUMBER_A_NODES","time_length","r","w","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    for NUMBER_A_NODES = 15:2:53
        try
            global_approx_params = [NUMBER_A_NODES,3,3,3,6,3]

            start_time = time()
            ss_star = steady_state(guess_R, guess_W, GLOBAL_PARAMS,global_approx_params, MODEL_PARAMS)

            time_length = round(time() - start_time; digits=0)
            r     = round(ss_star[2];     digits=7)
            w     = round(ss_star[3];     digits=7)
            len_r = round(ss_star[1][10]; digits=7)
            len_w = round(ss_star[1][11]; digits=7)
            cap_out=round(ss_star[1][12]; digits=7)
            cre_out=round(ss_star[1][13]; digits=7)
            occ_w = round(ss_star[1][14]; digits=7)
            occ_se= round(ss_star[1][15]; digits=7)
            occ_emp=round(ss_star[1][16]; digits=7)

            println_sameline((NUMBER_A_NODES,time_length,r,w,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))
            println(io,      (NUMBER_A_NODES,time_length,r,w,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))

            guess_R = r
            guess_W = w
        catch e
            println_sameline((NUMBER_A_NODES,e))
            println(io,(NUMBER_A_NODES,"failed to converge in time"))
        end
    end
end

open("Results/number_u_nodes_grid_results.txt","a") do io
    guess_R = -0.008174259548335908#(r_min+r_max)/2#
    guess_W = 0.23529682939411575#(w_min+w_max)/2#
    NUMBER_A_NODES = 49#15
    println_sameline(("NUMBER_U_NODES","time_length","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    println(io,("NUMBER_U_NODES","time_length","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    for NUMBER_U_NODES = 3:2:9
        try
            global_approx_params = [NUMBER_A_NODES,NUMBER_U_NODES,NUMBER_U_NODES,3,6,3]

            start_time = time()
            approx_object = build_skill_nodes(global_approx_params, MODEL_PARAMS)
            res = AllubErosa(guess_R,guess_W, GLOBAL_PARAMS, global_approx_params, MODEL_PARAMS, approx_object)

            time_length = round(time() - start_time; digits=0)
            len_r = round(res[10]; digits=7)
            len_w = round(res[11]; digits=7)
            cap_out=round(res[12]; digits=7)
            cre_out=round(res[13]; digits=7)
            occ_w = round(res[14]; digits=7)
            occ_se= round(res[15]; digits=7)
            occ_emp=round(res[16]; digits=7)

            println_sameline((NUMBER_A_NODES,time_length,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))
        catch e
            println_sameline((NUMBER_U_NODES,e))
            println(io,(NUMBER_U_NODES,"failed to converge in time"))
        end
    end
end

open("Results/number_u_nodes_ss_grid_results.txt","a") do io
    guess_R = -0.008174259548335908#(r_min+r_max)/2#
    guess_W = 0.23529682939411575#(w_min+w_max)/2#
    NUMBER_A_NODES = 49#15
    println_sameline(("NUMBER_U_NODES","time_length","r","w","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    println(io,("NUMBER_U_NODES","time_length","r","w","len_r","len_w","cap_out","cre_out","occ_w","occ_se","occ_emp"))
    for NUMBER_U_NODES = 3:2:9
        try
            global_approx_params = [NUMBER_A_NODES,NUMBER_U_NODES,NUMBER_U_NODES,3,6,3]

            start_time = time()
            ss_star = steady_state(guess_R, guess_W, GLOBAL_PARAMS,global_approx_params, MODEL_PARAMS)

            time_length = round(time() - start_time; digits=0)
            r     = round(ss_star[2];     digits=7)
            w     = round(ss_star[3];     digits=7)
            len_r = round(ss_star[1][10]; digits=7)
            len_w = round(ss_star[1][11]; digits=7)
            cap_out=round(ss_star[1][12]; digits=7)
            cre_out=round(ss_star[1][13]; digits=7)
            occ_w = round(ss_star[1][14]; digits=7)
            occ_se= round(ss_star[1][15]; digits=7)
            occ_emp=round(ss_star[1][16]; digits=7)

            println_sameline((NUMBER_U_NODES,time_length,r,w,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))
            println(io,      (NUMBER_U_NODES,time_length,r,w,len_r,len_w,cap_out,cre_out,occ_w,occ_se,occ_emp))

            guess_R = r
            guess_W = w
        catch e
            println_sameline((NUMBER_U_NODES,e))
            println(io,(NUMBER_U_NODES,"failed to converge in time"))
        end
    end
end
