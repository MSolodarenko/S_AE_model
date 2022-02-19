include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2
using XLSX

print_sameline("Loading functions for steady_state procedure")
include("Functions/steady_state.jl")

LOCAL_DIR = "$(@__DIR__)/Results/Calibration/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Calibration\\"
end
mkpath(LOCAL_DIR)


#                 1      2         3      4       5          6       7       8       9       10      11      12         13        14
#                 K/Y,   Credit/Y, occ_W, occ_SP, var_log_C, W-W,    W-SP,   SP-SP,  SP-EMP, EMP-SP, EMP-EMP,var_log_y, gini_y_w, gini_y_ent
TARGET_MOMENTS = [1.704, 0.6035,   0.702, 0.085,  0.26,      0.9853, 0.0089, 0.8366, 0.1349, 0.0842, 0.8873, 0.48,      0.43,     0.25]

XLSX.openxlsx("$(LOCAL_DIR)calibration_res.xlsx",mode="w") do xf
    sheet1 = xf[1]
    labels1=["CODE_NAME", "Weighted sum of deviations", "Sum of abs deviations", "Sum of squared deviations","K/Y","$(TARGET_MOMENTS[1])","Credit/Y","$(TARGET_MOMENTS[2])","occ_W","$(TARGET_MOMENTS[3])","occ_SP","$(TARGET_MOMENTS[4])","var_log_C","$(TARGET_MOMENTS[5])","W-W","$(TARGET_MOMENTS[6])","W-SP","$(TARGET_MOMENTS[7])","SP-SP","$(TARGET_MOMENTS[8])","SP-EMP","$(TARGET_MOMENTS[9])","EMP-SP","$(TARGET_MOMENTS[10])","EMP-EMP","$(TARGET_MOMENTS[11])","var_log_y","$(TARGET_MOMENTS[12])","gini_y_w","$(TARGET_MOMENTS[13])","gini_y_ent","$(TARGET_MOMENTS[14])","R","len_R","W","len_W"]
    XLSX.writetable!(sheet1, zeros(length(labels1)),labels1, anchor_cell=XLSX.CellRef("A1"))
    XLSX.rename!(sheet1,"Deviations")

    XLSX.addsheet!(xf,"Parameters")
    sheet2 = xf[2]
    labels2=["CODE_NAME", "lambda","beta","c_e","rho_m","sigma_eps_m","rho_eps_m_w","sigma_zeta","p_alpha","eta_alpha","mu_m_alpha","rho_alpha_m_w","sigma_eps_w"]
    XLSX.writetable!(sheet2, zeros(length(labels2)),labels2, anchor_cell=XLSX.CellRef("A1"))
end

# externally calibrated parameters of the model's economy (Italy)
#LAMBDA = 1.513028
#BETA = 0.917506
DELTA = 0.1
GAMMA = 0.16
ETA = 0.3
THETA = 0.54

#C_E = 0.007001
#RHO_M = 0.710266
RHO_W = 0.96
#SIGMA_EPS_M = 1.022857
#RHO_EPS_M_W = 0.361605
#SIGMA_ZETA = 0.091089
#P_ALPHA = 0.024445
#ETA_ALPHA = 5.896035
PROB_NODE1_ALPHA = 0.39
#MU_M_ALPHA = -4.817675
#RHO_ALPHA_M_W = 0.181454
SIGMA_ALPHA_W = 0.12
#SIGMA_EPS_W = 0.021947
CRRA = 1.0 # >0.0

W_MIN = 0.18
W_MAX = 0.3

# internally calibrated parameters
#                  1         2         3         4         5           6           7          8         9         10         11            12
#                  lambda    beta      c_e       rho_m     sigma_eps_m rho_eps_m_w sigma_zeta p_alpha   eta_alpha mu_m_alpha rho_alpha_m_w sigma_eps_w
init_int_params = [1.513028, 0.917506, 0.007001, 0.710266, 1.022857,   0.361605,   0.091089,  0.024445, 5.896035, -4.817675, 0.181454,     0.021947]

int_params_min  = [1.001,    0.86,     0.001,    0.70,     0.75,       0.215,      0.09,      0.02,     4.5,      -5.0,      0.04,         0.01]
int_params_max  = [2.000,    0.99,     0.070,    0.99,     1.33,       0.455,      0.28,      0.06,     6.0,      -2.0,      0.21,         0.30]

function params_from_int_params(i_ps)
    #         1        2        3      4      5    6      7        8        9      10           11           12          13          14      15                16          17              18              19         20
    #         lambda,  beta,    delta, gamma, eta, theta, c_e,     rho_m,   rho_w, sigma_eps_m, rho_eps_m_w, sigma_zeta, p_alpha, eta_alpha, prob_node1_alpha, mu_m_alpha, rho_alpha_m_w, sigma_alpha_w, sigma_eps_w, crra

    params = [i_ps[1], i_ps[2], DELTA, GAMMA, ETA, THETA, i_ps[3], i_ps[4], RHO_W, i_ps[5],     i_ps[6],     i_ps[7],    i_ps[8], i_ps[9],   PROB_NODE1_ALPHA, i_ps[10],   i_ps[11],      SIGMA_ALPHA_W, i_ps[12],    CRRA]
    return params
end
function int_params_from_params(ps)
    #                  1      2      3      4      5           6           7          8       9         10         11            12
    #                  lambda beta   c_e    rho_m  sigma_eps_m rho_eps_m_w sigma_zeta p_alpha eta_alpha mu_m_alpha rho_alpha_m_w sigma_eps_w
    init_int_params = [ps[1], ps[2], ps[7], ps[8], ps[10],     ps[11],     ps[12],    ps[13], ps[14],   ps[16],    ps[17],       ps[19]]
    return init_int_params
end
function calculate_deviations(res)
    #                 1        2         3        4        5          6             7             8             9             10            11            12         13        14
    #                 K/Y,     Credit/Y, occ_W,   occ_SP,  var_log_C, W-W,          W-SP,         SP-SP,        SP-EMP,       EMP-SP,       EMP-EMP,      var_log_y, gini_y_w, gini_y_ent
    model_moments  = [res[12], res[13],  res[14], res[15], res[17],   res[18][1,1], res[18][1,2], res[18][2,2], res[18][2,3], res[18][3,2], res[18][3,3], res[19],   res[20],  res[21]]

    deviations = (model_moments .- TARGET_MOMENTS)./TARGET_MOMENTS

    return deviations, model_moments
end
function calculate_aggregate_errors(deviations)
    agg_errs = zeros(3)
    calibration_weights = [1.0,1.0,1000.0,1000.0,0.0,100.0,100.0,100.0,100.0,100.0,100.0,10.0,10.0,10.0]
    calibration_weights ./= sum(calibration_weights)
    agg_errs[1] = sum(abs, deviations.*calibration_weights)*length(deviations)
    agg_errs[2] = sum(abs, deviations)
    agg_errs[3] = sum(deviations.^2)
    return agg_errs
end

# global parameters of the model's code
#                   1           2           3       4
#                gen_tol_x, gen_tol_f, distr_tol, val_tol
GLOBAL_PARAMS = [1e-6, 1e-4, 1e-9, 1e-7]

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [35,3,3,3,6,3]

text_output = false
fig_output = false
calc_add_results = false
CODE_NAME = 0
function func(cur_int_params)
    global CODE_NAME += 1
    println_sameline("Params - $(round.(cur_int_params; digits=6))")
    if !(int_params_min < cur_int_params < int_params_max)
        problematic_i = findfirst(!(int_params_min .< cur_int_params .< int_params_max))
        println_sameline("error: $problematic_i parameter out of bounds")
        return ones(length(TARGET_MOMENTS)).*(10000.0/length(TARGET_MOMENTS))
    end

    SS = steady_state(Inf,Inf, GLOBAL_PARAMS,GLOBAL_APPROX_PARAMS,params_from_int_params(cur_int_params))

    deviations, model_moments = calculate_deviations(SS[1])

    agg_errs = calculate_aggregate_errors(deviations)

    @save "$(LOCAL_DIR)SS_$(CODE_NAME).jld2" SS

    XLSX.openxlsx("$(LOCAL_DIR)calibration_res.xlsx",mode="rw") do xf
        sheet1 = xf[1]
        row1=[CODE_NAME, agg_errs[1], agg_errs[2], agg_errs[3], deviations[1], model_moments[1], deviations[2], model_moments[2], deviations[3], model_moments[3], deviations[4], model_moments[4], deviations[5], model_moments[5], deviations[6], model_moments[6], deviations[7], model_moments[7], deviations[8], model_moments[8], deviations[9], model_moments[9], deviations[10], model_moments[10], deviations[11], model_moments[11], deviations[12], model_moments[12], deviations[13], model_moments[13], deviations[14], model_moments[14], SS[2], SS[1][10], SS[3], SS[1][11]]
        sheet1["A$(CODE_NAME+1)"] = row1

        sheet2 = xf[2]
        row2 = hcat([CODE_NAME],cur_int_params)
        sheet2["A$(CODE_NAME+1)"] = row2
    end
end


candidates = []
push!(candidates, init_int_params)
#=
for i in 1:length(init_int_params)
    for k in [-0.1, -0.05, -0.01, 0.01, 0.05, 0.1]
        int_params = copy(init_int_params)
        percent_dev = k*int_params[i]
        int_params[i] += percent_dev
        int_params[i] = min(max(int_params_min[i], int_params[i]), int_params_max[i])
        push!(candidates, int_params)
    end
end
=#

#=
candidates = []
push!(candidates, init_int_params)
for i in 1:length(init_int_params)-1
    for ik in [-0.05,0.05]
        for j in i:length(init_int_params)
            for jk in [-0.05,0.05]
                int_params = copy(init_int_params)
                percent_dev_i = ik*int_params[i]
                int_params[i] += percent_dev_i
                int_params[i] = min(max(int_params_min[i], int_params[i]), int_params_max[i])
                percent_dev_j = jk*int_params[j]
                int_params[j] += percent_dev_j
                int_params[j] = min(max(int_params_min[j], int_params[j]), int_params_max[j])
                if int_params != init_int_params
                    push!(candidates, int_params)
                end
            end
        end
    end
end
println_sameline(length(candidates))
=#
println_sameline(length(candidates))

for i_ps in 1:length(candidates)
    println_sameline(i_ps)
    @time func(candidates[i_ps])
end
