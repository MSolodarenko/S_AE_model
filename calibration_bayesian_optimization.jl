include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2
using XLSX
using ProgressMeter


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

search_range = [(int_params_min[i], int_params_max[i]) for i in 1:length(init_int_params)]

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

    deviations = (model_moments .- TARGET_MOMENTS)#./TARGET_MOMENTS

    deviations[1:2]   ./= TARGET_MOMENTS[1:2]
    deviations[5]     /= TARGET_MOMENTS[5]
    deviations[12:14] ./= TARGET_MOMENTS[12:14]

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

function f(cur_int_params)

    println_sameline("Params - $(round.(cur_int_params; digits=6))")
    if !(int_params_min < cur_int_params < int_params_max)
        problematic_i = findfirst(!(int_params_min .< cur_int_params .< int_params_max))
        println_sameline("error: $problematic_i parameter out of bounds")
        return 10000.0+sum(cur_int_params.-((int_params_min.+int_params_max)./2))
    end

    is_read_success = true
    CODE_NAME = -1
    result = XLSX.openxlsx("$(LOCAL_DIR)calibration_res.xlsx",mode="rw") do xf
        sheet1 = xf[1]
        sheet2 = xf[2]

        CODE_NAME = XLSX.get_dimension(sheet1).stop.row_number

        for row_i in 2:CODE_NAME
            i_int_params = vec(Float64.(sheet2["B$row_i:M$row_i"]))
            if maximum(abs,(i_int_params - cur_int_params)./(i_int_params + cur_int_params)) < 0.01
                #=
                @load "$(LOCAL_DIR)SS_$(sheet2["A$(row_i)"]).jld2" SS
                deviations, model_moments = calculate_deviations(SS[1])
                agg_errs = calculate_aggregate_errors(deviations)
                =#
                agg_errs = [Float64(sheet1["B$row_i"])]
                if agg_errs[1] < 10000.0
                    print_sameline("error: cur_int_params have been previously calculated in SS_$(sheet1["A$row_i"])) (max diff is < 0.01)")
                    is_read_success = false
                    return agg_errs[1]
                end
            end
        end

        row1=[ CODE_NAME, 10000.0 ]
        sheet1["A$(CODE_NAME+1)"] = row1

        row2 = vcat([CODE_NAME],cur_int_params)
        sheet2["A$(CODE_NAME+1)"] = row2
    end

    if !is_read_success
        return result
    end

    SS = []
    try
        SS = steady_state(Inf,Inf, GLOBAL_PARAMS,GLOBAL_APPROX_PARAMS,params_from_int_params(cur_int_params))
    catch e
        println_sameline(e)
        return 10000.0+sum(cur_int_params.-((int_params_min.+int_params_max)./2))
    end

    deviations, model_moments = calculate_deviations(SS[1])

    agg_errs = calculate_aggregate_errors(deviations)

    println_sameline("Agg_error: $(agg_errs[1])")

    @save "$(LOCAL_DIR)SS_$(CODE_NAME).jld2" SS

    XLSX.openxlsx("$(LOCAL_DIR)calibration_res.xlsx",mode="rw") do xf
        sheet1 = xf[1]
        row1=[CODE_NAME, agg_errs[1], agg_errs[2], agg_errs[3], deviations[1], model_moments[1], deviations[2], model_moments[2], deviations[3], model_moments[3], deviations[4], model_moments[4], deviations[5], model_moments[5], deviations[6], model_moments[6], deviations[7], model_moments[7], deviations[8], model_moments[8], deviations[9], model_moments[9], deviations[10], model_moments[10], deviations[11], model_moments[11], deviations[12], model_moments[12], deviations[13], model_moments[13], deviations[14], model_moments[14], SS[2], SS[1][10], SS[3], SS[1][11]]
        sheet1["A$(CODE_NAME+1)"] = row1
    end
    return agg_errs[1]
end

#=
using BayesianOptimization, GaussianProcesses, Distributions
# Choose as a model an elastic GP with input dimensions 2.
# The GP is called elastic, because data can be appended efficiently.
model = ElasticGPE(12,                            # 2 input dimensions
                   mean = MeanConst(0.),
                   kernel =  Mat52Ard(zeros(12), 0.),
                   logNoise = -2.,
                   capacity = 3000)              # the initial capacity of the GP is 3000 samples.
#set_priors!(model.mean, [Normal(1, 2)])

# Optimize the hyperparameters of the GP using maximum a posteriori (MAP) estimates every 50 steps
modeloptimizer = MAPGPOptimizer(every = 20, noisebounds = [-4, 3],       # bounds of the logNoise
                                kernbounds = [[-3*ones(12); -3], [4*ones(12); 3]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
                                maxeval = 100)


x = []
y = []
XLSX.openxlsx("$(LOCAL_DIR)calibration_res.xlsx",mode="rw") do xf
    sheet1 = xf[1]
    sheet2 = xf[2]

    max_CODE_NAME = XLSX.get_dimension(sheet1).stop.row_number

    @showprogress for i in 2:max_CODE_NAME
        if !(sheet1["B$i"] === missing)
            append!(x, [vec(Float64.(sheet2["B$i:M$i"]))])
            append!(y, -Float64(sheet1["B$i"]))
        end
    end
end
append!(model, hcat(x...), y)

opt = BOpt(f,
           model,
           ExpectedImprovement(),
           modeloptimizer,
           int_params_min,
           int_params_max,
           maxiterations = 300,
           sense = Min,
           initializer_iterations = 0
          )

result = boptimize!(opt)
=#

using BlackBoxOptim
res = bboptimize(f; SearchRange=search_range, Method=:adaptive_de_rand_1_bin)
