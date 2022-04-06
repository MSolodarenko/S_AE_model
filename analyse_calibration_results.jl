using Plots
using JLD2
using XLSX
using ProgressMeter

LOCAL_DIR = "$(@__DIR__)/Results/Calibration/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Calibration\\"
end
mkpath(LOCAL_DIR)

DEVIATIONS, PARAMETERS = XLSX.openxlsx("$(LOCAL_DIR)calibration_res.xlsx",mode="rw") do xf
    sheet1 = xf[1]
    sheet2 = xf[2]

    max_CODE_NAME = XLSX.get_dimension(sheet1).stop.row_number

    devs = []
    pars = []

    @showprogress for row_i in 2:max_CODE_NAME
        if Float64(sheet1["B$row_i"]) < 10000
            cur_dev = vec(Float64.(sheet1["A$row_i:AJ$row_i"]))
            cur_par = vec(Float64.(sheet2["A$row_i:M$row_i"]))

            if abs(Float64(sheet1["AH$row_i"])) <= 1e-4 #interest rate
                if abs(Float64(sheet1["AJ$row_i"])) <= 1e-4 # wage

                    is_appropriate = true
                    #=
                    for i in [5,7,9,11,13,27,29,31]
                        if abs(cur_dev[i]) >= 0.5
                            is_appropriate = false
                        end
                    end
                    =#
                    if abs(cur_dev[5]) >= 0.05
                        is_appropriate = false
                    end
                    if abs(cur_dev[7]) >= 0.05
                        is_appropriate = false
                    end


                    if is_appropriate
                        append!(devs, [cur_dev] )
                        append!(pars, [cur_par] )
                    end
                end
            end
        end
    end

    return devs, pars
end

WEIGHTED_DEVS = [DEVIATIONS[i][2] for i in 1:length(DEVIATIONS)]
SUM_OF_ABS_DEVS = [DEVIATIONS[i][3] for i in 1:length(DEVIATIONS)]
SUM_OF_SQUARED_DEVS = [DEVIATIONS[i][4] for i in 1:length(DEVIATIONS)]

DEV_LIST = [WEIGHTED_DEVS, SUM_OF_ABS_DEVS, SUM_OF_SQUARED_DEVS]
DEV_LIST_NAMES = ["WEIGHTED_DEVS", "SUM_OF_ABS_DEVS", "SUM_OF_SQUARED_DEVS"]

K_YS = [DEVIATIONS[i][5] for i in 1:length(DEVIATIONS)]
CREDIT_YS = [DEVIATIONS[i][7] for i in 1:length(DEVIATIONS)]
OCC_WS = [DEVIATIONS[i][9] for i in 1:length(DEVIATIONS)]
OCC_SPS = [DEVIATIONS[i][11] for i in 1:length(DEVIATIONS)]
VAR_LOG_C = [DEVIATIONS[i][13] for i in 1:length(DEVIATIONS)]
VAR_LOG_YS = [DEVIATIONS[i][27] for i in 1:length(DEVIATIONS)]
GINI_Y_W = [DEVIATIONS[i][29] for i in 1:length(DEVIATIONS)]
GINI_Y_ENT = [DEVIATIONS[i][31] for i in 1:length(DEVIATIONS)]

IND_DEV_LIST = [K_YS, CREDIT_YS, OCC_WS, OCC_SPS, VAR_LOG_C, VAR_LOG_YS, GINI_Y_W, GINI_Y_ENT]
IND_DEV_LIST_NAMES = ["K_YS", "CREDIT_YS", "OCC_WS", "OCC_SPS", "VAR_LOG_C", "VAR_LOG_YS", "GINI_Y_W", "GINI_Y_ENT"]
#                       1       2           3           4           5           6           7               8

LAMBDAS = [PARAMETERS[i][2] for i in 1:length(DEVIATIONS)]
BETAS = [PARAMETERS[i][3] for i in 1:length(DEVIATIONS)]
C_ES = [PARAMETERS[i][4] for i in 1:length(DEVIATIONS)]
RHO_MS = [PARAMETERS[i][5] for i in 1:length(DEVIATIONS)]
SIGMA_EPS_MS = [PARAMETERS[i][6] for i in 1:length(DEVIATIONS)]
RHO_EPS_MS = [PARAMETERS[i][7] for i in 1:length(DEVIATIONS)]
SIGMA_ZETAS = [PARAMETERS[i][8] for i in 1:length(DEVIATIONS)]
P_ALPHAS = [PARAMETERS[i][9] for i in 1:length(DEVIATIONS)]
ETA_ALPHAS = [PARAMETERS[i][10] for i in 1:length(DEVIATIONS)]
MU_M_ALPHAS = [PARAMETERS[i][11] for i in 1:length(DEVIATIONS)]
RHO_ALPHA_MS = [PARAMETERS[i][12] for i in 1:length(DEVIATIONS)]
SIGMA_EPS_WS = [PARAMETERS[i][13] for i in 1:length(DEVIATIONS)]

PAR_LIST = [LAMBDAS, BETAS, C_ES, RHO_MS, SIGMA_EPS_MS, RHO_EPS_MS, SIGMA_ZETAS, P_ALPHAS, ETA_ALPHAS, MU_M_ALPHAS, RHO_ALPHA_MS, SIGMA_EPS_WS]
PAR_LIST_NAMES = ["LAMBDAS", "BETAS", "C_ES", "RHO_MS", "SIGMA_EPS_MS", "RHO_EPS_MS", "SIGMA_ZETAS", "P_ALPHAS", "ETA_ALPHAS", "MU_M_ALPHAS", "RHO_ALPHA_MS", "SIGMA_EPS_WS"]
#                   1           2       3       4           5               6           7               8           9               10              11
#throw(error)

plot_DEV = Array{Any}(undef,length(DEV_LIST))
for dev_i in 1:length(DEV_LIST)
    plot_DEV[dev_i] = scatter(DEV_LIST[dev_i], ylabel=DEV_LIST_NAMES[dev_i], legend=false)
    display(plot_DEV[dev_i])
    #sleep(10)
    min_iter = argmin(DEV_LIST[dev_i])
    max_iter = argmax(DEV_LIST[dev_i])
    display((DEV_LIST_NAMES[dev_i], "min:", DEVIATIONS[min_iter][1], DEV_LIST[dev_i][min_iter], "max:", DEVIATIONS[max_iter][1], DEV_LIST[dev_i][max_iter]))
end

plot_PAR_DEV = Array{Any}(undef, length(PAR_LIST), length(DEV_LIST))
for dev_i in 1:length(DEV_LIST)
    for par_i in 1:length(PAR_LIST)
        plot_PAR_DEV[par_i,dev_i] = scatter(PAR_LIST[par_i], DEV_LIST[dev_i], xlabel=PAR_LIST_NAMES[par_i], ylabel=DEV_LIST_NAMES[dev_i], legend=false)
        #display(plot_PAR_DEV[par_i,dev_i])
        #sleep(10)
    end
end

plot_PAR_IND_DEV = Array{Any}(undef, length(PAR_LIST), length(IND_DEV_LIST))
for dev_i in 1:length(IND_DEV_LIST)
    for par_i in 1:length(PAR_LIST)
        plot_PAR_IND_DEV[par_i,dev_i] = scatter(PAR_LIST[par_i], IND_DEV_LIST[dev_i], xlabel=PAR_LIST_NAMES[par_i], ylabel=IND_DEV_LIST_NAMES[dev_i], legend=false)
        #display(plot_PAR_IND_DEV[par_i,dev_i])
        #sleep(10)
    end
    min_iter = argmin(IND_DEV_LIST[dev_i])
    max_iter = argmax(IND_DEV_LIST[dev_i])
    zero_iter = argmin(abs.(IND_DEV_LIST[dev_i]))
    display((IND_DEV_LIST_NAMES[dev_i], "down:", DEVIATIONS[min_iter][1], IND_DEV_LIST[dev_i][min_iter], "up:", DEVIATIONS[max_iter][1], IND_DEV_LIST[dev_i][max_iter], "zero:", DEVIATIONS[zero_iter][1], IND_DEV_LIST[dev_i][zero_iter]))

end

for par_i in 1:length(PAR_LIST)
    display(plot(Tuple(plot_PAR_IND_DEV[par_i,:])..., layout=length(plot_PAR_IND_DEV[par_i,:])))
end
for dev_i in 1:length(IND_DEV_LIST)
    display(plot(Tuple(plot_PAR_IND_DEV[:,dev_i])..., layout=length(plot_PAR_IND_DEV[:,dev_i])))
end

# beta to Capital-to-Output
display(plot(plot_PAR_IND_DEV[2,1]))

# lambda to Credit-to-Output
display(plot(plot_PAR_IND_DEV[1,2]))
