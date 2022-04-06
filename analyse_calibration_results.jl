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

    error_pars = [] # total_dev >= 10000
    error_devs = []

    not_converged_pars = [] # len_r or len_w > 1e-4
    not_converged_devs = []

    good_converged_pars = [] # any deviation < 1%
    good_converged_devs = []

    converged_pars = [] # everything else
    converged_devs = []

    @showprogress for row_i in 2:max_CODE_NAME
        cur_par = vec(Float64.(sheet2["A$row_i:M$row_i"]))

        if Float64(sheet1["B$row_i"]) < 10000
            cur_dev = vec(Float64.(sheet1["A$row_i:AJ$row_i"]))

             #interest rate and # wage
            if abs(Float64(sheet1["AH$row_i"])) <= 1e-4 && abs(Float64(sheet1["AJ$row_i"])) <= 1e-4

                is_good = false

                for i in [5,7,9,11,13,27,29,31]
                    if abs(cur_dev[i]) < 0.01
                        is_good = true
                    end
                end

                if is_good
                    append!(good_converged_devs, [cur_dev] )
                    append!(good_converged_pars, [cur_par] )
                else
                    append!(converged_devs, [cur_dev] )
                    append!(converged_pars, [cur_par] )
                end
            else
                append!(not_converged_devs, [cur_dev] )
                append!(not_converged_pars, [cur_par] )
            end
        else
            append!(error_pars, [cur_par])
            append!(error_devs, [[cur_par[1], 125.0, 125.0, 125.0]])
        end
    end

    return [good_converged_devs, converged_devs, not_converged_devs, error_devs], [good_converged_pars, converged_pars, not_converged_pars, error_pars]
end

max_cats = 1
COLORS_LIST = ["green", "blue", "orange", "red"]

ID_LISTs = []

DEV_LISTs = []
DEV_LIST_NAMES = ["WEIGHTED_DEVS", "SUM_OF_ABS_DEVS", "SUM_OF_SQUARED_DEVS"]
#                   1                   2               3
IND_DEV_LISTs = []
IND_DEV_LIST_NAMES = ["K_YS", "CREDIT_YS", "OCC_WS", "OCC_SPS", "VAR_LOG_C", "VAR_LOG_YS", "GINI_Y_W", "GINI_Y_ENT"]
#                       1       2           3           4           5           6           7               8
PAR_LISTs = []
PAR_LIST_NAMES = ["LAMBDAS", "BETAS", "C_ES", "RHO_MS", "SIGMA_EPS_MS", "RHO_EPS_MS", "SIGMA_ZETAS", "P_ALPHAS", "ETA_ALPHAS", "MU_M_ALPHAS", "RHO_ALPHA_MS", "SIGMA_EPS_WS"]
#                   1           2       3       4           5               6           7               8           9               10              11

for cat = 1:max_cats
    ID = [PARAMETERS[cat][i][1] for i in 1:length(PARAMETERS[cat])]
    append!(ID_LISTs, [ID])

    WEIGHTED_DEVS = [DEVIATIONS[cat][i][2] for i in 1:length(DEVIATIONS[cat])]
    SUM_OF_ABS_DEVS = [DEVIATIONS[cat][i][3] for i in 1:length(DEVIATIONS[cat])]
    SUM_OF_SQUARED_DEVS = [DEVIATIONS[cat][i][4] for i in 1:length(DEVIATIONS[cat])]

    DEV_LIST = [WEIGHTED_DEVS, SUM_OF_ABS_DEVS, SUM_OF_SQUARED_DEVS]
    append!(DEV_LISTs, [DEV_LIST])

    if cat < 4
        K_YS = [DEVIATIONS[cat][i][5] for i in 1:length(DEVIATIONS[cat])]
        CREDIT_YS = [DEVIATIONS[cat][i][7] for i in 1:length(DEVIATIONS[cat])]
        OCC_WS = [DEVIATIONS[cat][i][9] for i in 1:length(DEVIATIONS[cat])]
        OCC_SPS = [DEVIATIONS[cat][i][11] for i in 1:length(DEVIATIONS[cat])]
        VAR_LOG_C = [DEVIATIONS[cat][i][13] for i in 1:length(DEVIATIONS[cat])]
        VAR_LOG_YS = [DEVIATIONS[cat][i][27] for i in 1:length(DEVIATIONS[cat])]
        GINI_Y_W = [DEVIATIONS[cat][i][29] for i in 1:length(DEVIATIONS[cat])]
        GINI_Y_ENT = [DEVIATIONS[cat][i][31] for i in 1:length(DEVIATIONS[cat])]

        IND_DEV_LIST = [K_YS, CREDIT_YS, OCC_WS, OCC_SPS, VAR_LOG_C, VAR_LOG_YS, GINI_Y_W, GINI_Y_ENT]
        append!(IND_DEV_LISTs, [IND_DEV_LIST])
    end

    LAMBDAS = [PARAMETERS[cat][i][2] for i in 1:length(PARAMETERS[cat])]
    BETAS = [PARAMETERS[cat][i][3] for i in 1:length(PARAMETERS[cat])]
    C_ES = [PARAMETERS[cat][i][4] for i in 1:length(PARAMETERS[cat])]
    RHO_MS = [PARAMETERS[cat][i][5] for i in 1:length(PARAMETERS[cat])]
    SIGMA_EPS_MS = [PARAMETERS[cat][i][6] for i in 1:length(PARAMETERS[cat])]
    RHO_EPS_MS = [PARAMETERS[cat][i][7] for i in 1:length(PARAMETERS[cat])]
    SIGMA_ZETAS = [PARAMETERS[cat][i][8] for i in 1:length(PARAMETERS[cat])]
    P_ALPHAS = [PARAMETERS[cat][i][9] for i in 1:length(PARAMETERS[cat])]
    ETA_ALPHAS = [PARAMETERS[cat][i][10] for i in 1:length(PARAMETERS[cat])]
    MU_M_ALPHAS = [PARAMETERS[cat][i][11] for i in 1:length(PARAMETERS[cat])]
    RHO_ALPHA_MS = [PARAMETERS[cat][i][12] for i in 1:length(PARAMETERS[cat])]
    SIGMA_EPS_WS = [PARAMETERS[cat][i][13] for i in 1:length(PARAMETERS[cat])]

    PAR_LIST = [LAMBDAS, BETAS, C_ES, RHO_MS, SIGMA_EPS_MS, RHO_EPS_MS, SIGMA_ZETAS, P_ALPHAS, ETA_ALPHAS, MU_M_ALPHAS, RHO_ALPHA_MS, SIGMA_EPS_WS]
    append!(PAR_LISTs, [PAR_LIST])
end
#throw(error)

plot_DEV = Array{Any}(undef,length(DEV_LISTs[1]))
for dev_i in 1#:length(DEV_LIST)
    plt = plot()
    for cat = 1:max_cats
        scatter!(plt, ID_LISTs[cat], DEV_LISTs[cat][dev_i], ylabel=DEV_LIST_NAMES[dev_i], color=COLORS_LIST[cat], legend=false)
        #display(plt)

        min_iter = argmin(DEV_LISTs[cat][dev_i])
        max_iter = argmax(DEV_LISTs[cat][dev_i])
        #display((DEV_LIST_NAMES[dev_i], "min:", ID_LISTs[cat][min_iter], DEV_LISTs[cat][dev_i][min_iter], "max:", ID_LISTs[cat][max_iter], DEV_LISTs[cat][dev_i][max_iter]))

    end
    plot_DEV[dev_i] = plt
    #sleep(10)
end
#throw(error)

plot_PAR_DEV = Array{Any}(undef, length(PAR_LISTs[1]), length(DEV_LISTs[1]))
for dev_i in 1#:length(DEV_LISTs[1])
    for par_i in 1:length(PAR_LISTs[1])
        plt = plot()
        for cat = 1:max_cats
            scatter!(plt, PAR_LISTs[cat][par_i], DEV_LISTs[cat][dev_i], color=COLORS_LIST[cat], xlabel=PAR_LIST_NAMES[par_i], ylabel=DEV_LIST_NAMES[dev_i], legend=false)
        end
        plot_PAR_DEV[par_i,dev_i] = plt
        #display(plot_PAR_DEV[par_i,dev_i])
        #sleep(10)
    end
end
#throw(error)

plot_PAR_IND_DEV = Array{Any}(undef, length(PAR_LISTs[1]), length(IND_DEV_LISTs[1]))
for dev_i in 1:length(IND_DEV_LISTs[1])
    for par_i in 1:length(PAR_LISTs[1])
        plt = plot()
        for cat = 1:min(max_cats,3)
            scatter!(plt, PAR_LISTs[cat][par_i], IND_DEV_LISTs[cat][dev_i], color=COLORS_LIST[cat], xlabel=PAR_LIST_NAMES[par_i], ylabel=IND_DEV_LIST_NAMES[dev_i], legend=false)

            min_iter = argmin(IND_DEV_LISTs[cat][dev_i])
            max_iter = argmax(IND_DEV_LISTs[cat][dev_i])
            zero_iter = argmin(abs.(IND_DEV_LISTs[cat][dev_i]))
            #display((IND_DEV_LIST_NAMES[dev_i], "down:", ID_LISTs[cat][min_iter], IND_DEV_LISTs[cat][dev_i][min_iter], "up:", ID_LISTs[cat][max_iter], IND_DEV_LISTs[cat][dev_i][max_iter], "zero:", ID_LISTs[cat][zero_iter], IND_DEV_LISTs[cat][dev_i][zero_iter]))

            min_iter = argmin(PAR_LISTs[cat][dev_i])
            max_iter = argmax(PAR_LISTs[cat][dev_i])
            #display((PAR_LIST_NAMES[dev_i], "down:", ID_LISTs[cat][min_iter], PAR_LISTs[cat][dev_i][min_iter], "up:", ID_LISTs[cat][max_iter], PAR_LISTs[cat][dev_i][max_iter]))

        end
        plot_PAR_IND_DEV[par_i,dev_i] = plt
        #display(plot_PAR_IND_DEV[par_i,dev_i])
        #sleep(10)

    end

end
#throw(error)

for par_i in 1:length(PAR_LISTs[1])
    display(plot(Tuple(plot_PAR_IND_DEV[par_i,:])..., layout=length(plot_PAR_IND_DEV[par_i,:])))
end
for dev_i in 1:length(IND_DEV_LISTs[1])
    display(plot(Tuple(plot_PAR_IND_DEV[:,dev_i])..., layout=length(plot_PAR_IND_DEV[:,dev_i])))
end

# beta to Capital-to-Output
#display(plot(plot_PAR_IND_DEV[2,1]))

# lambda to Credit-to-Output
#display(plot(plot_PAR_IND_DEV[1,2]))
