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
using SchumakerSpline

CODENAME = "SS_2092_69"
country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_big_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_big_grid\\$(country)_$(CODENAME)\\"
end

function create_plot(X,XLABEL::String,Y,YLABEL::String,Y1,Y2, IS_Y_PERCENTAGE::Bool=false, OCCUPATION::String="H", TICKSFONTSIZE::Int64=9)
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    elseif OCCUPATION=="ENT"
        COLOR="brown"
    end

    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=T,length=NUM_XTICKS))))
    if length(X) <= 10
        XTICKS = Int64.(round.(collect(range(1;stop=T,step=1))))
    end
    XTICKS = [1; XTICKS[2:end]]

    NUM_YTICKS = 10
    YLIMS1 = minimum([Y; Y1; Y2])
    YLIMS2 = maximum([Y; Y1; Y2])
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    Y12low = minimum([Y1; Y2])
    Y12high = maximum([Y1; Y2])
    SHARElow = (Y12low-YLIMS1)/(YLIMS2-YLIMS1)
    SHAREmid = (Y12high-Y12low)/(YLIMS2-YLIMS1)
    SHAREup = (YLIMS2-Y12high)/(YLIMS2-YLIMS1)
    #display([SHAREup; SHAREmid; SHARElow])

    LOWLENGTH = Int64(round(NUM_YTICKS*SHARElow))
    YTICKSlow = [YLIMS1]
    if LOWLENGTH>1
        YTICKSlow = collect(range(YLIMS1; stop=Y12low, length=LOWLENGTH))
        STEP = abs(YTICKSlow[2]-YTICKSlow[1])
        YTICKSlow = YTICKSlow[1:end-1]
    end

    MIDLENGTH = Int64(round(NUM_YTICKS*SHAREmid))
    YTICKSmid = []
    if MIDLENGTH>1
        YTICKSmid = collect(range(Y12low; stop=Y12high, length=MIDLENGTH))
        STEP = abs(YTICKSmid[2]-YTICKSmid[1])
        YTICKSmid = YTICKSmid[2:end-1]
        YLIMS1 = Y12low-LOWLENGTH*STEP
        YTICKSlow = collect(range(Y12low; stop=YLIMS1, step=-STEP))[end:-1:1]
    end

    HIGHLENGTH = Int64(round(NUM_YTICKS*SHAREup))
    YTICKShigh = []
    if LOWLENGTH+MIDLENGTH<3
        YTICKShigh = collect(range(Y12high; stop=YLIMS2, length=HIGHLENGTH))
    else
        if HIGHLENGTH==0 && YLIMS2!=Y12high
            YTICKShigh = [YLIMS2]
        else
            YLIMS2 = Y12high+HIGHLENGTH*STEP
            YTICKShigh = collect(range(Y12high; stop=YLIMS2, step=STEP))
        end
    end
    #YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))
    YTICKS = [YTICKSlow; YTICKSmid; YTICKShigh]

    if YTICKS[1] > minimum([Y; Y1; Y2])
        YTICKS[1] = minimum([Y; Y1; Y2])
        YLIMS1 = minimum([Y; Y1; Y2])
    end
    if YTICKS[end] < maximum([Y; Y1; Y2])
        YTICKS[end] = maximum([Y; Y1; Y2])
        YLIMS2 = maximum([Y; Y1; Y2])
    end

    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
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

    Y1line = ones(length(Y)).*Y1
    Y2line = ones(length(Y)).*Y2
    if length(Y) == 50
        Y1line[10+1:end] .= NaN
        Y2line[1:40-1] .= NaN
    elseif length(Y) == 10
        Y1line[3+1:end] .= NaN
        Y2line[1:8-1] .= NaN
    end
    plt = plot(collect(X), collect.([Y, Y1line,Y2line]),
                    #color=[COLOR "green" "red"],
                    color=[COLOR COLOR COLOR],
                    linestyle=[:solid :dot :dash],
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    xtickfontsize=TICKSFONTSIZE,
                    ytickfontsize=TICKSFONTSIZE,
                    xticks=XTICKS,
                    yticks = YTICKS,
                    ylims = YLIMS )
    # text annotation
    #=if Y1 != Y2
        if IS_Y_PERCENTAGE
            TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1*100;digits=DIGITS+1))%)"
        else
            TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1;digits=DIGITS+1)))"
        end
        POS = Y1+YLIMMARGIN*1.25
        if Y1 < Y2
            POS = Y1-YLIMMARGIN*1.25
            if POS < YLIMS1-YLIMMARGIN
                POS = Y1+YLIMMARGIN*1.25
            end
        end
        annotate!([X[end]], POS, text(TEXT1, COLOR, :right, 7))

        if IS_Y_PERCENTAGE
            TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2*100;digits=DIGITS+1))%)"
        else
            TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2;digits=DIGITS+1)))"
        end
        POS = Y2-YLIMMARGIN*1.25
        if POS < YLIMS1-YLIMMARGIN || Y1 < Y2
            POS = Y2+YLIMMARGIN*1.25
        end
        annotate!([X[end]], POS, text(TEXT2, COLOR, :right, 7))
    end=#
    return plt
end
function create_plot(X,XLABEL::String,Y,YLABEL::String, IS_Y_PERCENTAGE::Bool=false, OCCUPATION::String="H", TICKSFONTSIZE::Int64=9)
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    elseif OCCUPATION=="ENT"
        COLOR="brown"
    end

    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=T,length=NUM_XTICKS))))
    if length(X) <= 10
        XTICKS = Int64.(round.(collect(range(1;stop=T,step=1))))
    end
    XTICKS = [1; XTICKS[2:end]]

    NUM_YTICKS = 10
    YLIMS1 = minimum([Y; 0.0])
    YLIMS2 = maximum([Y; 0.0])
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YTICKSlow = []
    YTICKShigh = []
    numYTICKSlow = Int64(round(NUM_YTICKS*(0.0-YLIMS1)/(YLIMS2-YLIMS1)))
    YTICKSlow = []
    if numYTICKSlow > 0
        YTICKSlow = collect(range(YLIMS1; stop=0.0, length=numYTICKSlow+1))[1:end-1]
    end
    numYTICKShigh= NUM_YTICKS-numYTICKSlow#Int64(round(NUMYTICKS*(YLIMS2-0.0)/(YLIMS2-YLIMS1)))
    YTICKShigh = collect(range(0.0; stop=YLIMS2, length=numYTICKShigh+1))[2:end]
    YTICKS = [YTICKSlow; 0.0; YTICKShigh]
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
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

    plt = plot(collect(X), collect.([Y, zeros(length(Y))]),
                    #color=[COLOR "green" "red"],
                    color=[COLOR COLOR COLOR],
                    linestyle=[:solid :dot :dash],
                    legend=false,
                    xlabel=XLABEL,
                    ylabel=YLABEL,
                    xtickfontsize=TICKSFONTSIZE,
                    ytickfontsize=TICKSFONTSIZE,
                    xticks=XTICKS,
                    yticks = YTICKS,
                    ylims = YLIMS )
    # text annotation
    #=if Y1 != Y2
        if IS_Y_PERCENTAGE
            TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1*100;digits=DIGITS+1))%)"
        else
            TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1;digits=DIGITS+1)))"
        end
        POS = Y1+YLIMMARGIN*1.25
        if Y1 < Y2
            POS = Y1-YLIMMARGIN*1.25
            if POS < YLIMS1-YLIMMARGIN
                POS = Y1+YLIMMARGIN*1.25
            end
        end
        annotate!([X[end]], POS, text(TEXT1, COLOR, :right, 7))

        if IS_Y_PERCENTAGE
            TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2*100;digits=DIGITS+1))%)"
        else
            TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2;digits=DIGITS+1)))"
        end
        POS = Y2-YLIMMARGIN*1.25
        if POS < YLIMS1-YLIMMARGIN || Y1 < Y2
            POS = Y2+YLIMMARGIN*1.25
        end
        annotate!([X[end]], POS, text(TEXT2, COLOR, :right, 7))
    end=#
    return plt
end
function create_combined_plot(X,XLABEL::String,Ys,YLABELs,YLABEL,Y1s,Y2s, IS_Y_PERCENTAGE::Bool=false, OCCUPATION=["SP","EMP"], LEGENDPOS=false)
    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=T,length=NUM_XTICKS))))

    XTICKS = [1; XTICKS[2:end]]

    NUMYTICKS = 10
    YLIMS1 = minimum([minimum.(Ys); minimum.(Y1s); minimum.(Y2s)])
    YLIMS2 = maximum([maximum.(Ys); maximum.(Y1s); maximum.(Y2s)])
    YTICKSlow = []
    YTICKShigh = []
    # try
    #     YTICKSlow = collect(range(YLIMS1; stop=0.0, length=Int64(round(NUMYTICKS*(0.0-YLIMS1)/(YLIMS2-YLIMS1)))))
    #     STEP = YTICKSlow[2]-YTICKSlow[1]
    #     YTICKShigh = collect(range(0.0; stop=YLIMS2+STEP, step=STEP))[2:end]
    #     YLIMS2 = YTICKShigh[end]
    # catch e
    #     YTICKShigh = collect(range(0.0; stop=YLIMS2, length=Int64(round(NUMYTICKS*(YLIMS2-0.0)/(YLIMS2-YLIMS1)))))
    #     STEP = YTICKShigh[2]-YTICKShigh[1]
    #     YTICKSlow = collect(range(0.0; stop=YLIMS1-STEP, step=-STEP))[end:-1:2]
    # end
    numYTICKSlow = Int64(round(NUMYTICKS*(0.0-YLIMS1)/(YLIMS2-YLIMS1)))
    YTICKSlow = []
    if numYTICKSlow > 0
        YTICKSlow = collect(range(YLIMS1; stop=0.0, length=numYTICKSlow+1))[1:end-1]
    end
    numYTICKShigh= NUMYTICKS-numYTICKSlow#Int64(round(NUMYTICKS*(YLIMS2-0.0)/(YLIMS2-YLIMS1)))
    YTICKShigh = collect(range(0.0; stop=YLIMS2, length=numYTICKShigh+1))[2:end]
    YTICKS = [YTICKSlow; 0.0; YTICKShigh]
    #YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUMYTICKS))
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)

    COLORS=[]
    for y_i = 1:length(Ys)
        if OCCUPATION[y_i]=="W"
            push!(COLORS,"purple")
        elseif OCCUPATION[y_i]=="SP"
            push!(COLORS,"red")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"green")
        elseif OCCUPATION[y_i]=="ENT"
            push!(COLORS,"brown")
        else
            push!(COLORS,"blue")
        end

    end

    if IS_Y_PERCENTAGE
        DIGITS = Int(max(2, round(log10(0.01/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]))+0.5;digits=0) ))
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = plot()
    for y_i = 1:length(Ys)
        if YLABELs[y_i] == "W"
            YLABELs[y_i] = "Workers"
        elseif YLABELs[y_i] == "SP"
            YLABELs[y_i] = "Sole Prop."
        elseif YLABELs[y_i] == "EMP"
            YLABELs[y_i] = "Employers"
        elseif YLABELs[y_i] == "ENT"
            YLABELs[y_i] = "Entrepreneurs"
        end
        plot!(plt,collect(X), collect.([Ys[y_i], zeros(length(Ys[y_i]))]),
                        #color=[COLORS[y_i] "green" "red"],
                        color=[COLORS[y_i] "black"],
                        linestyle=[:solid :dot],
                        legend=LEGENDPOS,
                        xlabel=XLABEL,
                        label=[YLABELs[y_i] ""],
                        ylabel = YLABEL,
                        xticks = XTICKS,
                        yticks = YTICKS,
                        ylims = YLIMS )
        # text annotation
        #=if Y1s[y_i] != Y2s[y_i]
            if IS_Y_PERCENTAGE
                TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1s[y_i]*100;digits=DIGITS+1))%)"
            else
                TEXT1 = "Economy at λ=$(round(lambda_1;digits=2)) ($(round(Y1s[y_i];digits=DIGITS+1)))"
            end
            POS = Y1s[y_i]+YLIMMARGIN*1.25
            if Y1s[y_i] < Y2s[y_i]
                POS = Y1s[y_i]-YLIMMARGIN*1.25
                if POS < YLIMS1-YLIMMARGIN
                    POS = Y1s[y_i]+YLIMMARGIN*1.25
                end
            end
            annotate!([X[end]], POS, text(TEXT1, COLORS[y_i], :right, 7))

            if IS_Y_PERCENTAGE
                TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2s[y_i]*100;digits=DIGITS+1))%)"
            else
                TEXT2 = "Economy at λ=$(round(lambda_2;digits=2)) ($(round(Y2s[y_i];digits=DIGITS+1)))"
            end
            POS = Y2s[y_i]-YLIMMARGIN*1.25
            if POS < YLIMS1-YLIMMARGIN || Y1s[y_i] < Y2s[y_i]
                POS = Y2s[y_i]+YLIMMARGIN*1.25
            end
            annotate!([X[end]], POS, text(TEXT2, COLORS[y_i], :right, 7))
        end=#
    end
    return plt
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

        p = Progress((T*number_asset_grid*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes), dt=0.5,
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

    # compute cumulative conditional occupational choice probability after reform
    cum_density_distr_w = copy(init_density_distr_w)
    cum_density_distr_sp = copy(init_density_distr_sp)
    cum_density_distr_emp = copy(init_density_distr_emp)

    cum_occ_trans_matrix = zeros(T+1,3,3)

    # inequality measures for occupations fixed in time t=1 after reform
    function compute_household_measures(policy, r, w, lambda, a_grid)

        occ_choice, income, earnings, = compute_income_profile(a_grid,number_asset_grid,r, w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

        consumption = income .- policy
        wealth = ones(size(init_density_distr)).*a_grid
        capital_income = wealth.*r
        income = income .- wealth#earnings.+capital_income#
        savings = policy .- wealth
        #consumption = income.+wealth .- policy

        w_choice = Float64.(occ_choice.==1.0)
        sp_choice = Float64.(occ_choice.==2.0)
        emp_choice = Float64.(occ_choice.==3.0)

        return consumption, earnings, income, wealth, capital_income, savings, occ_choice,w_choice,sp_choice,emp_choice
    end
    function compute_inequality_measures(measure, distr)
        mean = sum(measure.*distr)

        #variance = sum(((measure.-mean).^2).*distr)/sum(distr)
        mean /= sum(distr)
        #=
        ddc = distr./sum(distr)
        ddc_vec = vec(ddc)
        measure_vec = vec(measure)
        index_non_zero = findall(x-> x>1e-5,ddc_vec)
        ddc_vec_non_zero = ddc_vec[index_non_zero]
        measure_vec_non_zero = measure_vec[index_non_zero]
        gini = sum([ ddc_vec_non_zero[y_i]*ddc_vec_non_zero[y_j]*abs(measure_vec_non_zero[y_i]-measure_vec_non_zero[y_j]) for y_i=1:length(ddc_vec_non_zero), y_j=1:length(ddc_vec_non_zero) ]) / (2*sum(measure_vec_non_zero.*ddc_vec_non_zero))
        =#
        return [mean, 0.0,0.0]#variance, gini]
    end
    # time periods, [consumption, earnings, income, wealth, capital_income, savings], [W,SP,EMP], [mean, var, gini]
    cum_household_measures = zeros(T+1,6,3, 3)
    # time periods, [consumption, earnings, income, wealth, capital_income, savings], [W,SP,EMP]*[W,SP,EMP], [mean, var, gini]
    cum_trans_household_measures = zeros(T+1,6,3,3, 3)

    # compute conditional occupational choice probability after reform
    density_distr_w = copy(init_density_distr_w)
    density_distr_sp = copy(init_density_distr_sp)
    density_distr_emp = copy(init_density_distr_emp)
    occ_trans_matrix = zeros(T+1,3,3)

    # Compute conditional occupational choice probability for 50 years in a steady state pre-reform
    ss_density_distr = ss_1[1][5]
    ss_occ_choice = ss_1[1][22]
    ss_w_choice = Float64.(ss_occ_choice.==1.0)
    ss_sp_choice = Float64.(ss_occ_choice.==2.0)
    ss_emp_choice = Float64.(ss_occ_choice.==3.0)
    ss_init_density_distr_w = ss_density_distr.*ss_w_choice
    ss_init_density_distr_sp = ss_density_distr.*ss_sp_choice
    ss_init_density_distr_emp = ss_density_distr.*ss_emp_choice
    ss_policy = ss_1[1][4]
    ss_a1_indices = ss_1[1][39]
    ss_lottery_prob_1 = ss_1[1][40]
    ss_lottery_prob_2 = ss_1[1][41]

    ss_cum_density_distr_w = copy(ss_init_density_distr_w)
    ss_cum_density_distr_sp = copy(ss_init_density_distr_sp)
    ss_cum_density_distr_emp = copy(ss_init_density_distr_emp)
    ss_cum_occ_trans_matrix = zeros(T+1,3,3)

    ss_consumption, ss_earnings, ss_income, ss_wealth, ss_capital_income, ss_savings,  = compute_household_measures(ss_policy, ss_1[1][44], ss_1[1][45], lambda_s[1], ss_1[1][3])
    # time periods, [consumption, earnings, income, wealth, savings], [W,SP,EMP], [mean, var, gini]
    ss_cum_household_measures = zeros(T+1,6,3, 3)
    # time periods, [consumption, earnings, income, wealth, savings], [W,SP,EMP]*[W,SP,EMP], [mean, var, gini]
    ss_cum_trans_household_measures = zeros(T+1,6,3,3, 3)

    # t = 0 & iterator = 1
    ss_cum_occ_trans_matrix[1,1,1] = 1.0
    ss_cum_occ_trans_matrix[1,2,2] = 1.0
    ss_cum_occ_trans_matrix[1,3,3] = 1.0

    Threads.@threads for (m,o) = collect(Iterators.product(1:6,1:3))
        dd = ss_cum_density_distr_w
        if o == 2
            dd = ss_cum_density_distr_sp
        elseif o == 3
            dd = ss_cum_density_distr_emp
        end
        measure = ss_consumption
        if m==2
            measure = ss_earnings
        elseif m==3
            measure = ss_income
        elseif m==4
            measure = ss_wealth
        elseif m==5
            measure = ss_capital_income
        elseif m==6
            measure = ss_savings
        end
        ss_cum_household_measures[1,m,o,:] .= compute_inequality_measures(measure, dd)
        for o2 = 1:3
            occ2 = ss_w_choice
            if o2 == 2
                occ2 = ss_sp_choice
            elseif o2 == 3
                occ2 = ss_emp_choice
            end
            if sum(dd.*occ2) == 0.0
                ss_cum_trans_household_measures[1,m,o,o2,:] .= 0.0
            else
                ss_cum_trans_household_measures[1,m,o,o2,:] .= compute_inequality_measures(measure.*occ2,dd.*occ2)
            end
            #println_sameline(sum(dd.*occ2))
        end
    end

    cum_occ_trans_matrix[1,1,1] = 1.0
    cum_occ_trans_matrix[1,2,2] = 1.0
    cum_occ_trans_matrix[1,3,3] = 1.0

    cum_household_measures[1,:,:,:,:] = copy(ss_cum_household_measures[1,:,:,:,:])
    cum_trans_household_measures[1,:,:,:,:] = copy(ss_cum_trans_household_measures[1,:,:,:,:])

    occ_trans_matrix[1,:,:] = copy(ss_1[1][18])

    # t = 1 & iterator = 2

    ss_cum_density_distr_w_pr = copy(ss_cum_density_distr_w).*0.0
    ss_cum_density_distr_sp_pr = copy(ss_cum_density_distr_sp).*0.0
    ss_cum_density_distr_emp_pr = copy(ss_cum_density_distr_emp).*0.0

    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        for a_i in 1:number_asset_grid
            ss_j_1 = ss_a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

            ss_cum_density_distr_w_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            ss_cum_density_distr_sp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            ss_cum_density_distr_emp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

            if ss_j_1 != number_asset_grid
                ss_cum_density_distr_w_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_sp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_emp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

            end
        end
    end

    ss_cum_density_distr_marginal_assets_w = sum(sum(sum(sum(ss_cum_density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
    ss_cum_density_distr_marginal_assets_sp = sum(sum(sum(sum(ss_cum_density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
    ss_cum_density_distr_marginal_assets_emp = sum(sum(sum(sum(ss_cum_density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

    Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
        ss_cum_temp_w = sum(ss_cum_density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
        ss_cum_temp_w = ss_cum_temp_w*P_u
        ss_cum_temp_sp = sum(ss_cum_density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
        ss_cum_temp_sp = ss_cum_temp_sp*P_u
        ss_cum_temp_emp = sum(ss_cum_density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
        ss_cum_temp_emp = ss_cum_temp_emp*P_u

        for zeta_prime_i in 1:number_zeta_nodes
            ss_cum_density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_w
            ss_cum_density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_sp
            ss_cum_density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_emp

        end
    end

    Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        ss_cum_density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_w
        ss_cum_density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_sp
        ss_cum_density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_emp
    end

    ss_cum_occ_trans_matrix[2, 1,1] = sum(ss_w_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
    ss_cum_occ_trans_matrix[2, 1,2] = sum(ss_sp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
    ss_cum_occ_trans_matrix[2, 1,3] = sum(ss_emp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)

    ss_cum_occ_trans_matrix[2, 2,1] = sum(ss_w_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
    ss_cum_occ_trans_matrix[2, 2,2] = sum(ss_sp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
    ss_cum_occ_trans_matrix[2, 2,3] = sum(ss_emp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)

    ss_cum_occ_trans_matrix[2, 3,1] = sum(ss_w_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
    ss_cum_occ_trans_matrix[2, 3,2] = sum(ss_sp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
    ss_cum_occ_trans_matrix[2, 3,3] = sum(ss_emp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)

    ss_cum_density_distr_w = copy(ss_cum_density_distr_w_pr)
    ss_cum_density_distr_sp = copy(ss_cum_density_distr_sp_pr)
    ss_cum_density_distr_emp = copy(ss_cum_density_distr_emp_pr)

    Threads.@threads for (m,o) = collect(Iterators.product(1:6,1:3))
        dd = ss_cum_density_distr_w
        if o == 2
            dd = ss_cum_density_distr_sp
        elseif o == 3
            dd = ss_cum_density_distr_emp
        end
        measure = ss_consumption
        if m==2
            measure = ss_earnings
        elseif m==3
            measure = ss_income
        elseif m==4
            measure = ss_wealth
        elseif m==5
            measure = ss_capital_income
        elseif m==6
            measure = ss_savings
        end
        ss_cum_household_measures[2,m,o,:] .= compute_inequality_measures(measure, dd)
        for o2 = 1:3
            occ2 = ss_w_choice
            if o2 == 2
                occ2 = ss_sp_choice
            elseif o2 == 3
                occ2 = ss_emp_choice
            end
            ss_cum_trans_household_measures[2,m,o,o2,:] .= compute_inequality_measures(measure.*occ2,dd.*occ2)
        end
    end

    cum_occ_trans_matrix[2,:,:] = copy(ss_cum_occ_trans_matrix[2,:,:])

    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        temp_d_w = Schumaker(ss_1[1][3], cumsum(ss_cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Constant,Constant))
        cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = diff([0.0;evaluate.(temp_d_w, asset_grid)])
        cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i])

        temp_d_sp = Schumaker(ss_1[1][3], cumsum(ss_cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Constant,Constant))
        cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = diff([0.0;evaluate.(temp_d_sp, asset_grid)])
        cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])

        temp_d_emp = Schumaker(ss_1[1][3], cumsum(ss_cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Constant,Constant))
        cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = diff([0.0;evaluate.(temp_d_emp, asset_grid)])
        cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])

    end

    cum_consumption, cum_earnings, cum_income, cum_wealth, cum_capital_income, cum_savings, cum_occ_choice,cum_w_choice,cum_sp_choice,cum_emp_choice = compute_household_measures(policy_s[1,:,:,:,:,:], r_s[1], w_s[1], lambda_s[1], asset_grid)
    Threads.@threads for (m,o) = collect(Iterators.product(1:6,1:3))
        dd = cum_density_distr_w
        if o == 2
            dd = cum_density_distr_sp
        elseif o == 3
            dd = cum_density_distr_emp
        end
        measure = cum_consumption
        if m==2
            measure = cum_earnings
        elseif m==3
            measure = cum_income
        elseif m==4
            measure = cum_wealth
        elseif m==5
            measure = cum_capital_income
        elseif m==6
            measure = cum_savings
        end
        cum_household_measures[2,m,o,:] .= compute_inequality_measures(measure,dd)
        for o2 = 1:3
            occ2 = cum_w_choice
            if o2 == 2
                occ2 = cum_sp_choice
            elseif o2 == 3
                occ2 = cum_emp_choice
            end
            cum_trans_household_measures[2,m,o,o2,:] .= compute_inequality_measures(measure.*occ2,dd.*occ2)
        end
    end

    occ_trans_matrix[2,:,:] = copy(ss_1[1][18])

    # t = 2:T
    p = Progress((T-1)*((number_asset_grid*number_u_nodes+2)*(number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes)), dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=100)
    for t=1:T-1 # & i = 3:T+1
        i = t+2

        cum_density_distr_w_pr = copy(cum_density_distr_w).*0.0
        cum_density_distr_sp_pr = copy(cum_density_distr_sp).*0.0
        cum_density_distr_emp_pr = copy(cum_density_distr_emp).*0.0

        ss_cum_density_distr_w_pr = copy(ss_cum_density_distr_w).*0.0
        ss_cum_density_distr_sp_pr = copy(ss_cum_density_distr_sp).*0.0
        ss_cum_density_distr_emp_pr = copy(ss_cum_density_distr_emp).*0.0

        density_distr_w_pr = copy(density_distr_w).*0.0
        density_distr_sp_pr = copy(density_distr_sp).*0.0
        density_distr_emp_pr = copy(density_distr_emp).*0.0

        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            for a_i in 1:number_asset_grid
                j_1 = a1_indices[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                cum_density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                cum_density_distr_sp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                cum_density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                density_distr_sp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                if j_1 != number_asset_grid
                    cum_density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    cum_density_distr_sp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    cum_density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                    density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_sp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                end

                ss_j_1 = ss_a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                ss_cum_density_distr_w_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_sp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_emp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                if ss_j_1 != number_asset_grid
                    ss_cum_density_distr_w_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    ss_cum_density_distr_sp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    ss_cum_density_distr_emp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                end
                next!(p)
            end
        end

        cum_density_distr_marginal_assets_w = sum(sum(sum(sum(cum_density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        cum_density_distr_marginal_assets_sp = sum(sum(sum(sum(cum_density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        cum_density_distr_marginal_assets_emp = sum(sum(sum(sum(cum_density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        ss_cum_density_distr_marginal_assets_w = sum(sum(sum(sum(ss_cum_density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        ss_cum_density_distr_marginal_assets_sp = sum(sum(sum(sum(ss_cum_density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        ss_cum_density_distr_marginal_assets_emp = sum(sum(sum(sum(ss_cum_density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        density_distr_marginal_assets_w = sum(sum(sum(sum(density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        density_distr_marginal_assets_sp = sum(sum(sum(sum(density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        density_distr_marginal_assets_emp = sum(sum(sum(sum(density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
            cum_temp_w = sum(cum_density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            cum_temp_w = cum_temp_w*P_u
            cum_temp_sp = sum(cum_density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            cum_temp_sp = cum_temp_sp*P_u
            cum_temp_emp = sum(cum_density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            cum_temp_emp = cum_temp_emp*P_u

            ss_cum_temp_w = sum(ss_cum_density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            ss_cum_temp_w = ss_cum_temp_w*P_u
            ss_cum_temp_sp = sum(ss_cum_density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            ss_cum_temp_sp = ss_cum_temp_sp*P_u
            ss_cum_temp_emp = sum(ss_cum_density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            ss_cum_temp_emp = ss_cum_temp_emp*P_u

            temp_w = sum(density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_w = temp_w*P_u
            temp_sp = sum(density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_sp = temp_sp*P_u
            temp_emp = sum(density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_emp = temp_emp*P_u

            for zeta_prime_i in 1:number_zeta_nodes
                cum_density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*cum_temp_w
                cum_density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*cum_temp_sp
                cum_density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*cum_temp_emp

                ss_cum_density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_w
                ss_cum_density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_sp
                ss_cum_density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_emp

                density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_w
                density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_sp
                density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_emp

                next!(p)
            end
        end

        Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            cum_density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*cum_density_distr_marginal_assets_w
            cum_density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*cum_density_distr_marginal_assets_sp
            cum_density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*cum_density_distr_marginal_assets_emp

            ss_cum_density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_w
            ss_cum_density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_sp
            ss_cum_density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_emp

            density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_w
            density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_sp
            density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_emp

            next!(p)
        end

        # compute inequality measures
        cum_consumption, cum_earnings, cum_income, cum_wealth, cum_capital_income, cum_savings, occ_choice,cum_w_choice,cum_sp_choice,cum_emp_choice = compute_household_measures(policy_s[t+1,:,:,:,:,:], r_s[t+1], w_s[t+1], lambda_s[t+1], asset_grid)
        Threads.@threads for (m,o) = collect(Iterators.product(1:6,1:3))
            dd = cum_density_distr_w_pr
            if o == 2
                dd = cum_density_distr_sp_pr
            elseif o == 3
                dd = cum_density_distr_emp_pr
            end
            measure = cum_consumption
            if m==2
                measure = cum_earnings
            elseif m==3
                measure = cum_income
            elseif m==4
                measure = cum_wealth
            elseif m==5
                measure = cum_capital_income
            elseif m==6
                measure = cum_savings
            end
            cum_household_measures[i,m,o,:] .= compute_inequality_measures(measure, dd)
            for o2 = 1:3
                occ2 = cum_w_choice
                if o2 == 2
                    occ2 = cum_sp_choice
                elseif o2 == 3
                    occ2 = cum_emp_choice
                end
                cum_trans_household_measures[i,m,o,o2,:] .= compute_inequality_measures(measure.*occ2,dd.*occ2)
            end
        end

        w_choice = Float64.(occ_choice.==1.0)
        sp_choice = Float64.(occ_choice.==2.0)
        emp_choice = Float64.(occ_choice.==3.0)

        # conditional cumulative probability of households - P( ot = OCCj | o1 = OCCi )
        # W->W
        cum_occ_trans_matrix[i, 1,1] = sum(w_choice.*cum_density_distr_w_pr)/sum(init_density_distr_w)
        # W->SP
        cum_occ_trans_matrix[i, 1,2] = sum(sp_choice.*cum_density_distr_w_pr)/sum(init_density_distr_w)
        # W->EMP
        cum_occ_trans_matrix[i, 1,3] = sum(emp_choice.*cum_density_distr_w_pr)/sum(init_density_distr_w)

        # SP->W
        cum_occ_trans_matrix[i, 2,1] = sum(w_choice.*cum_density_distr_sp_pr)/sum(init_density_distr_sp)
        # W->SP
        cum_occ_trans_matrix[i, 2,2] = sum(sp_choice.*cum_density_distr_sp_pr)/sum(init_density_distr_sp)
        # W->EMP
        cum_occ_trans_matrix[i, 2,3] = sum(emp_choice.*cum_density_distr_sp_pr)/sum(init_density_distr_sp)

        # EMP->W
        cum_occ_trans_matrix[i, 3,1] = sum(w_choice.*cum_density_distr_emp_pr)/sum(init_density_distr_emp)
        # EMP->SP
        cum_occ_trans_matrix[i, 3,2] = sum(sp_choice.*cum_density_distr_emp_pr)/sum(init_density_distr_emp)
        # EMP->EMP
        cum_occ_trans_matrix[i, 3,3] = sum(emp_choice.*cum_density_distr_emp_pr)/sum(init_density_distr_emp)

        # conditional cumulative occupational choice probability for steady state
        ss_cum_occ_trans_matrix[i, 1,1] = sum(ss_w_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
        ss_cum_occ_trans_matrix[i, 1,2] = sum(ss_sp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
        ss_cum_occ_trans_matrix[i, 1,3] = sum(ss_emp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)

        ss_cum_occ_trans_matrix[i, 2,1] = sum(ss_w_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
        ss_cum_occ_trans_matrix[i, 2,2] = sum(ss_sp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
        ss_cum_occ_trans_matrix[i, 2,3] = sum(ss_emp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)

        ss_cum_occ_trans_matrix[i, 3,1] = sum(ss_w_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
        ss_cum_occ_trans_matrix[i, 3,2] = sum(ss_sp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
        ss_cum_occ_trans_matrix[i, 3,3] = sum(ss_emp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)

        # inequality measures for steady state
        Threads.@threads for (m,o) = collect(Iterators.product(1:6,1:3))
            dd = ss_cum_density_distr_w_pr
            if o == 2
                dd = ss_cum_density_distr_sp_pr
            elseif o == 3
                dd = ss_cum_density_distr_emp_pr
            end
            measure = ss_consumption
            if m==2
                measure = ss_earnings
            elseif m==3
                measure = ss_income
            elseif m==4
                measure = ss_wealth
            elseif m==5
                measure = ss_capital_income
            elseif m==6
                measure = ss_savings
            end
            ss_cum_household_measures[i,m,o,:] .= compute_inequality_measures(measure, dd)
            for o2 = 1:3
                occ2 = ss_w_choice
                if o2 == 2
                    occ2 = ss_sp_choice
                elseif o2 == 3
                    occ2 = ss_emp_choice
                end
                ss_cum_trans_household_measures[i,m,o,o2,:] .= compute_inequality_measures(measure.*occ2,dd.*occ2)
            end
        end

        # conditional occupational choice probability for steady state
        occ_trans_matrix[i, 1,1] = sum(w_choice.*density_distr_w_pr)/sum(density_distr_w)
        occ_trans_matrix[i, 1,2] = sum(sp_choice.*density_distr_w_pr)/sum(density_distr_w)
        occ_trans_matrix[i, 1,3] = sum(emp_choice.*density_distr_w_pr)/sum(density_distr_w)

        occ_trans_matrix[i, 2,1] = sum(w_choice.*density_distr_sp_pr)/sum(density_distr_sp)
        occ_trans_matrix[i, 2,2] = sum(sp_choice.*density_distr_sp_pr)/sum(density_distr_sp)
        occ_trans_matrix[i, 2,3] = sum(emp_choice.*density_distr_sp_pr)/sum(density_distr_sp)

        occ_trans_matrix[i, 3,1] = sum(w_choice.*density_distr_emp_pr)/sum(density_distr_emp)
        occ_trans_matrix[i, 3,2] = sum(sp_choice.*density_distr_emp_pr)/sum(density_distr_emp)
        occ_trans_matrix[i, 3,3] = sum(emp_choice.*density_distr_emp_pr)/sum(density_distr_emp)


        # final in-loop step
        cum_density_distr_w = copy(cum_density_distr_w_pr)
        cum_density_distr_sp = copy(cum_density_distr_sp_pr)
        cum_density_distr_emp = copy(cum_density_distr_emp_pr)

        ss_cum_density_distr_w = copy(ss_cum_density_distr_w_pr)
        ss_cum_density_distr_sp = copy(ss_cum_density_distr_sp_pr)
        ss_cum_density_distr_emp = copy(ss_cum_density_distr_emp_pr)

        density_distr_w = w_choice.*(density_distr_w_pr.+density_distr_sp_pr.+density_distr_emp_pr)
        density_distr_sp = sp_choice.*(density_distr_w_pr.+density_distr_sp_pr.+density_distr_emp_pr)
        density_distr_emp = emp_choice.*(density_distr_w_pr.+density_distr_sp_pr.+density_distr_emp_pr)

        #temp_res = [sum(density_distr_w), sum(density_distr_sp), sum(density_distr_emp)]
        #println_sameline(temp_res)
    end

    println_sameline("Final cycle - Done")

    return cum_occ_trans_matrix,ss_cum_occ_trans_matrix, occ_trans_matrix, cum_household_measures,ss_cum_household_measures, cum_trans_household_measures,ss_cum_trans_household_measures
end
cum_occ_trans_matrix,ss_cum_occ_trans_matrix, occ_trans_matrix, cum_household_measures,ss_cum_household_measures, cum_trans_household_measures,ss_cum_trans_household_measures = f()

# graphing part
LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/Transition/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\Transition\\"
end
mkpath(LOCAL_DIR_OCCUPATION)

occ_text = ["W","SP","EMP"]

plts_cum_otm = Array{Any}(undef,3,3)
plts_cum_ssotm = Array{Any}(undef,3,3)
cum_diffotm = Array{Any}(undef,3,3)
plts_cum_diffotm = Array{Any}(undef,3,3)
plts_cum_wssotm = Array{Any}(undef,3,3)

plts_otm = Array{Any}(undef,3,3)
plts_otm_zoom = Array{Any}(undef,3,3)

for (i,j) = collect(Iterators.product(1:3,1:3))
    global T = length(cum_occ_trans_matrix[:,i,j])-1

    OCC_ = "W"
    if i==2
        OCC_="SP"
    elseif i==3
        OCC_="EMP"
    end

    cum_diffotm[i,j] = cum_occ_trans_matrix[2:T+1,i,j].-ss_cum_occ_trans_matrix[2:T+1,i,j]
    plts_cum_diffotm[i,j] = create_plot(1:T, "Time (Years)",cum_diffotm[i,j], ""#="$(occ_text[i])->$(occ_text[j])"=#, true,OCC_)
    plts_cum_wssotm[i,j] = plot(1:T, [cum_occ_trans_matrix[2:T+1,i,j] ss_cum_occ_trans_matrix[2:T+1,i,j]], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

    #plts_otm[i,j] = plot(0:T, occ_trans_matrix[1:T+1,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

    plts_otm[i,j] = create_plot(1:T, "Time (Years)",occ_trans_matrix[2:T+1,i,j],""#="$(occ_text[i])->$(occ_text[j])"=#,occ_trans_matrix[1,i,j],occ_trans_matrix[T+1,i,j], true,OCC_)
    plts_otm_zoom[i,j] = create_plot(1:10, "Time (Years)",occ_trans_matrix[2:10+1#=T+1=#,i,j],""#="$(occ_text[i])->$(occ_text[j])"=#,occ_trans_matrix[1,i,j],occ_trans_matrix[T+1,i,j], true,OCC_)

end

# it shows the difference between cumulative occupational mobility after reform and pre-reform
# plt_cum_diffotm = plot(plts_cum_diffotm[1,1],plts_cum_diffotm[1,2],plts_cum_diffotm[1,3], plts_cum_diffotm[2,1],plts_cum_diffotm[2,2],plts_cum_diffotm[2,3], plts_cum_diffotm[3,1],plts_cum_diffotm[3,2],plts_cum_diffotm[3,3], layout=(3,3))
# display(plt_cum_diffotm)
# savefig(plt_cum_diffotm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_diff_with_ss.png")

# TT = length(cum_occ_trans_matrix[:,1,1])-1
# ii=1
# LEGENDPOS = :bottomright
# pltW = create_combined_plot(collect(1:TT), "Time (Years)", [cum_diffotm[ii,1][2:end],cum_diffotm[ii,2][2:end],cum_diffotm[ii,3][2:end]], ["W","SP","EMP"],"Occupational shares for Workers from t=0",[cum_diffotm[ii,1][1],cum_diffotm[ii,2][1],cum_diffotm[ii,3][1]],[cum_diffotm[ii,1][end],cum_diffotm[ii,2][end],cum_diffotm[ii,3][end]], true,["W","SP","EMP"],LEGENDPOS)
# display(pltW)
# savefig(pltW,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_W_diff.png")
#
# ii=2
# LEGENDPOS = false
# pltSP = create_combined_plot(collect(1:TT), "Time (Years)", [cum_diffotm[ii,1][2:end],cum_diffotm[ii,2][2:end],cum_diffotm[ii,3][2:end]], ["W","SP","EMP"],"Occupational shares for Sole Prop. from t=0",[cum_diffotm[ii,1][1],cum_diffotm[ii,2][1],cum_diffotm[ii,3][1]],[cum_diffotm[ii,1][end],cum_diffotm[ii,2][end],cum_diffotm[ii,3][end]], true,["W","SP","EMP"],LEGENDPOS)
# display(pltSP)
# savefig(pltSP,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_SP_diff.png")
#
# ii=3
# LEGENDPOS = false
# pltEMP = create_combined_plot(collect(1:TT), "Time (Years)", [cum_diffotm[ii,1][2:end],cum_diffotm[ii,2][2:end],cum_diffotm[ii,3][2:end]], ["W","SP","EMP"],"Occupational shares for Employers from t=0",[cum_diffotm[ii,1][1],cum_diffotm[ii,2][1],cum_diffotm[ii,3][1]],[cum_diffotm[ii,1][end],cum_diffotm[ii,2][end],cum_diffotm[ii,3][end]], true,["W","SP","EMP"],LEGENDPOS)
# display(pltEMP)
# savefig(pltEMP,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_EMP_diff.png")

#plt = plot(pltW,pltSP,pltEMP, layout=(1,3))
#display(plt)
#savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_diff_with_ss.png")

# it shows cumulative occupational mobilities after reform and pre-reform
# plt_cum_wssotm = plot(plts_cum_wssotm[1,1],plts_cum_wssotm[1,2],plts_cum_wssotm[1,3], plts_cum_wssotm[2,1],plts_cum_wssotm[2,2],plts_cum_wssotm[2,3], plts_cum_wssotm[3,1],plts_cum_wssotm[3,2],plts_cum_wssotm[3,3], layout=(3,3))
# display(plt_cum_wssotm)
# savefig(plt_cum_wssotm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_with_ss.png")


# it shows the transitional share of occupation1 at time 1 who become occupation2 at time t
# plt_otm = plot(plts_otm[1,1],plts_otm[1,2],plts_otm[1,3], plts_otm[2,1],plts_otm[2,2],plts_otm[2,3], plts_otm[3,1],plts_otm[3,2],plts_otm[3,3], layout=(3,3))
# display(plt_otm)
# savefig(plt_otm,"$(LOCAL_DIR_OCCUPATION)time_occupational_mobility.png")
for (i,j) = collect(Iterators.product(1:3,1:3))
    plt = plts_otm_zoom[i,j]
    display(plt)
    savefig(plt,"$(LOCAL_DIR_OCCUPATION)time_occupational_mobility_zoomed_$(occ_text[i])_$(occ_text[j]).png")
end
# plt_otm_zoom = plot(plts_otm_zoom[1,1],plts_otm_zoom[1,2],plts_otm_zoom[1,3], plts_otm_zoom[2,1],plts_otm_zoom[2,2],plts_otm_zoom[2,3], plts_otm_zoom[3,1],plts_otm_zoom[3,2],plts_otm_zoom[3,3], layout=(3,3))
# display(plt_otm_zoom)
# savefig(plt_otm_zoom,"$(LOCAL_DIR_OCCUPATION)time_occupational_mobility_zoomed.png")

# it shows the inequality measures for fixed occupation in time t=1

LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/Transition/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\Transition\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#cum_household_measures
inequality_measures = ["mean"#=,"variance","gini"=#]
measures = ["consumption","earnings","income","wealth","capital income","savings"]
for m = 4#1:6
    for im = 1#:3
        # plts_cum_wsshm = Array{Any}(undef,3)
        # for o = 1:3
        #     T = length(cum_household_measures[:,m,o,im])-1
        #     plts_cum_wsshm[o] = plot(0:T, [cum_household_measures[1:T+1,m,o,im] ss_cum_household_measures[1:T+1,m,o,im]], legend=false,ylabel="$(occ_text[o])'s $(inequality_measures[im]) $(measures[m])")
        #
        # end
        # plt_cum_wsshm = plot(plts_cum_wsshm[1],plts_cum_wsshm[2],plts_cum_wsshm[3], layout=(1,3))
        # display(plt_cum_wsshm)
        #savefig(plt_cum_wsshm,"$(LOCAL_DIR_INEQUALITY)time_$(inequality_measures[im])_$(measures[m])_with_ss.png")

        plts_cum_wssmthm = Array{Any}(undef,3,3)
        for (i,j) = collect(Iterators.product(1:3,1:3))
            T = length(cum_trans_household_measures[:,m,i,j,im])-1
            plts_cum_wssmthm[i,j] = plot(1:T, cum_trans_household_measures[2:T+1,m,i,j,im].-ss_cum_trans_household_measures[2:T+1,m,i,j,im], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

        end

        plt_cum_wssmthm = plot(plts_cum_wssmthm[1,1],plts_cum_wssmthm[1,2],plts_cum_wssmthm[1,3], plts_cum_wssmthm[2,1],plts_cum_wssmthm[2,2],plts_cum_wssmthm[2,3], plts_cum_wssmthm[3,1],plts_cum_wssmthm[3,2],plts_cum_wssmthm[3,3], layout=(3,3))
        display(plt_cum_wssmthm)
        #savefig(plt_cum_wssmthm,"$(LOCAL_DIR_INEQUALITY)time_cumulative_$(inequality_measures[im])_$(measures[m])_with_ss.png")
    end
end


LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/Transition/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\Transition\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
occ_text = ["W","SP","EMP"]
inequality_measures = ["Mean"#=,"variance","gini"=#]
measures = ["Consumption","Earnings","Income","Wealth","Capital_Income","Savings"]
for m = 1:6
    for im = 1#:3
        #=plts_cum_diff_in_diff_hm = Array{Any}(undef,3)
        for o = 1:3
            T = length(cum_household_measures[:,m,o,im])-1

            plts_cum_diff_in_diff_hm[o] = plot(0:T, (cum_household_measures[1:T+1,m,o,im].-ss_cum_household_measures[1:T+1,m,o,im])./cum_household_measures[1,m,o,im], legend=false,ylabel="$(occ_text[o])'s $(inequality_measures[im]) $(measures[m]) (diff-in-diff%)")

        end
        plt_cum_diff_in_diff_hm = plot(plts_cum_diff_in_diff_hm[1],plts_cum_diff_in_diff_hm[2],plts_cum_diff_in_diff_hm[3], layout=(1,3))
        display(plt_cum_diff_in_diff_hm)
        savefig(plt_cum_diff_in_diff_hm,"$(LOCAL_DIR_INEQUALITY)time_$(inequality_measures[im])_$(measures[m])_diff_in_diff.png")
        =#
        T = 0
        cum_diff_percent_hm = Array{Any}(undef,3)
        cum_diff_hm = Array{Any}(undef,3)
        for o = 1:3
            T = length(cum_household_measures[:,m,o,im])-1
            cum_diff_hm[o] = (cum_household_measures[1:T+1,m,o,im].-ss_cum_household_measures[1:T+1,m,o,im])
            cum_diff_percent_hm[o] = cum_diff_hm[o]./ss_cum_household_measures[1:T+1,m,o,im]
            if m==6
                cum_diff_percent_hm[o] = cum_diff_hm[o]./ss_cum_household_measures[1:T+1,4,o,im]
            end
            # cum_diff_percent_hm[o] = cum_diff_hm[o]./((cum_household_measures[1:T+1,m,o,im].+ss_cum_household_measures[1:T+1,m,o,im])/2.0)
        end
        LEGENDPOS = false
        if m==6
            LEGENDPOS = :topright
        end
        plt = create_combined_plot(collect(1:T),"Time (Years)", [cum_diff_percent_hm[1][2:end],cum_diff_percent_hm[2][2:end],cum_diff_percent_hm[3][2:end]], ["W","SP","EMP"], ""#="$(inequality_measures[im]) $(measures[m]) (diff%)"=#, [cum_diff_percent_hm[1][1],cum_diff_percent_hm[2][1],cum_diff_percent_hm[3][1]],[cum_diff_percent_hm[1][end],cum_diff_percent_hm[2][end],cum_diff_percent_hm[3][end]],true, ["W","SP","EMP"],LEGENDPOS)
        display(plt)
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_$(inequality_measures[im])_$(measures[m])_diff_percent.png")

        # LEGENDPOS = false
        # if m==6
        #     LEGENDPOS = :bottomright
        # end
        # plt = create_combined_plot(collect(1:T),"Time (Years)", [cum_diff_hm[1][2:end],cum_diff_hm[2][2:end],cum_diff_hm[3][2:end]], ["W","SP","EMP"], "$(inequality_measures[im]) $(measures[m]) (diff)", [cum_diff_hm[1][1],cum_diff_hm[2][1],cum_diff_hm[3][1]],[cum_diff_hm[1][end],cum_diff_hm[2][end],cum_diff_hm[3][end]],false, ["W","SP","EMP"],LEGENDPOS)
        # display(plt)
        # savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_$(inequality_measures[im])_$(measures[m])_diff.png")
    end
end
