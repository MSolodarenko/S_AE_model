##############
# Script to track households throuout their transition paths
##############

include("Functions/print_sameline.jl")
include("Functions/profit.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots

using JLD2
using ProgressMeter
using SchumakerSpline

TIME_PERIODS = 50

function create_combined_plot(X,XLABEL::String,Ys,YLABELs,YLABEL,Y1s,Y2s, IS_Y_PERCENTAGE::Bool=false, OCCUPATION=["SP","EMP"], LEGENDPOS=false)
    TICKSFONTSIZE = 18

    NUM_XTICKS = 6
    XTICKS = Int64.(round.(collect(range(0;stop=length(X),length=NUM_XTICKS))))

    XTICKS = [1; XTICKS[2:end]]

    NUMYTICKS = 7#11
    YLIMS1 = minimum([minimum.(Ys); minimum.(Y1s); minimum.(Y2s)])
    YLIMS2 = maximum([maximum.(Ys); maximum.(Y1s); maximum.(Y2s)])
    YTICKSlow = []
    YTICKShigh = []

    if YLIMS1 >= 0.0
        YTICKSlow = []
        YTICKShigh = collect(range(0.0; stop=YLIMS2, length=NUMYTICKS))[2:end]
    elseif YLIMS2 <= 0.0
        YTICKSlow = collect(range(YLIMS1; stop=0.0, length=NUMYTICKS))[1:end-1]
        YTICKShigh = []
    elseif abs(YLIMS1) < (abs(YLIMS1)+abs(YLIMS2))/(NUMYTICKS-1)
        YTICKSlow = [0.0-(abs(YLIMS1)+abs(YLIMS2))/(NUMYTICKS-1)]
        YTICKShigh = collect(range(0.0; stop=YLIMS2, length=NUMYTICKS-1))[2:end]
    elseif abs(YLIMS2) < (abs(YLIMS1)+abs(YLIMS2))/(NUMYTICKS-1)
        YTICKSlow = collect(range(YLIMS1; stop=0.0, length=NUMYTICKS-1))[1:end-1]
        YTICKShigh = [0.0+(abs(YLIMS1)+abs(YLIMS2))/(NUMYTICKS-1)]
    # elseif abs(abs(YLIMS1) - abs(YLIMS2)) < (abs(YLIMS1)+abs(YLIMS2))/(NUMYTICKS-1)
    #     YLIMS2 = max(abs(YLIMS1),abs(YLIMS2))
    #     YLIMS1 = -YLIMS2
    #     stepYTICKS = (abs(YLIMS1)+abs(YLIMS2))/(NUMYTICKS-1)
    #     YTICKSlow = reverse(collect(range(0.0; stop=YLIMS1, step=-stepYTICKS))[2:end])
    #     YTICKShigh = collect(range(0.0; stop=YLIMS2, step=stepYTICKS))[2:end]
    elseif abs(YLIMS1) < abs(YLIMS2)
        numYTICKSlow = Int64(round((NUMYTICKS-1)*abs(YLIMS1)/(abs(YLIMS1)+abs(YLIMS2))))+1
        YTICKShigh = collect(range(0.0; stop=YLIMS2, length=NUMYTICKS-numYTICKSlow))[2:end]
        stepYTICKS = abs(YTICKShigh[1]-YTICKShigh[2])
        YTICKSlow = reverse(collect(range(0.0; step=-stepYTICKS, length=numYTICKSlow))[2:end])
    else#if abs(YLIMS1) >= abs(YLIMS2)
        numYTICKShigh = Int64(round((NUMYTICKS-1)*abs(YLIMS2)/(abs(YLIMS1)+abs(YLIMS2))))+1
        YTICKSlow = collect(range(YLIMS1; stop=0.0, length=NUMYTICKS-numYTICKShigh))[1:end-1]
        stepYTICKS = abs(YTICKSlow[1]-YTICKSlow[2])
        YTICKShigh = collect(range(0.0; step=stepYTICKS, length=numYTICKShigh))[2:end]
    end
    YTICKS = [YTICKSlow; 0.0; YTICKShigh]

    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015

    if maximum([maximum.(Ys); maximum.(Y1s); maximum.(Y2s)]) <= YTICKS[end-1] + YLIMMARGIN
        YTICKS = YTICKS[1:end-1]
        YLIMS2 = YTICKS[end]
    elseif minimum([minimum.(Ys); minimum.(Y1s); minimum.(Y2s)]) >= YTICKS[2] - YLIMMARGIN
        YTICKS = YTICKS[2:end]
        YLIMS1 = YTICKS[1]
    end

    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)

    YLIMS1 = YTICKS[1]
    YLIMS2 = YTICKS[end]
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
            # push!(COLORS,"yellow")
            push!(COLORS,RGB(255/255, 221/255, 0/255))
        else
            # push!(COLORS,"blue")
            push!(COLORS,RGB(0/255, 87/255, 183/255))
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
                        legendfontsize=ceil(Int64,TICKSFONTSIZE*0.72),
                        xlabel=XLABEL,
                        label=[YLABELs[y_i] ""],
                        linewidth=2.5, thickness_scaling = 1,
                        xtickfontsize=TICKSFONTSIZE,
                        ytickfontsize=TICKSFONTSIZE,
                        xguidefontsize=TICKSFONTSIZE,
                        yguidefontsize=TICKSFONTSIZE,
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

println_sameline("import data from the model with occupational mobility")
LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE\\"
end
@load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2
trans_res_w_occ = copy(trans_res)
ss_1_w_occ = copy(ss_1)

println_sameline("import data from the model without occupational mobility")
LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE_fixed_occ/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE_fixed_occ\\"
end
@load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2
trans_res_wo_occ = copy(trans_res)
ss_1_wo_occ = copy(ss_1)

println_sameline("extract parameters and variables connected to skills")
#These must be the same: ss_1_w_occ[4],ss_1_wo_occ[4]
approx_object = ss_1_w_occ[4]
global_approx_params = ss_1_w_occ[1][47]
model_params = ss_1_w_occ[5]

number_u_nodes = approx_object[7]
number_zeta_nodes = global_approx_params[4]
number_alpha_m_nodes = global_approx_params[5]
number_alpha_w_nodes = global_approx_params[6]

P_u = approx_object[4]
p_alpha = model_params[13]
P_zeta = approx_object[3]

stat_P_u = approx_object[5]
P_alpha = approx_object[6]

delta = model_params[3]
gamma = model_params[4]
eta = model_params[5]
theta = model_params[6]
c_e = model_params[7]
z_m_nodes = approx_object[1]
z_w_nodes = approx_object[2]

println_sameline("extract lambdas")
lambda_s_no_reform = ones(TIME_PERIODS).*ss_1_w_occ[5][1]
lambda_s_reform = copy(trans_res_w_occ[2])

println_sameline("extract factor prices paths")
r_s_no_reform_w_occ = ones(TIME_PERIODS).*ss_1_w_occ[2]
w_s_no_reform_w_occ = ones(TIME_PERIODS).*ss_1_w_occ[3]

r_s_reform_w_occ = copy(trans_res_w_occ[3])
w_s_reform_w_occ = copy(trans_res_w_occ[4])

r_s_no_reform_wo_occ = ones(TIME_PERIODS).*ss_1_wo_occ[2]
w_s_no_reform_wo_occ = ones(TIME_PERIODS).*ss_1_wo_occ[3]

r_s_reform_wo_occ = copy(trans_res_wo_occ[3])
w_s_reform_wo_occ = copy(trans_res_wo_occ[4])

println_sameline("create common asset_grid #ss_1[1][3]")
number_asset_grid = ss_1_w_occ[1][2]*5
a_max = max(ss_1_w_occ[1][1],ss_1_wo_occ[1][1])
asset_grid = exp.(collect(range(log(0.01+1); stop=log(a_max+1), length=number_asset_grid))).-1

println_sameline("intrapolate policy functions and initial distribution on a new asset_grid")
#   set the structure of policy functions and distribution - [1:TIME_PERIODS][occ-1:3][1:number_asset_grid,1:number_u_nodes,1:number_zeta_nodes,1:number_alpha_m_nodes,1:number_alpha_w_nodes]
#   policy_set = [policy, a1_indices, lottery_prob_1, lottery_prob_2]
println_sameline("policy_set_reform_w_occ...")
gen_policy_set_reform_w_occ = [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for t=1:TIME_PERIODS] for f=1:4]
p = Progress(TIME_PERIODS*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes*number_asset_grid, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
Threads.@threads for (t,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:TIME_PERIODS,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
    p_f = Schumaker(trans_res_w_occ[10], trans_res_w_occ[5][8][t,:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
    Threads.@threads for a_i in 1:number_asset_grid
        a1 = evaluate(p_f, asset_grid[a_i])
        gen_policy_set_reform_w_occ[1][t][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
        j_1 = sum(a1 .>= asset_grid)
        j = j_1 + 1

        gen_policy_set_reform_w_occ[2][t][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = max(1,min(Int64(round(j_1;digits=0)),number_asset_grid))
        if j <= number_asset_grid && j_1 >= 1
            gen_policy_set_reform_w_occ[3][t][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
            gen_policy_set_reform_w_occ[4][t][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
        elseif j_1 == number_asset_grid
            gen_policy_set_reform_w_occ[3][t][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
            gen_policy_set_reform_w_occ[4][t][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
        end
        next!(p)
    end
end
##########
policy_set_reform_w_occ = [ [ [ gen_policy_set_reform_w_occ[f][t] for o=1:3] for t=1:TIME_PERIODS] for f=1:4]
############

println_sameline("policy_set_reform_wo_occ...")
############
policy_set_reform_wo_occ = [ [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3] for t=1:TIME_PERIODS] for f=1:4]
############
p = Progress(TIME_PERIODS*3*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes*number_asset_grid, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
Threads.@threads for (t,(o,(u_i,(zeta_i,(alpha_m_i,alpha_w_i))))) in collect(Iterators.product(1:TIME_PERIODS,Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))))
    p_f = Schumaker(trans_res_wo_occ[10],trans_res_wo_occ[5][8][o][t,:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
    Threads.@threads for a_i in 1:number_asset_grid
        a1 = evaluate(p_f, asset_grid[a_i])
        policy_set_reform_wo_occ[1][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
        j_1 = sum(a1 .>= asset_grid)
        j = j_1 + 1

        policy_set_reform_wo_occ[2][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = max(1,min(Int64(round(j_1;digits=0)),number_asset_grid))
        if j <= number_asset_grid && j_1 >= 1
            policy_set_reform_wo_occ[3][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
            policy_set_reform_wo_occ[4][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
        elseif j_1 == number_asset_grid
            policy_set_reform_wo_occ[3][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
            policy_set_reform_wo_occ[4][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
        end
        next!(p)
    end
end

println_sameline("policy_set_no_reform_w_occ...")
gen_policy_set_no_reform_w_occ = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for f=1:4]
p = Progress(number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes*number_asset_grid, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
    p_f = Schumaker(ss_1_w_occ[1][3], ss_1_w_occ[1][4][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
    Threads.@threads for a_i in 1:number_asset_grid
        a1 = evaluate(p_f, asset_grid[a_i])
        gen_policy_set_no_reform_w_occ[1][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
        j_1 = sum(a1 .>= asset_grid)
        j = j_1 + 1

        gen_policy_set_no_reform_w_occ[2][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = max(1,min(Int64(round(j_1;digits=0)),number_asset_grid))
        if j <= number_asset_grid && j_1 >= 1
            gen_policy_set_no_reform_w_occ[3][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
            gen_policy_set_no_reform_w_occ[4][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
        elseif j_1 == number_asset_grid
            gen_policy_set_no_reform_w_occ[3][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
            gen_policy_set_no_reform_w_occ[4][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
        end
        next!(p)
    end
end
############
policy_set_no_reform_w_occ = [ [ [ gen_policy_set_no_reform_w_occ[f] for o=1:3] for t=1:TIME_PERIODS] for f=1:4]
############

println_sameline("policy_set_no_reform_wo_occ...")
init_policy_set_no_reform_wo_occ = [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3] for f=1:4]
p = Progress(3*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes*number_asset_grid, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
Threads.@threads for (o,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
    p_f = Schumaker(ss_1_wo_occ[1][3], ss_1_wo_occ[1][4][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
    Threads.@threads for a_i in 1:number_asset_grid
        a1 = evaluate(p_f, asset_grid[a_i])
        init_policy_set_no_reform_wo_occ[1][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
        j_1 = sum(a1 .>= asset_grid)
        j = j_1 + 1

        init_policy_set_no_reform_wo_occ[2][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = max(1,min(Int64(round(j_1;digits=0)),number_asset_grid))
        if j <= number_asset_grid && j_1 >= 1
            init_policy_set_no_reform_wo_occ[3][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
            init_policy_set_no_reform_wo_occ[4][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
        elseif j_1 == number_asset_grid
            init_policy_set_no_reform_wo_occ[3][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
            init_policy_set_no_reform_wo_occ[4][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
        end
        next!(p)
    end
end
############
policy_set_no_reform_wo_occ = [ [ [ init_policy_set_no_reform_wo_occ[f][o] for o=1:3] for t=1:TIME_PERIODS] for f=1:4]
############

# MODE_OF_EXPERIMENT = {1 - Take stationary distribution from the model_w_occ as starting distr,
#                       2 - Take stationary distribution from the model_wo_occ as starting distr,
#                       3 - Take an intersection of two stationary distributions}
function conduct_the_experiment(MODE_OF_EXPERIMENT, LOCAL_DIR)

    println_sameline("MODE_OF_EXPERIMENT=$MODE_OF_EXPERIMENT - initialise distribution paths - distr[t][o] t=1:TIME_PERIODS o=1:3 #W,SP,EMP")
    init_distr_w_occ = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    init_distr_wo_occ = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    # # extract initial distribution for the model with occupational mobility #ss[1][5]  - init_distr[o] o=1:3 #W,SP,EMP
    p = Progress(3*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
    Threads.@threads for (o,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
        if MODE_OF_EXPERIMENT==1
            d_f_w_occ = Schumaker(ss_1_w_occ[1][3], cumsum(ss_1_w_occ[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i].*Float64.(ss_1_w_occ[1][22][:,u_i,zeta_i,alpha_m_i,alpha_w_i].==Float64(o)) ); extrapolation=(Linear,Linear))
            init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(d_f_w_occ, asset_grid)])
            init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i]
        elseif MODE_OF_EXPERIMENT==2
            d_f_wo_occ = Schumaker(ss_1_wo_occ[1][3], cumsum(ss_1_wo_occ[1][5][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] ); extrapolation=(Linear,Linear))
            init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(d_f_wo_occ, asset_grid)])
            init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i]
        elseif MODE_OF_EXPERIMENT==3
            d_f_w_occ = Schumaker(ss_1_w_occ[1][3], cumsum(ss_1_w_occ[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i].*Float64.(ss_1_w_occ[1][22][:,u_i,zeta_i,alpha_m_i,alpha_w_i].==Float64(o)) ); extrapolation=(Linear,Linear))
            init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(d_f_w_occ, asset_grid)])
            d_f_wo_occ = Schumaker(ss_1_wo_occ[1][3], cumsum(ss_1_wo_occ[1][5][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] ); extrapolation=(Linear,Linear))
            init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(d_f_wo_occ, asset_grid)])

            intersection = min.(init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i], init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= intersection
            init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= intersection
        elseif MODE_OF_EXPERIMENT==4
            d_f_w_occ = Schumaker(ss_1_w_occ[1][3], cumsum(ss_1_w_occ[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i].*Float64.(ss_1_w_occ[1][22][:,u_i,zeta_i,alpha_m_i,alpha_w_i].==Float64(o)) ); extrapolation=(Linear,Linear))
            init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(d_f_w_occ, asset_grid)])
            d_f_wo_occ = Schumaker(ss_1_wo_occ[1][3], cumsum(ss_1_wo_occ[1][5][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] ); extrapolation=(Linear,Linear))
            init_distr_wo_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(d_f_wo_occ, asset_grid)])
        else
            throw("MODE_OF_EXPERIMENT indicator value is neither 1 nor 2")
        end


        # line below doesn't work because we have combinations of skills such that there are zero households of particular occupation who hold any assets
        # init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_1_w_occ[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i].*Float64.(ss_1_w_occ[1][22][:,u_i,zeta_i,alpha_m_i,alpha_w_i].==Float64(o)) )/sum(init_distr_w_occ[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
        next!(p)
    end

    # display(sum.(init_distr_w_occ))
    # throw(err)

    distr_no_reform_w_occ = [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3] for t=1:TIME_PERIODS]
    distr_reform_w_occ = [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3] for t=1:TIME_PERIODS]
    distr_no_reform_wo_occ = [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3] for t=1:TIME_PERIODS]
    distr_reform_wo_occ = [ [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3] for t=1:TIME_PERIODS]
    for o = 1:3
        distr_no_reform_w_occ[1][o] = copy(init_distr_w_occ[o])
        distr_reform_w_occ[1][o] = copy(init_distr_w_occ[o])
        distr_no_reform_wo_occ[1][o] = copy(init_distr_wo_occ[o])
        distr_reform_wo_occ[1][o] = copy(init_distr_wo_occ[o])
    end

    # inequality measures for occupations fixed in time t=1 after reform
    function compute_household_measures_w_occ(policy, r, w, lambda, a_grid)

        occ_choice, income, earnings, = compute_income_profile(a_grid,number_asset_grid,r, w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

        consumption = income .- policy
        wealth = ones(size(policy)).*a_grid
        capital_income = wealth.*r
        income = income .- wealth#earnings.+capital_income#
        savings = policy .- wealth
        #consumption = income.+wealth .- policy

        return consumption, earnings, income, wealth, capital_income, savings
    end
    function compute_household_measures_wo_occ(policy, r, w, lambda, a_grid)

        occ_choice, income, earnings, = compute_income_profile_fixed_occ(a_grid,number_asset_grid,r, w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

        consumption = Array{Any}(undef, 3)
        wealth = Array{Any}(undef, 3)
        capital_income = Array{Any}(undef, 3)
        savings = Array{Any}(undef, 3)

        for o = 1:3
            consumption[o] = income[o] .- policy[o]
            wealth[o] = ones(size(policy[o])).*a_grid
            capital_income[o] = wealth[o].*r
            income[o] = income[o] .- wealth[o]#earnings.+capital_income#
            savings[o] = policy[o] .- wealth[o]
            #consumption = income.+wealth .- policy
        end

        return consumption, earnings, income, wealth, capital_income, savings
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

    println_sameline("MODE_OF_EXPERIMENT=$MODE_OF_EXPERIMENT - initialise household measures - consumption, earnings, income, wealth, capital_income, savings")
    hm_nr_w_o = [ [ [ zeros(3) for m=1:6] for o=1:3] for t=1:TIME_PERIODS]
    hm_r_w_o = [ [ [ zeros(3) for m=1:6] for o=1:3] for t=1:TIME_PERIODS]
    hm_nr_wo_o = [ [ [ zeros(3) for m=1:6] for o=1:3] for t=1:TIME_PERIODS]
    hm_r_wo_o = [ [ [ zeros(3) for m=1:6] for o=1:3] for t=1:TIME_PERIODS]
    t=1
    p = Progress(2+3*2+3*6*4, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
    measures_nr_wo_o = compute_household_measures_wo_occ(policy_set_no_reform_wo_occ[1][t], r_s_no_reform_wo_occ[t], w_s_no_reform_wo_occ[t], lambda_s_no_reform[t], asset_grid)
    next!(p)
    measures_r_wo_o = compute_household_measures_wo_occ(policy_set_reform_wo_occ[1][t], r_s_reform_wo_occ[t], w_s_reform_wo_occ[t], lambda_s_reform[t], asset_grid)
    next!(p)
    Threads.@threads for o=1:3
        measures_nr_w_o = compute_household_measures_w_occ(policy_set_no_reform_w_occ[1][t][o], r_s_no_reform_w_occ[t], w_s_no_reform_w_occ[t], lambda_s_no_reform[t], asset_grid)
        next!(p)
        measures_r_w_o = compute_household_measures_w_occ(policy_set_reform_w_occ[1][t][o], r_s_reform_w_occ[t], w_s_reform_w_occ[t], lambda_s_reform[t], asset_grid)
        next!(p)
        Threads.@threads for m=1:6
            hm_nr_w_o[t][o][m] = compute_inequality_measures(measures_nr_w_o[m], distr_no_reform_w_occ[t][o])
            next!(p)
            hm_r_w_o[t][o][m] = compute_inequality_measures(measures_r_w_o[m], distr_reform_w_occ[t][o])
            next!(p)
            hm_nr_wo_o[t][o][m] = compute_inequality_measures(measures_nr_wo_o[m][o], distr_no_reform_wo_occ[t][o])
            next!(p)
            hm_r_wo_o[t][o][m] = compute_inequality_measures(measures_r_wo_o[m][o], distr_reform_wo_occ[t][o])
            next!(p)
        end
    end

    println_sameline("MODE_OF_EXPERIMENT=$MODE_OF_EXPERIMENT - intialise temporary arrays")
    d_a1_z0_nr_w_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    d_a1_z0_r_w_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    d_a1_z0_nr_wo_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    d_a1_z0_r_wo_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]

    new_d_nr_w_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    new_d_r_w_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    new_d_nr_wo_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]
    new_d_r_wo_o = [ zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes) for o=1:3]

    println_sameline("MODE_OF_EXPERIMENT=$MODE_OF_EXPERIMENT - start evolution of distributions...")
    for t = 1:TIME_PERIODS-1
        println_sameline("t=$t")

        # calculate distribution
        # over future capital and current level of skills
        p = Progress(3*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes*4, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
        Threads.@threads for (o,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
            d_a1_z0_nr_w_o[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= 0.0
            non_zero_asset_grid_iters = findall(!iszero,distr_no_reform_w_occ[t][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in non_zero_asset_grid_iters
                j_1 = Int32(policy_set_no_reform_w_occ[2][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i])
                d_a1_z0_nr_w_o[o][j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_no_reform_w_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_no_reform_w_occ[3][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                if j_1 != number_asset_grid
                    d_a1_z0_nr_w_o[o][j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_no_reform_w_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_no_reform_w_occ[4][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                end
            end
            next!(p)

            d_a1_z0_r_w_o[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= 0.0
            non_zero_asset_grid_iters = findall(!iszero,distr_reform_w_occ[t][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in non_zero_asset_grid_iters
                j_1 = Int32(policy_set_reform_w_occ[2][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i])
                d_a1_z0_r_w_o[o][j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_reform_w_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_reform_w_occ[3][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                if j_1 != number_asset_grid
                    d_a1_z0_r_w_o[o][j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_reform_w_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_reform_w_occ[4][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                end
            end
            next!(p)

            d_a1_z0_nr_wo_o[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= 0.0
            non_zero_asset_grid_iters = findall(!iszero,distr_no_reform_wo_occ[t][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in non_zero_asset_grid_iters
                j_1 = Int32(policy_set_no_reform_wo_occ[2][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i])
                d_a1_z0_nr_wo_o[o][j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_no_reform_wo_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_no_reform_wo_occ[3][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                if j_1 != number_asset_grid
                    d_a1_z0_nr_wo_o[o][j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_no_reform_wo_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_no_reform_wo_occ[4][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                end
            end
            next!(p)

            d_a1_z0_r_wo_o[o][:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= 0.0
            non_zero_asset_grid_iters = findall(!iszero,distr_reform_wo_occ[t][o][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in non_zero_asset_grid_iters
                j_1 = Int32(policy_set_reform_wo_occ[2][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i])
                d_a1_z0_r_wo_o[o][j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_reform_wo_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_reform_wo_occ[3][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                if j_1 != number_asset_grid
                    d_a1_z0_r_wo_o[o][j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr_reform_wo_occ[t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*policy_set_reform_wo_occ[4][t][o][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                end
            end
            next!(p)
        end

        # first calculate transition
        # for people that do not change alpha_m and alpha_w
        p = Progress(3*number_alpha_m_nodes*number_alpha_w_nodes*(4+number_zeta_nodes*4), dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
        Threads.@threads for (o,(alpha_m_i,alpha_w_i)) in collect(Iterators.product(1:3,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))
            # zeta is the transitory shock,
            # so add over all levels of zeta
            # and then draw new u_prime and new zeta_prime
            # temp_distr_sum_zeta
            temp_d_s_z_nr_w_o = sum(d_a1_z0_nr_w_o[o][:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_d_s_z_nr_w_o = temp_d_s_z_nr_w_o*P_u
            next!(p)

            temp_d_s_z_r_w_o = sum(d_a1_z0_r_w_o[o][:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_d_s_z_r_w_o = temp_d_s_z_r_w_o*P_u
            next!(p)

            temp_d_s_z_nr_wo_o = sum(d_a1_z0_nr_wo_o[o][:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_d_s_z_nr_wo_o = temp_d_s_z_nr_wo_o*P_u
            next!(p)

            temp_d_s_z_r_wo_o = sum(d_a1_z0_r_wo_o[o][:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_d_s_z_r_wo_o = temp_d_s_z_r_wo_o*P_u
            next!(p)

            for zeta_prime_i in 1:number_zeta_nodes
                new_d_nr_w_o[o][:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_d_s_z_nr_w_o
                next!(p)

                new_d_r_w_o[o][:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_d_s_z_r_w_o
                next!(p)

                new_d_nr_wo_o[o][:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_d_s_z_nr_wo_o
                next!(p)

                new_d_r_wo_o[o][:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_d_s_z_r_wo_o
                next!(p)
            end
        end

        # second calculate transition
        # for people that change alpha_m and alpha_w
        #   and therefore draw new alpha_m_prime, alpha_w_prime
        #       and zeta_prime, u_prime
        #           from stationary distributions for this shocks

        # calculate sum of capital of all people who change skills
        d_m_a_nr_w_o = Array{Any}(undef,3)
        d_m_a_r_w_o = Array{Any}(undef,3)
        d_m_a_nr_wo_o = Array{Any}(undef,3)
        d_m_a_r_wo_o = Array{Any}(undef,3)
        p = Progress(3*4, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
        Threads.@threads for o = 1:3
            d_m_a_nr_w_o[o] = sum(sum(sum(sum(d_a1_z0_nr_w_o[o],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            next!(p)
            d_m_a_r_w_o[o] = sum(sum(sum(sum(d_a1_z0_r_w_o[o],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            next!(p)
            d_m_a_nr_wo_o[o] = sum(sum(sum(sum(d_a1_z0_nr_wo_o[o],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            next!(p)
            d_m_a_r_wo_o[o] = sum(sum(sum(sum(d_a1_z0_r_wo_o[o],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            next!(p)
        end

        p = Progress(3*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes*4, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
        Threads.@threads for (o,(u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
            distr_no_reform_w_occ[t+1][o][:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] = (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*d_m_a_nr_w_o[o]
            next!(p)
            distr_reform_w_occ[t+1][o][:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] = (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*d_m_a_r_w_o[o]
            next!(p)
            distr_no_reform_wo_occ[t+1][o][:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] = (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*d_m_a_nr_wo_o[o]
            next!(p)
            distr_reform_wo_occ[t+1][o][:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] = (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*d_m_a_r_wo_o[o]
            next!(p)
        end

        p = Progress(3*4, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
        Threads.@threads for o = 1:3
            distr_no_reform_w_occ[t+1][o] .+= new_d_nr_w_o[o]
            next!(p)
            distr_reform_w_occ[t+1][o] .+= new_d_r_w_o[o]
            next!(p)
            distr_no_reform_wo_occ[t+1][o] .+= new_d_nr_wo_o[o]
            next!(p)
            distr_reform_wo_occ[t+1][o] .+= new_d_r_wo_o[o]
            next!(p)
        end

        p = Progress(2+3*2+3*6*4, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
        measures_nr_wo_o = compute_household_measures_wo_occ(policy_set_no_reform_wo_occ[1][t+1], r_s_no_reform_wo_occ[t+1], w_s_no_reform_wo_occ[t+1], lambda_s_no_reform[t+1], asset_grid)
        next!(p)
        measures_r_wo_o = compute_household_measures_wo_occ(policy_set_reform_wo_occ[1][t+1], r_s_reform_wo_occ[t+1], w_s_reform_wo_occ[t+1], lambda_s_reform[t+1], asset_grid)
        next!(p)
        Threads.@threads for o=1:3
            measures_nr_w_o = compute_household_measures_w_occ(policy_set_no_reform_w_occ[1][t+1][o], r_s_no_reform_w_occ[t+1], w_s_no_reform_w_occ[t+1], lambda_s_no_reform[t+1], asset_grid)
            next!(p)
            measures_r_w_o = compute_household_measures_w_occ(policy_set_reform_w_occ[1][t+1][o], r_s_reform_w_occ[t+1], w_s_reform_w_occ[t+1], lambda_s_reform[t+1], asset_grid)
            next!(p)
            Threads.@threads for m=1:6
                hm_nr_w_o[t+1][o][m] = compute_inequality_measures(measures_nr_w_o[m], distr_no_reform_w_occ[t+1][o])
                next!(p)
                hm_r_w_o[t+1][o][m] = compute_inequality_measures(measures_r_w_o[m], distr_reform_w_occ[t+1][o])
                next!(p)
                hm_nr_wo_o[t+1][o][m] = compute_inequality_measures(measures_nr_wo_o[m][o], distr_no_reform_wo_occ[t+1][o])
                next!(p)
                hm_r_wo_o[t+1][o][m] = compute_inequality_measures(measures_r_wo_o[m][o], distr_reform_wo_occ[t+1][o])
                next!(p)
            end
        end

    end

    println_sameline("MODE_OF_EXPERIMENT=$MODE_OF_EXPERIMENT - save the results in jld2")

    @save "$(LOCAL_DIR)transition_results.jld2" number_asset_grid a_max asset_grid policy_set_reform_w_occ policy_set_reform_wo_occ policy_set_no_reform_w_occ policy_set_no_reform_wo_occ distr_no_reform_w_occ distr_reform_w_occ distr_no_reform_wo_occ distr_reform_wo_occ hm_nr_w_o hm_r_w_o hm_nr_wo_o hm_r_wo_o

    println_sameline("MODE_OF_EXPERIMENT=$MODE_OF_EXPERIMENT - plot the results")
    # household measures - 6, graphs - 3: w_occ, wo_occ, diff(w_occ-wo_occ)
    #plots diff %
    ps_dp = Array{Any}(undef,6,3)
    #plots diff absolute
    ps_d = Array{Any}(undef,6,3)

    # household measures - 6, occupations - 3
    ps_occ_dp = Array{Any}(undef,6,3)

    measures_name = ["consumption", "earnings", "income", "wealth", "capital_income", "savings"]
    p = Progress(6*3, dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=25)
    for m = 1:6
        # w_occ
        w_occ_diff = Array{Any}(undef,3)
        w_occ_diff_percent = Array{Any}(undef,3)
        for o=1:3
            w_occ_diff[o] = zeros(TIME_PERIODS)
            w_occ_diff_percent[o] = zeros(TIME_PERIODS)
            for lt=1:TIME_PERIODS
                w_occ_diff[o][lt] = hm_r_w_o[lt][o][m][1]-hm_nr_w_o[lt][o][m][1]
                w_occ_diff_percent[o][lt] = w_occ_diff[o][lt]/hm_nr_w_o[lt][o][m][1]
                if m==6
                    w_occ_diff_percent[o][lt] = w_occ_diff[o][lt]/hm_nr_w_o[lt][o][4][1]
                end
            end
        end
        x_label = "Time (Years)"
        y_label = ""
        # if m==1
        #     y_label = "d occmob+finlib"
        # end
        LEGENDPOS = false
        if m==4
            LEGENDPOS = :topright
        end

        ps_dp[m,1] = create_combined_plot(collect(1:TIME_PERIODS),x_label, w_occ_diff_percent, ["W","SP","EMP"], y_label#="$(inequality_measures[im]) $(measures[m]) (diff%)"=#, [w_occ_diff_percent[1][1],w_occ_diff_percent[2][1],w_occ_diff_percent[3][1]],[w_occ_diff_percent[1][end],w_occ_diff_percent[2][end],w_occ_diff_percent[3][end]],true, ["W","SP","EMP"],LEGENDPOS)
        #display(ps_dp[m,1])
        savefig(ps_dp[m,1],"$(LOCAL_DIR)time_mean_$(measures_name[m])_diff_percent_effect_finlib_and_occmob.png")
        ps_d[m,1] = create_combined_plot(collect(1:TIME_PERIODS),x_label, w_occ_diff, ["W","SP","EMP"], y_label#="$(inequality_measures[im]) $(measures[m]) (diff%)"=#, [w_occ_diff[1][1],w_occ_diff[2][1],w_occ_diff[3][1]],[w_occ_diff[1][end],w_occ_diff[2][end],w_occ_diff[3][end]],false, ["W","SP","EMP"],LEGENDPOS)
        # display(ps_dp[m,1])
        next!(p)

        # wo_occ
        wo_occ_diff = Array{Any}(undef,3)
        wo_occ_diff_percent = Array{Any}(undef,3)
        for o=1:3
            wo_occ_diff[o] = zeros(TIME_PERIODS)
            wo_occ_diff_percent[o] = zeros(TIME_PERIODS)
            for lt=1:TIME_PERIODS
                wo_occ_diff[o][lt] = hm_r_wo_o[lt][o][m][1].-hm_nr_wo_o[lt][o][m][1]
                wo_occ_diff_percent[o][lt] = wo_occ_diff[o][lt]/hm_nr_w_o[lt][o][m][1]
                if m==6
                    wo_occ_diff_percent[o][lt] = wo_occ_diff[o][lt]/hm_nr_w_o[lt][o][4][1]
                end
            end
        end
        x_label = "Time (Years)"
        y_label = ""
        # if m==1
        #     y_label = "d finlib"
        # end
        LEGENDPOS = false
        # if m==6
        #     LEGENDPOS = :topright
        # end
        ps_dp[m,2] = create_combined_plot(collect(1:TIME_PERIODS),x_label, wo_occ_diff_percent, ["W","SP","EMP"], y_label#="$(inequality_measures[im]) $(measures[m]) (diff%)"=#, [wo_occ_diff_percent[1][1],wo_occ_diff_percent[2][1],wo_occ_diff_percent[3][1]],[wo_occ_diff_percent[1][end],wo_occ_diff_percent[2][end],wo_occ_diff_percent[3][end]],true, ["W","SP","EMP"],LEGENDPOS)
        #display(ps_dp[m,2])
        savefig(ps_dp[m,2],"$(LOCAL_DIR)time_mean_$(measures_name[m])_diff_percent_effect_finlib_part.png")
        ps_d[m,2] = create_combined_plot(collect(1:TIME_PERIODS),x_label, wo_occ_diff, ["W","SP","EMP"], y_label#="$(inequality_measures[im]) $(measures[m]) (diff%)"=#, [wo_occ_diff[1][1],wo_occ_diff[2][1],wo_occ_diff[3][1]],[wo_occ_diff[1][end],wo_occ_diff[2][end],wo_occ_diff[3][end]],false, ["W","SP","EMP"],LEGENDPOS)
        # display(ps_dp[m,2])
        next!(p)

        # diff between w_occ and wo_occ
        diff = Array{Any}(undef,3)
        diff_percent = Array{Any}(undef,3)
        for o=1:3
            diff[o] = w_occ_diff[o].-wo_occ_diff[o]
            diff_percent[o] = zeros(TIME_PERIODS)
            for lt=1:TIME_PERIODS
                diff_percent[o][lt] = diff[o][lt]/hm_nr_w_o[lt][o][m][1]
                if m==6
                    diff_percent[o][lt] = diff[o][lt]/hm_nr_w_o[lt][o][4][1]
                end
            end
        end
        x_label = "Time (Years)"#"Mean $(measures_name[m])"
        y_label = ""
        # if m==1
        #     y_label = "d occmob"
        # end
        LEGENDPOS = false
        # if m==6
        #     LEGENDPOS = :topright
        # end
        ps_dp[m,3] = create_combined_plot(collect(1:TIME_PERIODS),x_label, diff_percent, ["W","SP","EMP"], y_label, [diff_percent[1][1],diff_percent[2][1],diff_percent[3][1]],[diff_percent[1][end],diff_percent[2][end],diff_percent[3][end]],true, ["W","SP","EMP"],LEGENDPOS)
        #display(ps_dp[m,3])
        savefig(ps_dp[m,3],"$(LOCAL_DIR)time_mean_$(measures_name[m])_diff_percent_effect_occmob_part.png")
        ps_d[m,3] = create_combined_plot(collect(1:TIME_PERIODS),x_label, diff, ["W","SP","EMP"], y_label, [diff[1][1],diff[2][1],diff[3][1]],[diff[1][end],diff[2][end],diff[3][end]],false, ["W","SP","EMP"],LEGENDPOS)
        # display(ps_dp[m,3])
        next!(p)

        x_label = "Time (Years)"
        y_label = ""
        # if m==1
        #     y_label = "d occmob+finlib"
        # end
        LEGENDPOS = false
        if m==4
            LEGENDPOS = :topright
        end
        OCCS_NAME = ["W","SP","EMP"]
        for o=1:3
            # w_occ_diff_percent[o], wo_occ_diff_percent[o], diff_percent[o]
            ps_occ_dp[m,o] = create_combined_plot(collect(1:TIME_PERIODS),x_label, [w_occ_diff_percent[o], wo_occ_diff_percent[o], diff_percent[o]], ["Total","Direct","Indirect"], y_label, [w_occ_diff_percent[o][1], wo_occ_diff_percent[o][1], diff_percent[o][1]],[w_occ_diff_percent[o][end], wo_occ_diff_percent[o][end], diff_percent[o][end]],true, ["H","ENT","SP"],LEGENDPOS)
            display(ps_occ_dp[m,o])
            savefig(ps_occ_dp[m,o],"$(LOCAL_DIR)time_mean_$(measures_name[m])_diff_percent_effect_$(OCCS_NAME[o]).png")
        end
    end

    # display(plot(ps_dp[1,1],ps_dp[1,2],ps_dp[1,3], layout=(3,1)))
    # m = 1         2           3       4       5               6
    # consumption, earnings, income, wealth, capital_income, savings
    # display(plot(ps_dp[1,1],ps_dp[3,1],ps_dp[4,1],
    #              ps_dp[1,2],ps_dp[3,2],ps_dp[4,2],
    #              ps_dp[1,3],ps_dp[3,3],ps_dp[4,3],
    #              layout=(3,3)))
    # # W - purple; SP - red; EMP - green;
    # display(plot(ps_d[1,1],ps_d[3,1],ps_d[4,1],
    #              ps_d[1,2],ps_d[3,2],ps_d[4,2],
    #              ps_d[1,3],ps_d[3,3],ps_d[4,3],
    #              layout=(3,3)))

    # m = 1         2           3       4       5               6
    # consumption, earnings, income, wealth, capital_income, savings
    # o = 1    2                 3
    # workers, sole-proprietors, employers
    display(plot(ps_occ_dp[1,1],ps_occ_dp[1,2],ps_occ_dp[1,3],
                 ps_occ_dp[2,1],ps_occ_dp[2,2],ps_occ_dp[2,3],
                 ps_occ_dp[3,1],ps_occ_dp[3,2],ps_occ_dp[3,3],
                 layout=(3,3)))

    return [number_asset_grid, a_max, asset_grid, policy_set_reform_w_occ, policy_set_reform_wo_occ, policy_set_no_reform_w_occ, policy_set_no_reform_wo_occ, distr_no_reform_w_occ, distr_reform_w_occ, distr_no_reform_wo_occ, distr_reform_wo_occ, hm_nr_w_o, hm_r_w_o, hm_nr_wo_o, hm_r_wo_o, ps_d, ps_dp]
end

# MODE_OF_EXPERIMENT = {1 - Take stationary distribution from the model_w_occ as starting distr,
#                       2 - Take stationary distribution from the model_wo_occ as starting distr}

LOCAL_DIR_ALTERNATIVE = "$(@__DIR__)/Results/Transitionary/TDE_comparison/Inequality/Transition/Alternative_intersection/"
if Sys.iswindows()
    LOCAL_DIR_ALTERNATIVE = "$(@__DIR__)\\Results\\Transitionary\\TDE_comparison\\Inequality\\Transition\\Alternative_intersection\\"
end
mkpath(LOCAL_DIR_ALTERNATIVE)
results_of_the_experiment_3 = conduct_the_experiment(3, LOCAL_DIR_ALTERNATIVE)
# ps_dp3 = results_of_the_experiment_3[end]
# display(plot(ps_dp3[1,1],ps_dp3[3,1],ps_dp3[4,1],
#              ps_dp3[1,2],ps_dp3[3,2],ps_dp3[4,2],
#              ps_dp3[1,3],ps_dp3[3,3],ps_dp3[4,3],
#              layout=(3,3)))

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE_comparison/Inequality/Transition/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE_comparison\\Inequality\\Transition\\"
end
mkpath(LOCAL_DIR)
results_of_the_experiment_1 = conduct_the_experiment(1, LOCAL_DIR)
# ps_dp1 = results_of_the_experiment_1[end]
# display(plot(ps_dp1[1,1],ps_dp1[3,1],ps_dp1[4,1],
#              ps_dp1[1,2],ps_dp1[3,2],ps_dp1[4,2],
#              ps_dp1[1,3],ps_dp1[3,3],ps_dp1[4,3],
#              layout=(3,3)))

# LOCAL_DIR_ALTERNATIVE = "$(@__DIR__)/Results/Transitionary/TDE_comparison/Inequality/Transition/Alternative/"
# if Sys.iswindows()
#     LOCAL_DIR_ALTERNATIVE = "$(@__DIR__)\\Results\\Transitionary\\TDE_comparison\\Inequality\\Transition\\Alternative\\"
# end
# mkpath(LOCAL_DIR_ALTERNATIVE)
# results_of_the_experiment_2 = conduct_the_experiment(2, LOCAL_DIR_ALTERNATIVE)
# ps_dp2 = results_of_the_experiment_2[end]
# display(plot(ps_dp2[1,1],ps_dp2[3,1],ps_dp2[4,1],
#              ps_dp2[1,2],ps_dp2[3,2],ps_dp2[4,2],
#              ps_dp2[1,3],ps_dp2[3,3],ps_dp2[4,3],
#              layout=(3,3)))



# LOCAL_DIR_ALTERNATIVE = "$(@__DIR__)/Results/Transitionary/TDE_comparison/Inequality/Transition/Alternative_full_different_distr/"
# if Sys.iswindows()
#     LOCAL_DIR_ALTERNATIVE = "$(@__DIR__)\\Results\\Transitionary\\TDE_comparison\\Inequality\\Transition\\Alternative_full_different_distr\\"
# end
# mkpath(LOCAL_DIR_ALTERNATIVE)
# results_of_the_experiment_4 = conduct_the_experiment(4, LOCAL_DIR_ALTERNATIVE)
# ps_dp4 = results_of_the_experiment_4[end]
# display(plot(ps_dp4[1,1],ps_dp4[3,1],ps_dp4[4,1],
#              ps_dp4[1,2],ps_dp4[3,2],ps_dp4[4,2],
#              ps_dp4[1,3],ps_dp4[3,3],ps_dp4[4,3],
#              layout=(3,3)))
