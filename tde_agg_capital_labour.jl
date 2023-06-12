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

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE\\"
end

TIME_PERIODS = 50#100
@load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2

T = trans_res[1]
K_ds = trans_res[5][3]
K_ss = trans_res[5][4]
display(plot(1:T, [K_ds K_ss], legend=false))

K_s = (K_ds.+K_ss)./2
display(plot(1:T, K_s, legend=false))

L_ds = trans_res[5][5]
L_ss = trans_res[5][6]
display(plot(1:T, [L_ds L_ss], legend=false))

L_s = (L_ds.+L_ss)./2
display(plot(1:T, L_s, legend=false))
