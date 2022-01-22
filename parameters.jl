using QuantEcon
using JLD2
include("VAR_Markov.jl")

function build_skill_nodes(approx_ps, ps)
    ########################
    # persistant shocks to skills
    number_u_m_nodes = approx_ps[2]#7
    number_u_w_nodes = approx_ps[3]#7

    sigma_eps_m = ps[10]#1.14553769575612#sqrt(1.145) #Variance innovation managerial skill
    sigma_eps_w = ps[19]#0.07352027977526#sqrt(0.073) #Variance innovation working skill

    rho_m = ps[8]#0.78839752660496#0.788#*(1-0.041) #Autocorrelation managerial skill
    rho_w = ps[9]#0.96#*(1-0.041) #Autocorrelation of working-ability shocks

    rho_eps_m_w = ps[11]#0.3351447009#0.335 # Correlation between innovation

    # transitory shock to managerial skill
    number_zeta_nodes = approx_ps[4]#3

    sigma_zeta = ps[12]#0.211613239161#sqrt(0.212)

    # fixed effects on managerial skill
    number_alpha_m_nodes = approx_ps[5]#3#6 # n>=3

    p_alpha = ps[13]#1-0.958774210075#0.041   # probability of new skills draw
    eta_alpha = ps[14]#5.403585195245#5.40  # Pareto tail managerial skill
    prob_node1_alpha = ps[15]#0.22878505546#0.229  # probability of first fixed effect  on managerial skill
    prob_nodeend_alpha = 0.99999381691 #probability of last fixed effect on managerial skill

    mu_m_alpha = ps[16]#-2.975991405139#-2.976 # location of first fixed effect on managerial skill distribution

    number_alpha_w_nodes = approx_ps[6]#3

    rho_alpha_m_w = ps[17]#0.147638317002#0.145
    sigma_alpha_w = ps[18]#0.42#sqrt(0.42)

    ############################
    ## Approximation objects ###

    # transitory shock to managerial skill
    if number_zeta_nodes > 1
        mc_zeta = tauchen(number_zeta_nodes,0.0,sqrt(sigma_zeta),0.0,1)#rouwenhorst(number_zeta_nodes,0.0,sqrt(sigma_zeta))
        zeta_nodes = collect(mc_zeta.state_values)
        P_zeta = stationary_distributions(mc_zeta)[1]
    else
        zeta_nodes = [0.0]
        P_zeta = [1.0]
    end

    #display("Transitory shocks:")
    #display("P_zeta:"*string(P_zeta))
    #display("zeta_nodes:"*string(zeta_nodes))


    if number_u_m_nodes == 7 && number_u_w_nodes == 7
        ymat2 = open("$(@__DIR__)/ymat2.txt") do f
            read(f, String)
        end
        s1 = parse.(Float64,split(ymat2))
        u_mw_nodes = zeros(2,number_u_m_nodes*number_u_w_nodes)
        for i in 1:length(s1)
            if i%2 == 1
                u_mw_nodes[1,cld(i,2)] = s1[i]
            else
                u_mw_nodes[2,cld(i,2)] = s1[i]
            end
        end

        pimat = open("$(@__DIR__)/pimat.txt") do f
            read(f, String)
        end
        s2 = parse.(Float64,split(pimat))
        P_u_temp = zeros(number_u_m_nodes*number_u_w_nodes,number_u_m_nodes*number_u_w_nodes)
        for i in 1:length(s2)
            P_u_temp[i] = s2[i]
        end
        P_u = transpose(P_u_temp)

    elseif number_u_m_nodes > 1 && number_u_w_nodes > 1
        A2 = [rho_w 0.0; 0.0 rho_m]
        SIGMA = [sigma_eps_w rho_eps_m_w*sqrt(sigma_eps_w*sigma_eps_m); rho_eps_m_w*sqrt(sigma_eps_w*sigma_eps_m) sigma_eps_m]
        N = [number_u_w_nodes; number_u_m_nodes]
        P_u, u_mw_nodes, u_nodes = fn_var_to_markov(A2,SIGMA,N)
    elseif number_u_m_nodes > 1
        mc_u_m = tauchen(number_u_m_nodes,0.0,sqrt(sigma_eps_m))
        u_m_nodes = collect(mc_u_m.state_values)
        P_u_m = stationary_distributions(mc_u_m)[1]

        P_u = mc_u_m.p
        u_mw_nodes = zeros(2,number_u_m_nodes)
        u_mw_nodes[1,:] .= zeros(number_u_m_nodes)
        u_mw_nodes[2,:] .= u_m_nodes
    end

    #display(u_mw_nodes)

    number_u_nodes = size(u_mw_nodes,2)

    mc_u = MarkovChain(P_u)
    stat_P_u = stationary_distributions(mc_u)[1]

    mean_u_w = stat_P_u'*u_mw_nodes[1,:]
    mean_u_m = stat_P_u'*u_mw_nodes[2,:]
    mean_u_w2= stat_P_u'*(u_mw_nodes[1,:].^2)
    mean_u_m2= stat_P_u'*(u_mw_nodes[2,:].^2)
    var_u_w = mean_u_w2 - mean_u_w^2
    var_u_m = mean_u_m2 - mean_u_m^2

    #display("Persistent shocks:")
    #display("P_u"*string(P_u))
    #display("stat_P_u"*string(stat_P_u))
    #display("u_mw_nodes[1]"*string(u_mw_nodes[1,:]))
    #display("u_mw_nodes[2]"*string(u_mw_nodes[2,:]))


    # fixed effect to managerial skill
    if number_alpha_m_nodes > 1 && number_alpha_w_nodes > 1
        # fortran version
        cdf_m = collect(range(prob_node1_alpha; stop = prob_nodeend_alpha, length=number_alpha_m_nodes))
        stat_P_alpha_m_temp = zeros(number_alpha_m_nodes)
        for i in number_alpha_m_nodes-1:-1:2
            stat_P_alpha_m_temp[i] = cdf_m[i] - cdf_m[i-1]
        end
        stat_P_alpha_m_temp[1] = prob_node1_alpha
        stat_P_alpha_m_temp[end] = 1.0-prob_nodeend_alpha

        cut_off_m = zeros(number_alpha_m_nodes)
        cut_off_m[1] = stat_P_alpha_m_temp[1]/2
        cut_off_m[end] = 1 - stat_P_alpha_m_temp[end]/2
        for i in number_alpha_m_nodes-1:-1:2
            cut_off_m[i] = (cdf_m[i] + cdf_m[i-1])/2
        end

        stat_P_alpha_m_temp = stat_P_alpha_m_temp./sum(stat_P_alpha_m_temp)

        alpha_m_nodes_temp = zeros(number_alpha_m_nodes)
        alpha_m_nodes_temp = 1 ./ ((1 .- cut_off_m).^(1/eta_alpha))
        alpha_m_nodes_temp = (eta_alpha/(eta_alpha-1)).*alpha_m_nodes_temp
        alpha_m_nodes_temp[1] = 1.0

        alpha_m_nodes = alpha_m_nodes_temp #.+ (mu_m_alpha - alpha_m_nodes_temp[1])

        P_alpha_m = stat_P_alpha_m_temp

        #display(alpha_m_nodes')
        E_alpha_m = alpha_m_nodes'*P_alpha_m
        #std_alpha_m = sqrt(sum([(alpha_m_nodes[i] - E_alpha_m)^2*stat_P_alpha_m[i] for i = 1:number_alpha_m_nodes]))
        std_alpha_m = (alpha_m_nodes.^2)'*P_alpha_m - E_alpha_m^2

        var_alpha_w_cond_alpha_m = sigma_alpha_w - rho_alpha_m_w^2*std_alpha_m
        var_alpha_w_cond_alpha_m = max(1e-9,var_alpha_w_cond_alpha_m)

        # fixed effect to working skill
        alpha_w_nodes_temp = zeros(number_alpha_m_nodes,number_alpha_w_nodes)
        P_alpha_w_temp = zeros(number_alpha_m_nodes,number_alpha_w_nodes)

        mc_alpha_w = tauchen(number_alpha_w_nodes,0.0,sqrt(var_alpha_w_cond_alpha_m),0.0,1.5)#rouwenhorst(number_alpha_w_nodes,0.0,sigma_w_f)

        for alpha_m_i = 1:number_alpha_m_nodes

            alpha_w_nodes_temp[alpha_m_i,:] .= collect(mc_alpha_w.state_values) .+ ( rho_alpha_m_w * (alpha_m_nodes[alpha_m_i] - E_alpha_m) * sqrt(sigma_alpha_w / std_alpha_m) )
            P_alpha_w_temp[alpha_m_i,:] = stationary_distributions(mc_alpha_w)[1]
        end

        alpha_w_nodes = alpha_w_nodes_temp
        P_alpha_w = P_alpha_w_temp

        P_alpha = P_alpha_w.*P_alpha_m
    else
        alpha_m_nodes = [0.0]
        alpha_w_nodes = zeros(1,1)
        P_alpha = ones(1,1)
    end
    #display("Fixed effects")
    #display("alpha_m_nodes"*string(alpha_m_nodes))
    #display("alpha_w_nodes"*string(alpha_w_nodes))
    #display("P_alpha"*string(P_alpha))


    # construct nodes for combination of shocks
    z_m_nodes = zeros(number_u_nodes,number_alpha_m_nodes,number_zeta_nodes)
    z_w_nodes = zeros(number_u_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    for u_i in 1:number_u_nodes
        for alpha_m_i in 1:number_alpha_m_nodes
            for zeta_i in 1:number_zeta_nodes
                z_m_nodes[u_i,alpha_m_i,zeta_i] = (exp(u_mw_nodes[2,u_i])/exp(var_u_m/2)) * exp(alpha_m_nodes[alpha_m_i]) * exp(zeta_nodes[zeta_i])
            end
            for alpha_w_i in 1:number_alpha_w_nodes
                z_w_nodes[u_i,alpha_m_i,alpha_w_i] = (exp(u_mw_nodes[1,u_i])/exp(var_u_w/2)) * exp(alpha_w_nodes[alpha_m_i,alpha_w_i])
            end
        end
    end
    z_m_nodes .*= exp(mu_m_alpha)

    return [z_m_nodes, z_w_nodes, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes]
end
