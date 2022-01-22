include("print_sameline.jl")

include("parameters.jl")
include("AllubErosa.jl")

function steady_state(R, W, global_params, global_approx_params, model_params)

    rw_maxiters = 100

    #           1       2     3     4       5   6       7   8       9       10          11              12      13          14          15              16          17              18              19       20
    #         lambda, beta, delta, gamma, eta, theta, c_e, rho_m, rho_w, sigma_eps_m, rho_eps_m_w, sigma_zeta, p_alpha, eta_alpha, prob_node1_alpha, mu_m_alpha, rho_alpha_m_w, sigma_alpha_w, sigma_eps_w, crra
    #  model_params

    r_min = -model_params[3]
    r_max = 1/model_params[2]-1
    # approximation object for skills nodes
    #               1               2       3       4       5       6           7
    #               z_m_nodes, z_w_nodes, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes
    approx_object = build_skill_nodes(global_approx_params, model_params)

    println("Initiate the search for general equilibrium factor prices                                      \r")
    #R = (r_min+r_max)/2#R_
    #W = (w_min+w_max)/2#W_
    if R <= r_min || R >= r_max
        R = (r_min+r_max)/2
    end
    if W <= w_min || W >= w_max
        W = (w_min+w_max)/2
    end
    res = AllubErosa(R,W, global_params, global_approx_params, model_params, approx_object)

    len_r = res[10]
    len_w = res[11]
    total_len = abs(len_r)+abs(len_w)

    best_total_len = total_len
    best_len_r = len_r
    best_len_w = len_w
    best_res = copy(res)
    best_R = R
    best_W = W

    new_R = R*(1.0+len_r)/(1.0-len_r) + (1.0-2.0*r_min)*len_r/(1.0-len_r)#R*(1.0+sign(R)*100*len_r)/(1.0-sign(R)*100*len_r)#R + len_r*1e-1*2#
    new_R = min(max(r_min, new_R), r_max)
    new_W = W*(1.0+len_w)/(1.0-len_w)
    new_W = min(max(w_min, new_W), w_max)

    rw_iters = 1

    gen_tol_x = global_params[1]
    gen_tol_f = global_params[2]

    cand_r = [R]
    cand_w = [W]
    cand_len_r = [len_r]
    cand_len_w = [len_w]

    old_old_R = R
    old_R = R
    old_old_W = W
    old_W = W
    old_old_len_r = len_r
    old_len_r = len_r
    old_old_len_w = len_w
    old_len_w = len_w
    old_total_len = total_len
    old_new_R = new_R
    old_new_W = new_W
    new_R = new_R*0.1 + R*(1.0-0.9)
    new_W = new_W*0.1 + W*(1.0-0.9)

    println_sameline("#0 - r:$(round(R*100;digits=4))%, w:$(round(W;digits=6)) - total_len:$(round(total_len;digits=6)) - len_r:$(round(len_r;digits=6)) - len_w:$(round(len_w;digits=6)) - new_r:$(round(new_R*100;digits=4))%, new_w:$(round(new_W;digits=6))")

    while rw_iters <= rw_maxiters && (abs(len_r) > gen_tol_f || abs(len_w) > gen_tol_f)

        old_old_R = old_R
        old_R = R
        old_old_W = old_W
        old_W = W
        old_old_len_r = old_len_r
        old_len_r = len_r
        old_old_len_w = old_len_w
        old_len_w = len_w
        old_total_len = total_len
        old_new_R = new_R
        old_new_W = new_W

        R = new_R
        W = new_W
        try
            res = AllubErosa(R,W, global_params, global_approx_params, model_params, approx_object)

            len_r = res[10]
            len_w = res[11]
            total_len = abs(len_r)+abs(len_w)

            new_R = R - len_r*(R-old_R)/(len_r-old_len_r)
            best_new_R = R - len_r*(R-best_R)/(len_r-best_len_r)
            new_R = (new_R+best_new_R)/2.0
            if new_R <= r_min
                new_R = 0.9*r_min + (1.0-0.9)*R
            elseif new_R >= r_max
                new_R = 0.9*r_max + (1.0-0.9)*R
            end
            new_W = W - len_w*(W-old_W)/(len_w-old_len_w)
            best_new_W = W - len_w*(W-best_W)/(len_w-best_len_w)
            new_W = (new_W+best_new_W)/2.0
            if new_W <= w_min
                new_W = 0.9*w_min + (1.0-0.9)*W
            elseif new_W >= w_max
                new_W = 0.9*w_max + (1.0-0.9)*W
            end

            if total_len < best_total_len
                best_total_len = total_len
                best_len_r = len_r
                best_len_w = len_w
                best_res = copy(res)
                best_R = R
                best_W = W
            end

        catch e
            println(e)
            R = old_R
            old_R = old_old_R
            W = old_W
            old_W = old_old_W
            old_len_r = old_old_len_r
            old_len_w = old_old_len_w
            total_len = old_total_len

            new_R = (R + new_R)/2.0
            new_W = (W + new_W)/2.0

            rw_iters -= 1

        end


        if max(abs(new_R-R), abs(new_W-W)) < gen_tol_x
            new_R = R + len_r
            new_R = min(max(r_min, new_R), r_max)
            new_W = W + len_w
            new_W = min(max(w_min, new_W), w_max)
        end
        println_sameline("#$(rw_iters) - r:$(round(R*100;digits=4))%, w:$(round(W;digits=6)) - total_len:$(round(total_len;digits=6)) - len_r:$(round(len_r;digits=6)) - len_w:$(round(len_w;digits=6)) - new_r:$(round(new_R*100;digits=4))%, new_w:$(round(new_W;digits=6))")

        # visualisation
        #=
        append!(cand_r, R)
        append!(cand_w, W)
        append!(cand_len_r, len_r)
        append!(cand_len_w, len_w)
        p1 = plot(cand_r,cand_w, legend=false, title = "$(rw_iters)")
        p21 = plot(cand_r,[cand_len_r, zeros(length(cand_len_r))], legend=false)
        p22 = plot(cand_w,[cand_len_r, zeros(length(cand_len_r))], legend=false)
        p2 = plot(p21,p22, layout=(2,1), title = "len_r")
        p31 = plot(cand_r,[cand_len_w, zeros(length(cand_len_w))], legend=false)
        p32 = plot(cand_w,[cand_len_w, zeros(length(cand_len_w))], legend=false)
        p3 = plot(p31,p32, layout=(2,1), title = "len_w")
        #p4 = plot(cand_r,cand_w,abs.(cand_len_r).+abs.(cand_len_w), legend=false)
        p41 = plot(cand_r,[abs.(cand_len_r).+abs.(cand_len_w), zeros(length(cand_r))], legend=false)
        p42 = plot(cand_w,[abs.(cand_len_r).+abs.(cand_len_w), zeros(length(cand_w))], legend=false)
        p4 = plot(p41,p42, layout=(2,1), title = "len")
        display(plot(p1,p2,p3,p4, layout=(2,2)))
        =#

        rw_iters += 1

    end

    r_min = -DELTA#-delta
    r_max = 1/BETA-1#1/beta-1

    if total_len < best_total_len #rw_iters < rw_maxiters
        println("The interest rate - $(R*100)% and the wage - $(W) with total_len - $(total_len)")
        return [res, R, W, approx_object, model_params]
    else
        println("The interest rate - $(best_R*100)% and the wage - $(best_W) with total_len - $(best_total_len)")
        return [best_res, best_R, best_W, approx_object, model_params]
    end
end
