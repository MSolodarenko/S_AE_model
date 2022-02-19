include("print_sameline.jl")

include("parameters.jl")
include("AllubErosa.jl")

function draw_best_grid(best_total_len_grid, best_len_r_grid,best_len_w_grid, best_R_grid,best_W_grid, R,W, len_r,len_w, new_R,new_W)
    len_r_grid = []
    len_w_grid = []
    R_grid = []
    W_grid = []
    for r_i=1:2
        for w_i=1:2
            if best_total_len_grid[r_i,w_i] != Inf
                push!(len_r_grid, best_len_r_grid[r_i,w_i])
                push!(len_w_grid, best_len_w_grid[r_i,w_i])
                push!(R_grid, best_R_grid[r_i,w_i])
                push!(W_grid, best_W_grid[r_i,w_i])
            end
        end
    end

    p1 = scatter(len_r_grid, len_w_grid, legend=false)
    for i1 = 1:length(len_r_grid)
        for i2 = i1+1:length(len_w_grid)
            plot!(p1, [len_r_grid[i1],len_r_grid[i2]], [len_w_grid[i1],len_w_grid[i2]], legend=false)
        end
    end
    hline!(p1, [0.0])
    vline!(p1, [0.0])
    scatter!(p1, [len_r], [len_w], legend=false)

    p2 = scatter(R_grid, W_grid, legend=false)
    for i1 = 1:length(R_grid)
        for i2 = i1+1:length(W_grid)
            plot!(p2, [R_grid[i1],R_grid[i2]], [W_grid[i1],W_grid[i2]], legend=false)
        end
    end
    scatter!(p2, [R,new_R], [W,new_W], legend=false)
    plot!(p2, [R,new_R], [W,new_W], arrow=true, legend=false)

    display(plot(p1,p2, layout=(1,2)))
end

function steady_state(R, W, global_params, global_approx_params, model_params)

    rw_maxiters = 60#100

    #           1       2     3     4       5   6       7   8       9       10          11              12      13          14          15              16          17              18              19       20
    #         lambda, beta, delta, gamma, eta, theta, c_e, rho_m, rho_w, sigma_eps_m, rho_eps_m_w, sigma_zeta, p_alpha, eta_alpha, prob_node1_alpha, mu_m_alpha, rho_alpha_m_w, sigma_alpha_w, sigma_eps_w, crra
    #  model_params

    r_min = -model_params[3]
    r_max = 1/model_params[2]-1

    w_min = 0.18#0.01#0.25
    w_max = 0.3#0.99#0.5
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
    res = []
    try
        println_sameline("Initial guess - R:$R - W:$W")
        res = AllubErosa(R,W, global_params, global_approx_params, model_params, approx_object)
    catch e
        R = (R+(r_min+r_max)/2)/2
        W = (W+(w_min+w_max)/2)/2
        println_sameline(("Move guess to the middle - R:$R - W:$W - ",e))
        try
            res = AllubErosa(R,W, global_params, global_approx_params, model_params, approx_object)
        catch e
            R = (r_min+r_max)/2
            W = (w_min+w_max)/2
            println_sameline(("Set guess as the middle - R:$R - W:$W - ",e))
            try
                res = AllubErosa(R,W, global_params, global_approx_params, model_params, approx_object)
            catch e
                println_sameline("Error in iteration #0")
                throw(error(e))
            end
        end
    end

    len_r = res[10]
    len_w = res[11]
    total_len = abs(len_r)+abs(len_w)

    new_R = R*(1.0+len_r)/(1.0-len_r) + (1.0-2.0*r_min)*len_r/(1.0-len_r)#R*(1.0+sign(R)*100*len_r)/(1.0-sign(R)*100*len_r)#R + len_r*1e-1*2#
    new_R = new_R*0.1 + R*(1.0-0.1)
    if new_R <= r_min
        new_R = 0.9*r_min + (1.0-0.9)*R
    elseif new_R >= r_max
        new_R = 0.9*r_max + (1.0-0.9)*R
    end

    new_W = W*(1.0+len_w)/(1.0-len_w)
    new_W = new_W*0.1 + W*(1.0-0.1)
    if new_W <= w_min
        new_W = 0.9*w_min + (1.0-0.9)*W
    elseif new_W >= w_max
        new_W = 0.9*w_max + (1.0-0.9)*W
    end

    println_sameline("#0 - r:$(round(R*100;digits=4))%, w:$(round(W;digits=6)) - total_len:$(round(total_len;digits=6)) - len_r:$(round(len_r;digits=6)) - len_w:$(round(len_w;digits=6)) - new_r:$(round(new_R*100;digits=4))%, new_w:$(round(new_W;digits=6))")

    best_total_len_grid = ones(2,2).*Inf
    best_len_r_grid = zeros(2,2)
    best_len_w_grid = zeros(2,2)
    best_R_grid = zeros(2,2)
    best_W_grid = zeros(2,2)
    if len_r > 0.0
        r_i = 1
    else
        r_i = 2
    end
    if len_w > 0.0
        w_i = 1
    else
        w_i = 2
    end
    if best_total_len_grid[r_i,w_i] > total_len
        best_total_len_grid[r_i,w_i] = total_len
        best_len_r_grid[r_i,w_i] = len_r
        best_len_w_grid[r_i,w_i] = len_w
        best_R_grid[r_i,w_i] = R
        best_W_grid[r_i,w_i] = W
        println_sameline("#0 - New best_R_grid[$(r_i)] and best_W_grid[$(w_i)]")
    end
    draw_best_grid(best_total_len_grid, best_len_r_grid,best_len_w_grid, best_R_grid,best_W_grid, R,W, len_r,len_w, new_R,new_W)

    best_total_len = total_len
    best_len_r = len_r
    best_len_w = len_w
    best_res = copy(res)
    best_R = R
    best_W = W
    println_sameline("#0 - New best R and W")

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

    while rw_iters <= rw_maxiters && (abs(len_r) > gen_tol_f || abs(len_w) > gen_tol_f)

        old_old_old_R = old_old_R
        old_old_R = old_R
        old_R = R
        old_old_old_W = old_old_W
        old_old_W = old_W
        old_W = W
        old_old_old_len_r = old_old_len_r
        old_old_len_r = old_len_r
        old_len_r = len_r
        old_old_old_len_w = old_old_len_w
        old_old_len_w = old_len_w
        old_len_w = len_w
        old_total_len = total_len
        old_new_R = new_R
        old_new_W = new_W

        R = new_R
        W = new_W
        try
            if abs(R - best_R) < gen_tol_x && abs(W - best_W) < gen_tol_x
                throw(error("same as best_R and best_W"))
            end
            if abs(R - old_R) < gen_tol_x && abs(W - old_W) < gen_tol_x
                throw(error("same as old_R and old_W"))
            end
            if abs(R - old_old_R) < gen_tol_x && abs(W - old_old_W) < gen_tol_x
                throw(error("same as old_old_R and old_old_W"))
            end
            if abs(R - old_old_old_R) < gen_tol_x && abs(W - old_old_old_W) < gen_tol_x
                throw(error("same as old_old_old_R and old_old_old_W"))
            end
            res = AllubErosa(R,W, global_params, global_approx_params, model_params, approx_object)

            len_r = res[10]
            len_w = res[11]
            total_len = abs(len_r)+abs(len_w)

            if len_r > 0.0
                r_i = 1
                r_i_= 2
            else
                r_i = 2
                r_i_= 1
            end
            if len_w > 0.0
                w_i = 1
                w_i_= 2
            else
                w_i = 2
                w_i_= 1
            end
            if best_total_len_grid[r_i_,w_i_] != Inf && best_total_len_grid[r_i,w_i_] != Inf && (len_w - len_r*(len_w-best_len_w_grid[r_i_,w_i_])/(len_r-best_len_r_grid[r_i_,w_i_])) < 0.0
                println_sameline("#$(rw_iters) - new_R and new_W from the triangle in the two opposite quadrants best_R_W [$r_i_,$w_i_] and [$r_i,$w_i_]")
                new_R1 = R - len_r*(R-best_R_grid[r_i_,w_i_]) /(len_r-best_len_r_grid[r_i_,w_i_])
                new_R2 = R - len_r*(R-best_R_grid[r_i,w_i_]) /(len_r-best_len_r_grid[r_i,w_i_])

                new_W1 = W - len_w*(W-best_W_grid[r_i_,w_i_]) /(len_w-best_len_w_grid[r_i_,w_i_])
                new_W2 = W - len_w*(W-best_W_grid[r_i,w_i_]) /(len_w-best_len_w_grid[r_i,w_i_])

                k_1 = (len_w-best_len_w_grid[r_i_,w_i_])/(len_r-best_len_r_grid[r_i_,w_i_])
                b_1 = len_w - len_r*k_1

                best_len_r_star = best_len_r_grid[r_i,w_i_]
                best_len_w_star = best_len_w_grid[r_i,w_i_]
                best_new_R_star = new_R2
                best_new_W_star = new_W2

                k_2 = best_len_w_star/best_len_r_star
                b_2 = 0.0
                x_starstar = (b_2-b_1)/(k_1-k_2)
                y_starstar = k_1*x_starstar + b_1
                d_1 = sqrt(x_starstar^2+y_starstar^2)
                d_2 = sqrt(best_len_r_star^2+best_len_w_star^2)
                relax = 1.0-d_1/(d_1+d_2)

                new_R = new_R1*relax + best_new_R_star*(1-relax)
                new_W = new_W1*relax + best_new_W_star*(1-relax)

            elseif best_total_len_grid[r_i_,w_i_] != Inf && best_total_len_grid[r_i_,w_i] != Inf && (len_w - len_r*(len_w-best_len_w_grid[r_i_,w_i_])/(len_r-best_len_r_grid[r_i_,w_i_])) > 0.0
                println_sameline("#$(rw_iters) - new_R and new_W from the triangle in the two opposite quadrants best_R_W [$r_i_,$w_i_] and [$r_i_,$w_i]")
                new_R1 = R - len_r*(R-best_R_grid[r_i_,w_i_]) /(len_r-best_len_r_grid[r_i_,w_i_])
                new_R3 = R - len_r*(R-best_R_grid[r_i_,w_i]) /(len_r-best_len_r_grid[r_i_,w_i])

                new_W1 = W - len_w*(W-best_W_grid[r_i_,w_i_]) /(len_w-best_len_w_grid[r_i_,w_i_])
                new_W3 = W - len_w*(W-best_W_grid[r_i_,w_i]) /(len_w-best_len_w_grid[r_i_,w_i])

                k_1 = (len_w-best_len_w_grid[r_i_,w_i_])/(len_r-best_len_r_grid[r_i_,w_i_])
                b_1 = len_w - len_r*k_1

                best_len_r_star = best_len_r_grid[r_i_,w_i]
                best_len_w_star = best_len_w_grid[r_i_,w_i]
                best_new_R_star = new_R3
                best_new_W_star = new_W3

                k_2 = best_len_w_star/best_len_r_star
                b_2 = 0.0
                x_starstar = (b_2-b_1)/(k_1-k_2)
                y_starstar = k_1*x_starstar + b_1
                d_1 = sqrt(x_starstar^2+y_starstar^2)
                d_2 = sqrt(best_len_r_star^2+best_len_w_star^2)
                relax = 1.0-d_1/(d_1+d_2)

                new_R = new_R1*relax + best_new_R_star*(1-relax)
                new_W = new_W1*relax + best_new_W_star*(1-relax)

            elseif best_total_len_grid[r_i_,w_i_] != Inf && best_total_len_grid[r_i,w_i] != Inf && (len_w - len_r*(len_w-best_len_w_grid[r_i_,w_i_])/(len_r-best_len_r_grid[r_i_,w_i_])) * (best_len_w_grid[r_i_,w_i_] - best_len_r_grid[r_i_,w_i_]*(best_len_w_grid[r_i_,w_i_]-best_len_w_grid[r_i,w_i])/(best_len_r_grid[r_i_,w_i_]-best_len_r_grid[r_i,w_i])) < 0.0
                println_sameline("#$(rw_iters) - new_R and new_W from the triangle in two quadrants best_R_W [$r_i_,$w_i_] and [$r_i,$w_i]")
                new_R1 = R - len_r*(R-best_R_grid[r_i_,w_i_]) /(len_r-best_len_r_grid[r_i_,w_i_])
                new_R4 = R - len_r*(R-best_R_grid[r_i,w_i]) /(len_r-best_len_r_grid[r_i,w_i])

                new_W1 = W - len_w*(W-best_W_grid[r_i_,w_i_]) /(len_w-best_len_w_grid[r_i_,w_i_])
                new_W4 = W - len_w*(W-best_W_grid[r_i,w_i]) /(len_w-best_len_w_grid[r_i,w_i])

                k_1 = (len_w-best_len_w_grid[r_i_,w_i_])/(len_r-best_len_r_grid[r_i_,w_i_])
                b_1 = len_w - len_r*k_1

                best_len_r_star = best_len_r_grid[r_i,w_i]
                best_len_w_star = best_len_w_grid[r_i,w_i]
                best_new_R_star = new_R4
                best_new_W_star = new_W4

                k_2 = best_len_w_star/best_len_r_star
                b_2 = 0.0
                x_starstar = (b_2-b_1)/(k_1-k_2)
                y_starstar = k_1*x_starstar + b_1
                d_1 = sqrt(x_starstar^2+y_starstar^2)
                d_2 = sqrt(best_len_r_star^2+best_len_w_star^2)
                relax = 1.0-d_1/(d_1+d_2)

                new_R = new_R1*relax + best_new_R_star*(1-relax)
                new_W = new_W1*relax + best_new_W_star*(1-relax)

            elseif best_total_len_grid[r_i_,w_i_] != Inf
                println_sameline("#$(rw_iters) - new_R and new_W from the opposite quadrant best_R_W_grid[$(r_i_),$(w_i_)]")
                new_R = R - len_r*(R-best_R_grid[r_i_,w_i_]) /(len_r-best_len_r_grid[r_i_,w_i_])
                #best_new_R = R - len_r*(R-best_R)/(len_r-best_len_r)
                #new_R = (new_R+best_new_R)/2.0

                new_W = W - len_w*(W-best_W_grid[r_i_,w_i_]) /(len_w-best_len_w_grid[r_i_,w_i_])
                #best_new_W = W - len_w*(W-best_W)/(len_w-best_len_w)
                #new_W = (new_W+best_new_W)/2.0
            elseif best_total_len_grid[r_i,w_i_] != Inf && best_total_len_grid[r_i_,w_i] != Inf
                println_sameline("#$(rw_iters) - new_R and new_W from the combination of adjacent quadrants best_R_W_grid [$r_i,$w_i_] and [$r_i_,$w_i]")
                new_R1 = R - len_r*(R-best_R_grid[r_i,w_i_]) /(len_r-best_len_r_grid[r_i,w_i_])
                new_R2 = R - len_r*(R-best_R_grid[r_i_,w_i]) /(len_r-best_len_r_grid[r_i_,w_i])
                #new_R = (new_R1+new_R2)/2.0
                #best_new_R = R - len_r*(R-best_R)/(len_r-best_len_r)
                #new_R = (new_R+best_new_R)/2.0

                new_W1 = W - len_w*(W-best_W_grid[r_i,w_i_]) /(len_w-best_len_w_grid[r_i,w_i_])
                new_W2 = W - len_w*(W-best_W_grid[r_i_,w_i]) /(len_w-best_len_w_grid[r_i_,w_i])
                #new_W = (new_W1+new_W2)/2.0
                #best_new_W = W - len_w*(W-best_W)/(len_w-best_len_w)
                #new_W = (new_W+best_new_W)/2.0

                d1 = sqrt(best_R_grid[r_i,w_i_]^2 + best_W_grid[r_i,w_i_]^2)
                d2 = sqrt(best_R_grid[r_i_,w_i]^2 + best_W_grid[r_i_,w_i]^2)
                relax = 1.0-d1/(d1+d2)
                new_R = relax*new_R1 + (1.0-relax)*new_R2
                new_W = relax*new_W1 + (1.0-relax)*new_W2
            elseif best_total_len_grid[r_i,w_i] != Inf
                println_sameline("#$(rw_iters) - new_R and new_W from the same quadrant best_R_W_grid[$(r_i),$(w_i)]")
                _new_R = R - len_r*(R-best_R_grid[r_i,w_i]) /(len_r-best_len_r_grid[r_i,w_i])
                best_new_R = R - len_r*(R-best_R)/(len_r-best_len_r)
                new_R = (_new_R+best_new_R)/2.0

                _new_W = W - len_w*(W-best_W_grid[r_i,w_i]) /(len_w-best_len_w_grid[r_i,w_i])
                best_new_W = W - len_w*(W-best_W)/(len_w-best_len_w)
                new_W = (_new_W+best_new_W)/2.0
            else
                println_sameline("#$(rw_iters) - new_R and new_W from the combination of best_R_W and old_R_W")
                _new_R     = R - len_r*(R-old_R) /(len_r-old_len_r)
                best_new_R = R - len_r*(R-best_R)/(len_r-best_len_r)
                new_R = (_new_R+best_new_R)/2.0

                _new_W     = W - len_w*(W-old_W) /(len_w-old_len_w)
                best_new_W = W - len_w*(W-best_W)/(len_w-best_len_w)
                new_W = (_new_W+best_new_W)/2.0
            end
            _new_R     = R - len_r*(R-old_R) /(len_r-old_len_r)
            new_R = _new_R*0.25+new_R*(1.0-0.25)

            _new_W     = W - len_w*(W-old_W) /(len_w-old_len_w)
            new_W = _new_W*0.25+new_W*(1.0-0.25)

            if new_R <= r_min
                new_R = 0.9*r_min + (1.0-0.9)*best_R
            elseif new_R >= r_max
                new_R = 0.9*r_max + (1.0-0.9)*best_R
            end

            if new_W <= w_min
                new_W = 0.9*w_min + (1.0-0.9)*best_W
            elseif new_W >= w_max
                new_W = 0.9*w_max + (1.0-0.9)*best_W
            end

            println_sameline("#$(rw_iters) - r:$(round(R*100;digits=4))%, w:$(round(W;digits=6)) - total_len:$(round(total_len;digits=6)) - len_r:$(round(len_r;digits=6)) - len_w:$(round(len_w;digits=6)) - new_r:$(round(new_R*100;digits=4))%, new_w:$(round(new_W;digits=6))")
            draw_best_grid(best_total_len_grid, best_len_r_grid,best_len_w_grid, best_R_grid,best_W_grid, R,W, len_r,len_w, new_R,new_W)
            if best_total_len_grid[r_i,w_i] > total_len
                best_total_len_grid[r_i,w_i] = total_len
                best_len_r_grid[r_i,w_i] = len_r
                best_len_w_grid[r_i,w_i] = len_w
                best_R_grid[r_i,w_i] = R
                best_W_grid[r_i,w_i] = W
                println_sameline("#$(rw_iters) - New best_R_grid[$(r_i)] and best_W_grid[$(w_i)]")
            end

            if total_len < best_total_len
                best_total_len = total_len
                best_len_r = len_r
                best_len_w = len_w
                best_res = copy(res)
                best_R = R
                best_W = W
                println_sameline("#$(rw_iters) - New best R and W")
            end

        catch e
            println_sameline(("#$(rw_iters) - half the R and W",e))

            #new_R = (old_R + new_R)/2.0
            new_R = (old_R + new_R + best_R)/3.0
            new_R = min(max(r_min, new_R), r_max)
            #new_W = (old_W + new_W)/2.0
            new_W = (old_W + new_W + best_W)/3.0
            new_W = min(max(w_min, new_W), w_max)

            total_len = old_total_len

            println_sameline("#$(rw_iters) - r:$(round(R*100;digits=4))%, w:$(round(W;digits=6)) - total_len:-.------ - len_r:-.------ - len_w:-.------ - new_r:$(round(new_R*100;digits=4))%, new_w:$(round(new_W;digits=6))")
            draw_best_grid(best_total_len_grid, best_len_r_grid,best_len_w_grid, best_R_grid,best_W_grid, R,W, len_r,len_w, new_R,new_W)
            R = old_R
            old_R = old_old_R
            old_len_r = old_old_len_r
            W = old_W
            old_W = old_old_W
            old_len_w = old_old_len_w

            #rw_iters -= 1

        end

        if max(abs(new_R-R), abs(new_W-W)) < gen_tol_x
            #=
            new_R = R + len_r
            new_R = min(max(r_min, new_R), r_max)
            new_W = W + len_w
            new_W = min(max(w_min, new_W), w_max)

            println_sameline("#$(rw_iters) - r:$(round(R*100;digits=4))%, w:$(round(W;digits=6)) - total_len:$(round(total_len;digits=6)) - len_r:$(round(len_r;digits=6)) - len_w:$(round(len_w;digits=6)) - new_r:$(round(new_R*100;digits=4))%, new_w:$(round(new_W;digits=6))")
            =#
            println_sameline("Factor prices have convereged, but markets are not clear")
            rw_iters = rw_maxiters
        end

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

    if total_len < best_total_len #rw_iters < rw_maxiters
        println("The interest rate - $(R*100)% and the wage - $(W) with total_len - $(total_len)")
        return [res, R, W, approx_object, model_params]
    else
        println("The interest rate - $(best_R*100)% and the wage - $(best_W) with total_len - $(best_total_len)")
        return [best_res, best_R, best_W, approx_object, model_params]
    end
end
