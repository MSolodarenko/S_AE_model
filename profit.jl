#########################
# production function
function production(tm::Float64,zm::Float64,k::Float64,zw::Float64,nd::Float64, gamma, eta, theta)
    return (tm*zm)^gamma * k^eta * ((1-tm)*zw + nd)^theta
end
# function that returns profit of entrepreneurial endeavour
function profit(zm::Float64,zw::Float64,k::Float64,nd::Float64,tm::Float64,r::Float64,w::Float64, delta, gamma, eta, theta)
    return production(tm,zm,k,zw,nd, gamma, eta, theta) - w*nd - (r+delta)*k
end

# unconstrained optimal choice of
# capital demand, total labor demand, outside labour demand and managerial time and profit
#       for self-employed
function self_optimal(zm::Float64,zw::Float64,r::Float64, w::Float64, delta, gamma, eta, theta)
    tm = gamma/(theta+gamma)
    k_ent = ((zm*tm)^gamma)*eta/(r+delta)
    l = (1-tm)*zw
    k = (k_ent * l^theta)^(1/(1-eta))
    y = profit(zm,zw,k,0.0,tm,r,w, delta, gamma, eta, theta)
    return [k, l, 0.0, tm, y]
end

# unconstrained optimal choice of
# capital demand, total labor demand, outside labour demand and managerial time and profit
#       for employers
function entrep_optimal(zm::Float64,zw::Float64,r::Float64, w::Float64, delta, gamma, eta, theta, c_e)
    k_ent= zm^gamma * eta/(r+delta)
    k = ((r+delta)*theta/(w*eta))^theta
    k = (k_ent*k)^(1.0/(1.0-theta-eta))
    nd = k * (r+delta)*theta/(w*eta)
    y = profit(zm,zw,k,nd,1.0,r,w, delta, gamma, eta, theta) - c_e
    return [k, nd, nd, 1.0, y]
end

# constrained by collateral optimal choice of
# capital demand, total labor demand, outside labour demand and managerial time and profit
#       for self-employed
function self_constrained(a::Float64,zm::Float64,zw::Float64,r::Float64,w::Float64, lambda, delta, gamma, eta, theta)
    tm = gamma/(theta+gamma)
    k = lambda*a
    l = (1-tm)*zw
    y = profit(zm,zw,k,0.0,tm,r,w, delta, gamma, eta, theta)
    return [k, l, 0.0, tm, y]
end

# constrained by collateral optimal choice of
# capital demand, total labor demand, outside labour demand and managerial time and profit
#       for employers
function entrep_constrained(a::Float64,zm::Float64,zw::Float64,r::Float64,w::Float64, lambda, delta, gamma, eta, theta, c_e)
    k = lambda*a
    tm = ( (zm^gamma * (k^eta) * theta/w) * ( (theta*zw/gamma)^(theta-1.0) ) )^(1.0/(1.0-theta-gamma))
	tm = min(max(0.0,tm),1.0)
	nd = ( (zm*tm)^gamma * (k^eta) * theta/w )^(1.0/(1.0-theta)) - ((1.0-tm)*zw)
	nd = max(nd,0.0)
    y = profit(zm,zw,k,nd,tm,r,w, delta, gamma, eta, theta) - c_e
    return [k, nd+(1.0-tm)*zw, nd, tm, y]
end

##########
# function that returns
#   occupation choice
#   income
#   earnings
#   capital excess demand from market
#   capital demand
#   credit
#   labour excess demand from market
#   labour demand
#   labour supply
#   output
#   cost_of_employing
function compute_income_profile(a::Float64,zm::Float64,zw::Float64,r::Float64,w::Float64, lambda, delta, gamma, eta, theta, c_e)

    se_opt = self_optimal(zm,zw,r,w, delta, gamma, eta, theta)
    se_const = self_constrained(a,zm,zw,r,w, lambda, delta, gamma, eta, theta)
    emp_opt = entrep_optimal(zm,zw,r,w, delta, gamma, eta, theta, c_e)
    emp_const = entrep_constrained(a,zm,zw,r,w, lambda, delta, gamma, eta, theta, c_e)

    income_w = w*zw + (1.0+r)*a #### !!! change (1+r)*a -> r*a !!! ####

    income_se = 0.0
    if se_opt[1] >= se_const[1] # k >= lambda*a
        income_se = se_const
    else# k < lambda*a
        income_se = se_opt
    end

    income_emp = 0.0
    if emp_opt[1] >= emp_const[1] # k >= lambda*a
        income_emp = emp_const
    else# k < lambda*a
        income_emp = emp_opt
    end

	occ_choice = 1
	if income_w >= max(income_se[5],income_emp[5])+(1+r)*a #### !!! change (1+r)*a -> r*a !!! ####
		occ_choice = 1
	elseif income_se[5] >= income_emp[5]
		occ_choice = 2
	else
		occ_choice = 3
	end

    if occ_choice == 1 #income_w > max(income_se[5],income_emp[5])+(1+r)*a
        # Workers
        #occ_choice = 1
        income = income_w
        credit = 0.0
        deposit = a
        earnings = w*zw
        capital_excess = -a
        capital_d = 0.0
        labour_excess = -zw
        labour_d = 0.0
        labour_s = zw
        output = 0.0

        cost_of_employing = 0.0
		managerial_input = 0.0
    elseif occ_choice == 2 #income_se[5] > income_emp[5]
        # Self-employed
        #occ_choice = 2
        income = income_se[5] + (1+r)*a #### !!! change (1+r)*a -> r*a !!! ####
        credit = max(0.0,income_se[1]-a)
        deposit = max(a-income_se[1],0.0)
        earnings = income_se[5]#income - a - r*deposit
        capital_excess = -a+income_se[1]
        capital_d = income_se[1]
        labour_excess = 0.0
        labour_d = income_se[2]
        labour_s = income_se[2]
        # income = k,l,nd,tm,y
        # production = tm,zm,k,zw,nd
        output = production(income_se[4],zm,income_se[1],zw,income_se[3], gamma, eta, theta)

        cost_of_employing = 0.0
		managerial_input = income_se[4]*zm
    else #occ_choice == 3
        # Employers
        #occ_choice = 3
        income = income_emp[5] + (1+r)*a #### !!! change (1+r)*a -> r*a !!! ####
        credit = max(0.0,income_emp[1]-a)
        deposit = max(a-income_emp[1],0.0)
        earnings = income_emp[5]#income - a - r*deposit#income_emp[5] + r*(a-deposit)#
        capital_excess = -a + income_emp[1]
        capital_d = income_emp[1]
        labour_excess = income_emp[3]
        labour_d = income_emp[2]
        labour_s = income_emp[2] - income_emp[3]
        # income = k,l,nd,tm,y
        # production = tm,zm,k,zw,nd
        output = production(income_emp[4],zm,income_emp[1],zw,income_emp[3], gamma, eta, theta)

        cost_of_employing = c_e
		managerial_input = income_emp[4]*zm
    end
    return occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input
end

function compute_income_profile(a_nodes::Array{Float64,1},number_a_nodes::Int64,r::Float64, w::Float64, number_zeta_nodes::Int64, number_alpha_m_nodes::Int64, number_alpha_w_nodes::Int64, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
    occ_choice = Array{Int32}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    income = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    earnings = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    capital_excess = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    capital_d = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    credit = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_excess = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_d = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_s = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    deposit = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    output = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    cost_of_employing = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	managerial_input = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    for a_i in 1:number_a_nodes
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            res = compute_income_profile(a_nodes[a_i],z_m_nodes[u_i,alpha_m_i,zeta_i],z_w_nodes[u_i,alpha_m_i,alpha_w_i],r,w, lambda, delta, gamma, eta, theta, c_e)

            occ_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[1]
            income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[2]
            earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[3]
            capital_excess[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[4]
            capital_d[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[5]
            credit[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[6]
            labour_excess[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[7]
            labour_d[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[8]
            labour_s[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[9]
            deposit[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[10]
            output[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[11]

            cost_of_employing[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[12]
			managerial_input[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[13]
        end
    end

    return occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input
end

function compute_income_and_earnings(a::Float64,zm::Float64,zw::Float64,r::Float64,w::Float64, lambda, delta, gamma, eta, theta, c_e)
	return compute_income_profile(a,zm,zw,r,w, lambda, delta, gamma, eta, theta, c_e)[2:3]
end

function compute_income_and_earnings(a_nodes::Array{Float64,1},number_a_nodes::Int64,r::Float64, w::Float64, number_zeta_nodes::Int64, number_alpha_m_nodes::Int64, number_alpha_w_nodes::Int64, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
    income = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    earnings = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    for a_i in 1:number_a_nodes
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            res = compute_income_and_earnings(a_nodes[a_i],z_m_nodes[u_i,alpha_m_i,zeta_i],z_w_nodes[u_i,alpha_m_i,alpha_w_i],r,w, lambda, delta, gamma, eta, theta, c_e)

            income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[1]
            earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[2]

        end
    end

    return income, earnings
end

# computation of income prodile for fixed occupation
function compute_income_profile(occ_choice::Int64,a::Float64,zm::Float64,zw::Float64,r::Float64,w::Float64, lambda, delta, gamma, eta, theta, c_e)

    se_opt = self_optimal(zm,zw,r,w, delta, gamma, eta, theta)
    se_const = self_constrained(a,zm,zw,r,w, lambda, delta, gamma, eta, theta)
    emp_opt = entrep_optimal(zm,zw,r,w, delta, gamma, eta, theta, c_e)
    emp_const = entrep_constrained(a,zm,zw,r,w, lambda, delta, gamma, eta, theta, c_e)

    income_w = w*zw + (1.0+r)*a

    income_se = 0.0
    if se_opt[1] >= se_const[1] # k >= lambda*a
        income_se = se_const
    else# k < lambda*a
        income_se = se_opt
    end

    income_emp = 0.0
    if emp_opt[1] >= emp_const[1] # k >= lambda*a
        income_emp = emp_const
    else# k < lambda*a
        income_emp = emp_opt
    end

    if occ_choice == 1 #income_w > max(income_se[5],income_emp[5])+(1+r)*a
        # Workers
        #occ_choice = 1
        income = income_w
        credit = 0.0
        deposit = a
        earnings = w*zw
        capital_excess = -a
        capital_d = 0.0
        labour_excess = -zw
        labour_d = 0.0
        labour_s = zw
        output = 0.0

        cost_of_employing = 0.0
		managerial_input = 0.0
    elseif occ_choice == 2 #income_se[5] > income_emp[5]
        # Self-employed
        #occ_choice = 2
        income = income_se[5] + (1+r)*a
        credit = max(0.0,income_se[1]-a)
        deposit = max(a-income_se[1],0.0)
        earnings = income_se[5]#income - a - r*deposit
        capital_excess = -a+income_se[1]
        capital_d = income_se[1]
        labour_excess = 0.0
        labour_d = income_se[2]
        labour_s = income_se[2]
        # income = k,l,nd,tm,y
        # production = tm,zm,k,zw,nd
        output = production(income_se[4],zm,income_se[1],zw,income_se[3], gamma, eta, theta)

        cost_of_employing = 0.0
		managerial_input = income_se[4]*zm
    elseif occ_choice == 3
        # Employers
        #occ_choice = 3
        income = income_emp[5] + (1+r)*a
        credit = max(0.0,income_emp[1]-a)
        deposit = max(a-income_emp[1],0.0)
        earnings = income_emp[5]#income - a - r*deposit#income_emp[5] + r*(a-deposit)#
        capital_excess = -a + income_emp[1]
        capital_d = income_emp[1]
        labour_excess = income_emp[3]
        labour_d = income_emp[2]
        labour_s = income_emp[2] - income_emp[3]
        # income = k,l,nd,tm,y
        # production = tm,zm,k,zw,nd
        output = production(income_emp[4],zm,income_emp[1],zw,income_emp[3], gamma, eta, theta)

        cost_of_employing = c_e
		managerial_input = income_emp[4]*zm
	else
		throw("occ_choice is out of bound")
    end
    return occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input
end

function compute_income_profile(occ::Int64,a_nodes::Array{Float64,1},number_a_nodes::Int64,r::Float64, w::Float64, number_zeta_nodes::Int64, number_alpha_m_nodes::Int64, number_alpha_w_nodes::Int64, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
    occ_choice = Array{Int32}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    income = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    earnings = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    capital_excess = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    capital_d = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    credit = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_excess = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_d = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_s = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    deposit = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    output = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    cost_of_employing = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	managerial_input = Array{Float64}(undef,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    for a_i in 1:number_a_nodes
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            res = compute_income_profile(occ,a_nodes[a_i],z_m_nodes[u_i,alpha_m_i,zeta_i],z_w_nodes[u_i,alpha_m_i,alpha_w_i],r,w, lambda, delta, gamma, eta, theta, c_e)

            occ_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[1]
            income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[2]
            earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[3]
            capital_excess[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[4]
            capital_d[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[5]
            credit[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[6]
            labour_excess[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[7]
            labour_d[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[8]
            labour_s[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[9]
            deposit[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[10]
            output[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[11]

            cost_of_employing[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[12]
			managerial_input[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = res[13]
        end
    end

    return occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input
end

function compute_income_profile_fixed_occ(a_nodes::Array{Float64,1},number_a_nodes::Int64,r::Float64, w::Float64, number_zeta_nodes::Int64, number_alpha_m_nodes::Int64, number_alpha_w_nodes::Int64, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
	occ_choice = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	income = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	earnings = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	capital_excess = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	capital_d = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	credit = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	labour_excess = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	labour_d = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	labour_s = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	deposit = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	output = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	cost_of_employing = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
	managerial_input = Array{Any}(undef, 3,number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

	Threads.@threads for occ = 1:3
		occ_choice[occ,:,:,:,:,:], income[occ,:,:,:,:,:], earnings[occ,:,:,:,:,:], capital_excess[occ,:,:,:,:,:], capital_d[occ,:,:,:,:,:], credit[occ,:,:,:,:,:], labour_excess[occ,:,:,:,:,:], labour_d[occ,:,:,:,:,:], labour_s[occ,:,:,:,:,:], deposit[occ,:,:,:,:,:], output[occ,:,:,:,:,:], cost_of_employing[occ,:,:,:,:,:], managerial_input[occ,:,:,:,:,:] = compute_income_profile(occ,a_nodes,number_a_nodes,r,w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
	end

	return occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input
end
