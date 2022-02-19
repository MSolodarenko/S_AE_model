using ProgressMeter
#=
filename = "Results/number_a_nodes_grid_results.txt"
# read the file to lines
file_lines = open(filename) do f
    lines = readlines(f)
    (lines)
end

num_of_records = length(file_lines)-1
num_a_nodes = zeros(num_of_records)
time_length = zeros(num_of_records)
len_r = zeros(num_of_records)
len_w = zeros(num_of_records)
cap_out = zeros(num_of_records)
cre_out = zeros(num_of_records)
occ_w = zeros(num_of_records)
occ_se = zeros(num_of_records)
ooc_emp = zeros(num_of_records)
@showprogress for i = 1:num_of_records
    line_i = i+1
    #display(file_lines[7])
    line = split(file_lines[line_i], ",")
    try
        num_a_nodes[i] = parse(Int64,split(line[1], "(")[2])
        time_length[i] = parse(Float64,line[2])
        len_r[i] = parse(Float64,line[3])
        len_w[i] = parse(Float64,line[4])
        cap_out[i] = round(parse(Float64,line[5]); digits=2)
        cre_out[i] = round(parse(Float64,line[6]); digits=2)
        occ_w[i] = round(parse(Float64,line[7]); digits=2)
        occ_se[i] = round(parse(Float64,line[8]); digits=2)
        ooc_emp[i] = round(parse(Float64,split(line[9],")")[1]); digits=2)
    catch e
        display(i)
        display(line)
        throw(e)
    end
end

dif_num_a_nodes = diff(num_a_nodes)
dif_time_length = diff(time_length)
dif_len_r = diff(len_r)
dif_len_w = diff(len_w)
dif_cap_out = diff(cap_out)
dif_cre_out = diff(cre_out)
dif_occ_w = diff(occ_w)
dif_occ_se = diff(occ_se)
dif_ooc_emp = diff(ooc_emp)

num_a_nodes[findmin(abs,dif_cap_out)[2]]
num_a_nodes[findmin(abs,dif_cre_out)[2]]
num_a_nodes[findmin(abs,dif_occ_w)[2]]
num_a_nodes[findmin(abs,dif_occ_se)[2]]
num_a_nodes[findmin(abs,dif_ooc_emp)[2]]

ind = findmin(abs,abs.(dif_cap_out).+abs.(dif_cre_out).+abs.(dif_occ_w).+abs.(dif_occ_se).+abs.(dif_ooc_emp))[2]
display([ind,num_a_nodes[ind],time_length[ind],len_r[ind],len_w[ind],cap_out[ind],cre_out[ind],occ_w[ind],occ_se[ind],ooc_emp[ind]])
display([ind,dif_num_a_nodes[ind],dif_time_length[ind],dif_len_r[ind],dif_len_w[ind],dif_cap_out[ind],dif_cre_out[ind],dif_occ_w[ind],dif_occ_se[ind],dif_ooc_emp[ind]])
=#

filename = "Results/number_a_nodes_ss_grid_results.txt"
# read the file to lines
file_lines = open(filename) do f
    lines = readlines(f)
    (lines)
end

num_of_records = length(file_lines)-1
num_a_nodes = zeros(num_of_records)
time_length = zeros(num_of_records)
r = zeros(num_of_records)
w = zeros(num_of_records)
len_r = zeros(num_of_records)
len_w = zeros(num_of_records)
cap_out = zeros(num_of_records)
cre_out = zeros(num_of_records)
occ_w = zeros(num_of_records)
occ_se = zeros(num_of_records)
ooc_emp = zeros(num_of_records)
@showprogress for i = 1:num_of_records
    line_i = i+1
    #display(file_lines[7])
    line = split(file_lines[line_i], ",")
    try
        num_a_nodes[i] = parse(Int64,split(line[1], "(")[2])
        time_length[i] = parse(Float64,line[2])
        r[i] = round(parse(Float64,line[3]); digits=6)
        w[i] = round(parse(Float64,line[4]); digits=6)
        len_r[i] = parse(Float64,line[5])
        len_w[i] = parse(Float64,line[6])
        cap_out[i] = round(parse(Float64,line[7]); digits=2)
        cre_out[i] = round(parse(Float64,line[8]); digits=2)
        occ_w[i] = round(parse(Float64,line[9]); digits=2)
        occ_se[i] = round(parse(Float64,line[10]); digits=2)
        ooc_emp[i] = round(parse(Float64,split(line[11],")")[1]); digits=2)
    catch e
        display(i)
        display(line)
        throw(e)
    end
end

dif_num_a_nodes = diff(num_a_nodes)
dif_time_length = diff(time_length)
dif_r = diff(r)
dif_w = diff(w)
dif_len_r = diff(len_r)
dif_len_w = diff(len_w)
dif_cap_out = diff(cap_out)
dif_cre_out = diff(cre_out)
dif_occ_w = diff(occ_w)
dif_occ_se = diff(occ_se)
dif_ooc_emp = diff(ooc_emp)

num_a_nodes[findmin(abs,dif_r)[2]]
num_a_nodes[findmin(abs,dif_w)[2]]
num_a_nodes[findmin(abs,dif_cap_out)[2]]
num_a_nodes[findmin(abs,dif_cre_out)[2]]
num_a_nodes[findmin(abs,dif_occ_w)[2]]
num_a_nodes[findmin(abs,dif_occ_se)[2]]
num_a_nodes[findmin(abs,dif_ooc_emp)[2]]

ind = findmin(abs,abs.(dif_cap_out).+abs.(dif_cre_out).+abs.(dif_occ_w).+abs.(dif_occ_se).+abs.(dif_ooc_emp))[2]
display([ind,num_a_nodes[ind],time_length[ind],len_r[ind],len_w[ind],cap_out[ind],cre_out[ind],occ_w[ind],occ_se[ind],ooc_emp[ind]])
display([ind,dif_num_a_nodes[ind],dif_time_length[ind],dif_len_r[ind],dif_len_w[ind],dif_cap_out[ind],dif_cre_out[ind],dif_occ_w[ind],dif_occ_se[ind],dif_ooc_emp[ind]])

temp_res = abs.(dif_cap_out).+abs.(dif_cre_out).+abs.(dif_occ_w).+abs.(dif_occ_se).+abs.(dif_ooc_emp)
