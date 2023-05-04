using Plots 

# primatives ----------------------------------------------------------
u1(c) = log(c);
u2(c) = sqrt(c);
u3(c) = (c^(-1))/(-1);
u4(c) = (c^(-4))/(-4);
u5(c) = (c^(-9))/(-9);

num_data = 5;
num_fin_grid = 100;

dom = range(0.05,2,length = num_fin_grid);
dat_dom = range(0.05, 2, length = num_data);

u1_dat = map(u1, dat_dom); 
u2_dat = map(u2, dat_dom); 
u3_dat = map(u3, dat_dom); 
u4_dat = map(u4, dat_dom); 
u5_dat = map(u5, dat_dom); 

u1_true = map(u1, dom); 
u2_true = map(u2, dom); 
u3_true = map(u3, dom); 
u4_true = map(u4, dom); 
u5_true = map(u5, dom); 

# inputs: nodes, desired grid

# Linear Interpolation 
# note that the nodes include the top and bottom of the domain
function lin_int(nodes_x, nodes_y, grid1)
    n = length(grid1);
    int_fun = zeros(n, 1);
    count = 1;
    for i in grid1
        ind = findmax(sign.(nodes_x .- i))[2] - 1;
        Ax = (nodes_x[ind + 1] - i)/(nodes_x[ind + 1] - nodes_x[ind]);
        int_fun[count] = Ax*nodes_y[ind] + (1 - Ax)*nodes_y[ind + 1];
        count += 1;
    end
    return int_fun
end

# Natural spline
function cubic_spline(nodes_x, nodes_y, grid1)
    n_y = length(nodes_y)
    # find second derivatives 
    # create matrix
    # go row by row 
    c_mat = zeros(n_y-2, n_y-2)
    s_vec = zeros(n_y-2, 1)

    s_aux = (nodes_y[2] - nodes_y[1])/(nodes_x[2] - nodes_x[1])

    for i = 1:(n_y-3)
        c_mat[i, i] = (nodes_x[i+2] - nodes_x[i])/3
        c_mat[i, i+1] = (nodes_x[i+2] - nodes_x[i+1])/6
        c_mat[i+1, i] = c_mat[i, i+1]
        s_vec[i] = (nodes_y[i+2] - nodes_y[i+1])/(nodes_x[i+2] - nodes_x[i+1]) - s_aux
        s_aux = s_vec[i] + s_aux
    end

    c_mat[n_y-2, n_y-2] = (nodes_x[n_y] - nodes_x[n_y-2])/3
    s_vec[n_y-2] = (nodes_y[n_y] - nodes_y[n_y-1])/(nodes_x[n_y] - nodes_x[n_y-1]) - s_aux

    # solve system 
    y_pp = c_mat\s_vec
    y_pp = vcat(0, y_pp)
    y_pp = vcat(y_pp, 0)
    n = length(grid1);
    int_fun = zeros(n, 1);
    count1 = 1

    for j in grid1

        # locate closest on x grid 
        ind = findmax(sign.(nodes_x .- j))[2] - 1

        Ax = (nodes_x[ind + 1] - j)/(nodes_x[ind + 1] - nodes_x[ind])
        Bx = 1 - Ax
        Cx = (1/6)*(Ax^3 - Ax)*(nodes_x[ind+1] - nodes_x[ind])^2
        Dx = (1/6)*(Bx^3 - Bx)*(nodes_x[ind+1] - nodes_x[ind])^2

        int_fun[count1] = Ax*nodes_y[ind] + Bx*nodes_y[ind+1] + Cx*y_pp[ind] + Dx*y_pp[ind+1] 
        count1 = count1 + 1 
    end
    
    return int_fun
end

# global poly with a newton basis 


# Iterpolate 
inter1 = lin_int(dat_dom, u1_dat, dom)
inter12 = cubic_spline(dat_dom, u1_dat, dom)
p1 = plot(dom, [u1_true inter1 inter12])


inter2 = lin_int(dat_dom, u2_dat, dom)
inter22 = cubic_spline(dat_dom, u2_dat, dom)
p2 = plot(dom, [u2_true inter2 inter22])
