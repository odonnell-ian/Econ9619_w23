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

# inputs: nodes, desired grid

# Linear Interpolation 
# note that the nodes include the top and bottom of the domain
function lin_int(nodes_x, nodes_y, grid1)
    n = length(grid1);
    int_fun = zeros(n, 1);
    count = 1;
    for i in grid1
        ind = findmax(sign.(nodes_x .- i))[2] - 1;
        Ax = (nodes_x[ind + 1] - grid1[ind])/(nodes_x[ind + 1] - nodes_x[ind]);
        int_fun[count] = Ax*nodes_y[ind] + (1 - Ax)*nodes_y[ind + 1];
        count += 1;
    end
    return int_fun
end

# Cubic spline


# Other spline 