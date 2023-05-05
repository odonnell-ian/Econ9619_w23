using Plots 
using Optim 
using Interpolations
using ForwardDiff
using Alert
using Roots

# Set the initial parameters 
alpha = 1/3;
z = 1;
sig = 2; 
η = 1;
β = 0.98;
δ = 0.05;
χ = 8;

# Steady states 
r_ss = (1/β) - 1 + δ;
D = (r_ss^((η+alpha)/(1-alpha)))*(1/(z*alpha))^((η+1)/(1-alpha))*(alpha/(1 - alpha));
k_ss = (χ^(-1/sig)*D^(-1/sig)*((r_ss/alpha) - δ)^(-1))^(sig/(sig+η));
l_ss = (r_ss/(z*alpha))^(1/(1-alpha))*k_ss;
y_ss = z*k_ss^alpha*l_ss^(1-alpha);
c_ss = y_ss - k_ss + (1-δ)*k_ss;
w_ss = z*(1-alpha)*k_ss^alpha*l_ss^(-alpha);

# Problem 1, Part a --------------------------------------------------------
k_low = 10^(-5);
k_high = 2*k_ss;
grid_size = 50; 
max_iter = 200; 
#l_len = 30;

# = range(10^(-6), 1, length = l_len);
k_grid = range(k_low, k_high, length = grid_size);


tol = 10^-6;

# lhs = l -> (z*(k_ss^alpha)*(l^(1-alpha)) + (1-δ)*k_ss - k_ss)^(-sig)*(1-alpha)*z*(k_ss^alpha)*l^(-alpha) - χ*l^η

# lhs(l_ss)
# error("hi")
# define utility function
function util(x)

    lhs = l -> ((z*(x[1]^alpha)*(l^(1-alpha)) + (1-δ)*x[1] - x[2])^(-sig)*(1-alpha)*z*(x[1]^alpha)*l^(-alpha) - χ*l^η)^2

    opt_aux = optimize(lhs, 0, 1, Brent())

    l_curr = Optim.minimizer(opt_aux)

    c = z*(x[1]^alpha)*(l_curr^(1-alpha)) - x[2] + x[1]*(1-δ)

    u_val = ((c^(1-sig))/(1-sig)) - χ*((l_curr^(1+η))/(1+η)) 

    if c < 10^(-6)
        u_val = u_val - 100*(c - 10^(-6))^2 
    end

    return u_val
end


inner_optimizer = NelderMead()

function Tv(v_old1, k_grid1)

    v_sitp = linear_interpolation(k_grid1, v_old1)

    v_new = zeros(grid_size);
    g_k = zeros(grid_size);
    g_l = zeros(grid_size);
   
    for i in 1:grid_size
        #print(i,"\n")
        u_curr = x2 -> -util([k_grid1[i], x2]) - β * v_sitp(x2)


        k_max = max(min(z*(k_grid[i]^alpha) + (1-δ)*k_grid[i] - 10^(-6), k_high), k_ss) 
        #k_max = k_high

        sol_curr = optimize(u_curr, k_low, k_max, Brent())
        
        
        v_new[i] = -Optim.minimum(sol_curr)
        g_k[i] = Optim.minimizer(sol_curr)


    end
    return v_new, g_k, g_l
end

# define function to iterate over the value function 
function vfi1()
    #v_old = zeros(grid_size); # try chainging this
    v_old = ones(grid_size)/50 
    v_curr = zeros(grid_size);
    g_k1 = ones(grid_size)/20;
    g_l1 = ones(grid_size)/20;
    v_dist = 1;
    count = 0;
    while count <= max_iter && v_dist > tol 
        v_curr, g_k1, g_l1 = Tv(v_old, k_grid);
        v_dist = maximum(abs.(v_curr./v_old.-1));
        count += 1; 
        v_old = v_curr; 
    end

    return v_curr, g_k1, g_l1, count, v_dist
end

# record time and number of iterations 
@time x = vfi1()

p6 = plot(k_grid, z*(k_grid.^alpha).*(x[3].^(1-alpha)) - x[2] + (1-δ)*k_grid)