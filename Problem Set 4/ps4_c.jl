using Plots 
using Optim 
using Interpolations
using ForwardDiff
using Alert
using Roots
using IntervalRootFinding

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

# define utility function
# function util(x)

#     lhs = l -> (z*x[1]^alpha*l^(1-alpha) + (1-δ)*x[1] - x[2])*(1-alpha)*z*(x[1]^alpha)*l^(-alpha) - χ*l^η

#     l_curr = find_zero(lhs, (0, 1))

#     c = z*(x[1]^alpha)*(l_curr^(1-alpha)) - x[2] + x[1]*(1-δ)

#     u_val = ((c^(1-sig))/(1-sig)) - χ*((l_curr^(1+η))/(1+η))

#     if c <= 10^(-6)
#         u_val = u_val - 4*(c - 10^(-6))^2 
#     end

#     return u_val
# end

# define utility function
# function util(x)

#     lhs = l -> (z*x[1]^alpha*l^(1-alpha) + (1-δ)*x[1] - x[2])*(1-alpha)*z*(x[1]^alpha)*l^(-alpha) - χ*l^η # wrong 

#     l_curr = find_zero(lhs, (0, 1))

#     c = z*(x[1]^alpha)*(l_curr^(1-alpha)) - x[2] + x[1]*(1-δ)

#     u_val = ((c^(1-sig))/(1-sig)) - χ*((l_curr^(1+η))/(1+η))

#     # if c <= 10^(-6)
#     #     u_val = u_val - 4*(c - 10^(-6))^2 
#     # end

#     return u_val
# end

# v_itp = interpolate(ones(grid_size), BSpline(Cubic(Line(OnGrid()))))
# v_sitp = Interpolations.scale(v_itp, k_grid)

# error("hi")

function FOC(x)

    F = ones(2)

    # implied c 
    c = z*(x[1]^alpha)*(x[2]^(1-alpha)) - x[3] + x[1]*(1-δ)
    # captial 
    F[1] = c^(-sig)*(-1) + β*x[4]^2
    # labour 
    F[2] = (c^(-sig)*(1-alpha)*z*(x[1]^alpha)*x[2]^(-alpha) - χ*x[2]^η)^2

    pen = 0
    if c < 10^(-5)
        pen = 400*(c-10^(-5))^2
    end

    return F[1] + F[2] + pen 

end

inner_optimizer = NelderMead()

function Tv(v_old1, k_grid1)

    #v_sitp = linear_interpolation(k_grid1, v_old1)
    v_itp = interpolate(v_old1, BSpline(Cubic(Line(OnGrid()))))
    v_sitp = Interpolations.scale(v_itp, k_grid)

    v_new = zeros(grid_size);
    g_k = zeros(grid_size);
    g_l = zeros(grid_size);
   
    for i in 1:grid_size
        #print(i,"\n")
        foc_curr = x2 -> FOC([k_grid1[i], x2[1], x2[2], Interpolations.gradient(v_sitp, x2[2])[1]])
    

        k_max = k_max = max(min(z*(k_grid[i]^alpha) + (1-δ)*k_grid[i] - 10^(-6), k_high), k_ss) 

        sol_curr = optimize(foc_curr, [k_low, 0], [k_max, 1], [0.5*k_low+0.5*k_max, 0.4], Fminbox(inner_optimizer))
        
        v_new[i] = -Optim.minimum(sol_curr)
        g_k[i] = Optim.minimizer(sol_curr)[1]
        g_l[i] = Optim.minimizer(sol_curr)[2]

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