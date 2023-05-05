using Plots 
using Optim 
using Interpolations
using ForwardDiff
using Alert

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
grid_size = 30; 
max_iter = 500; 
#l_len = 30;



# = range(10^(-6), 1, length = l_len);
k_grid = range(k_low, k_high, length = grid_size);
new_grid = 

tol = 10^-6;

# define utility function
function util(x)

    c = z*(x[1]^alpha)*(x[2]^(1-alpha)) - x[3] + x[1]*(1-δ)

    u_val = ((c^(1-sig))/(1-sig)) - χ*((x[2]^(1+η))/(1+η))

    if c < 10^(-6)
        u_val = u_val - (10^10)*(c - 10^(-6))^2 
    end

    return u_val
end



# A_x = 1.:2.:40.
# A = [log(x) for x in A_x]
# itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
# sitp = scale(itp, A_x)
# sitp(3.) # exactly log(3.)
# sitp(3.5) # approximately log(3.5)

# u_curr = x2 -> -util([k_grid[1], x2[2], x2[1]]) - β * sitp(x2[1])
# t3 = u_curr([3,0.1])
# t6 = ForwardDiff.gradient(u_curr,[5,0.1])
# t1 = optimize(u_curr, [3,0], [10, 1], [5,0.1])

#v_itp = interpolate(ones(grid_size), BSpline(Cubic(Line(OnGrid()))))
#v_sitp = Interpolations.scale(v_itp, k_grid)

#v_lin = linear_interpolation(k_grid, ones(grid_size))
#v_lin(2.3)
# u_curr = x2 -> -util([k_grid[1], x2[2], x2[1]]) - β * v_sitp(x2[1])
#u_curr = x2 -> -util([k_grid[1], x2[2], x2[1]]) - β * v_lin(x2[1])

#k_max = max(min(z*(k_grid[1]^alpha) + (1-δ)*k_grid[1] - 10^(-6), k_high), k_ss)
inner_optimizer = NelderMead()
#sol_curr = optimize(u_curr, [k_low, 0], [k_max, 1], [0.5*k_low+0.5*k_max, l_ss], Fminbox(inner_optimizer)) 
#sol_curr = optimize(u_curr, [k_low, 0], [k_max, 1], [0.5*k_low+0.5*k_max, l_ss]) 
#alert("check here\n")
#error("hello")
# define function get new value function 
function Tv(v_old1, k_grid1)

    #v_itp = interpolate(v_old1, BSpline(Cubic(Line(OnGrid()))))
    v_sitp = linear_interpolation(k_grid1, v_old1)
    #v_sitp = Interpolations.scale(v_itp, k_grid1)

    v_new = zeros(grid_size);
    g_k = zeros(grid_size);
    g_l = zeros(grid_size);

    for i in 1:grid_size

        u_curr = x2 -> (-1)*(util([k_grid1[i], x2[2], x2[1]]) + β * v_sitp(x2[1]))

        #k_min = k_low
        #k_max = max(min(z*(k_grid[i]^alpha) + (1-δ)*k_grid[i] - 10^(-6), k_high), k_ss) 
        k_max = k_high

        # not sure how to check for derivative over all of the boundary area 
        #t2 = u_curr([k_ss, 0.5])
        #print(t2, "\n")

        sol_curr = optimize(u_curr, [k_low, 0], [k_max, 1], [0.5*k_low+0.5*k_max, l_ss], Fminbox(inner_optimizer))
        
        
        v_new[i] = -Optim.minimum(sol_curr)
        g_k[i] = Optim.minimizer(sol_curr)[1]
        g_l[i] = Optim.minimizer(sol_curr)[2]

        # v_aux = zeros(grid_size, l_len);
        # for j in 1:grid_size
        #     for k in 1:l_len 
        #         c_curr = max(0, z*(k_grid1[i]^alpha)*(l_grid[k]^(1-alpha)) - k_grid1[j] + k_grid1[i]*(1-δ));
        #         #c_curr = z*(k_grid1[i]^alpha)*(l_grid[k]^(1-alpha)) - k_grid1[j] + k_grid1[i]*(1-δ);
        #         v_aux[j,k] = util(c_curr, l_grid[k]) + β*v_old1[j];
        #     end
        # end
        # v_new[i], g1 = findmax(v_aux);
        # g_k[i] = g1[1];
        # g_l[i] = g1[2];
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
        v_old = 0.95*v_old + 0.05*v_curr; 
    end

    return v_curr, g_k1, g_l1, count, v_dist
end

# record time and number of iterations 
@time x = vfi1()

# plot value function 
p1 = plot(k_grid, x[1]);
# plot policy functions 
p2 = plot(k_grid,  [x[2] k_grid]);
p3 = plot(k_grid, x[3]);

# Calculate euler errors and plot 
# capital first 
k_for = x[2];
k_itp = linear_interpolation(k_grid, x[2])
k_for2 = k_itp(x[2])
l_curr = x[3]
l_itp = linear_interpolation(k_grid, x[3])
l_for = l_itp(x[3])

c_plan = z.*k_grid.^alpha.*(l_curr.^(1-alpha)) - k_for + (1-δ).*k_grid;
rhs_k = (z.*k_grid.^alpha.*(l_curr.^(1-alpha)) - k_for + (1-δ).*k_grid).^(-sig);
lhs1_k = (z.*k_for.^alpha.*(l_for.^(1-alpha)) - k_for2 + (1-δ).*k_for).^(-sig); # check here if things work
lhs2_k = ((alpha*z).*k_for.^(alpha - 1).*(l_for.^(1-alpha)) .+ (1 - δ));

l_err = (rhs_k.*((z*(1-alpha))*k_grid.^(alpha).*l_curr.^(-alpha)))./(χ * l_curr.^η).-1;

resid = β*(lhs1_k .* lhs2_k)./rhs_k .- 1;

p4 = scatter(k_grid, resid);
hline!([0.01 -0.01]);
savefig(p4, "resid.png")

p6 = plot(k_grid, z*(k_grid.^alpha).*(x[3].^(1-alpha)) - x[2] + (1-δ)*k_grid)
