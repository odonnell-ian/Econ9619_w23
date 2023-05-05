using Plots 

# Set the initial parameters 
alpha = 1/3;
z = 1;
sig = 2; 
η = 1;
β = 0.95;
δ = 0.75;
χ = 34;

# Steady states 
r_ss = (1/β) - 1 + δ;
D = (r_ss^((η+alpha)/(1-alpha)))*(1/(z*alpha))^((η+1)/(1-alpha))*(alpha/(1 - alpha));
k_ss = (χ^(-1/sig)*D^(-1/sig)*((r_ss/alpha) - δ)^(-1))^(sig/(sig+η));
l_ss = (r_ss/(z*alpha))^(1/(1-alpha))*k_ss;
y_ss = z*k_ss^alpha*l_ss^(1-alpha);
c_ss = y_ss - k_ss + (1-δ)*k_ss;
w_ss = z*(1-alpha)*k_ss^alpha*l_ss^(-alpha);

# Problem 5, Part a --------------------------------------------------------
k_low = 10^(-5);
k_high = 2*k_ss;
grid_size = 50; 
max_iter = 300; 
l_len = 30;

l_grid = range(10^(-6), 1, length = l_len);
k_grid = range(k_low, k_high, length = grid_size);


tol = 10^-6;

# define utility function
function util(c1, l1)
    return ((c1^(1-sig))/(1-sig)) - χ*((l1^(1+η))/(1+η))
end
#plot(error("hello")
# define function get new value function 
function Tv(v_old1, k_grid1)
    v_new = zeros(grid_size);
    g_k = zeros(grid_size);
    g_l = zeros(grid_size);
    g1 = [0 0];

    for i in 1:grid_size
        v_aux = zeros(grid_size, l_len);
        for j in 1:grid_size
            for k in 1:l_len 
                c_curr = max(0, z*(k_grid1[i]^alpha)*(l_grid[k]^(1-alpha)) - k_grid1[j] + k_grid1[i]*(1-δ));
                #c_curr = z*(k_grid1[i]^alpha)*(l_grid[k]^(1-alpha)) - k_grid1[j] + k_grid1[i]*(1-δ);
                v_aux[j,k] = util(c_curr, l_grid[k]) + β*v_old1[j];
            end
        end
        v_new[i], g1 = findmax(v_aux);
        g_k[i] = g1[1];
        g_l[i] = g1[2];
    end
    return v_new, g_k, g_l
end

# define function to iterate over the value function 
function vfi1()
    v_old = zeros(grid_size);
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

# plot value function 
p1 = plot(k_grid, x[1]);
# plot policy functions 
p2 = plot(k_grid, [k_grid[Int.(x[2])] k_grid]);
p3 = plot(k_grid, l_grid[Int.(x[3])]);

# Calculate euler errors and plot 
# capital first 
k_for = k_grid[Int.(x[2])];
k_for2 = k_grid[Int.(x[2][Int.(x[2])])];
l_curr = l_grid[Int.(x[3])];
l_for = l_grid[Int.(x[3][Int.(x[2])])];

c_plan = z.*k_grid.^alpha.*(l_curr.^(1-alpha)) - k_for + (1-δ).*k_grid;
rhs_k = (z.*k_grid.^alpha.*(l_curr.^(1-alpha)) - k_for + (1-δ).*k_grid).^(-sig);
lhs1_k = (z.*k_for.^alpha.*(l_for.^(1-alpha)) - k_for2 + (1-δ).*k_for).^(-sig); # check here if things work
lhs2_k = ((alpha*z).*k_for.^(alpha - 1).*(l_for.^(1-alpha)) .+ (1 - δ));

l_err = (rhs_k.*((z*(1-alpha))*k_grid.^(alpha).*l_curr.^(-alpha)))./(χ * l_curr.^η).-1;

resid = β*(lhs1_k .* lhs2_k)./rhs_k .- 1;

p4 = scatter(k_grid, resid);
hline!([0.01 -0.01]);
savefig(p4, "resid.png")
savefig(p1, "value.png")
savefig(p2, "capital.png")
savefig(p3, "labour.png")



