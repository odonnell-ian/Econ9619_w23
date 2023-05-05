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
grid_size = 30; 
max_iter = 400; 
zs = 5
#l_len = 30;

# = range(10^(-6), 1, length = l_len);
k_grid = range(k_low, k_high, length = grid_size);

function rouwenhurst(rho, sd, N)
    mu_eps = 0

    q = (rho+1)/2
    nu = ((N-1)/(1-rho^2))^(1/2) * sd

    P = [q 1-q;1-q q]


    for i=2:N-1
    P = q*[P zeros(i,1); zeros(1,i+1)] + (1-q)*[zeros(i,1) P;zeros(1,i+1)] + (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P]
    P[2:i,:] = P[2:i,:]/2
    end

    return P, range(mu_eps/(1-rho)-nu,mu_eps/(1-rho)+nu,N)

end





tol = 10^-6;

# define utility function
function util(x)

    c = x[4]*(x[1]^alpha)*(x[2]^(1-alpha)) - x[3] + x[1]*(1-δ)

    u_val = ((c^(1-sig))/(1-sig)) - χ*((x[2]^(1+η))/(1+η))

    if c < 10^(-6)
        u_val = u_val - (10^10)*(c - 10^(-6))^2 
    end

    return u_val
end


inner_optimizer = NelderMead()

function Tv(v_old1, k_grid1, p, zg)
    z_len = length(zg)
    
    v_sitp = [linear_interpolation(k_grid1, v_old1[:,n]) for n in 1:z_len]

    v_new = zeros(grid_size, z_len);
    g_k = zeros(grid_size, z_len);
    g_l = zeros(grid_size, z_len);
    
    for j in 1:z_len
        
        for i in 1:grid_size

            u_curr = x2 -> (-1)*(util([k_grid1[i], x2[2], x2[1], zg[j]]) + sum(β * p[j,:] .* [v_sitp[l](x2[1]) for l in 1:z_len]))

            #k_min = k_low
            k_max = max(min(z*(k_grid[i]^alpha) + (1-δ)*k_grid[i] - 10^(-6), k_high), k_ss) 
            #k_max = k_high

            sol_curr = optimize(u_curr, [k_low, 0], [k_max, 1], [0.5*k_low+0.5*k_max, l_ss], Fminbox(inner_optimizer))
        
            v_new[i, j] = -Optim.minimum(sol_curr)
            g_k[i, j] = Optim.minimizer(sol_curr)[1]
            g_l[i, j] = Optim.minimizer(sol_curr)[2]

        end
    end
    return v_new, g_k, g_l
end

# define function to iterate over the value function 
function vfi1()
    #v_old = zeros(grid_size); # try chainging this
    p1, z_grid = rouwenhurst(0.9, 0.1, zs)
    z_grid = exp.(z_grid)
    v_old = ones(grid_size, zs)/50 
    v_curr = zeros(grid_size, zs);
    g_k1 = ones(grid_size), zs/20;
    g_l1 = ones(grid_size, zs)/20;
    v_dist = 1;
    count = 0;
    while count <= max_iter && v_dist > tol 
        v_curr, g_k1, g_l1 = Tv(v_old, k_grid, p1, z_grid);
        v_dist = maximum(abs.(v_curr./v_old.-1));
        count += 1; 
        v_old = 0.95*v_old + 0.05*v_curr; 
    end

    return v_curr, g_k1, g_l1, count, v_dist, z_grid
end

# record time and number of iterations 
@time x = vfi1()

p7 = contour(x[6], k_grid, x[1])