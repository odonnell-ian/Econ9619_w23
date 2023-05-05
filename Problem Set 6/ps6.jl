using Plots 
using Optim 
using Interpolations
using ForwardDiff
using Alert
using Roots
using LinearAlgebra

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


function Tv(v_old1, k_grid1, p, zg)
    z_len = length(zg)
    
    v_new = zeros(grid_size, z_len);
   

    k_end = repeat(k_grid, 1, z_len)
    l_end = zeros(grid_size, z_len)
    
    for j in 1:z_len
        
        for i in 1:grid_size
            
            ev = sum(β *  p[j,:] .* [v_old1[l](k_grid1[i]) for l in 1:z_len])

            # if ev < 0 
            #     error("Bad news")
            # end

            v_diff = ForwardDiff.derivative(v_old1[j], k_grid[i])
            #v_diff = Interpolations.gradient(v_old1[j], k_grid[i])[1]
            
            if v_diff < 0 # not sure what to here?????
                v_diff = 10^(-6)
                c_curr = v_diff
            else
                c_curr = v_diff^(-1/sig)
            end

            # if c_curr < 0 
            #     error("Bad news")
            # end

            k_hat = l1 -> ((χ*l1^(η+alpha))/((1-alpha)*zg[j]*v_diff))^(1/alpha) 

            l_fun = l3 -> (c_curr + k_grid1[i] - zg[j]*(k_hat(l3)^alpha)*(l3^(1-alpha)) - (1-δ)*k_hat(l3))^2

            l_prob = optimize(l_fun, 0, 1, Brent())

            l_est = Optim.minimizer(l_prob)

            k_curr = k_hat(l_est)

            v_new[i, j] = ((c_curr^(1-sig))/(1-sig)) - χ*((l_est^(1+η))/(1+η)) + ev

            k_end[i, j] = k_curr

            l_end[i, j] = l_est 

        end
    end
    #v_itp = interpolate(v_old1, BSpline(Cubic(Line(OnGrid()))))
    #v_sitp = Interpolations.scale(v_itp, k_grid1)
    # t1 = interpolate(v_new[:,1], BSpline(Cubic(Line(OnGrid()))))
    # t2 = Interpolations.scale(t1, k_grid1)
    # print(typeof(k_end[:,1]))
    # print(typeof(k_grid1))
    #t1 = interpolate(k_end[:,1], v_new[:,1], Gridded(Linear()))
    
    #v_new = [interpolate((k_end[:,n],), v_new[:,n], Gridded(Linear())) for n in 1:z_len]
    v_new = [interpolate((k_grid1,), v_new[:,n], Gridded(Linear())) for n in 1:z_len] # this is not correct, but the above line gives bounds errors 
    return v_new, k_end, l_end
end
# A = rand(20)
# A_x = 1.0:2.0:40.0
# nodes = (A_x,)
# itp = interpolate(nodes, A, Gridded(Linear()))
# itp(2.0)
# p1, z_grid = rouwenhurst(0.9, 0.1, zs)
#     z_grid = exp.(z_grid)
#     v_old = ones(grid_size, zs)/50 
#     t2 = interpolate((k_grid,), rand(30), Gridded(Linear()))
#     v_old = [interpolate(k_grid, v_old[:,1], Gridded(Linear())) for n in 1:5]
#     v_diff = Interpolations.gradient(v_old[1], k_grid[1])[1]
#     error("hi")

# define function to iterate over the value function 
function vfi1()
    #v_old = zeros(grid_size); # try chainging this
    p1, z_grid = rouwenhurst(0.9, 0.1, zs)
    z_grid = exp.(z_grid)
    v_old = repeat(k_grid, 1, zs) 
    v_old = [interpolate((k_grid,), v_old[:,n], Gridded(Linear())) for n in 1:zs]
    #v_old = [Interpolations.scale(interpolate(v_old[:,n], BSpline(Cubic(Line(OnGrid())))), k_grid) for n in 1:zs]
    # v_diff = Interpolations.gradient(v_old, k_grid[1])[1]
    # error("hi")
    v_curr = zeros(grid_size, zs);
    g_k1 = ones(grid_size), zs/20;
    g_l1 = ones(grid_size, zs)/20;
    v_dist = 1;
    count = 0;
    #print(typeof(v_old))
    while count <= max_iter && v_dist > tol 
        v_curr, g_k1, g_l1 = Tv(v_old, k_grid, p1, z_grid);
        
        v_dist = abs(maximum(maximum([abs.(v_curr[kk]./v_old[kk]).-1 for kk in 1:zs])))
       #print(typeof(v_dist), v_dist)
        count += 1; 
        v_old = v_curr; 
    end

    return v_curr, g_k1, g_l1, count, v_dist, z_grid
end

# record time and number of iterations 
@time x = vfi1()
 x1 = zeros(30, 5)
 for i = 1:30
    for j = 1:5 
        x1[i, j] = x[1][j](k_grid[i])
    end
end
p7 = contour(x[6], k_grid, x1)
p8 = contour(x[6], k_grid, x[2])
p9 = contour(x[6], k_grid, x[3])

savefig(p7, "p7.png")
savefig(p8, "p8.png")
savefig(p9, "p9.png")