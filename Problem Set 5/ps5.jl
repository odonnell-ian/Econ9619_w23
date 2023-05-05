using Plots 
using Random 
using Distributions
using Statistics
using StatsBase

Random.seed!(0406)

# Problem 1 ---------------------------------------------
sam = 10000
sig1 = 1
rho1 = 0.9

# simulate markov chain 
d = Normal(0, sig1)
d2 = Normal()
η = rand(d, sam-1)

z = zeros(sam)

for i = 2:sam 
    z[i] = rho1*z[i-1] + η[i-1]
end

# tauchen 
function tauchen(rho, sd, N, w)
    d = 2 * w * sd / ( (1-rho^2)^0.5 * (N-1) )  # distance between two points
    z_1 = - w * sd / ( (1-rho^2)^0.5 )
    z_N = + w * sd / ( (1-rho^2)^0.5 )

    z_curr = z_1:d:z_N;
    pi=zeros(N,N);

    for i = 1 : N
        for j = 2 : N-1  
            pi[i,j] = cdf(d2, (z_curr[j]+d/2-rho*z_curr[i])/sd) - cdf(d2, (z_curr[j]-d/2-rho*z_curr[i])/sd);        
        end
        pi[i,1] = cdf(d2, (z_curr[1]+d/2-rho*z_curr[i])/sd) ; 
        pi[i,N] = 1 - cdf(d2, (z_curr[N]-d/2-rho*z_curr[i])/sd);
    end
    return pi, z_curr
end

mc_tauchen5, grid_t5 = tauchen(rho1, sig1, 5, 3)
mc_tauchen15, grid_t15 = tauchen(rho1, sig1, 15, 3)

# Rouwenhorst 
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

mc_row5, grid_r5 = rouwenhurst(rho1, sig1, 5)
mc_row15, grid_r15 = rouwenhurst(rho1, sig1, 15) 

# simulate
function mc_sample_path(P, init, sample_size)
    # function from quant econ 
   @assert size(P)[1] == size(P)[2] # square required
    N = size(P)[1] # should be square

    # create vector of discrete RVs for each row
    dists = [Categorical(P[i, :]) for i in 1:N]

    # setup the simulation
    X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
    X[1] = init # set the initial state

    for t in 2:sample_size
        dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
        X[t] = rand(dist) # draw new value
    end
    return X
end

sim_t_5 = grid_t5[mc_sample_path(mc_tauchen5, 3, sam)]
sim_t_15 = grid_t15[mc_sample_path(mc_tauchen15, 8, sam)]

sim_r_5 = grid_r5[mc_sample_path(mc_row5, 3, sam)]
sim_r_15 = grid_r15[mc_sample_path(mc_row15, 8, sam)]

# find moments 
mz = mean(z)
mt5 = mean(sim_t_5)
mt15 = mean(sim_t_15)

mr5 = mean(sim_r_5)
mr15 = mean(sim_r_15)

vz = var(z)
vt5 = var(sim_t_5)
vt15 = var(sim_t_15)

vr5 = var(sim_r_5)
vr15 = var(sim_r_15)

skz = skewness(z)
skt5 = skewness(sim_t_5)
skt15 = skewness(sim_t_15)

skr5 = skewness(sim_r_5)
skr15 = skewness(sim_r_15)

kz = kurtosis(z)
kt5 = kurtosis(sim_t_5)
kt15 = kurtosis(sim_t_15)

kr5 = kurtosis(sim_r_5)
kr15 = kurtosis(sim_r_15)

# autocorrlation 
az = autocor(z, [0,1,2,3,4], demean = false)
at5 = autocor(sim_t_5, [0,1,2,3,4], demean = false)
at15 = autocor(sim_t_15, [0,1,2,3,4], demean = false)
ar5 = autocor(sim_r_5, [0,1,2,3,4], demean = false)
ar15 = autocor(sim_r_15, [0,1,2,3,4], demean = false)

p1 = plot(az)
p2 = plot(at5)
p3 = plot(at15)
p4 = plot(ar5)
p5 = plot(ar15)


# histograms 
p6 = histogram(z)
p7 = histogram(sim_t_5)
p8 = histogram(sim_t_15)
p9 = histogram(sim_r_5)
p11 = histogram(sim_r_15)



savefig(p1, "p1.png")
savefig(p2, "p2.png")
savefig(p3, "p3.png")
savefig(p4, "p4.png")
savefig(p5, "p5.png")
savefig(p6, "p6.png")
savefig(p7, "p7.png")
savefig(p8, "p8.png")
savefig(p9, "p9.png")
savefig(p11, "p11.png")

