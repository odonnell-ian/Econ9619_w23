using Plots

# Initial Parameters - assume that beta is 0.95
a = 1/3;
z = 1;
β = 0.95;
n = 100;

# initial steady state
k_ss0 = (β*a*z)^(1/(1-a));
y_ss0 = z*k_ss0^a;
c_ss0 = y_ss0 - k_ss0;
r_ss0 = z*a*k_ss0^(a-1);
w_ss0 = z*(k_ss0^a)*(1-a);

# calculate period 1/after shock values 
k1 = 0.8*k_ss0;
z1 = 1.05*z;
y1 = z1*k1^a;
c1= y1 - k1;
r1 = z1*a*k1^(a-1);
w1 = z1*(k1^a)*(1-a);

# calculate end steady state as a check 
k_ss2 = (β*a*z1)^(1/(1-a));
y_ss2 = z1*k_ss2^a;
c_ss2 = y_ss2 - k_ss2;
r_ss2 = z1*a*k_ss2^(a-1);
w_ss2 = z1*(k_ss2^a)*(1-a);

# initalize paths
k_path = zeros(n);
k_path[1] = k1;

y_path = zeros(n);
y_path[1] = y1;

c_path = zeros(n);
c_path[1] = c1;

r_path = zeros(n);
r_path[1] = r1;

w_path = zeros(n);
w_path[1] = w1;

# define a function for g 

# write a loop to fill in the k path and others 
for i = 2:n
    k_path[i] = β*a*z1*k_path[i-1]^a;
    y_path[i] = z1*k_path[i]^a;
    c_path[i] = y_path[i] - k_path[i];
    r_path[i] = z1*a*k_path[i]^(a-1);
    w_path[i] = z1*(k_path[i]^a)*(1-a);
end

# make plots 
pk = plot(1:n, k_path, lw=3, title="Capital", label = "Path", xlabel = "Time", ylabel = "k")
hline!([k_ss0 k_ss2], ls=:dot, lw = 3, label = ["Initial SS" "Final SS"])
savefig(pk, "fig1.png")

py = plot(1:n, y_path, lw=3, title="Output", label = "Path", xlabel = "Time", ylabel = "y")
hline!([y_ss0 y_ss2], ls=:dot, lw = 3, label = ["Initial SS" "Final SS"])
savefig(py, "fig2.png")

pc = plot(1:n, c_path, lw=3, title="Consumption", label = "Path", xlabel = "Time", ylabel = "c")
hline!([c_ss0 c_ss2], ls=:dot, lw = 3, label = ["Initial SS" "Final SS"])
savefig(pc, "fig3.png")

pr = plot(1:n, r_path, lw=3, title="Return on Capital", label = "Path", xlabel = "Time", ylabel = "r")
hline!([r_ss0 r_ss2], ls=:dot, lw = 3, label = ["Initial SS" "Final SS"])
savefig(pr, "fig4.png")

pw = plot(1:n, w_path, lw=3, title="Wage", label = "Path", xlabel = "Time", ylabel = "w")
hline!([w_ss0 w_ss2], ls=:dot, lw = 3, label = ["Initial SS" "Final SS"])
savefig(pw, "fig5.png")


