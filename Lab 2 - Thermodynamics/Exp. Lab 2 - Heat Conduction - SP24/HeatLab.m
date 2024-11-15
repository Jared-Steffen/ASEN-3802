%% ASEN 3802 Lab 2 -- Heat Conduction

% Jared Steffen
% Reece Fountain
% Nathan Malyszek

% Housekeeping
clc
clear
close all


%% Material Properties

% Aluminum 7075-T651
rho_al = 2810; %[kg/m^3]
cp_al = 960; %[J/(kg*K)]
k_al = 130; %[W/(m*K)]
alpha_al = k_al / (rho_al*cp_al); %[m^2/s]

% Brass C360
rho_br = 8500; %[kg/m^3]
cp_br = 380; %[J/(kg*K)]
k_br = 115; %[W/(m*K)]
alpha_br = k_br / (rho_br * cp_br); %[m^2/s]

% Stainless Steel T-303 Annealed
rho_ss = 8000; %[kg/m^3]
cp_ss = 500; %[J/(kg*K)]
k_ss = 16.2; %[W/(m*K)]
alpha_ss = k_ss / (rho_ss * cp_ss); %[m^2/s]

% Cross Sectional Area for All Rods
d = 1 * 0.0254; %[in] -> [m]
A = pi * (d/2)^2; %[m^2]

%% Question 1

% Exmaple Thermocouple Temps and Locations
Th1 = [17.6;21.61;25.13;29.22;34.92;38.10;45.21;47.01]; %[K]
Th_loc1 = [(11/8):0.5:(39/8)]' .* 0.0254; %[m]

% Line of Best Fit
poly = polyfit(Th_loc1,Th1,1);

% T0 (y intercept) and H (slope)
T0 = poly(2);
H = poly(1);

% Equaiton of Best Fit Line
best_fit = H .* [0;Th_loc1] + T0;

% Plot
figure()
plot(Th_loc1 .* 100,Th1,'o')
hold on
plot([0;Th_loc1] .* 100,best_fit)
xlabel('Thermocouple Location [cm]')
ylabel('Thermocouple Temperature [°C]')
grid on; grid minor
hold off
title('Thermocouple Reading Along Bar')
legend('Data','Best Fit','Location','northwest')

%% Question 2

% Hand-written derivation achieved results that bn = (8*H*L*(-1)^(n))/(pi^2*(2*n-1)^2)
%% Question 3

% Variables
n = 0:10;
L = 5.875 * 0.0254; %[m]
x = 4.875 * 0.0254; %[m]

% Fourier Coefficients
t1 = 1; %[s]
t2 = 1000; %[s]
Fo1 = (alpha_al * t1)/L^2;
Fo2 = (alpha_al * t2)/L^2;

% Solve u(x,t) for t=1s and t=1000s
for i = 1:length(n)
    [u_xt1(i),~] = heat_eq_sol(x,t1,T0,H,alpha_al,L,n(i));
    [u_xt1000(i),~] = heat_eq_sol(x,t2,T0,H,alpha_al,L,n(i));
end

% Note: 1 term is sufficient for t = 1000s (steady state)

% Plot to Compare
figure()
plot(n,u_xt1,'o-')
hold on
plot(n,u_xt1000,'o-')
grid on; grid minor
xlabel('Index Number (n)')
ylabel('Thermocouple 8 Temperature [°C]')
hold off
title('Fourier Convergence Analysis')
legend('TC 8 Temperature @ t = 1s','TC 8 Temperature @ t = 1000s','Location','east')


%% Question 4

% Solving u(x,t) w/ Various Thermal Diffusivity
alpha_range = alpha_al/2:alpha_al/2:2*alpha_al;
t = 1:1000; %[s]
n = 1;
for j = 1:length(alpha_range)
    for i = 1:length(t)
        [u_xt2(i,j),~] = heat_eq_sol(x,t(i),T0,H,alpha_range(j),L,n);
    end
end

% Plot
figure();
plot(t,u_xt2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Thermocouple 8 Temperature [°C]')
title('Thermal Diffusivity Sensitivity Analysis')
legend('0.5*\alpha','\alpha','1.5*\alpha','2*\alpha','Location','northwest')

% Note: Increasing alpha leads to faster steady state of thermocouples

%% Plotting all test cases

% use dir and strsplit to make extracting material, voltage, and amperage
% easy to use for plotting titles
a=dir('*mA');
datafiles = dir("*mA");

for i=1:length(a)
    load(a(i).name)
% how to get voltage and amperage from file names?
% - options include strsplit, regex, etc.
% ultimately, we need to use the format of each file name
% 'material'_'volts'V_'amps'mA
b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3
% {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
%mat = strsplit(b{1});
v = strsplit(b{2},'V'); % volts are always in the second portion
ampval= strsplit(b{3},'mA'); % amps are always in the third portion
materials(i) = b(1);
volts(i) = str2num(v{1}); % convert string to number (vector)
amps(i) = str2num(ampval{1});
end

for i = 1:size(a,1)
data = load(a(i).name);
data_Th_steady(i,:) = data(end,3:10)';
material = materials(i);
voltage = volts(i);
current = amps(i);

time_end(i) = start_location_plot(data,material,voltage,current);
end

%% Find Analytical and and Experimental Slopes and Plot to Compare

data_Th_steady = data_Th_steady';

% Slopes/Graphs For Aluminum, Brass, Steel
 for i = 1:length(find(strcmp('Aluminum',materials) == 1))
      [H_an_al(i),H_exp_al(i),T0_al(i),best_fit_exp_al(:,i),best_fit_an_al(:,i)] = slopes(Th_loc1,data_Th_steady(:,i),A,amps(i),volts(i),k_al);
      figure();
      plot([0;Th_loc1]*100,best_fit_exp_al(:,i))
      hold on
      plot([0;Th_loc1]*100,best_fit_an_al(:,i))
      plot(Th_loc1.*100,data_Th_steady(:,i),'o')
      xlabel('Thermocouple Location [cm]')
      ylabel('Thermocouple Temperature [°C]')
      grid on; grid minor
      title('Aluminum Steady State Temparture Distribution for Analytical and Experimental Models')
      legend('Experimental','Analytical','Location','northwest')
 end
 hold off


 for i = 1:length(find(strcmp('Brass',materials) == 1))
      [H_an_br(i),H_exp_br(i),T0_br(i),best_fit_exp_br(:,i),best_fit_an_br(:,i)] = slopes(Th_loc1,data_Th_steady(:,i+2),A,amps(i+2),volts(i+2),k_br);
      figure();
      plot([0;Th_loc1]*100,best_fit_exp_br(:,i))
      hold on
      plot([0;Th_loc1]*100,best_fit_an_br(:,i))
      plot(Th_loc1.*100,data_Th_steady(:,i+2),'o')
      xlabel('Thermocouple Location [cm]')
      ylabel('Thermocouple Temperature [°C]')
      grid on; grid minor
      title('Brass Steady State Temparture Distribution for Analytical and Experimental Models')
      legend('Experimental','Analytical','Location','northwest')
 end
 hold off

  for i = 1:length(find(strcmp('Steel',materials) == 1))
      [H_an_ss(i),H_exp_ss(i),T0_ss(i),best_fit_exp_ss(:,i),best_fit_an_ss(:,i)] = slopes(Th_loc1,data_Th_steady(:,i+4),A,amps(i+4),volts(i+4),k_ss);
      figure();
      plot([0;Th_loc1]*100,best_fit_exp_ss(:,i))
      hold on
      plot([0;Th_loc1]*100,best_fit_an_ss(:,i))
      plot(Th_loc1.*100,data_Th_steady(:,i+4),'o')
      xlabel('Thermocouple Location [cm]')
      ylabel('Thermocouple Temperature [°C]')
      grid on; grid minor
      title('Steel Steady State Temparture Distribution for Analytical and Experimental Models')
      legend('Experimental','Analytical','Location','northwest')
  end
  hold off


%% Model 1A

% Time Vectors
time_al1 = 0:10:time_end(1);
time_al2 = 0:10:time_end(2);
time_br1 = 0:10:time_end(3);
time_br2 = 0:10:time_end(4);
time_steel = 0:10:time_end(5);

% Solve 1D Heat Equation For All Cases
for j = 1:length(Th_loc1)
    for i = 1:length(time_al1)
        [u_xt_al1(i,j),~] = heat_eq_sol(Th_loc1(j),time_al1(i),T0_al(1),H_an_al(1),alpha_al,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_al2)
        [u_xt_al2(i,j),~] = heat_eq_sol(Th_loc1(j),time_al2(i),T0_al(2),H_an_al(2),alpha_al,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br1)
        [u_xt_br1(i,j),~] = heat_eq_sol(Th_loc1(j),time_br1(i),T0_br(1),H_an_br(1),alpha_br,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br2)
        [u_xt_br2(i,j),~] = heat_eq_sol(Th_loc1(j),time_br2(i),T0_br(2),H_an_br(2),alpha_br,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_steel)
        [u_xt_steel(i,j),~] = heat_eq_sol(Th_loc1(j),time_steel(i),T0_ss,H_an_ss,alpha_ss,L,1);
    end
end

% Load Data
data_al1 = load(a(1).name);
data_al2 = load(a(2).name);
data_br1 = load(a(3).name);
data_br2 = load(a(4).name);
data_ss = load(a(5).name);

% Plot Model 1A
plot_model(u_xt_al1,time_al1,data_al1,materials(1),volts(1),amps(1),'1A');
plot_model(u_xt_al2,time_al2,data_al2,materials(2),volts(2),amps(2),'1A');
plot_model(u_xt_br1,time_br1,data_br1,materials(3),volts(3),amps(3),'1A');
plot_model(u_xt_br2,time_br2,data_br2,materials(4),volts(4),amps(4),'1A');
plot_model(u_xt_steel,time_steel,data_ss,materials(5),volts(5),amps(5),'1A');

%% Model 1B

% Solve 1D Heat Equation For All Cases w/ Experimental Slope
for j = 1:length(Th_loc1)
    for i = 1:length(time_al1)
        [u_xt_al1_B(i,j),~] = heat_eq_sol(Th_loc1(j),time_al1(i),T0_al(1),H_exp_al(1),alpha_al,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_al2)
        [u_xt_al2_B(i,j),~] = heat_eq_sol(Th_loc1(j),time_al2(i),T0_al(2),H_exp_al(2),alpha_al,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br1)
        [u_xt_br1_B(i,j),~] = heat_eq_sol(Th_loc1(j),time_br1(i),T0_br(1),H_exp_br(1),alpha_br,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br2)
        [u_xt_br2_B(i,j),~] = heat_eq_sol(Th_loc1(j),time_br2(i),T0_br(2),H_exp_br(2),alpha_br,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_steel)
        [u_xt_steel_B(i,j),~] = heat_eq_sol(Th_loc1(j),time_steel(i),T0_ss,H_exp_ss,alpha_ss,L,1);
    end
end

% Plot Model 1B
plot_model(u_xt_al1_B,time_al1,data_al1,materials(1),volts(1),amps(1),'1B');
plot_model(u_xt_al2_B,time_al2,data_al2,materials(2),volts(2),amps(2),'1B');
plot_model(u_xt_br1_B,time_br1,data_br1,materials(3),volts(3),amps(3),'1B');
plot_model(u_xt_br2_B,time_br2,data_br2,materials(4),volts(4),amps(4),'1B');
plot_model(u_xt_steel_B,time_steel,data_ss,materials(5),volts(5),amps(5),'1B');

%% Model 2

% Initital Temp Distribution Slope
M_exp_al1 = init_dist_plot(data_al1,Th_loc1,materials(1));
M_exp_al2 = init_dist_plot(data_al2,Th_loc1,materials(2));
M_exp_br1 = init_dist_plot(data_br1,Th_loc1,materials(3));
M_exp_br2 = init_dist_plot(data_br2,Th_loc1,materials(4));
M_exp_ss = init_dist_plot(data_ss,Th_loc1,materials(5));

% Solve 1D Heat Equation For All Cases w/ Experimental Slope and Initial Temp Distributions
for j = 1:length(Th_loc1)
    for i = 1:length(time_al1)
        u_xt_al1_2(i,j) = heat_eq_sol2(Th_loc1(j),time_al1(i),T0_al(1),H_exp_al(1),alpha_al,L,1,M_exp_al1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_al2)
        u_xt_al2_2(i,j) = heat_eq_sol2(Th_loc1(j),time_al2(i),T0_al(2),H_exp_al(2),alpha_al,L,1,M_exp_al2);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br1)
        u_xt_br1_2(i,j) = heat_eq_sol2(Th_loc1(j),time_br1(i),T0_br(1),H_exp_br(1),alpha_br,L,1,M_exp_br1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br2)
        u_xt_br2_2(i,j) = heat_eq_sol2(Th_loc1(j),time_br2(i),T0_br(2),H_exp_br(2),alpha_br,L,1,M_exp_br2);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_steel)
        u_xt_steel_2(i,j) = heat_eq_sol2(Th_loc1(j),time_steel(i),T0_ss,H_exp_ss,alpha_ss,L,1,M_exp_ss);
    end
end

% Plot vs Experimental Data
plot_model(u_xt_al1_2,time_al1,data_al1,materials(1),volts(1),amps(1),'2')
plot_model(u_xt_al2_2,time_al2,data_al2,materials(2),volts(2),amps(2),'2');
plot_model(u_xt_br1_2,time_br1,data_br1,materials(3),volts(3),amps(3),'2');
plot_model(u_xt_br2_2,time_br2,data_br2,materials(4),volts(4),amps(4),'2');
plot_model(u_xt_steel_2,time_steel,data_ss,materials(5),volts(5),amps(5),'2');

%% Model 3

% Find Adjusted Thermal Diffusivity
alpha_adj_al = alpha_var(u_xt_al1_B(end-50:end,8),alpha_al,x,time_al1,T0_al(1),H_exp_al(1),L,1);
alpha_adj_br = alpha_var(u_xt_br1_B(end-50:end,8),alpha_br,x,time_br1,T0_br(1),H_exp_br(1),L,1);
alpha_adj_ss = alpha_var(u_xt_steel_B(end-50:end,8),alpha_ss,x,time_steel,T0_ss,H_exp_ss,L,1);

% Solve 1D Heat Equation For All Cases w/ Experimental Slope and Adjusted Thermal Diffusivity
for j = 1:length(Th_loc1)
    for i = 1:length(time_al1)
        [u_xt_al1_3(i,j),Fo_al1(i)] = heat_eq_sol(Th_loc1(j),time_al1(i),T0_al(1),H_exp_al(1),alpha_adj_al,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_al2)
        [u_xt_al2_3(i,j),Fo_al2(i)] = heat_eq_sol(Th_loc1(j),time_al2(i),T0_al(2),H_exp_al(2),alpha_adj_al,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br1)
        [u_xt_br1_3(i,j),Fo_br1(i)] = heat_eq_sol(Th_loc1(j),time_br1(i),T0_br(1),H_exp_br(1),alpha_adj_br,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_br2)
        [u_xt_br2_3(i,j),Fo_br2(i)] = heat_eq_sol(Th_loc1(j),time_br2(i),T0_br(1),H_exp_br(2),alpha_adj_br,L,1);
    end
end

for j = 1:length(Th_loc1)
    for i = 1:length(time_steel)
        [u_xt_steel_3(i,j),Fo_ss(i)] = heat_eq_sol(Th_loc1(j),time_steel(i),T0_ss,H_exp_ss,alpha_adj_ss,L,1);
    end
end

% Plot Results
plot_model(u_xt_al1_3,time_al1,data_al1,materials(1),volts(1),amps(1),'3');
plot_model(u_xt_al2_3,time_al2,data_al2,materials(2),volts(2),amps(2),'3');
plot_model(u_xt_br1_3,time_br1,data_br1,materials(3),volts(3),amps(3),'3');
plot_model(u_xt_br2_3,time_br2,data_br2,materials(4),volts(4),amps(4),'3');
plot_model(u_xt_steel_3,time_steel,data_ss,materials(5),volts(5),amps(5),'3');

% Fourier Number and Time to Steady State
[t_al_ss,Fo_al_new] = FoNum(Fo_al1,alpha_adj_al,time_al1,L);
[t_br_ss,Fo_br_new] = FoNum(Fo_br1,alpha_adj_br,time_br1,L);
[t_steel_ss,Fo_steel_new] = FoNum(Fo_ss,alpha_adj_ss,time_steel,L);

%% Comparison of Models

% Gathering Thermocouple 8 Data for Models 2 and 3
% Aluminum 25V
al1_mod2 = u_xt_al1_2(:,end);
al1_mod3 = u_xt_al1_3(:,end);

% Aluminum 30V
al2_mod2 = u_xt_al2_2(:,end);
al2_mod3 = u_xt_al2_3(:,end);

% Brass 25V
br1_mod2 = u_xt_br1_2(:,end);
br1_mod3 = u_xt_br1_3(:,end);

% Brass 30V
br2_mod2 = u_xt_br2_2(:,end);
br2_mod3 = u_xt_br2_3(:,end);

% Steel 22V
ss_mod2 = u_xt_steel_2(:,end);
ss_mod3 = u_xt_steel_3(:,end);

% Plot
plot_errorbars(data_al1,al1_mod2,al1_mod3,time_al1,materials(1),volts(1))
plot_errorbars(data_al2,al2_mod2,al2_mod3,time_al2,materials(2),volts(2))
plot_errorbars(data_br1,br1_mod2,br1_mod3,time_br1,materials(3),volts(3))
plot_errorbars(data_br2,br2_mod2,br2_mod3,time_br2,materials(4),volts(4))
plot_errorbars(data_ss,ss_mod2,ss_mod3,time_steel,materials(5),volts(5))

%% Functions
function [u_xt,Fo] = heat_eq_sol(x,t,T0,H,alpha,L,n)
%----------------------------------------------------------------
% Function to solve u(x,t) for a given n and material properties
% Inputs: 
% x - thermocouple location(s) [m]
% t - time vector [s]
% T0 - initial temperature [°C]
% H - steady state temperature distrubtion slope [°C/m]
% alpha - thermal diffusivity [m^2/s]
% L - rod length (not true length) [m]
% n - number of iterations [unitless]
%
% Outputs:
% u_xt - 1D heat equation solution [°C]
%----------------------------------------------------------------

% n Variable Functions
bn = @(n) (8*H*L*(-1)^(n))/(pi^2*(2*n-1)^2);
lambda_n = @(n) (2*n-1)*pi/(2*L);

% Run Summation for n Terms
gx_sum = 0;
for i = 1:n
    gx_sum = gx_sum + bn(i)*sin(lambda_n(i)*x)*exp(-(lambda_n(i)^2*alpha*t));
end
    
% Calculate u(x,t)
u_xt = T0 + H*x + gx_sum;


FoFunc = @(a,t,l)(a*t)/((L)^2); %anon fn, a "alpha", t "time", L "length"

% Fourier calculation
Fo = FoFunc(alpha,t,L);


end

function time_end = start_location_plot(DataSet, Material, Voltage, Current)
%----------------------------------------------------------------
% Function to solve from end of time vector and plot experimental data
% Inputs: 
% DataSet - data for a specific test case [several units]
% Material - material type
% Voltage - test voltage value [V]
% Current - test current value [mA]
%
% Outputs:
% time_end - end of time vector [s]
%----------------------------------------------------------------

% extract time vector from first column
time = (DataSet(:,1));

% look at last thermocouple data in order to find where to start plot
val = DataSet(:,end);

% find where to start plot using zeroth, first, and second derivatives of
% line
x_start = find(val == val(1),1,"last");

dx = diff(val);
dx_start = find(dx > 0.1,1);

d2x = diff(diff(val));
d2x_start = find(d2x > 0.2,1);

starting_vals = [x_start dx_start d2x_start];

% find the smallest starting index out of the three methods to use for plotting
start_val = min(starting_vals);

% change time vector to reflect new starting point
time = time(start_val:end);
time = time - time(1);
time_end = time(end);

% plot, ignoring first and second columns to only plot relative data
figure()
for i = 3:size(DataSet,2)
    plot(time, DataSet(start_val:end,i))
    hold on
end
grid on; grid minor
xlabel("Time [s]")
ylabel("Temperature [°C]")
title("Test Case for " + Material + " at " + Voltage + "V and " + Current + "mA")

end

function [H_an,H_exp,T0,best_fit_exp,best_fit_an] = slopes(Th_loc,Th,A,I,V,k)
%----------------------------------------------------------------
% Function to solve for steady temp distribution slope, initial temp, and
% best fit lines
% Inputs: 
% Th_loc - thermocouple location(s) [m]
% Th - test case temperature [°C]
% A - cross-sectional area of rod [m^2]
% I - test case current [mA]
% V - test case voltage [V]
% k - thermal conductivity of material [W/(m*K)]
%
% Outputs:
% H_an - analytical steady state temperature distribution slope [°C/m]
% H_exp - experimental steady state temperature distribution slope [°C/m]
% T0 - initial temperature [°C]
% best_fit_exp = line of best fit equation for experimental data
% best_fit_an = line of best fit equation for analytical data
%----------------------------------------------------------------

% Polyfit Data
poly = polyfit(Th_loc,Th,1);

% Experimental Slope and Y-Intercept
H_exp = poly(1);
T0 = poly(2);

% Analytical Slope
I = I / 1e3;
Q_dot = V*I;
H_an = Q_dot / (k*A);

% Best Fit Line
best_fit_an = H_an .* [0;Th_loc] + T0;
best_fit_exp = H_exp .* [0;Th_loc] + T0;

end

function time_end = plot_model(u_xt,time_mat,DataSet, Material, Voltage, Current, Model)
%----------------------------------------------------------------
% Function to solve from end of time vector and plot models
% Inputs: 
% u_xt - 1D heat equation solution [°C]
% time_mat - original time vector for a specific material [s]
% DataSet - data for a specific test case [several units]
% Material - material type
% Voltage - test voltage value [V]
% Current - test current value [mA]
% Model - model number [1A,1B,2,3]
%
% Outputs:
% time_end - end of time vector [s]
%----------------------------------------------------------------

% extract time vector from first column
time = (DataSet(:,1));

% look at last thermocouple data in order to find where to start plot
val = DataSet(:,end);

% find where to start plot using zeroth, first, and second derivatives of
% line
x_start = find(val == val(1),1,"last");

dx = diff(val);
dx_start = find(dx > 0.1,1);

d2x = diff(diff(val));
d2x_start = find(d2x > 0.2,1);

starting_vals = [x_start dx_start d2x_start];

% find the smallest starting index out of the three methods to use for plotting
start_val = min(starting_vals);

% change time vector to reflect new starting point
time = time(start_val:end);
time = time - time(1);
time_end = time(end);

% plot, ignoring first and second columns to only plot relative data
figure()
for i = 3:size(DataSet,2)
    p1 = plot(time, DataSet(start_val:end,i),'-r');
    hold on
end
p2 = plot(time_mat, u_xt,'-b');
hold off
grid on; grid minor
xlabel("Time [s]")
ylabel("Temperature [°C]")
title("Model " + Model + " " + Material + " at " + Voltage + "V and " + Current + "mA")
legend([p1(1);p2(1)],{'Experimental';'Model'},'Location','southeast')

end

function M_exp = init_dist_plot(data,Th_loc,Material)
%----------------------------------------------------------------
% Function to solve for initial temp distribution slope
% Inputs: 
% data - data for specific test case [°C]
% Th_loc - thermocouple location(s) [m]
% Material - material type
%
% Outputs:
% M_exp -initial temperature distribution slope [°C/m]
%----------------------------------------------------------------

% Find and Polyfit Initial Temps
init_temp = data(1,3:end);
polyfit_data = polyfit(Th_loc,init_temp,1);

% Slope/Y-Intercept
M_exp = polyfit_data(1);
y_int = polyfit_data(2);

% Best Fit Line and Average (for y-scale)
best_fit = M_exp .* Th_loc + y_int;
avg = mean(init_temp);

% Plot
figure()
plot(Th_loc .* 100,init_temp,'o')
hold on
plot(Th_loc .* 100,best_fit)
grid on; grid minor
xlabel('Thermocouple Location [cm]')
ylabel('Thermocouple Temperature [°C]')
ylim([avg-1 avg+1])
title(Material + " Initial Temperature Distribution")
legend('Experimental Data','Best Fit')

end

function u_xt = heat_eq_sol2(x,t,T0,H,alpha,L,n,M)
%----------------------------------------------------------------
% Function to solve u(x,t) for a given n and material properties
% Inputs: 
% x - thermocouple location(s) [m]
% t - time vector [s]
% T0 - initial temperature [°C]
% H - steady state temperature distrubtion slope [°C/m]
% alpha - thermal diffusivity [m^2/s]
% L - rod length (not true length) [m]
% n - number of iterations [unitless]
% M - initial temperature distribution slope [°C/m]
%
% Outputs:
% u_xt - 1D heat equation solution [°C]
%----------------------------------------------------------------

% n Variable Functions
bn = @(n) (-8*(M-H)*L*(-1)^(n))/(pi^2*(2*n-1)^2);
lambda_n = @(n) (2*n-1)*pi/(2*L);

% Run Summation for n Terms
gx_sum = 0;
for i = 1:n
    gx_sum = gx_sum + bn(i)*sin(lambda_n(i)*x)*exp(-(lambda_n(i)^2*alpha*t));
end
    
% Calculate u(x,t)
u_xt = T0 + H*x + gx_sum;

end


function alpha_adj = alpha_var(u_xt,alpha,x,t,T0,H,L,n)
%----------------------------------------------------------------
% Function to solve for adjusted thermal diffusivity
% Inputs: 
% u_xt - 1D heat equation solution [°C]
% alpha - thermal diffusivity [m^2/s]
% x - thermocouple location(s) [m]
% t - time vector [s]
% T0 - initial temperature [°C]
% H - steady state temperature distrubtion slope [°C/m]
% alpha - thermal diffusivity [m^2/s]
% L - rod length (not true length) [m]
% n - number of iterations [unitless]
%
% Outputs:
% alpha_adj - adjusted thermal diffusivity [m^2/s]
%----------------------------------------------------------------

% Find Average of Steady State Variations
d = mean(u_xt);

% Range of Thermal Diffusivity to Test
alpha_range = alpha/1000:alpha/1000:1.5*alpha;

% Test Last Thermocouple with Alpha Range
for j = 1:length(alpha_range)
    for i = 1:length(t)
        [u_xt_alpha_var(i,j),~] = heat_eq_sol(x,t(i),T0(1),H(1),alpha_range(j),L,n);
    end
end

% Find Which Alpha Leads to Closest Steady State Average
e = find(u_xt_alpha_var(end,:) >= d, 1, 'first');
alpha_adj = alpha_range(e);

end

function [t_ss,Fo] = FoNum(Fo_exp,alpha_adj,time,L)
%----------------------------------------------------------------
% Function to solve for time to steady state and fourier number
% Inputs: 
% Fo_exp - fourier number from experimental data [unitless]
% alpha_adj - adjusted thermal diffusivity [m^2/s]
% time - time vector [s]
% L - rod length (not true length) [m]
%
% Outputs:
% t_ss - time to steady state [s]
% Fo - corresponding fourier number [unitless]
%----------------------------------------------------------------

Fo_pos = find(Fo_exp > 2,1,'first');

t_ss = time(Fo_pos);

Fo = (alpha_adj * t_ss)/L^2;

end

function [] = plot_errorbars(DataSet,mod2,mod3,tmod,Material,Voltage)
%----------------------------------------------------------------
% Function to plot thermocouple 8 final temp w/ error bars
% Inputs: 
% DataSet - data for a specific test case [several units]
% mod2 - model 2 final thermocouple 8 data for a specific test case [°C]
% mod3 - model 3 final thermocouple 8 data for a specific test case [°C]
% t_exp - experimental time vector [s]
% tmod - model time vector [s]
% Material - material type
% Voltage - test voltage value [V]
%
% Outputs:
% Thermocouple 8 final temperature comparison bar plot w/ error bars
%----------------------------------------------------------------

% extract time vector from first column
time = (DataSet(:,1));

% look at last thermocouple data in order to find where to start plot
val = DataSet(:,end);

% find where to start plot using zeroth, first, and second derivatives of
% line
x_start = find(val == val(1),1,"last");

dx = diff(val);
dx_start = find(dx > 0.1,1);

d2x = diff(diff(val));
d2x_start = find(d2x > 0.2,1);

starting_vals = [x_start dx_start d2x_start];

% find the smallest starting index out of the three methods to use for plotting
start_val = min(starting_vals);

% change time vector to reflect new starting point
time = time(start_val:end);
time = time - time(1);


% Error Vectors
err = 2 * zeros(length(DataSet(start_val:end,end)),1);
for i = 1:length(err)
    if rem(i,20) == 0 % if remainder or i/20 is zero sets value to 2
        err(i) = 2; 
    else
        err(i) = NaN; % otherwise removes data point
    end
end

% Plot Comparison
figure();
errorbar(time,DataSet(start_val:end,end),err,"LineStyle","none")
hold on
plot(time,DataSet(start_val:end,end))
plot(tmod,mod2)
plot(tmod,mod3)
grid on; grid minor
xlabel('Time [s]')
ylabel('Thermocouple Temperature [°C]')
legend('Experimental Error','Experimental Data','Model 2','Model 3','Location','southeast')
title(Material + " at " + Voltage + " V Final Model Comparisons")
hold off

end
