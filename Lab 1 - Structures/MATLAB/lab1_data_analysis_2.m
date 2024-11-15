%% House Keeping
clear,clc,close all

%% Read In Data
case_1_data = readmatrix("Case 1.txt");
case_2_data = readmatrix("Case 2.txt");
case_3_data = readmatrix("Case 3.txt");

case_1_data = case_1_data * 4.4482216153;
case_2_data = case_2_data * 4.4482216153;
case_3_data = case_3_data * 4.4482216153;


%% Sort Data

% Applied force for each case
case_1_applied_force = case_1_data(:,1);
case_2_applied_force = case_2_data(:,1);
case_3_applied_force = case_3_data(:,1);

% Reaction forces for each case
case_1_load_sensor_1 = case_1_data(:,2); % Case 1
case_1_load_sensor_2 = case_1_data(:,3);
case_1_load_sensor_3 = case_1_data(:,4);

case_2_load_sensor_1 = case_2_data(:,2); % Case 2
case_2_load_sensor_2 = case_2_data(:,3);
case_2_load_sensor_3 = case_2_data(:,4);

case_3_load_sensor_1 = case_3_data(:,2); % Case 3
case_3_load_sensor_2 = case_3_data(:,3);
case_3_load_sensor_3 = case_3_data(:,4);

% Internal force for each case
case_1_inline_load_cells = case_1_data(:,5) - mean(case_1_data(1:10,5));
case_2_inline_load_cells = case_2_data(:,5) - mean(case_2_data(1:10,5));
case_3_inline_load_cells = case_3_data(:,5) - mean(case_3_data(1:10,5));

% Displacement for each case
case_1_LVDT = (case_1_data(:,6)./4.4482216153)*25.4;
case_2_LVDT = (case_2_data(:,6)./4.4482216153)*25.4;
case_3_LVDT = (case_3_data(:,6)./4.4482216153)*25.4;

%% Finding Loading Location

% Both cases assume that loading sensor 3 = Rb
% Case 2
L = 4;
a1 = (L.*case_2_load_sensor_3)./case_2_applied_force;

% Case 3
a2 = (L.*case_3_load_sensor_3)./case_3_applied_force;

E = 70E9;
I = 2.4761E-6;
RA = max(case_3_load_sensor_1 + case_3_load_sensor_2);
RB = max(case_3_load_sensor_3);
P = (RA+RB)/2;
L = 4;
Lp = 2;
a = 0:.25:2;
% a3 = 0.08:0.01:2.04;
% a4 = 4:-0.01:2.04;
% 
% v_left = 1/(E*I) .*((RA/6)*x^3 - (RA/2) .* a3.^2.*x);
% v_right = 1/(E*I) .*((RA/6 - P/12)*x^3+(P/4).*a3.*x^2+((P/4).*a4.^2-(P/2).*a3.*a4-(RA/2).*a4.^2).*x);
% v_left = v_left .* 1000;
% v_right = v_right .* 1000;
% 
% v_tot = v_left + v_right;
% 
v_c = max(case_3_LVDT);

C3 = ((P/2).*a -(1/6)*P*(Lp^2)+(1/Lp).*((-P/2).*a.*(L.^2) - (1/6)*RA*(L^3)+(P/6)*(L^3)-P.*a.*(Lp.^2)+P.*a.*L.*Lp+(P/2)*Lp^3-(P/2)*(Lp^2)*L))./(1-Lp+L);
C1 = C3 -P.*a.*Lp+(P/2)*Lp;
C2 = (-P/2).*a.*L.^2 - (1/6)*RA*L^3+(P/6)*L^3-C1.*L;

v_r = (1/(E*I)).*((P/2).*a.*Lp^2 + (1/6)*(RA-P)*(Lp^3)+C1*Lp+C2);
v_l = (1/(E*I)).*((1/6)*RA*(Lp^3)+C3*Lp);

v_tot = v_l + v_r;

location1 = a(3); %m
location2 = (L - location1); %m

errorv = abs((-v_tot(3)*1000 - v_c)/v_c) *100;


% %% Linear Regression

% Case 1
P_R1_1 = polyfit(case_1_applied_force, case_1_load_sensor_1, 1);
P_R2_1 = polyfit(case_1_applied_force, case_1_load_sensor_2, 1);
P_R3_1 = polyfit(case_1_applied_force, case_1_load_sensor_3, 1);
P_I_1  = polyfit(case_1_applied_force, case_1_inline_load_cells, 1);
P_D_1 = polyfit(case_1_applied_force, case_1_LVDT, 1);

fit_R1_1 = polyval(P_R1_1, case_1_applied_force);
fit_R2_1 = polyval(P_R2_1, case_1_applied_force);
fit_R3_1 = polyval(P_R3_1, case_1_applied_force);
fit_I_1 = polyval(P_I_1, case_1_applied_force);
fit_D_1 = polyval(P_D_1, case_1_applied_force);

% Case 2
P_R1_2 = polyfit(case_2_applied_force, case_2_load_sensor_1, 1);
P_R2_2 = polyfit(case_2_applied_force, case_2_load_sensor_2, 1);
P_R3_2 = polyfit(case_2_applied_force, case_2_load_sensor_3, 1);
P_I_2  = polyfit(case_2_applied_force, case_2_inline_load_cells, 1);
P_D_2 = polyfit(case_2_applied_force, case_2_LVDT, 1);

fit_R1_2 = polyval(P_R1_2, case_2_applied_force);
fit_R2_2 = polyval(P_R2_2, case_2_applied_force);
fit_R3_2 = polyval(P_R3_2, case_2_applied_force);
fit_I_2 = polyval(P_I_2, case_2_applied_force);
fit_D_2 = polyval(P_D_2, case_2_applied_force);

% Case 3
P_R1_3 = polyfit(case_3_applied_force, case_3_load_sensor_1, 1);
P_R2_3 = polyfit(case_3_applied_force, case_3_load_sensor_2, 1);
P_R3_3 = polyfit(case_3_applied_force, case_3_load_sensor_3, 1);
P_I_3  = polyfit(case_3_applied_force, case_3_inline_load_cells, 1);
P_D_3 = polyfit(case_3_applied_force, case_3_LVDT, 1);

fit_R1_3 = polyval(P_R1_3, case_3_applied_force);
fit_R2_3 = polyval(P_R2_3, case_3_applied_force);
fit_R3_3 = polyval(P_R3_3, case_3_applied_force);
fit_I_3 = polyval(P_I_3, case_3_applied_force);
fit_D_3 = polyval(P_D_3, case_3_applied_force);

%% Standard Deviation of the forces 
N = 110; %The number of measurements taken 

%Reaction Force #1 
sigma_force_R1_1 = ((1)/(N-1)) .* (mean(case_1_load_sensor_1) - case_1_load_sensor_1).^2;
sigma_force_R1_C1 = sqrt(sigma_force_R1_1);

sigma_force_R1_2 = ((1)/(N-1)) .* (mean(case_1_load_sensor_2) - case_1_load_sensor_2).^2;
sigma_force_R1_C2 = sqrt(sigma_force_R1_2);

sigma_force_R1_3 = ((1)/(N-1)) .* (mean(case_1_load_sensor_3) - case_1_load_sensor_3).^2;
sigma_force_R1_C3 = sqrt(sigma_force_R1_3);


%Reaction Force #2 
sigma_force_R2_1 = ((1)/(N-1)) .* (mean(case_2_load_sensor_1) - case_2_load_sensor_1).^2;
sigma_force_R2_C1 = sqrt(sigma_force_R2_1);

sigma_force_R2_2 = ((1)/(N-1)) .* (mean(case_2_load_sensor_2) - case_2_load_sensor_2).^2;
sigma_force_R2_C2 = sqrt(sigma_force_R2_2);

sigma_force_R2_3 = ((1)/(N-1)) .* (mean(case_2_load_sensor_3) - case_2_load_sensor_3).^2;
sigma_force_R2_C3 = sqrt(sigma_force_R2_3);

%Reaction Force #3
sigma_force_R3_1 = ((1)/(N-1)) .* (mean(case_3_load_sensor_1) - case_3_load_sensor_1).^2;
sigma_force_R3_C1 = sqrt(sigma_force_R3_1);

N3 = 70;
sigma_force_R3_2 = ((1)/(N3-1)) .* (mean(case_3_load_sensor_2) - case_3_load_sensor_2).^2;
sigma_force_R3_C2 = sqrt(sigma_force_R3_2);

sigma_force_R3_3 = ((1)/(N3-1)) .* (mean(case_3_load_sensor_3) - case_3_load_sensor_3).^2;
sigma_force_R3_C3 = sqrt(sigma_force_R3_3);

%Internal Forces 
Sigma_internal_force_1 = ((1)/(N-1)) .* (mean(case_1_inline_load_cells) - case_1_inline_load_cells).^2;
Sigma_internal_force_C1 = sqrt(Sigma_internal_force_1);

Sigma_internal_force_2 = ((1)/(N-1)) .* (mean(case_2_inline_load_cells) - case_2_inline_load_cells).^2;
Sigma_internal_force_C2 = sqrt(Sigma_internal_force_2);

Sigma_internal_force_3 = ((1)/(N3-1)) .* (mean(case_3_inline_load_cells) - case_3_inline_load_cells).^2;
Sigma_internal_force_C3 = sqrt(Sigma_internal_force_3);

%Displacement
Sigma_displacement_1 = ((1)/(N-1)) .* (mean(case_1_LVDT) - case_1_LVDT).^2;
Sigma_displacement_C1 = sqrt(Sigma_displacement_1);

Sigma_displacement_2 = ((1)/(N-1)) .* (mean(case_2_LVDT) - case_2_LVDT).^2;
Sigma_displacement_C2 = sqrt(Sigma_displacement_2);

Sigma_displacement_3 = ((1)/(N3-1)) .* (mean(case_3_LVDT) - case_3_LVDT).^2;
Sigma_displacement_C3 = sqrt(Sigma_displacement_3);
 

%% Uncertainty Analysis


%% R^2

% SSR Values
SSR_R1_1 = 0; % Case 1
SSR_R2_1 = 0;
SSR_R3_1 = 0;
SSR_I_1 = 0;
SSR_D_1 = 0;

SSR_R1_2 = 0; % Case 2
SSR_R2_2 = 0;
SSR_R3_2 = 0;
SSR_I_2 = 0;
SSR_D_2 = 0;

SSR_R1_3 = 0; % Case 3
SSR_R2_3 = 0;
SSR_R3_3 = 0;
SSR_I_3 = 0;
SSR_D_3 = 0;

% Find sum of the residuals
for i = 1:length(case_1_load_sensor_1) % Case 1 and 2

    SSR_R1_1 = SSR_R1_1 + ((case_1_load_sensor_1(i) - fit_R1_1(i)).^2);
    SSR_R2_1 = SSR_R2_1 + ((case_1_load_sensor_2(i) - fit_R2_1(i)).^2);
    SSR_R3_1 = SSR_R3_1 + ((case_1_load_sensor_3(i) - fit_R3_1(i)).^2);
    SSR_I_1 = SSR_I_1 + ((case_1_inline_load_cells(i) - fit_I_1(i)).^2);
    SSR_D_1 = SSR_D_1 + ((case_1_LVDT(i) - fit_D_1(i)).^2);

    SSR_R1_2 = SSR_R1_2 + ((case_2_load_sensor_1(i) - fit_R1_2(i)).^2);
    SSR_R2_2 = SSR_R2_2 + ((case_2_load_sensor_2(i) - fit_R2_2(i)).^2);
    SSR_R3_2 = SSR_R3_2 + ((case_2_load_sensor_3(i) - fit_R3_2(i)).^2);
    SSR_I_2 = SSR_I_2 + ((case_2_inline_load_cells(i) - fit_I_2(i)).^2);
    SSR_D_2 = SSR_D_2 + ((case_2_LVDT(i) - fit_D_2(i)).^2);

end

for i = 1:length(case_3_load_sensor_1) % Case 3

    SSR_R1_3 = SSR_R1_3 + ((case_3_load_sensor_1(i) - fit_R1_3(i)).^2);
    SSR_R2_3 = SSR_R2_3 + ((case_3_load_sensor_2(i) - fit_R2_3(i)).^2);
    SSR_R3_3 = SSR_R3_3 + ((case_3_load_sensor_3(i) - fit_R3_3(i)).^2);
    SSR_I_3 = SSR_I_3 + ((case_3_inline_load_cells(i) - fit_I_3(i)).^2);
    SSR_D_3 = SSR_D_3 + ((case_3_LVDT(i) - fit_D_3(i)).^2);
end

% SST values
SST_R1_1 = 0; % Case 1
SST_R2_1 = 0;
SST_R3_1 = 0;
SST_I_1 = 0;
SST_D_1 = 0;

SST_R1_2 = 0; % Case 2
SST_R2_2 = 0;
SST_R3_2 = 0;
SST_I_2 = 0;
SST_D_2 = 0;

SST_R1_3 = 0; % Case 3
SST_R2_3 = 0;
SST_R3_3 = 0;
SST_I_3 = 0;
SST_D_3 = 0;

% Mean values
M_R1_1 = mean(case_1_load_sensor_1); % Case 1
M_R2_1 = mean(case_1_load_sensor_2);
M_R3_1 = mean(case_1_load_sensor_3);
M_I_1 = mean(case_1_inline_load_cells);
M_D_1 = mean(case_1_LVDT);

M_R1_2 = mean(case_2_load_sensor_1); % Case 2
M_R2_2 = mean(case_2_load_sensor_2);
M_R3_2 = mean(case_2_load_sensor_3);
M_I_2 = mean(case_2_inline_load_cells);
M_D_2 = mean(case_2_LVDT);

M_R1_3 = mean(case_3_load_sensor_1); % Case 3
M_R2_3 = mean(case_3_load_sensor_2);
M_R3_3 = mean(case_3_load_sensor_3);
M_I_3 = mean(case_3_inline_load_cells);
M_D_3 = mean(case_3_LVDT);

% Finding SST Values
for i = 1:length(case_1_load_sensor_1) % Case 1 and 2

    SST_R1_1 = SST_R1_1 + ((case_1_load_sensor_1(i) - M_R1_1).^2);
    SST_R2_1 = SST_R2_1 + ((case_1_load_sensor_2(i) - M_R2_1).^2);
    SST_R3_1 = SST_R3_1 + ((case_1_load_sensor_3(i) - M_R3_1).^2);
    SST_I_1 = SST_I_1 + ((case_1_inline_load_cells(i) - M_I_1).^2);
    SST_D_1 = SST_D_1 + ((case_1_LVDT(i) - M_D_1).^2);

    SST_R1_2 = SST_R1_2 + ((case_2_load_sensor_1(i) - M_R1_2).^2);
    SST_R2_2 = SST_R2_2 + ((case_2_load_sensor_2(i) - M_R2_2).^2);
    SST_R3_2 = SST_R3_2 + ((case_2_load_sensor_3(i) - M_R3_2).^2);
    SST_I_2 = SST_I_2 + ((case_2_inline_load_cells(i) - M_I_2).^2);
    SST_D_2 = SST_D_2 + ((case_2_LVDT(i) - M_D_2).^2);

end

for i = 1:length(case_3_load_sensor_1) % Case 3

    SST_R1_3 = SST_R1_3 + ((case_3_load_sensor_1(i) - M_R1_3).^2);
    SST_R2_3 = SST_R2_3 + ((case_3_load_sensor_2(i) - M_R2_3).^2);
    SST_R3_3 = SST_R3_3 + ((case_3_load_sensor_3(i) - M_R3_3).^2);
    SST_I_3 = SST_I_3 + ((case_3_inline_load_cells(i) - M_I_3).^2);
    SST_D_3 = SST_D_3 + ((case_3_LVDT(i) - M_D_3).^2);

end

% Finding R^2 Values
case_1_R_Squared(1) = 1 - (SSR_R1_1./SST_R1_1); % Case 1
case_1_R_Squared(2) = 1 - (SSR_R2_1./SST_R2_1);
case_1_R_Squared(3) = 1 - (SSR_R3_1./SST_R3_1);
case_1_R_Squared(4) = 1 - (SSR_I_1./SST_I_1);
case_1_R_Squared(5) = 1 - (SSR_D_1./SST_D_1);

case_2_R_Squared(1) = 1 - (SSR_R1_2./SST_R1_2); % Case 2
case_2_R_Squared(2) = 1 - (SSR_R2_2./SST_R2_2);
case_2_R_Squared(3) = 1 - (SSR_R3_2./SST_R3_2);
case_2_R_Squared(4) = 1 - (SSR_I_2./SST_I_2);
case_2_R_Squared(5) = 1 - (SSR_D_2./SST_D_2);

case_3_R_Squared(1) = 1 - (SSR_R1_3./SST_R1_3); % Case 3
case_3_R_Squared(2) = 1 - (SSR_R2_3./SST_R2_3);
case_3_R_Squared(3) = 1 - (SSR_R3_3./SST_R3_3);
case_3_R_Squared(4) = 1 - (SSR_I_3./SST_I_3);
case_3_R_Squared(5) = 1 - (SSR_D_3./SST_D_3);


%% Plots

figure(1) % Reaction Force 1
hold on

%scatter(case_1_applied_force, case_1_load_sensor_1,"LineWidth", 1.2) %Blue 
errorbar(case_1_applied_force,case_1_load_sensor_1,ones(size(sigma_force_R1_C1)))
%scatter(case_2_applied_force, case_2_load_sensor_1,"LineWidth", 1.2) %Red 
errorbar(case_2_applied_force,case_2_load_sensor_1,sigma_force_R2_C1)
%scatter(case_3_applied_force, case_3_load_sensor_1,"LineWidth", 1.2) %Yellow 
errorbar(case_3_applied_force,case_3_load_sensor_1,sigma_force_R3_C1)

plot(case_1_applied_force, fit_R1_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
plot(case_2_applied_force, fit_R1_2,"LineWidth", 1.2, "Color", [0.8500 0.3250 0.0980])
plot(case_3_applied_force, fit_R1_3,"LineWidth", 1.2, "Color", [0.9290 0.6940 0.1250])



xlabel("Applied Force (N)")
ylabel("Measured Reaction Force (N)")
title("Reaction Force 1 vs Applied Force")
legend('Case 1', 'Case 2', 'Case 3')
hold off


figure(2) % Reaction Force 2
hold on

%scatter(case_1_applied_force, case_1_load_sensor_2,"LineWidth", 1.2)
errorbar(case_1_applied_force,case_1_load_sensor_2,ones(size(sigma_force_R1_C2)))
%errorbar(case_1_applied_force,case_1_load_sensor_2,sigma_force_R3_C3)
%scatter(case_2_applied_force, case_2_load_sensor_2,"LineWidth", 1.2)
errorbar(case_2_applied_force,case_2_load_sensor_2,ones(size(sigma_force_R2_C2)))
%scatter(case_3_applied_force, case_3_load_sensor_2,"LineWidth", 1.2)
errorbar(case_3_applied_force,case_3_load_sensor_2,ones(size(sigma_force_R3_C2)))

plot(case_1_applied_force, fit_R2_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
plot(case_2_applied_force, fit_R2_2,"LineWidth", 1.2, "Color", [0.8500 0.3250 0.0980])
plot(case_3_applied_force, fit_R2_3,"LineWidth", 1.2, "Color", [0.9290 0.6940 0.1250])

xlabel("Applied Force (N)")
ylabel("Measured Reaction Force (N)")
title("Reaction Force 2 vs Applied Force")
legend('Case 1', 'Case 2', 'Case 3')
hold off

figure(3) % Reaction Force 3

hold on
%scatter(case_1_applied_force, case_1_load_sensor_3,"LineWidth", 1.2)
errorbar(case_1_applied_force,case_1_load_sensor_3,ones(size(sigma_force_R1_C3)))
%scatter(case_2_applied_force, case_2_load_sensor_3,"LineWidth", 1.2)
errorbar(case_2_applied_force,case_2_load_sensor_3,ones(size(sigma_force_R2_C3)))
%scatter(case_3_applied_force, case_3_load_sensor_3,"LineWidth", 1.2)
errorbar(case_3_applied_force,case_3_load_sensor_3,ones(size(sigma_force_R3_C3)))

plot(case_1_applied_force, fit_R3_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
plot(case_2_applied_force, fit_R3_2,"LineWidth", 1.2, "Color", [0.8500 0.3250 0.0980])
plot(case_3_applied_force, fit_R3_3,"LineWidth", 1.2, "Color", [0.9290 0.6940 0.1250])

xlabel("Applied Force (N)")
ylabel("Measured Reaction Force (N)")
title("Reaction Force 3 vs Applied Force")
legend('Case 1', 'Case 2', 'Case 3')
hold off

figure(4) % Internal Force

hold on
%scatter(case_1_applied_force, case_1_inline_load_cells,"LineWidth", 1.2)
errorbar(case_1_applied_force,case_1_inline_load_cells,ones(size(Sigma_internal_force_C1)))
%scatter(case_2_applied_force, case_2_inline_load_cells,"LineWidth", 1.2)
errorbar(case_2_applied_force,case_2_inline_load_cells,ones(size(Sigma_internal_force_C2)))
%scatter(case_3_applied_force, case_3_inline_load_cells,"LineWidth", 1.2)
errorbar(case_3_applied_force,case_3_inline_load_cells,ones(size(Sigma_internal_force_C3)))

plot(case_1_applied_force, fit_I_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
plot(case_2_applied_force, fit_I_2,"LineWidth", 1.2, "Color", [0.8500 0.3250 0.0980])
plot(case_3_applied_force, fit_I_3,"LineWidth", 1.2, "Color", [0.9290 0.6940 0.1250])

xlabel("Applied Force (N)")
ylabel("Measured Internal Force (N)")
title("Internal Force vs Applied Force")
legend('Case 1', 'Case 2', 'Case 3')
hold off

figure(5) % Displacement

hold on
% scatter(case_1_applied_force, case_1_LVDT,"LineWidth", 1.2)
errorbar(case_1_applied_force,case_1_LVDT,ones(size(Sigma_displacement_C1)))
% scatter(case_2_applied_force, case_2_LVDT,"LineWidth", 1.2)
errorbar(case_2_applied_force,case_2_LVDT,ones(size(Sigma_displacement_C2)))
% scatter(case_3_applied_force, case_3_LVDT,"LineWidth", 1.2)
errorbar(case_3_applied_force,case_3_LVDT,ones(size(Sigma_displacement_C3)))
% xlim([0 300])
% ylim([-0.2 1.6])

plot(case_1_applied_force, fit_D_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
plot(case_2_applied_force, fit_D_2,"LineWidth", 1.2, "Color", [0.8500 0.3250 0.0980])
plot(case_3_applied_force, fit_D_3,"LineWidth", 1.2, "Color", [0.9290 0.6940 0.1250])

xlabel("Applied Force (N)")
ylabel("Measured Mid-Point Displacement (mm)")
title("Measured Mid-Point Displacement vs Applied Force")
legend('Case 1', 'Case 2', 'Case 3')
hold off

%% Plots Question 2

P1 = 0:44.48:222.4;
x = 2;
% Reaction Force 1

figure(6)
subplot(3,1,1)
title('Reactions Forces at F0')
hold on
plot(case_1_applied_force, fit_R1_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
grid on; grid minor
plot(P1,P1/4)
plot(222.4,55.56,'ok')
xlabel('Applied Force [N]')
ylabel('F0 Reaction Force [N]')
hold off
subplot(3,1,2)
title('Reactions Forces at F1')
hold on
plot(case_1_applied_force, fit_R2_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
grid on; grid minor
plot(P1,P1/4)
plot(222.4,55.63,'ok')
xlabel('Applied Force [N]')
ylabel('F1 Reaction Force [N]')
hold off
subplot(3,1,3)
title('Reactions Forces at F2')
hold on
plot(case_1_applied_force, fit_R3_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
grid on; grid minor
plot(P1,P1/2)
plot(222.4,111.19,'ok')
xlabel('Applied Force [N]')
ylabel('F2 Reaction Force [N]')
Lgnd = legend('Experimental','Analytical','ANSYS','Location','bestoutside');
Lgnd.Position(1) = 0.75;
Lgnd.Position(2) = 0.45;
hold off
print('rf comparison','-r300','-dpng')

% Internal Forces

figure(7)
hold on
grid on; grid minor
plot(case_1_applied_force, fit_I_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])
plot(P1,P1)
plot(222.4,450.14,'ok')
xlabel('Applied Load [N]')
ylabel('Displacement [mm]')
title('Comparison of Internal Forces For All 3 Models ')
legend('Experimental','Analytical','ANSYS','Location','northwest')
hold off
print('internal comparison','-r300','-dpng')

% Displacement

figure(8)
hold on
plot(case_1_applied_force, fit_D_1,"LineWidth", 1.2, "Color", [0 0.4470 0.7410])

v_x_case1 = abs((1/(E*I) * (P1/12 * x.^3 - P1.*x)) .* 1000);
grid on; grid minor
plot(P1,v_x_case1)
plot(222.4,1.85,'ok')
xlabel('Applied Load [N]')
ylabel('Displacement [mm]')
title('Comparison of Displacement For All 3 Models ')
legend('Experimental','Analytical','ANSYS','Location','northwest')
hold off
print('displacement comparison','-r300','-dpng')

%% Plots Question 3

% First Case
figure(9)
subplot(2,1,1)
plot(case_2_applied_force(1:60,:), case_2_inline_load_cells(1:60,:),"LineWidth", 1.2)
title("Internal Force vs Applied Load")
xlabel("Applied Load [N]")
ylabel("Internal Force [N]")
grid on; grid minor
subplot(2,1,2)
plot(case_2_applied_force(1:60,:), case_2_LVDT(1:60,:),"LineWidth", 1.2)
title("Displacement vs Applied Load")
xlabel("Applied Load [N]")
ylabel("Displacement [mm]")
grid on; grid minor
print('case 2','-r300','-dpng')
% Second Case
figure(10)
subplot(2,1,1)
plot(case_3_applied_force(1:40,:), case_3_inline_load_cells(1:40,:),"LineWidth", 1.2)
title("Internal Force vs Applied Load")
xlabel("Applied Load [N]")
grid on; grid minor
ylabel("Internal Force [N]")
subplot(2,1,2)
plot(case_3_applied_force(1:40,:), case_3_LVDT(1:40,:),"LineWidth", 1.2)
title("Displacement vs Applied Load")
xlabel("Applied Load [N]")
ylabel("Displacement [mm]")
grid on; grid minor
print('case 3','-r300','-dpng')




