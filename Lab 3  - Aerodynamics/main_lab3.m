%% HK
clc;
close all;
clear;

%% Code

alpha_0006 = 10;    % degrees


[m_0006, p_0006, t_0006] = NACAdata('0006');
[m_0012, p_0012, t_0012] = NACAdata('0012');
[m_0018, p_0018, t_0018] = NACAdata('0018');
[m_2412, p_2412, t_2412] = NACAdata('2412');
[m_4412, p_4412, t_4412] = NACAdata('4412');

N = 25:35:500;    % Iteration number
c = 1;
allow = 0.01;
error = 1;
c_l_0006_actual = 1.145;        % Number retreived when we input N=5000 --> we assume this is the asymptotic c_l
cl_tolerance = 0;
c_l_0006 = ones(1,length(N));
TAT_slope = 2*pi^2 / 180;

for i=1:length(N)

[x_b_0006, y_b_0006] = NACA_Airfoils(m_0006,p_0006,t_0006,c,N(i));          % Finding x and y of airfoil for a given # of Panels (Problem 1)
[c_l_0006(i)] = Vortex_Panel(x_b_0006,y_b_0006,alpha_0006);                 % Finding Cl given certain array of panels

error = 1 - (c_l_0006(i) / c_l_0006_actual);

if(error < allow) && cl_tolerance == 0                      % Checking if error is below allowable
    cl_tolerance = c_l_0006(i);                             % taking the values of n and cl that give us tolerable error
    n_tolerance = N(i);
end
    
end

alpha_vary = linspace(-5,10,30);  % varying alpha for problem 2
c_l_0006_2 = ones(1,length(alpha_vary));
c_l_0012 = ones(1,length(alpha_vary));
c_l_0018 = ones(1,length(alpha_vary));      % Preallocating for speed
c_l_2412 = ones(1,length(alpha_vary));
c_l_4412 = ones(1,length(alpha_vary));
TAT_y = zeros(1,length(alpha_vary));

for i=1:length(alpha_vary)

[x_b_0006_2, y_b_0006_2] = NACA_Airfoils(m_0006,p_0006,t_0006,c,n_tolerance);        % Finding x and y of airfoil using pervious n panels (Problem 2)
[c_l_0006_2(i)] = Vortex_Panel(x_b_0006_2,y_b_0006_2,alpha_vary(i));                 % Finding Cl given certain alpha

[x_b_0012, y_b_0012] = NACA_Airfoils(m_0012,p_0012,t_0012,c,n_tolerance);        % Finding x and y of airfoil using pervious n panels (Problem 2)
[c_l_0012(i)] = Vortex_Panel(x_b_0012,y_b_0012,alpha_vary(i));

[x_b_0018, y_b_0018] = NACA_Airfoils(m_0018,p_0018,t_0018,c,n_tolerance);        % Finding x and y of airfoil using pervious n panels (Problem 2)
[c_l_0018(i)] = Vortex_Panel(x_b_0018,y_b_0018,alpha_vary(i));

[x_b_2412, y_b_2412] = NACA_Airfoils(m_2412,p_2412,t_2412,c,n_tolerance);        % Finding x and y of airfoil using pervious n panels (Problem 3)
[c_l_2412(i)] = Vortex_Panel(x_b_2412,y_b_2412,alpha_vary(i));

[x_b_4412, y_b_4412] = NACA_Airfoils(m_4412,p_4412,t_4412,c,n_tolerance);        % Finding x and y of airfoil using pervious n panels (Problem 3)
[c_l_4412(i)] = Vortex_Panel(x_b_4412,y_b_4412,alpha_vary(i));

TAT_y(i) = TAT_slope * alpha_vary(i);

% Symmetric Airfoils (Problem 2 and 3)
if(abs(c_l_0006_2(i)) < 0.05)
   alpha_L0_0006_2 = alpha_vary(i);
end

if(abs(c_l_0012(i)) < 0.05)
   alpha_L0_0012 = alpha_vary(i);       % finding zero lift angle of attack
end

if(abs(c_l_0018(i)) < 0.05)
   alpha_L0_0018 = alpha_vary(i);
end
% Cambered airfoils (Problem 3)
if(abs(c_l_2412(i)) < 0.05)
   alpha_L0_2412 = alpha_vary(i);       % finding zero lift angle of attack
end

if(abs(c_l_4412(i)) < 0.05)
   alpha_L0_4412 = alpha_vary(i);
end

end

% Lift Slops (Problems 2 and 3)
p_cl0006_2 = polyfit(alpha_vary,c_l_0006_2,1);
p_cl0012 = polyfit(alpha_vary,c_l_0012,1);
p_cl0018 = polyfit(alpha_vary,c_l_0018,1);
p_cl2412 = polyfit(alpha_vary,c_l_2412,1);
p_cl4412 = polyfit(alpha_vary,c_l_4412,1);
p_TAT = polyfit(alpha_vary,TAT_y,1);
a_0_0006 = p_cl0006_2(1);
a_0_0012 = p_cl0012(1);
a_0_0018 = p_cl0018(1);
a_0_2412 = p_cl2412(1);
a_0_4412 = p_cl4412(1);
a_0_TAT = p_TAT(1);

% TAT For Cambered Airfoils
dzdx1_2412 = @(x) (m_2412 * (2*p_2412 -1 + cos(x)))/p_2412^2;
dzdx2_2412 = @(x) (m_2412 * (2*p_2412 -1 + cos(x)))/(1-p_2412)^2;
dzdx1_4412 = @(x) (m_4412 * (2*p_4412 -1 + cos(x)))/p_4412^2;
dzdx2_4412 = @(x) (m_4412 * (2*p_4412 -1 + cos(x)))/(1-p_4412)^2;

x_theta = @(x) (cos(x) -1);

int1_2412 = @(x) dzdx1_2412(x).*x_theta(x);
int2_2412 = @(x) dzdx2_2412(x).*x_theta(x);
int1_4412 = @(x) dzdx1_4412(x).*x_theta(x);
int2_4412 = @(x) dzdx2_4412(x).*x_theta(x);

integ_2412 = integral(int1_2412,0,acos(1-2*p_2412)) + integral(int2_2412,acos(1-2*p_2412),pi);
integ_4412 = integral(int1_4412,0,acos(1-2*p_4412)) + integral(int2_4412,acos(1-2*p_2412),pi);

alpha_L0_TAT_2412 = -(1/pi) * integ_2412*(180/pi);
alpha_L0_TAT_4412 = -(1/pi) * integ_4412*(180/pi);

TAT_2412 = (2*pi^2)/180 * (alpha_vary - alpha_L0_TAT_2412);
TAT_4412 = (2*pi^2)/180 * (alpha_vary - alpha_L0_TAT_4412);

p_cl2412_TAT = polyfit(alpha_vary,TAT_2412,1);
p_cl4412_TAT = polyfit(alpha_vary,TAT_4412,1);

a_0_2412_TAT = p_cl2412_TAT(1);
a_0_4412_TAT = p_cl4412_TAT(1);

% Theory of Wings Sections from Abbot and von Doenhoff Data (TWS)
TWS_0006x = linspace(-4.033,10.04,10);
TWS_0006y = linspace(-0.4165,0.919,10);
TWS_0012x = linspace(-4.05,9.97,10);
TWS_0012y = linspace(-0.4089,0.9355,10);
TWS_2412x = linspace(-4.0625,10,10);
TWS_2412y = linspace(-0.1875,1.1625,10);
TWS_4412x = linspace(-3.9626,10.093,10);
TWS_4412y = linspace(0,1.45,10);

% Problem 4 functions calls and variables
b=1;
a0_t = 2*pi;
a0_r = 2*pi;
c_t = 1;
c_r = 5;
aero_t = 0;
aero_r = 0;
geo_t = 4;
geo_r = 4;
N2 = 50;

% [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N2);


%% CMD Line Prints
% Problem 1
fprintf('Problem 1: \n')
fprintf('\n')
fprintf("NACA 0006 Cl at 99 percent = %d \n", cl_tolerance)
fprintf("NACA 0006 - Number of Panels at 99 percent Cl = %d \n", n_tolerance)
fprintf('\n')

% Problem 2
fprintf('Problem 2: \n')
fprintf('\n')
fprintf("NACA 0006 lift slope = %d \n",a_0_0006)
fprintf("NACA 0006 zero lift AOA = %d \n",alpha_L0_0006_2)
fprintf('\n')
fprintf("NACA 0012 lift slope = %d \n",a_0_0012)
fprintf("NACA 0012 zero lift AOA = %d \n",alpha_L0_0012)
fprintf('\n')
fprintf("NACA 0018 lift slope = %d \n",a_0_0018)
fprintf("NACA 0018 zero lift AOA = %d \n",alpha_L0_0018)
fprintf('\n')
fprintf("Symmetric TAT lift slope = %d \n",a_0_TAT)
fprintf("Symmetric TAT zero lift AOA = %d \n",0)
fprintf('\n')

% Problem 3
fprintf('Problem 3: \n')
fprintf('\n')
fprintf("NACA 0012 lift slope = %d \n",a_0_0012)
fprintf("NACA 0012 zero lift AOA = %d \n",alpha_L0_0012)
fprintf("TAT NACA 0012 lift slope = %d \n",a_0_TAT)
fprintf("TAT NACA 0012 zero lift AOA = %d \n",0)
fprintf('\n')
fprintf("NACA 2412 lift slope = %d \n",a_0_2412)
fprintf("NACA 2412 zero lift AOA = %d \n",alpha_L0_2412)
fprintf("TAT NACA 2412 lift slope = %d \n",a_0_2412_TAT)
fprintf("TAT NACA 2412 zero lift AOA = %d \n",alpha_L0_TAT_2412)
fprintf('\n')
fprintf("NACA 4412 lift slope = %d \n",a_0_4412)
fprintf("NACA 4412 zero lift AOA = %d \n",alpha_L0_4412)
fprintf("TAT NACA 4412 lift slope = %d \n",a_0_4412_TAT)
fprintf("TAT NACA 4412 zero lift AOA = %d \n",alpha_L0_TAT_4412)
fprintf('\n')

%% Plotting

% Problem 1
figure(1)          
plot(N,c_l_0006)
hold on
yline(c_l_0006_actual,'k')
yline(cl_tolerance,'--r')
hold off
xlabel("Number of Panels")
ylabel("C_l")
title("C_l vs. Number of Panels (Problem 1)")
legend('C_l','Actual C_l','C_l Tolerance','Location','southeast')
ylim([1.105 1.15])
grid on

% Problem 2
figure(2)
plot(alpha_vary,c_l_0006_2,'r')
hold on
plot(alpha_vary, c_l_0012,'b')
plot(alpha_vary, c_l_0018,'g')

plot(TWS_0006x,TWS_0006y,'o--r')
plot(TWS_0012x,TWS_0012y,'o--b')

plot(alpha_vary, TAT_y,'k')

ylabel("C_l of Airfoil")
xlabel("\alpha (degrees)")
title("C_l vs \alpha (Problem 2)")
legend('0006', '0012', '0018','NACA 0006 Data','NACA 0012 Data', 'TAT','Location','northwest')
grid on

% Problem 3
figure(3)
plot(alpha_vary, c_l_0012,'r')
hold on
plot(alpha_vary,c_l_2412,'b')
plot(alpha_vary, c_l_4412,'g')

plot(TWS_0012x,TWS_0012y,'o--r')
plot(TWS_2412x,TWS_2412y,'o--b')
plot(TWS_4412x,TWS_4412y,'o--g')

plot(alpha_vary,TAT_y,'-.r')
plot(alpha_vary,TAT_2412,'-.b')
plot(alpha_vary,TAT_4412,'-.g')

ylabel("C_l of Airfoil")
xlabel("\alpha (degrees)")
title("C_l vs \alpha (Problem 3)")
legend('0012', '2412', '4412','NACA 0012 Data','NACA 2412 Data','NACA 4412 Data' ,'TAT 0012','TAT 2412','TAT 4412','Location','northwest')
grid on


function [x_b, y_b] = NACA_Airfoils(m,p,t,c,N)

x = linspace(c,0,N);  % Starting at TE going to LE

y_t = (t*c / 0.2) * (0.2969.*sqrt(x/c) - 0.126.*(x/c) - 0.3516.*((x/c).^2) + 0.2843.*((x/c).^3) - 0.1036.*((x/c).^4)); 

% preallocate
y_c = zeros(1,length(x));
dy_c = zeros(1,length(x));
x_U = zeros(1,length(x));
x_L = zeros(1,length(x));
y_U = zeros(1,length(x));
y_L = zeros(1,length(x));

    for i=1:length(x)
    if x(i) <= p*c

    y_c(i) = m*(x(i)/p^2)*(2*p - x(i)/c);
    dy_c(i) = -2*m * (x(i) - c*p) / (c*p^2);

    elseif x(i) > p*c
    
    y_c(i) = m*((c - x(i)) / (1-p)^2) * (1 + x(i)/c - 2*p); 
    dy_c(i) = -2*m*(x(i) - c*p) / (c * (p-1)^2);

    end

    squiggly = atan(dy_c);
    x_U(i) = x(i) - y_t(i)*sin(squiggly(i));
    x_L(i) = x(i) + y_t(i)*sin(squiggly(i));
    y_U(i) = y_c(i) + y_t(i)*cos(squiggly(i));
    y_L(i) = y_c(i) - y_t(i)*cos(squiggly(i));
    end

x_U = fliplr(x_U);  % Flips vectors to make sure we are going clockwise
y_U = fliplr(y_U);
 
x_b = [x_L , x_U(2:end)];       % (2:end) so that we do not duplicated leading edge value
y_b = [y_L , y_U(2:end)];

x_b(isnan(x_b)) = 0;
y_b(isnan(y_b)) = 0;

end

function [m, p, t] = NACAdata(str)

m = str(1);
m = str2double(m) / 100;

p = str(2);
p = str2double(p) / 10;

t1 = str(3);
t2 = str(4);

t = strcat(t1, t2);
t = str2double(t) / 100;

end

function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

% Wing geometry
S = 0.5 * (c_t + c_r) * b;
AR = b^2/S;

top = pi/(2*N);
theta = linspace(top,pi/2,N);        % Making theta vec (left half of span)
y = -b/2 * cos(theta);              % Getting y vec

% preallocate
c = zeros(length(theta),1);
a0 = zeros(length(theta),1);
a_L0 = zeros(length(theta),1);
a_geo = zeros(length(theta),1);
alpha_vec = zeros(N,1);
A_vec = zeros(N,N);

for i = 1:length(theta)
    c(i) = c_r + y(i) * (c_t - c_r) / (y(1) - y(end)); %linear interpolation
    a0(i) = a0_r + y(i) * (a0_t - a0_r) / (y(1) - y(end));
    a_L0(i) = aero_r + y(i) * (aero_t - aero_r) / (y(1) - y(end));
    a_geo(i) = geo_r + y(i) * (geo_t - geo_r) / (y(1) - y(end));
    alpha_vec(i) = a_geo(i) - a_L0(i); % Vector math stuff start
    for j = 1:length(theta)
        A_vec(i,j) = 4*b / (a0(i)*c(i)) * sin((2*j-1)*theta(i)) + (2*j-1) * sin((2*j-1)*theta(i))/sin(theta(i));
    end    
end

A_n = A_vec\alpha_vec;



end


function [CL] = Vortex_Panel(XB,YB,ALPHA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                           %
%                                  %
% XB  = Boundary Points x-location %
% YB  = Boundary Points y-location %
% ALPHA = AOA in degrees           %
%                                  %
% Output:                          %
%                                  %
% CL = Sectional Lift Coefficient  %
% improves efficiency by preallocating matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Convert to Radians %
%%%%%%%%%%%%%%%%%%%%%%

ALPHA = ALPHA*pi/180;

%%%%%%%%%%%%%%%%%%%%%
% Compute the Chord %
%%%%%%%%%%%%%%%%%%%%%

CHORD = max(XB)-min(XB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Number of Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = max(size(XB,1),size(XB,2))-1;
MP1 = M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate Matrices for Efficiency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(1,M);
Y = zeros(1,M);
S = zeros(1,M);
THETA = zeros(1,M);
SINE = zeros(1,M);
COSINE = zeros(1,M);
RHS = zeros(1,M);
CN1 = zeros(M);
CN2 = zeros(M);
CT1 = zeros(M);
CT2 = zeros(M);
AN = zeros(M);
AT = zeros(M);
V = zeros(1,M);
CP = zeros(1,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-Panel Relationships:                                  %
%                                                             %
% Determine the Control Points, Panel Sizes, and Panel Angles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    IP1 = I+1;
    X(I) = 0.5*(XB(I)+XB(IP1));
    Y(I) = 0.5*(YB(I)+YB(IP1));
    S(I) = sqrt( (XB(IP1)-XB(I))^2 +( YB(IP1)-YB(I))^2 );
    THETA(I) = atan2( YB(IP1)-YB(I), XB(IP1)-XB(I) );
    SINE(I) = sin( THETA(I) );
    COSINE(I) = cos( THETA(I) );
    RHS(I) = sin( THETA(I)-ALPHA );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:             %
%                                        %
% Determine the Integrals between Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    for J = 1:M
        if I == J
            CN1(I,J) = -1.0;
            CN2(I,J) = 1.0;
            CT1(I,J) = 0.5*pi;
            CT2(I,J) = 0.5*pi;
        else
            A = -(X(I)-XB(J))*COSINE(J) - (Y(I)-YB(J))*SINE(J);
            B = (X(I)-XB(J))^2 + (Y(I)-YB(J))^2;
            C = sin( THETA(I)-THETA(J) );
            D = cos( THETA(I)-THETA(J) );
            E = (X(I)-XB(J))*SINE(J) - (Y(I)-YB(J))*COSINE(J);
            F = log( 1.0 + S(J)*(S(J)+2*A)/B );
            G = atan2( E*S(J), B+A*S(J) );
            P = (X(I)-XB(J)) * sin( THETA(I) - 2*THETA(J) ) ...
              + (Y(I)-YB(J)) * cos( THETA(I) - 2*THETA(J) );
            Q = (X(I)-XB(J)) * cos( THETA(I) - 2*THETA(J) ) ...
              - (Y(I)-YB(J)) * sin( THETA(I) - 2*THETA(J) );
            CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C+D*E)*G/S(J);
            CN1(I,J) = 0.5*D*F + C*G - CN2(I,J);
            CT2(I,J) = C + 0.5*P*F/S(J) + (A*D-C*E)*G/S(J);
            CT1(I,J) = 0.5*C*F - D*G - CT2(I,J);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:           %
%                                      %
% Determine the Influence Coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    AN(I,1) = CN1(I,1);
    AN(I,MP1) = CN2(I,M);
    AT(I,1) = CT1(I,1);
    AT(I,MP1) = CT2(I,M);
    for J = 2:M
        AN(I,J) = CN1(I,J) + CN2(I,J-1);
        AT(I,J) = CT1(I,J) + CT2(I,J-1);
    end
end
AN(MP1,1) = 1.0;
AN(MP1,MP1) = 1.0;
for J = 2:M
    AN(MP1,J) = 0.0;
end
RHS(MP1) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the gammas %
%%%%%%%%%%%%%%%%%%%%%%%%

GAMA = AN\RHS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Tangential Veloity and Coefficient of Pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    V(I) = cos( THETA(I)-ALPHA );
    for J = 1:MP1
        V(I) = V(I) + AT(I,J)*GAMA(J);
    end
    CP(I) = 1.0 - V(I)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Sectional Coefficient of Lift %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIRCULATION = sum(S.*V);
CL = 2*CIRCULATION/CHORD;

end











