
%% Housekeeping
clc; clear; close all;
%Housekeeping

%% TESTING




%% Question 2

%Material Properties Sourced from Onlinemetals.com and Wikipedia

%Alpha Parameters
% Al 7075, Brass C360, SS T303
Rho = 1E-6*[2810 8500 8000]; %[kg/cm^3]
Cp = [960 380 500]; %[J/kg-K]
K = 1E-2*[130 115 16.2]; %[W/cm-K]


%Calculating alpha systematically
alphaFn = @(rho,cp,k) k/(rho*cp); %[cm^2/s]

alpha = zeros(length(K),1);

for i = 1:length(K)
alpha(i,1) = alphaFn(Rho(i),Cp(i),K(i));
end

%alpha = .4*alpha;
%% Question 3  - Experimental Steady State Solution [ESSS] (HESSS)

%{
Find the temperature 𝑇0 at position 𝑥0 where 𝑥0 is located 1
3/8
inches to the left of the first thermocouple (Th1). To do this, use the 
following experimental data (taken from a sample at steady state) and
reference Figure 2. Also determine 𝐻, the slope of the steady state
distribution. Note the thermocouple measurements are located 0.5 in apart
as in Figure 2. Remember that a steady state distribution is quasi-linear
so a least-squares fit of the temperature data can be used to determine 
the temperature slope, as it relates to the distance between the
thermocouples (Th1-Th8) and 𝑥0. Compare the temperature readings with the
calculated best-fit line.
%}

%%Importing Data
Q3_Temperature = [17.60,21.61,25.13,29.22,34.92,38.10,45.21,47.01];

%% Formatting Imported Data

%Format: Aluminum 25V -> AL for Aluminum, 25 for 25 Volts <AL><28>
AL25 = load('Aluminum_25V_240mA');
AL28 = load('Aluminum_28V_270mA');
BR26 = load('Brass_26V_245mA');
BR29 = load('Brass_29V_273mA');
ST21 = load('Steel_21V_192mA');

%Retrieving end temperature row
AL25_Temperature = AL25(end,2:end);
AL28_Temperature = AL28(end,2:end);
BR26_Temperature = BR26(end,2:end);
BR29_Temperature = BR29(end,2:end);
ST21_Temperature = ST21(end,2:end);

%Metadata TITLES
Q3_MD.P = struct;
AL25_MD.P = struct;
AL28_MD.P = struct;
BR26_MD.P = struct;
BR29_MD.P = struct;
ST21_MD.P = struct;

Q3_MD.T = 'Question 3';
AL25_MD.T = 'Aluminum 25V';
AL28_MD.T = 'Aluminum 28V';
BR26_MD.T = 'Brass 26V';
BR29_MD.T = 'Brass 29V';
ST21_MD.T = 'Steel 21V';

Q3_MD.P = 1;
AL25_MD.P = 2;
AL28_MD.P = 3;
BR26_MD.P = 4;
BR29_MD.P = 5;
ST21_MD.P = 6;


% Plugging in Data for ESSS "Experimental Method"

ESSS_Q3 = ESSSFN(Q3_Temperature,Q3_MD);
ESSS_AL25 = ESSSFN(AL25_Temperature,AL25_MD);
ESSS_AL28 = ESSSFN(AL28_Temperature,AL28_MD);
ESSS_BR26 = ESSSFN(BR26_Temperature,BR26_MD);
ESSS_BR29 = ESSSFN(BR29_Temperature,BR29_MD);
ESSS_ST21 = ESSSFN(ST21_Temperature,ST21_MD);

%% Analytic Method (HANA)

%Convergence Paramters
tsmall = 1; %[s]
tlong = 1000; %[s]

% Rod Parameters
T0Q3 = ESSS_Q3(1); % Base T0 temperature
HQ3 = ESSS_Q3(2); % H thermal diffusivity
x0Heater =  8.89+3.4925+1.27; % [cm] distance between x0 and left end of heater
xTh8 = 8.89+3.4925; % [cm] distance to thermocouple 8 

%Punching in the functions

[HANA_TS,Fo_TS] = uFN(T0Q3,HQ3,alpha(1),tsmall); %small time calculation
[HANA_TL,Fo_TL] = uFN(T0Q3,HQ3,alpha(1),tlong); %big time calculation


%% Question 6 Sensitivity Study on Diffusivity

time = 1:1:1000;

Q6 = zeros(length(K),length(time));
for i = 1 : length(K)
    [Q6(i,:),Fo] = uFN(T0Q3,HQ3,alpha(i),time);
end

%% Results 1: Steady State Distributions

%index:
% AL25, AL28, BR26, BR29, ST21

% Parameters
R1_A = 5.06707; %[cm^2]
R1_K = [K(1),K(1),K(2),K(2),K(3)]; %[W/(cm*k)]
R1_V = [25,28,26,29,21]; %[V]
R1_I = [.240,.270,.245,.273,.192]; %[A]

% Getting Hana based on H = Qdot/(k*A)
QdotFn = @(V,I) (V*I);

Qdot = zeros(1,length(R1_K));
for i = 1 : length(R1_K)
    Qdot(i) = QdotFn(R1_V(i),R1_I(i));
end


HANAFN = @(QDOT,K,A) QDOT/(K*A);

Hana = zeros(1,length(R1_K));
for i = 1 : length(R1_K)
    Hana(1,i) = HANAFN(Qdot(i),R1_K(i),R1_A); %[degK/cm^2-K]
end


%% Results 2: Time-Depedent Temperature Profiles (Model IA & IB)

%index:
% AL25, AL28, BR26, BR29, ST21


% Results 2a: Model 1A
%==============================================================
indexAL28 = 206:length(AL28);
AL28 = AL28(indexAL28,:);
AL28(:,1) = AL28(:,1) - AL28(1,1);


AL25_Time = AL25(:,1)';
AL28_Time = AL28(:,1)';
BR26_Time = BR26(:,1)';
BR29_Time = BR29(:,1)';
ST21_Time = ST21(:,1)';


% Computing Analytical Temperatyre over Time

%function reference:
% function [uVec,Fo] = uFN(T0,H,alpha,t) % base temp "T0", H "H", time span "t"
[R2HAT_AL25,R2Fo_AL25] = uFN(ESSS_AL25(1),Hana(1),alpha(1),AL25_Time);
[R2HAT_AL25,R2Fo_AL25] = uFN(ESSS_AL25(1),Hana(1),alpha(1),AL25_Time);
[R2HAT_AL28,R2Fo_AL28] = uFN(ESSS_AL28(1),Hana(2),alpha(1),AL28_Time);
[R2HAT_BR26,R2Fo_BR26] = uFN(ESSS_BR26(1),Hana(3),alpha(2),BR26_Time);
[R2HAT_BR29,R2Fo_BR29] = uFN(ESSS_BR29(1),Hana(4),alpha(2),BR29_Time);
[R2HAT_ST21,R2Fo_ST21] = uFN(ESSS_ST21(1),Hana(5),alpha(3),ST21_Time);


% Results 2b: Model IB
%===================================

Hexp = [ESSS_AL25(2) ESSS_AL28(2) ESSS_BR26(2) ESSS_BR29(2) ESSS_ST21(2)];

[R2HYB_AL25,R2FoHYB_AL25] = uFN(ESSS_AL25(1),Hexp(1),alpha(1),AL25_Time);
[R2HYB_AL28,R2FoHYB_AL28] = uFN(ESSS_AL28(1),Hexp(2),alpha(1),AL28_Time);
[R2HYB_BR26,R2FoHYB_BR26] = uFN(ESSS_BR26(1),Hexp(3),alpha(2),BR26_Time);
[R2HYB_BR29,R2FoHYB_BR29] = uFN(ESSS_BR29(1),Hexp(4),alpha(2),BR29_Time);
[R2HYB_ST21,R2FoHYB_ST21] = uFN(ESSS_ST21(1),Hexp(5),alpha(3),ST21_Time);


%% Results 3 - Initial State Distributions (Model II)

%fn ref
%Function Handle - [Returns T0, Hexp, Hana] = fn name(Temperature Data)
%function [ESSSCALC] = ESSS_FN(ESSS_Temperature,metadata)


%metadata
AL25_MD.T = 'Initial Temperature Distribution AL 25V';
AL28_MD.T = 'Initial Temperature Distribution AL 28V';
BR26_MD.T = 'Initial Temperature Distribution BR 26V';
BR29_MD.T = 'Initial Temperature Distribution BR 29V';
ST21_MD.T = 'Initial Temperature Distribution ST 21V';

AL25_MD.P = 7;
AL28_MD.P = 8;
BR26_MD.P = 9;
BR29_MD.P = 10;
ST21_MD.P = 11;


% "MIS" M_exp Initial State Distribution"
% "IST" Initial State Temperatures
IST_AL25 = AL25(1,2:end);
IST_AL28 = AL28(1,2:end);
IST_BR26 = BR26(1,2:end);
IST_BR29 = BR29(1,2:end);
IST_ST21 = ST21(1,2:end);

MIS_AL25 = ESSSFN(IST_AL25,AL25_MD);
MIS_AL28 = ESSSFN(IST_AL28,AL28_MD);
MIS_BR26 = ESSSFN(IST_BR26,BR26_MD);
MIS_BR29 = ESSSFN(IST_BR29,BR29_MD);
MIS_ST21 = ESSSFN(IST_ST21,ST21_MD);

%calculating gx for inital slope
GXFN = @(m,h) m-h;

gAL25 = GXFN(MIS_AL25(2),ESSS_AL25(2));
gAL28 = GXFN(MIS_AL28(2),ESSS_AL28(2));
gBR26 = GXFN(MIS_BR26(2),ESSS_BR26(2));
gBR29 = GXFN(MIS_BR29(2),ESSS_BR29(2));
gST21 = GXFN(MIS_ST21(2),ESSS_ST21(2));

%function reference:
%function [uVec,Fo] = MIIFN(T0,H,alpha,t,g) % base temp "T0", H "H", time span "t", g inital temperature span


[R3_AL25,R3Fo_AL25] = MIIFN(ESSS_AL25(1),Hexp(1),alpha(1),AL25_Time,gAL25);
[R3_AL28,R3Fo_AL28] = MIIFN(ESSS_AL28(1),Hexp(2),alpha(1),AL28_Time,gAL28);
[R3_BR26,R3Fo_BR26] = MIIFN(ESSS_BR26(1),Hexp(3),alpha(2),BR26_Time,gBR26);
[R3_BR29,R3Fo_BR29] = MIIFN(ESSS_BR29(1),Hexp(4),alpha(2),BR29_Time,gBR29);
[R3_ST21,R3Fo_ST21] = MIIFN(ESSS_ST21(1),Hexp(5),alpha(3),ST21_Time,gST21);


%% Results 4/5 - Variance in Thermal Diffusivities (Model III)
% 5: Time to Steady State

[tssAL25,foAL25] = TSSFN(R2FoHYB_AL25);
[tssAL28,foAL28] = TSSFN(R2FoHYB_AL28);
[tssBR26,foBR26] = TSSFN(R2FoHYB_BR26);
[tssBR29,foBR29] = TSSFN(R2FoHYB_BR29);
[tssST21,foST21] = TSSFN(R2FoHYB_ST21);


alpTss = @(FO,TSS) FO*14.9225^2/TSS;


aAdjAL25 = alpTss(R2FoHYB_AL25(tssAL25/10),tssAL25);
aAdjAL28 = alpTss(R2FoHYB_AL28(tssAL28/10),tssAL28);
aAdjBR26 = alpTss(R2FoHYB_BR26(tssBR26/10),tssBR26);
aAdjBR29 = alpTss(R2FoHYB_BR29(tssBR29/10),tssBR29);
aAdjST21 = alpTss(R2FoHYB_ST21(tssST21/10),tssST21);

[R4_AL25,R4Fo_AL25] = MIIFN(ESSS_AL25(1),Hexp(1),aAdjAL25,AL25_Time,gAL25);
[R4_AL28,R4Fo_AL28] = MIIFN(ESSS_AL28(1),Hexp(2),aAdjAL28,AL28_Time,gAL28);
[R4_BR26,R4Fo_BR26] = MIIFN(ESSS_BR26(1),Hexp(3),aAdjBR26,BR26_Time,gBR26);
[R4_BR29,R4Fo_BR29] = MIIFN(ESSS_BR29(1),Hexp(4),aAdjBR29,BR29_Time,gBR29);
[R4_ST21,R4Fo_ST21] = MIIFN(ESSS_ST21(1),Hexp(5),aAdjST21,ST21_Time,gST21);


%% PlottingIST_BR26 = BR26(1,2:end);

%▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
%▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
%▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

%% Question 3 Plotting

% IN ESSS Function to be more factorable 

%% Question 5 Plotting
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting the analysis

%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<|
if plotswitch == 1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|

figure(12)
hold on


%Plots convergence
plot(HANA_TL(1:(end-1),9),'c.-',MarkerSize=25);
plot(HANA_TS(1:(end-1),9),'bo-','LineWidth',1.4); %could be 2:end, or 1:(end-1)


%Annotation
title('HANA using Fourier Convergence')
xlabel('Index [n]')
ylabel('TC 8 Temperature [°C]')
legend('TC 8 Temperature @t=1s','TC 8 Temperature @t=1000s',Location='east')
xticks(1:10)
axis padded

xl = xlim; %why is this x/ylim function so convuluted to use????
yl = ylim;

% text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
text(1.5,HANA_TS(1,9)-3, sprintf('n=1,t=1s:\n %.3f [°C]',HANA_TS(1,9)),'FontSize' ...
    ,10)
text(9,HANA_TS(end,9)+4, sprintf('n=10,t=1s:\n %.3f [°C]',HANA_TS(end,9)),'FontSize' ...
    ,10)
text(1.5,HANA_TL(1,9)-9, sprintf('n=1,t=1000s:\n %.3f [°C]',HANA_TL(1,9)),'FontSize' ...
    ,10)
text(8.5,HANA_TL(end,9)-4, sprintf('n=10,t=1000s:\n %.3f [°C]',HANA_TL(end,9)),'FontSize' ...
    ,10)

text(.35*xl(2),.4*yl(2), sprintf('Fourier Number (1s)= %.3f',Fo_TS))
if Fo_TS > 0.2
    text(.35*xl(2),.4*yl(2)-2,'One term sufficient')
else
    text(.35*xl(2),.4*yl(2)-2,'One term insufficient')
end

text(.35*xl(2),.4*yl(2)-6, sprintf('Fourier Number (1000s)= %.3f',Fo_TL))
if Fo_TL > 0.2
    text(.35*xl(2),.4*yl(2)-8,'One term sufficient')
else
    text(.35*xl(2),.4*yl(2)-8,'One term insufficient')
end

print('-f7','PLOT ASEN 3802 - Q5 HANA','-dpng')

hold off
else
end

%% Question 6 Plotting

%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotswitch == 1 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(8)
hold on

%Plots linear fit line
plot(Q6(1,:),'Color',"#F83810",'LineWidth',1.4)
plot(Q6(2,:),'Color',"#F88830",'LineWidth',1.4)
plot(Q6(3,:),'Color',"#F8C030",'LineWidth',1.4)

%Annotation
title('Sensitivity Study on Diffusivity')
xlabel('Time [s]')
ylabel('Temperature [°C]')
legend('Aluminum 7075','Brass C360','SS T303','Location','best')

%text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
%text(ESSS_Distance(1)+.4,ESSS_LSM(1)-.2, sprintf('T0 = %.2f [°C]',ESSS_LSM(1)))
%text(ESSS_Distance(5)+.5,ESSS_LSM(5)-.5, sprintf('H = %.2f [°C/cm]',ESSS_P(1) ))

axis([0 1000 Q6(1,1)*(.8) Q6(1,end)*1.04])
print('-f8','PLOT ASEN 3802 - Q6 TEMP V TIME','-dpng')



hold off
else
end

%% Results 2 : Time-Depedent Temperature Profiles

%Model IA


%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotswitch == 1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(13)
sgtitle('Results 2a - Time Depedent Time Profiles (Model IA)')
hold on

%Plots linear fit line
subplot(2,3,1)
hold on
plot(AL25_Time,R2HAT_AL25,'r')
plot(AL25_Time,AL25(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')



subplot(2,3,2)
hold on
plot(AL28_Time,R2HAT_AL28,'r')
plot(AL28_Time,AL28(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,3)
hold on
plot(BR26_Time,R2HAT_BR26,'r')
plot(BR26_Time,BR26(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')

subplot(2,3,4)
hold on
plot(BR29_Time,R2HAT_BR29,'r')
plot(BR29_Time,BR29(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,5)
hold on
plot(ST21_Time,R2HAT_ST21,'r')
plot(ST21_Time,ST21(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


%Annotation

lgnd = legend('Model Temperature','','','','','','','','','Data Temperature');
lgnd.Position(1) = 0.75;
lgnd.Position(2) = 0.15;

%text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
%text(ESSS_Distance(1)+.4,ESSS_LSM(1)-.2, sprintf('T0 = %.2f [°C]',ESSS_LSM(1)))
%text(ESSS_Distance(5)+.5,ESSS_LSM(5)-.5, sprintf('H = %.2f [°C/cm]',ESSS_P(1) ))

%axis([0 1000 Q6(1,1)*(.8) Q6(1,end)*1.04])
%print('-f1','PLOT ASEN 3802 - R2 TEMP V TIME','-dpng')



hold off
else
end

%% Results 2b : Time-Depedent Temperature Profiles

%Model IB


%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotswitch == 1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(14)
sgtitle('Results 2b - Time Depedent Time Profiles (Model IB)')
hold on

%Plots linear fit line
subplot(2,3,1)
hold on
plot(AL25_Time,R2HYB_AL25,'r')
plot(AL25_Time,AL25(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')



subplot(2,3,2)
hold on
plot(AL28_Time,R2HYB_AL28,'r')
plot(AL28_Time,AL28(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,3)
hold on
plot(BR26_Time,R2HYB_BR26,'r')
plot(BR26_Time,BR26(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')

subplot(2,3,4)
hold on
plot(BR29_Time,R2HYB_BR29,'r')
plot(BR29_Time,BR29(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,5)
hold on
plot(ST21_Time,R2HYB_ST21,'r')
plot(ST21_Time,ST21(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


%Annotation

lgnd = legend('Model Temperature','','','','','','','','','Data Temperature');
lgnd.Position(1) = 0.75;
lgnd.Position(2) = 0.15;

%text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
%text(ESSS_Distance(1)+.4,ESSS_LSM(1)-.2, sprintf('T0 = %.2f [°C]',ESSS_LSM(1)))
%text(ESSS_Distance(5)+.5,ESSS_LSM(5)-.5, sprintf('H = %.2f [°C/cm]',ESSS_P(1) ))

%axis([0 1000 Q6(1,1)*(.8) Q6(1,end)*1.04])
%print('-f1','PLOT ASEN 3802 - R2 TEMP V TIME','-dpng')



hold off
else
end


%% Results 3 : Initial State Distributions

%Model II


%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotswitch == 1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(15)
sgtitle('Results 3 - Initial State Distributions (Model II)')
hold on

%Plots linear fit line
subplot(2,3,1)
hold on
plot(AL25_Time,R3_AL25,'r')
plot(AL25_Time,AL25(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')



subplot(2,3,2)
hold on
plot(AL28_Time,R3_AL28,'r')
plot(AL28_Time,AL28(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,3)
hold on
plot(BR26_Time,R3_BR26,'r')
plot(BR26_Time,BR26(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')

subplot(2,3,4)
hold on
plot(BR29_Time,R3_BR29,'r')
plot(BR29_Time,BR29(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,5)
hold on
plot(ST21_Time,R3_ST21,'r')
plot(ST21_Time,ST21(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


%Annotation

lgnd = legend('Model Temperature','','','','','','','','','Data Temperature');
lgnd.Position(1) = 0.75;
lgnd.Position(2) = 0.15;

%text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
%text(ESSS_Distance(1)+.4,ESSS_LSM(1)-.2, sprintf('T0 = %.2f [°C]',ESSS_LSM(1)))
%text(ESSS_Distance(5)+.5,ESSS_LSM(5)-.5, sprintf('H = %.2f [°C/cm]',ESSS_P(1) ))

%axis([0 1000 Q6(1,1)*(.8) Q6(1,end)*1.04])
%print('-f1','PLOT ASEN 3802 - R2 TEMP V TIME','-dpng')



hold off
else
end

%% Results 4 : Variance in Thermal Diffusivities (Model III)

%Model III


%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotswitch == 1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(16)
sgtitle('Results 4 - Variance in Thermal Diffusivities (Model III)')
hold on

%Plots linear fit line
subplot(2,3,1)
hold on
plot(AL25_Time,R4_AL25,'r')
plot(AL25_Time,AL25(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')



subplot(2,3,2)
hold on
plot(AL28_Time,R4_AL28,'r')
plot(AL28_Time,AL28(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,3)
hold on
plot(BR26_Time,R4_BR26,'r')
plot(BR26_Time,BR26(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')

subplot(2,3,4)
hold on
plot(BR29_Time,R4_BR29,'r')
plot(BR29_Time,BR29(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


subplot(2,3,5)
hold on
plot(ST21_Time,R4_ST21,'r')
plot(ST21_Time,ST21(:,2:end),'b')
hold off
axis padded
xlabel('Time [s]')
ylabel('Temperature [°C]')


%Annotation

lgnd = legend('Model Temperature','','','','','','','','','Data Temperature');
lgnd.Position(1) = 0.7;
lgnd.Position(2) = 0.25;

%text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
%text(ESSS_Distance(1)+.4,ESSS_LSM(1)-.2, sprintf('T0 = %.2f [°C]',ESSS_LSM(1)))
%text(ESSS_Distance(5)+.5,ESSS_LSM(5)-.5, sprintf('H = %.2f [°C/cm]',ESSS_P(1) ))

%axis([0 1000 Q6(1,1)*(.8) Q6(1,end)*1.04])
%print('-f1','PLOT ASEN 3802 - R2 TEMP V TIME','-dpng')



hold off
else
end

%% Report - Hypothesize 


figure(17)
sgtitle('Hypothesis - TC8 Final Temperatures per Model')
hold on

barlabel = categorical({'Data','Model IA','Model IB','Model II','Model III'});
xB = reordercats(barlabel,{'Data','Model IA','Model IB','Model II','Model III'});

HAL25 = [AL25(end,end) R2HAT_AL25(end,end) R2HYB_AL25(end,end) R3_AL25(end,end) R4_AL25(end,end)];
HAL28 = [AL28(end,end) R2HAT_AL28(end,end) R2HYB_AL28(end,end) R3_AL28(end,end) R4_AL28(end,end)];
HBR26 = [BR26(end,end) R2HAT_BR26(end,end) R2HYB_BR26(end,end) R3_BR26(end,end) R4_BR26(end,end)];
HBR29 = [BR29(end,end) R2HAT_BR29(end,end) R2HYB_BR29(end,end) R3_BR29(end,end) R4_BR29(end,end)];
HST21 = [ST21(end,end) R2HAT_ST21(end,end) R2HYB_ST21(end,end) R3_ST21(end,end) R4_ST21(end,end)];

err = 2+zeros(length(HAL25));

%Plots linear fit line
subplot(2,3,1)
hold on
bar(xB,HAL25,'FaceColor',[.88 .24 .05],'EdgeColor',[.9 .05 .05])
errorbar(xB,HAL25,err,err,'Color',[.1 .0002 .04],'LineStyle','none')
hold off
axis padded
ylabel('Temperature [°C]')



subplot(2,3,2)
hold on
bar(xB,HAL28,'FaceColor',[.88 .24 .05],'EdgeColor',[.9 .05 .05])
errorbar(xB,HAL28,err,err,'Color',[.1 .0002 .04],'LineStyle','none')
hold off
axis padded
ylabel('Temperature [°C]')


subplot(2,3,3)
hold on
bar(xB,HBR26,'FaceColor',[.88 .24 .05],'EdgeColor',[.9 .05 .05])
errorbar(xB,HBR26,err,err,'Color',[.1 .0002 .04],'LineStyle','none')
hold off
axis padded
ylabel('Temperature [°C]')

subplot(2,3,4)
hold on
bar(xB,HBR29,'FaceColor',[.88 .24 .05],'EdgeColor',[.9 .05 .05])
errorbar(xB,HBR29,err,err,'Color',[.1 .0002 .04],'LineStyle','none')
hold off
axis padded
ylabel('Temperature [°C]')


subplot(2,3,5)
hold on
bar(xB,HST21,'FaceColor',[.88 .24 .05],'EdgeColor',[.9 .05 .05])
errorbar(xB,HST21,err,err,'Color',[.1 .0002 .04],'LineStyle','none')
hold off
axis padded
ylabel('Temperature [°C]')

hold off


%% Reporting Findings Tabulated

%Steady State Distributions "SSD"

%{
Report your findings by (a) generating plots of the steady state
temperatures at each location according to the experimental and analytical
steady state slopes found above, overlayed onto the experimental data, and
(b) creating a table containing 𝑇0 [°C], 𝐻𝑒𝑥 𝑝 [°C/m], and 𝐻𝑎𝑛 [°C/m] for
each data set.
%}



%rowtitles = {Q3_MD.T;AL25_MD.T;AL28_MD.T;BR26_MD.T;BR29_MD.T;ST21_MD.T};
%T0 = subsref(ESSS,struct('type','()','subs',{{1,1}}))

% I have no idea on how to make a struct and use it correctly, perhaps I
% need to actually just not use a struct, or make it so that it pulls the
% T0s into one field and the HESSs into another..., but that probably takes
% as long as writing it simply and absent mindly, making things modular /
% elegant in coding srsly is sometimes just no worth it, this is the most i
% got too..: [T0,~] = subsref(ESSS,struct('type','()','subs',{{1,1}}));

T0 = [ESSS_Q3(1) ESSS_AL25(1) ESSS_AL28(1) ESSS_BR26(1) ESSS_BR29(1) ESSS_ST21(1)];
HESSS = [ESSS_Q3(2) ESSS_AL25(2) ESSS_AL28(2) ESSS_BR26(2) ESSS_BR29(2) ESSS_ST21(2)];
%Change later, T0, ... needs to be actual variables
%TABLE NEEDS TO HAVE: T0, HESS, HANA
%SSD_Table = table(T0,HESS,HANA,'RowNames',rowtitles);

%% Functions

%▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
%▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
%▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

%% ESSS Function

%Function Handle - [Returns T0, Hexp, Hana] = fn name(Temperature Data)
%function [ESSSCALC] = ESSS_FN(ESSS_Temperature,metadata)
function ESSSCALC = ESSSFN(ESSS_Temperature,metadata)


%Distance Vectors
ESSS_Dis_0 = [0 .5413];
ESSS_Dis_TC = linspace(ESSS_Dis_0(end),(ESSS_Dis_0(end)+8.89),8);

%Combine vectors for RSM calculationg, NOT true distance vector
ESSS_Dis_LSF = [ESSS_Dis_0(2),ESSS_Dis_TC(2:end)];

%Recreating the true Distance vector
ESSS_Distance = [0 ESSS_Dis_LSF];
ESSS_Dis_TC = ESSS_Distance(2:end);

%%Calculating the linear regression through polyfit
[ESSS_P,ESSS_Err] = polyfit(ESSS_Dis_LSF,ESSS_Temperature,1);

%Calculating the vector from linear regression, including T0 at 0
[ESSS_LSM,ESSS_Delta] = polyval(ESSS_P,ESSS_Distance,ESSS_Err);

%Calculating transduced temp with best fit line
%unsure what to do with this data...
ESSS_Diff_Err = ESSS_Temperature-(ESSS_LSM(2:end));


%Returning the called values
%EXPCALC [ T0, H_EXP/ESS "HESS" , H_ANA "HANA"]
ESSSCALC = [ESSS_P(2) ESSS_P(1)];


%=========================================================================

%%% Plotting the analysis


%%% PLOTTING SWITCH %%%%
plotswitch = true; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotswitch == 1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


figure(metadata.P)
hold on

%Plots linear fit line
plot(ESSS_Distance,ESSS_LSM,'b','LineWidth',1.4)

%Plots prediction interval lines
plot(ESSS_Distance,ESSS_LSM+2*ESSS_Delta,'c--')
plot(ESSS_Distance,ESSS_LSM-2*ESSS_Delta,'c--')

%Plots original raw data
plot(ESSS_Dis_TC,ESSS_Temperature,'ok','LineWidth',1.4)

%Plots desired anwser
plot(ESSS_Distance(1),ESSS_LSM(1),'*r','MarkerSize',12,'LineWidth',2)

%Annotation
title(metadata.T,'ESSS using LSF')
xlabel('Distance along Rod [cm]')
ylabel('Temperature [°C]')
legend('Linear Fit','Prediction Interval','','Acquired Data','T0','Location','southeast')

% text(ESSS_Distance(1)+.2,ESSS_LSM(1)-.2,)
text(ESSS_Distance(1)+.4,ESSS_LSM(1)-.2, sprintf('T0 = %.3f [°C]',ESSS_LSM(1)))
text(ESSS_Distance(5)+.5,ESSS_LSM(5)-.5, sprintf('H = %.3f [°C/cm]',ESSS_P(1) ))
 
axis padded
%axis([-.5 12.05 ESSS_LSM(1)-2 ESSS_LSM(end)+2])
print('-f1','PLOT ASEN 3802 - Q3 ESSS','-dpng')

hold off
else
end

end


%% uFn Function Begin

%function [u(x) function , Fourier Number] =
% HANAFN(base temp @ x0, heat trans coef,therm diffusivity, distance vector
% time span <value>, fourier terms)

function [uVec,Fo] = uFN(T0,H,alpha,t) % base temp "T0", H "H", time span "t"

L = 14.9225;
n = 10; %iteration count

%Distance Vectors
% I could make the distance vectors for HANA with more resolution, then
% just denote where TH8 is imatln the plot
HANA_Dis_0 = [0 .5413];
HANA_Dis_TC = linspace(HANA_Dis_0(end),(HANA_Dis_0(end)+8.89),8);

%Combine vectors for RSM calculationg, NOT true distance vector
HANA_Dis_LSF = [HANA_Dis_0(2),HANA_Dis_TC(2:end)];

%Recreating the true Distance vector
HANA_Distance = [0 HANA_Dis_LSF];
%HANA_Distance = [HANA_Distance HANA_Distance(end)+1.27];

X = HANA_Distance;


% Declaring fourier term transient solution with zeroes
if size(t,2) == 1000
    transientSumm = zeros(n,length(X)); 
elseif size(t,2) ~= 1
    transientSumm = zeros(length(t),length(X));
else
    transientSumm = zeros(n,length(X)); 
end

%Switch for Question 5 to Question 6
if size(t,2) == 1000
    %Question 6 parameters
    i = 1;
    %alpha = .1:.1:.8; NEED TO PASS THROUGH ALPHA WHEN WANTED
    X = X(9);

    bn_num = 8*H*L*((-1)^i); % bn nominator
    bn_den = pi*2*((2*i)-1)^2; % bn denominator
    bn = bn_num/bn_den; % bn evaluation
    
    %lambda n calculation, done by fractions as well
    lambdan_num = ((2*i)-1)*pi; %lambda nominator
    lambdan_den = 2*L; % lambda denominator
    lambdan = lambdan_num/lambdan_den; % lambda evaluation

    %summ func inner value
    summ_exp_qty = -(lambdan^2)*alpha.*t;
    summ_exp = exp(summ_exp_qty);
    summ_inner = sin(lambdan*X).*summ_exp; % outer summ expression
    summ = bn * summ_inner; % entire summ expression

    transientSumm = summ;
elseif size(t,2) ~= 1
%need to make sure that the x axis is the x position, then that the y is
%the time

    for j = 1:length(X)    
        for i = 1:n
    %bn calculation, done by fractions to ensure expression
    bn_num = 8*H*L*((-1)^i); % bn nominator
    bn_den = pi*2*((2*i)-1)^2; % bn denominator
    bn = bn_num/bn_den; % bn evaluation
    
    %lambda n calculation, done by fractions as well
    lambdan_num = ((2*i)-1)*pi; %lambda nominator
    lambdan_den = 2*L; % lambda denominator
    lambdan = lambdan_num/lambdan_den; % lambda evaluation

    %summ func inner value
    summ_exp_qty = -(lambdan^2)*alpha.*t;
    summ_exp = exp(summ_exp_qty);
    summ_inner = sin(lambdan*X(j)).*summ_exp; % outer summ expression
    summ = bn * summ_inner; % entire summ expression

    %add up fourier series
    transientSumm(:,j) = transientSumm(:,j) + summ';
        end
    end
else
    for i = 1:(n) %i being each iteration, n+1 to get to max

    %bn calculation, done by fractions to ensure expression
    bn_num = 8*H*L*((-1)^i); % bn nominator
    bn_den = pi*2*((2*i)-1)^2; % bn denominator
    bn = bn_num/bn_den; % bn evaluation
    
    %lambda n calculation, done by fractions as well
    lambdan_num = ((2*i)-1)*pi; %lambda nominator
    lambdan_den = 2*L; % lambda denominator
    lambdan = lambdan_num/lambdan_den; % lambda evaluation

    %summ func inner value
    summ_exp_qty = -(lambdan^2)*alpha.*t;
    summ_exp = exp(summ_exp_qty);
    summ_inner = sin(lambdan*X).*summ_exp; % outer summ expression
    summ = bn * summ_inner; % entire summ expression

    %add up fourier series
    transientSumm(i+1,:) = transientSumm(i,:) + summ;
    end
end

%rant below
%{
figuring out this entire for loop and size thing was about the worst
thing ever, it could have been way easier if I just drew out the loop as
onto a whiteboard so then I could have just visualized it instead of trying
every kind of solution, it was way too much really, also im not really sure
if this is right beause should the first index be atleast somewhat
influenced by the transient solution because n = 1 means that the first
series should influence the calculation, which means that row 1 is only t0
+ hx, so should the plot instead plot 2:end???
%}





% U(x) function combines steady and transient soln
uFn = @(x,t0,h,SUMM) t0 + h.*x + SUMM;

% Returning vector, calculating the uFN Function
uVec = uFn(X,T0,H,transientSumm); %returns the number 

% Fourier number function
FoFn = @(a,t,l)(a*t)/((L)^2); %anon fn, a "alpha", t "time", L "length"

% Fourier calculation
Fo = FoFn(alpha,t,L);

end


%% MIS MODEL II ufn



%function [u(x) function , Fourier Number] =
% HANAFN(base temp @ x0, heat trans coef,therm diffusivity, distance vector
% time span <value>, fourier terms)

function [uVec,Fo] = MIIFN(T0,H,alpha,t,g) % base temp "T0", H "H", time span "t", g inital temperature span

L = 14.9225; %test??
n = 10; %iteration count

%Distance Vectors
% I could make the distance vectors for HANA with more resolution, then
% just denote where TH8 is imatln the plot
HANA_Dis_0 = [0 .5413];
HANA_Dis_TC = linspace(HANA_Dis_0(end),(HANA_Dis_0(end)+8.89),8);

%Combine vectors for RSM calculationg, NOT true distance vector
HANA_Dis_LSF = [HANA_Dis_0(2),HANA_Dis_TC(2:end)];

%Recreating the true Distance vector
HANA_Distance = [0 HANA_Dis_LSF];
%HANA_Distance = [HANA_Distance HANA_Distance(end)+1.27];

X = HANA_Distance;


% Declaring fourier term transient solution with zeroes

transientSumm = zeros(length(t),length(X));

    for j = 1:length(X)    
        for i = 1:n
    
    %lambda n calculation, done by fractions as well
    lambdan_num = ((2*i)-1)*pi; %lambda nominator
    lambdan_den = 2*L; % lambda denominator
    lambdan = lambdan_num/lambdan_den; % lambda evaluation


    bnnum = (2*g)*(sin(L*lambdan)-(L*lambdan*cos(L*lambdan)));
    bnden = L*lambdan^2;
    bn = bnnum / bnden;


    %summ func inner value
    summ_exp_qty = -(lambdan^2)*alpha.*t;
    summ_exp = exp(summ_exp_qty);
    summ_inner = sin(lambdan*X(j)).*summ_exp; % outer summ expression
    summ = bn * summ_inner; % entire summ expression

    %add up fourier series
    transientSumm(:,j) = transientSumm(:,j) + summ';
        end
    end




% U(x) function combines steady and transient soln
uFn = @(x,t0,h,SUMM) t0 + h.*x + SUMM;

% Returning vector, calculating the uFN Function
uVec = uFn(X,T0,H,transientSumm); %returns the number 

% Fourier number function
FoFn = @(a,t,l)(a*t)/((L)^2); %anon fn, a "alpha", t "time", L "length"

% Fourier calculation
Fo = FoFn(alpha,t,L);


end

%% Results 4 & 5, Time to Steady State / alpha variance

% Calculating alpha



function [tss,fo] = TSSFN(FO)

FOlog = FO > 2;

FoTssIdx = find(FOlog,1,'first');

fo = FO(FoTssIdx);
tss = FoTssIdx * 10;



end
