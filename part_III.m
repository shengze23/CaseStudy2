close all;
load COVID_STL.mat;
%% 
%================Delta Wave Interactions================
%Initialize population for two regions

POP_1 = 2876487;
POP_2 = 1263619;

%Travel matrix 1
S2to1 = 0.09;  %travel rate of Susceptible pop2 to region 1
I2to1 = 0.11;
R2to1 = 0.08;

T_1 = [S2to1,0    ,0     ,0;
      0     ,I2to1 ,0    ,0;
      0     ,0     ,R2to1,0;    
      0     ,0     ,0    ,0];


%Transition Matrix 1
r_infec1 = unknowns_opt_D(1);
r_recover1 = unknowns_opt_D(3);
r_reinfec1 = unknowns_opt_D(2);
r_death1 = unknowns_opt_D(4);

A_1 = [1-r_infec1 - S1to2, r_recover1 ,r_reinfec1 ,0;
     r_infec1, 1-r_recover1-r_death1 - I1to2,0    ,0;
     0  ,0            ,1-r_reinfec1 - R1to2 ,0;    
     0     ,r_death1, 0 ,1];

%Travel matrix 2
S1to2 = 0.15;  %travel rate of Susceptible pop1 to region 2
I1to2 = 0.11;
R1to2 = 0.13;

T_2 = [S1to2, 0    ,0    ,0;
      0     ,I1to2 ,0    ,0;
      0     ,0     ,R1to2,0;    
      0     ,0     ,0    ,0];


%Transition Matrix 2
r_infec2 = 0.15;
r_recover2 = 0.02;
r_reinfec2 = 0.04;
r_death2 = 0.03;

A_2 = [1-r_infec2 - S2to1, r_recover2  ,r_reinfec2 ,0;
     r_infec2, 1-r_recover2-r_death2 - I2to1,0    ,0;
     0  ,0            ,1-r_reinfec2 - R2to1,0;    
     0     ,r_death2, 0 ,1];

%concatenating 4 matrices 

A_travel = [A_1,T_1;T_2,A_2];

%Initial SIRD conditions of region 1 and 2

sir_ini = [POP_1,0,0,0,POP_2,0,0,0]';

%days of interactions
t = 1000;

sir_tot = zeros(8,t);
sir_tot(:,1) = sir_ini;

for i = 2:t
    sir_ini = A_travel * sir_ini;
    sir_tot(:,i) = sir_ini;
end

%extract data from interactions model
I_1 = sir_tot(2,:);
I_2 = sir_tot(6,:);

D_1 = sir_tot(4,:);
D_2 = sir_tot(8,:);

cum_cases_1 = cumsum(I_1)*POP_1;
cum_cases_2 = cumsum(I_2)*POP_2;

%plotting the data

figure;
hold on;
plot(cum_cases_1);
plot(cum_cases_2);
legend('Cases 1','Cases 2');
xlabel('Days')
ylabel('Cases')
hold off;

%plotting the combined cum cases and deaths


%-----------------------Restricted travel Policy Delta ------------------
%policy - infected population in two regions cannot travel, others are
%encouraged to not travel







%% 
%=======================Omicron Wave Interactions========================
%Initialize population for two regions

POP_1 = 2876487;
POP_2 = 1263619;

%Transition Matrix 1
r_infec1 = unknowns_opt_O(1);
r_recover1 = unknowns_opt_O(3);
r_reinfec1 = unknowns_opt_O(2);
r_death1 = unknowns_opt_O(4);

A_1 = [1-r_infec1 - S1to2, r_recover1 ,r_reinfec1 ,0;
     r_infec1, 1-r_recover1-r_death1 - I1to2,0    ,0;
     0  ,0            ,1-r_reinfec1 - R1to2 ,0;    
     0     ,r_death1, 0 ,1];

%Transition Matrix 2
r_infec2 = 0.15;
r_recover2 = 0.02;
r_reinfec2 = 0.04;
r_death2 = 0.03;

A_2 = [1-r_infec2 - S2to1, r_recover2  ,r_reinfec2 ,0;
     r_infec2, 1-r_recover2-r_death2 - I2to1,0    ,0;
     0  ,0            ,1-r_reinfec2 - R2to1,0;    
     0     ,r_death2, 0 ,1];

%Travel matrix 1
S2to1 = 0.09;  %travel rate of Susceptible pop2 to region 1
I2to1 = 0.11;
R2to1 = 0.08;

T_1 = [S2to1,0    ,0     ,0;
      0     ,I2to1 ,0    ,0;
      0     ,0     ,R2to1,0;    
      0     ,0     ,0    ,0];

%Travel matrix 2
S1to2 = 0.15;  %travel rate of Susceptible pop1 to region 2
I1to2 = 0.11;
R1to2 = 0.13;

T_2 = [S1to2, 0    ,0    ,0;
      0     ,I1to2 ,0    ,0;
      0     ,0     ,R1to2,0;    
      0     ,0     ,0    ,0];

%concatenating 4 matrices 

A_travel = [A_1,T_1;T_2,A_2];

%Initial SIRD conditions of region 1 and 2

sir_ini = [POP_1,0,0,0,POP_2,0,0,0]';

%days of interactions
t = 500;

sir_tot = zeros(8,t);
sir_tot(:,1) = sir_ini;


for i = 2:t
    sir_ini = A_travel * sir_ini;
    sir_tot(:,i) = sir_ini;
end

%extract data from interactions model
I_1 = sir_tot(2,:);
I_2 = sir_tot(6,:);

D_1 = sir_tot(4,:);
D_2 = sir_tot(8,:);

cum_cases_1 = cumsum(I_1)*POP_1;
cum_cases_2 = cumsum(I_2)*POP_2;

%plotting the data

figure;
hold on;
plot(cum_cases_1);
plot(cum_cases_2);
legend('Cases 1 OMI','Cases 2 OMI');
xlabel('Days')
ylabel('Cases')
hold off;


%plotting the combined cum cases and deaths


%-----------------------Restricted travel Policy Omicron ------------------
%policy infected population in two 




%% 
