clear;
close all;
load mockdata2023.mat;
%% 
%mock cases and deaths

t = 400;

figure;
hold on;
plot(newInfections);
plot(cumulativeDeaths);
legend('Infections','Deaths')
hold off;


%6 stages S I R D vaccinated Breakthrough
vax_date = 120;

%
r_infec = 0.1;
r_reinfec = 0.1;
r_recover = 0.1;
r_death = 0.05;
r_vax = 0.3;
r_vaxBreak = 0.04;

%stages before vax
S0 = 0.9;
I0 = 0.1;
R0 = 0;
D0 = 0;


%stages after vax
S1 = 0.6;
I1 = 0.1;
R1 = 0.1;
D1 = 0.05;
V1 = 0.1;
VB1 = 0.05;

unknowns = [r_infec,r_reinfec,r_recover,r_death,r_vax,r_vaxBreak ...
    ,S0,I0,R0,D0,S1,I1,R1,D1,V1,VB1]';


A = [1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0];
b = 1;
Aeq = [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1];
beq = 1;

lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

fun = @(unknowns)ModelCompare_mock(unknowns,t,newInfections,cumulativeDeaths);

unknowns_opt_V = fmincon(fun, unknowns, A, b ,Aeq ,beq, lb, ub);

A_vax=[1-unknowns_opt_V(1)-unknowns_opt_V(5) ,0 ,unknowns_opt_V(2),0 ,0  ,0;
unknowns_opt_V(2)   ,1-unknowns_opt_V(3)-unknowns_opt_V(4) ,0 ,0  ,0  ,0;
0  ,unknowns_opt_V(3),1-unknowns_opt_V(2) - unknowns_opt_V(5),0 ,1-unknowns_opt_V(6),unknowns_opt_V(3);    
0     ,unknowns_opt_V(4), 0 ,1,  0,  0;
unknowns_opt_V(5)   ,0 , unknowns_opt_V(5)  ,0 ,0 ,0;
0,  0,  0,  0,  unknowns_opt_V(6) ,1-unknowns_opt_V(3)
];

A_noVax = [1-unknowns_opt_V(1) ,0 ,unknowns_opt_V(2),0 ,0  ,0;
unknowns_opt_V(2)   ,1-unknowns_opt_V(3)-unknowns_opt_V(4) ,0 ,0  ,0  ,0;
0  ,unknowns_opt_V(3),1-unknowns_opt_V(2),0 ,0,0;    
0     ,0, 0 ,1,  0,  0;
0   ,0 , 0  ,0 ,1 ,0;
0,  0,  0,  0,  0 ,1
];

x1 = unknowns_opt_V(11:16);

x0 = [unknowns_opt_V(7:10);0;0];

%simulate model
t = 400;

x = zeros(6,t);

x(:,1) = x0;

for i = 2:vax_date-1
    % Model equations
    x0 = A_noVax * x0;
    x(:,i) = x0;
    x(5,:) = 0;
    x(6,:) = 0;
end

x(:,vax_date) = x1;

for i = vax_date+1:t
    x1 = A_vax * x1;
    x(:,i) = x1;
end

% Extract active counts
%S = x(1,:);
I = x(2,:);
%R0 = x(3,:);
D = x(4,:);
V = x(5,:);
VB = x(6,:);

%convert active data of SIRD into cumulative data
new_cases = I + VB;
cum_deaths = D;




%graph the results
figure;
hold on;
plot(new_cases);
plot(newInfections);
legend('Model','newInfec')
hold off;

figure;
hold on;
plot(cum_deaths);
plot(cumulativeDeaths);
legend('Model','Deaths')
hold off;





%% 
%defining function for fmincon before vax


%defining function for fmincon during vax
function err = ModelCompare_mock(unknowns0, t, newInfections, cumulativeDeaths)

vax_date = 120;

%
r_infec0 = unknowns0(1,:);
r_reinfec0 = unknowns0(2,:);
r_recover0 = unknowns0(3,:);
r_death0 = unknowns0(4,:);
r_vax0 = unknowns0(5,:);
r_vaxBreak0 = unknowns0(6,:);

%stages before vax
S0 = unknowns0(7,:);
I0 = unknowns0(8,:);
R0 = unknowns0(9,:);
D0 = unknowns0(10,:);


%stages after vax
S1 = unknowns0(11,:);
I1 = unknowns0(12,:);
R1 = unknowns0(13,:);
D1 = unknowns0(14,:);
V1 = unknowns0(15,:);
VB1 = unknowns0(16,:);


x1 = [S1, I1, R1, D1, V1, VB1]';

x0 = [S0, I0, R0, D0, 0, 0]';

A_vax0=[1-r_infec0-r_vax0    ,0          ,r_reinfec0      ,0      ,0      ,0;
r_infec0               ,1-r_recover0-r_death0 ,0     ,0     ,0      ,0;
0  ,r_recover0    ,1-r_reinfec0 - r_vax0 ,0 ,1-r_vaxBreak0,r_recover0;    
0     ,r_death0, 0 ,1,  0,  0;
r_vax0   ,0 , r_vax0  ,0 ,0 ,0;
0,  0,  0,  0,  r_vaxBreak0 ,1-r_recover0
];

A_noVax0 = [1-r_infec0    ,0          ,r_reinfec0      ,0      ,0      ,0;
r_infec0               ,1-r_recover0-r_death0 ,0     ,0     ,0      ,0;
0  ,r_recover0    ,1-r_reinfec0 ,0 ,0,0;    
0     ,r_death0, 0 ,1,  0,  0;
0   ,0 , 0  ,0 ,1 ,0;
0,  0,  0,  0,   0,1
];

%simulate model
t = 400;

x = zeros(6,t);

x(:,1) = x0;

for i = 2:vax_date-1
    % Model equations
    x0 = A_noVax0*x0;
    x(:,i) = x0;
    x(5,:) = 0;
    x(6,:) = 0;
end

x(:,vax_date) = x1;

for i = vax_date+1:t
    x1 = A_vax0 * x1;
    x(:,i) = x1;
end

% Extract active counts
%S = x(1,:);
I = x(2,:);
%R0 = x(3,:);
D = x(4,:);
V = x(5,:);
VB = x(6,:);

%convert active data of SIRD into cumulative data
new_cases = I + VB;
cum_deaths = D;

% Calculate error
cases_err = sum((newInfections - new_cases).^2);
deaths_err = sum((cumulativeDeaths - cum_deaths).^2);
err = cases_err + deaths_err;
end
