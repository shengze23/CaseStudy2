clear;
close all;
load COVID_STL.mat;

%find the indices for Delta period
delta_start_end = find(dates >= datetime(2021,6,30) & dates <= datetime(2021,11,1));

delta_start = find(dates == datetime(2021,6,30));

delta_dates = dates(delta_start_end);
delta_cases = cases_STL(delta_start_end);
delta_deaths = deaths_STL(delta_start_end);

%plotting delta period
figure;
hold on;
plot(delta_dates,delta_cases);
plot(delta_dates,delta_deaths);
title('delta wave')
hold off;


%find the indices for Omicron wave
omicron_start_end = find(dates >= datetime(2021,10,27) &  dates <= datetime(2022,3,29));

omicron_start = find(dates == datetime(2021,10,27));

omicron_dates = dates(omicron_start_end);
omicron_cases = cases_STL(omicron_start_end);
omicron_deaths = deaths_STL(omicron_start_end);

%plot omicron
figure;
hold on;
plot(omicron_dates,omicron_cases);
plot(omicron_dates,omicron_deaths);
title('omicron wave')
hold off;

%% 

%fmincon for delta wave=================================================
% Initial assumption
r_infec = 0.1;
r_reinfec = 0.1;
r_recover = 0.1;
r_death = 0.1;

S = 0.75;
I = 0.1;
R = 0.1;
D = 0.05;

x0 = [S,I,R,D]';

unknowns = [r_infec,r_reinfec,r_recover,r_death,S,I,R,D]';

A = [1,1,1,1,0,0,0,0];
b = 1;
Aeq = [0,0,0,0,1,1,1,1];
beq = 1;

%lb = [0,0,0,0,0.7,0.2,0,0.05]'; %lower bound
%ub = [1,1,1,1,0.8,0.8,0.8,0.5]'; %upper bound

lb = [0, 0, 0, 0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1];

options = optimoptions('fmincon','Algorithm','interior-point');

fun = @(unknowns)ModelCompare(unknowns,delta_dates,delta_cases,delta_deaths);

unknowns_opt_D = fmincon(fun, unknowns, A, b ,Aeq ,beq, lb, ub);%,[],options);

%graph the results
x_opt0 = unknowns_opt_D(5:8);
xtot = zeros(4,length(delta_dates));

A_opt = [1-unknowns_opt_D(1) ,unknowns_opt_D(3) ,unknowns_opt_D(2) ,0;
     unknowns_opt_D(1) ,1-unknowns_opt_D(3)-unknowns_opt_D(4) ,0 ,0;
     0 ,0 ,1-unknowns_opt_D(2) ,0;
     0 ,unknowns_opt_D(4) ,0 ,1
    ];

xtot(:,1) = x_opt0;

for t = 2:length(delta_dates)
    x_opt0 = A_opt * x_opt0;
    xtot(:,t) = x_opt0;
end
S_opt = xtot(1,:);
I_opt = xtot(2,:);
R_opt = xtot(3,:);
D_opt = xtot(4,:);

oneVector = ones(1,length(delta_dates));

cum_cases = cumsum(I_opt)*POP_STL + delta_cases(1);
cum_deaths = D_opt*POP_STL + delta_deaths(1);

figure;
tiledlayout(2,2);

nexttile
hold on;
plot(delta_dates,cum_cases);
plot(delta_dates,delta_cases);
title('Model Cases vs. Real Cases')
legend('Model Cases','Real Cases')
hold off;

nexttile
hold on;
plot(delta_dates,cum_deaths);
plot(delta_dates,delta_deaths);
title('Model Deaths vs. Real Deaths')
legend('Model Deaths','Real Deaths')
hold off;

nexttile
hold on;
plot(delta_dates,delta_cases);
plot(delta_dates,delta_deaths);
title('Delta Wave')
legend('Cases','Deaths')
hold off;

nexttile
hold on;
plot(delta_dates,cum_cases);
plot(delta_dates,delta_cases);
plot(delta_dates,cum_deaths);
plot(delta_dates,delta_deaths);
title('Model vs. Real Data')
legend('Model Cases','Real Cases','Model Deaths','Real Deaths')
hold off;


%print out the matrix A and initial conditions
disp(A_opt);
disp(unknowns_opt_D(5:8));


%% ------------------------------------------------------------
%Omicron Wave
fun = @(unknowns)ModelCompare(unknowns,omicron_dates,omicron_cases,omicron_deaths);

unknowns_opt_O = fmincon(fun, unknowns, A, b ,Aeq ,beq, lb, ub);%,[],options);

%graph the results
x_opt0 = unknowns_opt_O(5:8);
xtot = zeros(4,length(omicron_dates));

A_opt = [1-unknowns_opt_O(1) ,0 ,unknowns_opt_O(2) ,0;
     unknowns_opt_O(1) ,1-unknowns_opt_O(3)-unknowns_opt_O(4) ,0 ,0;
     0 ,unknowns_opt_O(3) ,1-unknowns_opt_O(2) ,0;
     0 ,unknowns_opt_O(4) ,0 ,1
    ];

xtot(:,1) = x_opt0;

for t = 2:length(omicron_dates)
    x_opt0 = A_opt * x_opt0;
    xtot(:,t) = x_opt0;
end
S_opt = xtot(1,:);
I_opt = xtot(2,:);
R_opt = xtot(3,:);
D_opt = xtot(4,:);


cum_cases = cumsum(I_opt)*POP_STL + omicron_cases(1);
cum_deaths = D_opt*POP_STL + omicron_deaths(1);


%plotting model and real data and compare them together 
figure;
tiledlayout(2,2);

nexttile
hold on;
plot(omicron_dates,cum_cases);
plot(omicron_dates,omicron_cases);
title('Model Cases vs. Real Cases')
legend('Model Cases','Real Cases')
hold off;

nexttile
hold on;
plot(omicron_dates,cum_deaths);
plot(omicron_dates,omicron_deaths);
title('Model Deaths vs. Real Deaths')
legend('Model Deaths','Real Deaths')
hold off;

nexttile
hold on;
plot(omicron_dates,omicron_cases);
plot(omicron_dates,omicron_deaths);
title('omicron Wave')
legend('Cases','Deaths')
hold off;

nexttile
hold on;
plot(omicron_dates,cum_cases);
plot(omicron_dates,omicron_cases);
plot(omicron_dates,cum_deaths);
plot(omicron_dates,omicron_deaths);
title('Model vs. Real Data')
legend('Model Cases','Real Cases','Model Deaths','Real Deaths')
hold off;

%print out the matrix A and initial conditions
disp(A_opt);
disp(unknowns_opt_O(5:8));



%% 
%Policy that reduces 25% of the cases and deaths
%Addding a state Q - quarantine

%50% of population in each state except for deceased will remain at home 
%so no more interactions will happen between them

unknowns_Q = [unknowns_opt_O;0];

x_opt0 = unknowns_Q(5:9);
xtot = zeros(5,length(omicron_dates));
Q = 0.75; % percentage of the population is quarantined

A_opt = [1 - unknowns_opt_O(1)*(1-Q), 0, unknowns_opt_O(2), 0,  0;
         unknowns_opt_O(1)*(1-Q), 1 - unknowns_opt_O(3) - unknowns_opt_O(4), 0, 0,  0;
         0, unknowns_opt_O(3), 1 - unknowns_opt_O(2), 0,  0;
         0, unknowns_opt_O(4)*(1-Q), 0, 1,  0;
         0, 0, 0, 0, 1;
         ];

xtot(:,1) = x_opt0;

for t = 2:length(omicron_dates)
    x_opt0 = A_opt * x_opt0;
    xtot(:,t) = x_opt0;
end
S_opt = xtot(1,:);
I_opt = xtot(2,:);
R_opt = xtot(3,:);
D_opt = xtot(4,:);
Q = xtot(5,:);
%
%50 of the Susceptible population will be quarantined, assuming that
%quarantined population will not die or be infected
cum_cases = cumsum(I_opt)*POP_STL + omicron_cases(1);
cum_deaths = D_opt*POP_STL + omicron_deaths(1);
%
%plotting the quarantined model
figure;
tiledlayout(2,2);

nexttile
hold on;
plot(omicron_dates,cum_cases);
plot(omicron_dates,omicron_cases);
title('Model Cases vs. Real Cases')
legend('Model Cases','Real Cases')
hold off;

nexttile
hold on;
plot(omicron_dates,cum_deaths);
plot(omicron_dates,omicron_deaths);
title('Model Deaths vs. Real Deaths')
legend('Model Deaths','Real Deaths')
hold off;

nexttile
hold on;
plot(omicron_dates,omicron_cases);
plot(omicron_dates,omicron_deaths);
title('omicron Wave')
legend('Cases','Deaths')
hold off;

nexttile
hold on;
plot(omicron_dates,cum_cases);
plot(omicron_dates,omicron_cases);
plot(omicron_dates,cum_deaths);
plot(omicron_dates,omicron_deaths);
title('Model vs. Real Data')
legend('Model Cases','Real Cases','Model Deaths','Real Deaths')
hold off;

%calculate the percentage of death and cases in the quarantine model 

OMI_deaths = omicron_deaths(length(omicron_dates));
OMI_cases = omicron_cases(length(omicron_dates));

OMI_t = length(omicron_dates); %durantion of omicron

diff_cases = omicron_cases(OMI_t) - cum_cases(OMI_t);
diff_cases_percen = ((diff_cases)/OMI_cases)*100;
diff_deaths = omicron_deaths(OMI_t) - cum_deaths(OMI_t);
diff_deaths_percen = ((diff_deaths)/OMI_deaths)*100;

fprintf('The reduced infected population percentage is: %f\n',diff_cases_percen);
fprintf('The reduced deceased population percentage is: %f\n',diff_deaths_percen);


%% 
%defining function for fmincon
function err = ModelCompare(unknowns0, dates, cases_STL, deaths_STL)

%
r_infec0 = unknowns0(1,:);
r_reinfec0 = unknowns0(2,:);
r_recover0 = unknowns0(3,:);
r_death0 = unknowns0(4,:);

S0 = unknowns0(5,:);
I0 = unknowns0(6,:);
R0 = unknowns0(7,:);
D0 = unknowns0(8,:);

x0 = [S0, I0, R0, D0]';


A = [1-r_infec0, r_recover0  ,r_reinfec0 ,0;
     r_infec0, 1-r_recover0-r_death0,0    ,0;
     0  ,0  ,1-r_reinfec0 ,0;    
     0     ,r_death0, 0 ,1];


%simulate model
t = length(dates);

x = zeros(4,t);

x(:,1) = x0;

for i = 2:t
    % Model equations
    x0 = A*x0;
    x(:,i) = x0;
end

% Extract active counts
%S = x(1,:);
I = x(2,:);
%R0 = x(3,:);
D = x(4,:);

%convert active data of SIRD into cumulative data
cum_cases = cumsum(I)*1290497 + cases_STL(1);
cum_deaths = D*1290497 + deaths_STL(1);

% Calculate error
cases_err = sum((cases_STL - cum_cases).^2);
deaths_err = sum((deaths_STL - cum_deaths).^2);
err = cases_err + deaths_err;

end
%% 