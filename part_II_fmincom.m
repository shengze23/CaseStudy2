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

%draft
r = 0.9914;
y = 0.001;
z = 0.0017;
m = 0.0017;
j = 0.0017;
S = 0.49;
I = 0.013;
R = 0.49;
D = 0.0009;

x0 = [S,I,R,D]';
x = zeros(4,length(delta_dates));
A = [r ,y ,m ,0;
     1-r ,z ,0 ,0;
     0 ,j ,1-m ,0;
     0 ,1-y-z-j ,0 ,1
    ];
x(:,1) = x0;
for t = 2:length(delta_cases)
   x(:,t) = A * x(:,t-1);
end
I = x(2,:);
D = x(4,:);

cum_cases = I*POP_STL + delta_cases(1);
cum_deaths = D + delta_deaths(1);

figure;
hold on;
plot(delta_dates,delta_cases);
plot(delta_dates,cum_cases);
legend('original','wode')
hold off;


%% 

%fmincon for delta
% Initial assumption

r = 0.9725;
y = 0.0505;
z = 0.05;
m = 0.0125;
j = 0.4905;
S = 0.63;
I = 0.138;
R = 0.1613;
D = 0;

x0 = [S,I,R,D]';

unknowns = [r,y,z,m,j,S,I,R,D]';

A = [1,1,1,1,1,0,0,0,0];
b = 1;
Aeq = [0,0,0,0,0,1,1,1,1];
beq = 1;

%lb = [0.97,0.001,0.3,0.005,0,0.5,0.05,0.04,0.001]'; % Lower bounds
%ub = [0.999,0.3,0.7,0.2,0.5,1,0.2,0.5,0.1]'; % Upper bounds
lb = [0.99,0,0,0,0,0,0,0,0]';
ub = [1,1,1,1,1,1,1,1,1]';

options = optimoptions('fmincon','Algorithm','interior-point');

fun = @(unknowns)ModelCompare(unknowns,delta_dates,delta_cases,delta_deaths);

unknowns_opt = fmincon(fun, unknowns, A, b ,Aeq ,beq, lb, ub,[],options);

%defining function for fmincon
function err = ModelCompare(unknowns0, dates, cases_STL, deaths_STL)

r0 = 0.95;
y0 = 0.01;
z0 = 0.02;
m0 = 0.03;
j0 = 0.01;
S0 = 0.9;
I0 = 0.1;
R0 = 0;
D0 = 0;
%
%unknowns0 = [r0,y0,z0,m0,j0,S0,I0]';
r0 = unknowns0(1,:);
y0 = unknowns0(2,:);
z0 = unknowns0(3,:);
m0 = unknowns0(4,:);
j0 = unknowns0(5,:);
S0 = unknowns0(6,:);
I0 = unknowns0(7,:);
R0 = unknowns0(8,:);
D0 = unknowns0(9,:);

x0 = [S0, I0, R0, D0]';

%transition matrix
A = [r0 ,y0 ,m0 ,0;
     1-r0 ,z0 ,0 ,0;
     0 ,j0 ,1-m0 ,0;
     0 ,1-y0-z0-j0 ,0 ,1
    ];
%simulate model
t = length(dates);

x = zeros(4,t);

x(:,1) = x0;

for i = 2:t
    % Model equations  
    x(:,t) = A*x(:,t-1);
end

% Extract active counts
%S0 = x(1,:);
I = x(2,:);
%R0 = x(3,:);
D = x(4,:);

%convert active data of SIRD into cumulative data
cum_cases = cumsum(I).*1290497 + cases_STL(1);
cum_deaths = D.*1290497 + deaths_STL(1);

% Calculate error
cases_err = norm(cum_cases - cases_STL);
deaths_err = norm(cum_deaths - deaths_STL);
err = cases_err + deaths_err;
end
%% 


