clear;
close all;
load COVID_STL.mat;
%figure;
%hold on;
%plot(dates,cases_STL);
%plot(dates,deaths_STL);
%legend('cases','deaths')  
%hold off;
%find the indices for Delta period
delta_start_end = find(dates > datetime(2021,6,30) & dates < datetime(2021,10,26));
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

SS=0.99;
SI=1-SS;

IS=0.68;
II=0.2;
IR=0.01;
ID=1-IS-II-IR;

RS=0.01;
RR=1-RS;

a=0.5;
b=0.01;
c=0.48;
d=0.01;

A=[SS,IS,RS,0;SI,II,0,0;0,IR,RR,0;0,ID,0,1];
x0=[a b c d]';

xtotal(:,1) = x0;
time=length(delta_dates);
for t=2:time
    x0 = A * x0;
    xtotal(:,t) = x0;
end
S=xtotal(1,:);
I=xtotal(2,:);
R=xtotal(3,:);
D=xtotal(4,:);

cum_i = cumsum(I)*POP_STL;
cum_d = cumsum(D)*POP_STL;

cum_i_delta = cum_i + delta_cases(1);
cum_d_delta = cum_d + delta_deaths(1);



figure;
hold on;
plot(delta_dates,cum_i_delta);
plot(delta_dates,cum_d_delta);
legend('cum delta infection', 'cum delta death');
hold off;
title('hand tune data');





