clear;
close all;

%basic model without reinfection

%initializing a matrix x that contains the percentage of SIRD from 0 - time
x = zeros(4,500);

x_0 = [0.9,0.1,0,0]';

%transition matrix A
A = [0.95 0.04 0 0;
    0.05 0.85 0 0; 
    0 0.1 1 0; 
    0 0.01 0 1];

x(:,1) = x_0;

for t = 2:500
    x_0 = A * x_0;  % Update the state using the transition matrix
    x(:, t) = x_0;
end

%get the S I R D data from the matrix x
S = x(1,:);
I = x(2,:);
R = x(3,:);
D = x(4,:);

figure;
hold on;
plot(S);
plot(I);
plot(R);
plot(D);
title('Basic SIRD Model')
ylabel('propotion')
xlabel('time')
legend('S','I','R','D')
hold off;


%reinfection possible

% 0.05 of the recovered will become suscpetible again
A_reinfec = [0.95 0.04 0.05 0;
     0.05 0.85 0 0;
     0    0.10 0.95 0;
     0    0.01 0 1];

x_reinfec = zeros(4,3000);

x_reinfec(:,1) = x_0;

for t = 2:3000
    x_0 = A_reinfec * x_0; 
    x_reinfec(:, t) = x_0;
end 

S = x_reinfec(1,:);
I = x_reinfec(2,:);
R = x_reinfec(3,:);
D = x_reinfec(4,:);

figure;
hold on;
plot(S);
plot(I);
plot(R);
plot(D);
title('Reinfection-possible SIRD Model')
ylabel('propotion')
xlabel('time')
legend('S','I','R','D')
hold off;





