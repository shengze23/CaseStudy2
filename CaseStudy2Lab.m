clear all;
close all;
time=1000;
xtotal=zeros(4,1000);
xtotal(:,1)=1;
x=[0.9 0.1 0 0]';
A=[0.95 0.04 0 0; 0.05 0.85 0 0; 0 0.10 1 0; 0 0.01 0 1];
for t=2:time
    x=A*x;
    xtotal(:,t)=x;
end
S=xtotal(1,:);
I=xtotal(2,:);
R=xtotal(3,:);
D=xtotal(4,:);
figure;

plot(S);
hold on;
plot(I);
plot(R);
plot(D);
hold off;

