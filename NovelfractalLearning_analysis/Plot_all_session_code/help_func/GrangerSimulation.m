close all; clear all; clc;
%Granger causality values from bi-variate time series generated from a 1-st order AR model in Roebroeck et al. (NeuroImage 2005)
%[x(t+1), y(t+1)]^T=[-0.9, 0; I, -0.9][x(t), y(t)]^T+[ex(t), ey(t)], where [ex(t) ey(t)]^T are bi-variate noise with a covariance matrix S=[1 0; 0 1];

I=1; %strong influence would be something like 0.5
n=10000;
A=[-0.9 0; I -0.9];
S=[1 0; 0 1];
D=10; %additional time points to reach "steady state"

%calling AR model. make time series
v=arsim(zeros(1,2),A,S,n+2000+D);
v=v(1:n,:);

subplot(211); hold on; plot(v(:,1),'r'); %time series "X";
legend({'X'});
subplot(212); hold on; plot(v(:,2),'b'); %time series "Y";
legend({'Y'});
xlabel('time (sample)');

%calculate the granger
[granger, granger_F, granger_p]=etc_granger(v,1);
% granger
% granger_F
% granger_p

close all; clear all; clc
% Lin et al (Human Brain Mapp 2009).
n=10000;
I=0.5; %strong influence
A=[-0.9 0; I -0.9];
S=[1 0; 0 1];
D=10; %additional time points to reach "steady state"

%calling AR model to simulate the time series;
v=arsim(zeros(1,2),A,S,n+2000+D);
v=v(1:n,:);

subplot(211); hold on; plot(v(:,1),'r'); %time series "X";
legend({'X'});
subplot(212); hold on; plot(v(:,2),'b'); %time series "Y";
legend({'Y'});
xlabel('time (sample)');

%calculate the granger
[granger, granger_F, granger_p]=etc_granger(v,1);