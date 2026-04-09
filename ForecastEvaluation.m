% Reset everything, clear memory and screen
clear; close all; clc;
% Start clock
tic;

addpath('Functions')

[data,names] = xlsread(['Out-of-Sample Forecasts.xlsx'],'1850-1900 anomaly'); % load US data

date = data(1:7,1);
TempActual = data(1:7,2);

Forecasts = data(1:7,3:end);
ModelNames = names(3:end);
[H,N] = size(Forecasts);



for n=1:N
    for t=1:H
SquaredErrors(t,n) = (Forecasts(t,n)-TempActual(t))^2;
AbsoluteErrors(t,n) = abs(Forecasts(t,n)-TempActual(t));
    end
end
MSE = mean(SquaredErrors,1);
MAE = mean(AbsoluteErrors,1);


figure
for n=1:N
plot(date,SquaredErrors(:,n));
hold on
legend(ModelNames)
end

figure
for n=1:N
plot(date,AbsoluteErrors(:,n));
hold on
legend(ModelNames)
end