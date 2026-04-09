clear all; clc;


p = 2;   
data = xlsread('Final Dataset Forecast SSP.xlsx','hard_constraintapp1'); % load US data
[soft dat] = xlsread('Final Dataset Forecast SSP.xlsx','Soft_constraint I'); % load US data

hind = 0; % 0 - baseline & 1 - adverse
idx_ns = find(data(1,:)'==1);

data= data(2:end,:);
%data = data(100:end,:);

Y0 = data(1:p,:);  % save the first 8 obs as the initial conditions
Y = data(p+1:end,:);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end

%% Set the conditions
T2=27; % Forecast horizon
n_con=1;% number of conditioning hard variables 
n_soft=2;%number of conditioning soft variables
n_t = n_con + n_soft;


if hind == 0
    r_lower = soft(2:end,[1,3])';
    r_upper = soft(2:end,[2,4])';
    y_c = soft(2:end,5)';
else
    r_lower = soft(2:end,[6,8])';
    r_upper = soft(2:end,[7,9])';
   y_c = soft(2:end,10)';
end

%%
data = xlsread('Final Dataset Forecast SSP.xlsx','hard_constraintapp1'); % load US data
[soft dat] = xlsread('Final Dataset Forecast SSP.xlsx','Soft_constraint I'); % load US data

hind = 0; % 0 - baseline & 1 - adverse
idx_ns = find(data(1,:)'==1);

data= data(2:end,:);
%data = data(100:end,:);

Y0 = data(1:p,:);  % save the first 8 obs as the initial conditions
Y = data(p+1:end,:);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end


PathsCO2 = [vertcat(Y(109:end,7),r_lower(1,:)'),vertcat(Y(109:end,7),r_upper(1,:)')];
PathsCH4 = [vertcat(Y(109:end,8),r_lower(2,:)'),vertcat(Y(109:end,8),r_upper(2,:)')];
PathsNO2 = [vertcat(Y(109:end,9),y_c')];


%% Set the conditions
hind = 1; % 0 - baseline & 1 - adverse
idx_ns = find(data(1,:)'==1);


if hind == 0
    r_lower = soft(2:end,[1,3])';
    r_upper = soft(2:end,[2,4])';
    y_c = soft(2:end,5)';
else
    r_lower = soft(2:end,[6,8])';
    r_upper = soft(2:end,[7,9])';
    y_c = soft(2:end,10)';
end



PathsCO2 = [PathsCO2, vertcat(Y(109:end,7),r_lower(1,:)'),vertcat(Y(109:end,7),r_upper(1,:)')];
PathsCH4 = [PathsCH4, vertcat(Y(109:end,8),r_lower(2,:)'),vertcat(Y(109:end,8),r_upper(2,:)')];
PathsNO2 = [PathsNO2, vertcat(Y(109:end,9),y_c')];


T3 = length(PathsCO2);
T4 = length(r_lower(1,:));



figure
subplot(3,1,1)
plot(1:T3-T4,PathsCO2(1:T3-T4,1),'k','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,1),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,2),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,3),'r-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,4),'r-.','LineWidth',3);
hold on
xticks([1 4 9 14 19 24 29 34 39 44 49 54 59 64 71 76 81 86 91])
xticklabels({'1960', '1963','1968','1973','1978','1983','1988','1993','1998','2003','2008','2013','2018','2023','2030','2035','2040','2045','2050'})
xlim([1 91])
grid on
ylabel('ppm','FontSize',12)
legend('Actual','Soft Constraints Optimistic','','Soft Constraints Adverse','Location', 'northwest')
title('SSP Scenario Pathway for CO_2','FontSize',14)
subplot(3,1,2)
plot(1:T3-T4,PathsCH4(1:T3-T4,1),'k','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,1),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,2),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,3),'r-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,4),'r-.','LineWidth',3);
hold on
xticks([1 4 9 14 19 24 29 34 39 44 49 54 59 64 71 76 81 86 91])
xticklabels({'1960', '1963','1968','1973','1978','1983','1988','1993','1998','2003','2008','2013','2018','2023','2030','2035','2040','2045','2050'})
xlim([1 91])
grid on
ylabel('ppb','FontSize',12)
legend('Actual','Soft Constraints Optimistic','','Soft Constraints Adverse','Location', 'northwest')
title('SSP Scenario Pathway for CH_4','FontSize',14)
subplot(3,1,3)
plot(1:T3-T4,PathsNO2(1:T3-T4,1),'k','LineWidth',3);
hold on
plot(T3-T4:T3,PathsNO2(T3-T4:end,1),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsNO2(T3-T4:end,2),'r-.','LineWidth',3);
hold on
xticks([1 4 9 14 19 24 29 34 39 44 49 54 59 64 71 76 81 86 91])
xticklabels({'1960', '1963','1968','1973','1978','1983','1988','1993','1998','2003','2008','2013','2018','2023','2030','2035','2040','2045','2050'})
xlim([1 91])
grid on
ylabel('ppb','FontSize',12)
legend('Actual','Hard Constraint Optimistic','Hard Constraint Adverse','Location', 'northwest')
title('SSP Scenario Pathway for NO_2','FontSize',14)


%% Short Range


PathsCO2 = PathsCO2(38:end,:)
PathsCH4 = PathsCH4(38:end,:)
PathsNO2 = PathsNO2(38:end,:)


T3 = length(PathsCO2);
T4 = length(r_lower(1,:));

figure
subplot(3,1,1)
plot(1:T3-T4,PathsCO2(1:T3-T4,1),'k','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,1),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,2),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,3),'r-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCO2(T3-T4:end,4),'r-.','LineWidth',3);
hold on
xticks([1 7 12 17 22 27 34 39 44 49 54])
xticklabels({'1999','2003','2008','2013','2018','2023','2030','2035','2040','2045','2050'})
xlim([1 54])
grid on
ylabel('ppm','FontSize',12)
legend('Actual','Soft Constraints Optimistic','','Soft Constraints Adverse','Location', 'northwest')
title('SSP Scenario Pathway for CO_2','FontSize',14)
subplot(3,1,2)
plot(1:T3-T4,PathsCH4(1:T3-T4,1),'k','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,1),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,2),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,3),'r-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsCH4(T3-T4:end,4),'r-.','LineWidth',3);
hold on
xticks([1 7 12 17 22 27 34 39 44 49 54])
xticklabels({'1999','2003','2008','2013','2018','2023','2030','2035','2040','2045','2050'})
xlim([1 54])
grid on
ylabel('ppb','FontSize',12)
legend('Actual','Soft Constraints Optimistic','','Soft Constraints Adverse','Location', 'northwest')
title('SSP Scenario Pathway for CH_4','FontSize',14)
subplot(3,1,3)
plot(1:T3-T4,PathsNO2(1:T3-T4,1),'k','LineWidth',3);
hold on
plot(T3-T4:T3,PathsNO2(T3-T4:end,1),'b-.','LineWidth',3);
hold on
plot(T3-T4:T3,PathsNO2(T3-T4:end,2),'r-.','LineWidth',3);
hold on
xticks([1 7 12 17 22 27 34 39 44 49 54])
xticklabels({'1999','2003','2008','2013','2018','2023','2030','2035','2040','2045','2050'})
xlim([1 54])
grid on
ylabel('ppb','FontSize',12)
legend('Actual','Hard Constraint Optimistic','Hard Constraint Adverse','Location', 'northwest')
title('SSP Scenario Pathway for NO_2','FontSize',14)
