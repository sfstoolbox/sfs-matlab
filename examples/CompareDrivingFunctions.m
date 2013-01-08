% compare different Driving Functions
clc 
clear all
close all

%% Use debug = 1 to plot the results, otherwise use debug = 0
debug = 0;
%% Load Driving Function Results
DrivingFunctionWithoutWeights = load('DrivingFunctionWithoutWeights.mat');
DrivingFunctionWithPointWeights = load('DrivingFunctionWithPointWeights.mat');
DrivingFunctionWithAllWeights = load('DrivingFunctionWithAllWeights.mat');
phi = load('phi.mat');
phi = phi.phi;

DrivingFunctionWithoutWeights = 20*log10(abs(DrivingFunctionWithoutWeights.D_plot(1,:))); 
DrivingFunctionWithPointWeights = 20*log10(abs(DrivingFunctionWithPointWeights.D_plot(1,:)));
DrivingFunctionWithAllWeights = 20*log10(abs(DrivingFunctionWithAllWeights.D_plot(1,:)));

%% Plot results with 'point plot'
if debug
figure
subplot(1,2,1)
plot(phi,DrivingFunctionWithoutWeights,'r.')
hold on
plot(phi,DrivingFunctionWithPointWeights,'b.')
hold on
plot(phi,DrivingFunctionWithAllWeights,'g.')
grid on
legend('D without weights','D with weights from equally spaced points','D with weights from equally spaced points and surface weights','location','South')
ylabel('amplitude / db')
xlabel('$\varphi$ / rad','interpreter','latex')

subplot(1,2,2)
plot(phi,DrivingFunctionWithoutWeights-DrivingFunctionWithPointWeights,'r.')
hold on
plot(phi,DrivingFunctionWithoutWeights-DrivingFunctionWithAllWeights,'b.')
hold on
plot(phi,DrivingFunctionWithPointWeights-DrivingFunctionWithAllWeights,'g.')
grid on
legend('D_{error} between no weighting and point weights','D_{error} between no weights and point+surface weights','D_{error} between point weights and point+surface weights')
ylabel('amplitude / db')
xlabel('$\varphi$ / rad','interpreter','latex')
end
%% sort Driving Function values ascending
[DrivingFunctionWithoutWeights,idx] = sort(DrivingFunctionWithoutWeights,'ascend');
DrivingFunctionWithPointWeights = sort(DrivingFunctionWithPointWeights,'ascend');
DrivingFunctionWithAllWeights = sort(DrivingFunctionWithAllWeights,'ascend');
phi = phi(1,idx);
%% plot results with 'line plot'
if debug
    
figure
subplot(1,2,1)
plot(DrivingFunctionWithoutWeights,phi,'r')
hold on
plot(DrivingFunctionWithPointWeights,phi,'b')
hold on
plot(DrivingFunctionWithAllWeights,phi,'g')
grid on
legend('D without weights','D with weights from equally spaced points','D with weights from equally spaced points and surface weights')
xlabel('amplitude / db')
ylabel('$\varphi$ / rad','interpreter','latex')

subplot(1,2,2)
plot(DrivingFunctionWithoutWeights-DrivingFunctionWithPointWeights,phi,'r')
hold on
plot(DrivingFunctionWithoutWeights-DrivingFunctionWithAllWeights,phi,'b')
hold on
plot(DrivingFunctionWithPointWeights-DrivingFunctionWithAllWeights,phi,'g')
grid on
legend('D_{error} between no weighting and point weights','D_{error} between no weights and point+surface weights','D_{error} between point weights and point+surface weights')
xlabel('amplitude / db')
ylabel('$\varphi$ / rad','interpreter','latex')
end
