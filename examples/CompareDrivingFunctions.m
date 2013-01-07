% compare different Driving Functions
clc 
clear all
close all

DrivingFunctionWithoutWeights = load('DrivingFunctionWithoutWeights.mat');
DrivingFunctionWithPointWeights = load('DrivingFunctionWithPointWeights.mat');
DrivingFunctionWithAllWeights = load('DrivingFunctionWithAllWeights.mat');

DrivingFunctionWithoutWeights = 20*log10(abs(DrivingFunctionWithoutWeights.D_plot(1,:)))-40; 
DrivingFunctionWithPointWeights = 20*log10(abs(DrivingFunctionWithPointWeights.D_plot(1,:)));
DrivingFunctionWithAllWeights = 20*log10(abs(DrivingFunctionWithAllWeights.D_plot(1,:)));

figure
plot(DrivingFunctionWithoutWeights,'r')
hold on
plot(DrivingFunctionWithPointWeights,'b')
hold on
plot(DrivingFunctionWithAllWeights,'g')
grid on
legend('D without weights','D with weights from equally spaced points',' D with weights from equally spaced points and surface weights')
ylabel('amplitude / db')

figure
plot(DrivingFunctionWithoutWeights-DrivingFunctionWithPointWeights,'r')
hold on
plot(DrivingFunctionWithoutWeights-DrivingFunctionWithAllWeights,'b')
hold on
plot(DrivingFunctionWithPointWeights-DrivingFunctionWithAllWeights,'g')
grid on
legend('D_{error} between no weighting and point weights','D_{error} between no weights and point+surface weights','D_{error} between point weights and point+surface weights')
ylabel('amplitude / db')