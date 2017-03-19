%% Test PCA
%% Main File
close all
clear
set(0,'DefaultFigureWindowStyle','docked')

%% DataSet
load('Data')
Tuning = DataSet(trial,50);

%% Preprocessing
%Firing Rate
x = -5:0.1:5;
w = gaussmf(x,[2.5 0]);
Tuning = Convolution(Tuning,w);
%Plot
PlotTuning(Tuning)

%% PCA
F = []
for d = 1:Tuning.Nd
    f = Tuning.Train{1,d}.FiringRate;
    F = [F,f];
end
figure
imagesc(F)
Sf = cov(F');
figure
imagesc(Sf);
[P,L] = eig(Sf);
L = diag(L);
figure
plot(L);
