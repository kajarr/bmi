%% Main File
close all
clear
set(0,'DefaultFigureWindowStyle','docked')

%% DataSet
load('Data')
Tuning = DataSet(trial,50);

%% Preprocessing
%Eliminate Neuron
Tuning = EliminateUnit(Tuning,[38,49,52,76],'Spikes');
%Firing Rate
x = -5:0.5:5;
w = gaussmf(x,[2.5 0]);
Tuning = Convolution(Tuning,w);
%Plot
PlotTuning(Tuning)

%% Regression
%Preferred Direction
B = GetPreDirection(Tuning,0);
%AutoRegressive Model
%A = AutoRegression(Tuning,0);
%Residual Errors
%[Q,R] = Residual(Tuning,A,B);
%save('Noise','Q','R');
load('Noise')

%% Kalman Filter
%Prepare
% R = diag(diag(R));
T = [1 0 -1 0;0 1 0 -1];
H = B(:,2:3)*T;
zo = B(:,1);
A = eye(4);
Q = eye(4);
R = diag(diag(R));
%Kalman Filter
K = KalmanDecoder(A,Q,H,R);


%% Test
h1 = figure;
hold on
axis equal
X = [];
for d = 1:Tuning.Nd
F = Tuning.Train{1,d}.FiringRate;
x = [0,0,0,0]';
    %Iterate
    for t = 2:Tuning.Nt
        z = F(:,t) - zo;
        [K,v] = KalmanUpdate(K,x,z); 
        X = [X,x];
    end
    figure(h1)
    hold on
    plot(X(1,:),X(2,:))
end







