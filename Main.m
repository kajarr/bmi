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
w = gaussmf(-5:0.5:5,[5 0]);
Tuning = Convolution(Tuning,w);
%Plot
PlotTuning(Tuning)

%% Regression
%Preferred Direction
B = GetPreDirection(Tuning,0);
%AutoRegressive Model
A = AutoRegression(Tuning,0);
%Residual Errors
%[Q,R] = Residual(Tuning,A,B);
%save('Noise','Q','R');
load('Noise')

%% Kalman Filter
%Prepare
[Na,~] = size(A);
A = [zeros(Na,1),A];
A = [zeros(1,Na+1);A];
A(1,1) = 1;
Q = [zeros(Na,1),Q];
Q = [zeros(1,Na+1);Q];
R = diag(diag(R));
H = B';
%Kalman Filter
K = KalmanDecoder(A,Q,H,R);
A = eye(3);
Q = eye(3);

%% Test
h1 = figure;
h2 = figure;
hold on
axis equal
for d = 1:Tuning.Nd
F = Tuning.Test{1,d}.FiringRate;
v = [0,0,0]';
    %Iterate
    V = zeros(3,Tuning.Nt);
    X = zeros(3,Tuning.Nt+1);
    for t = 1:Tuning.Nt
        z = F(:,t);
        [K,v] = KalmanUpdate(K,v,z); 
        V(:,t) = v;
        X(:,t+1)  = X(:,t) + v;
    end
    figure(h1)
    hold on
    plot(X(2,:),X(3,:))
    figure(h2)
    subplot(2,4,d)
    plot(V')
end







