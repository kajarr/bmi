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
%PCA
Tuning = PCA(Tuning);


%% Regression
%Preferred Direction
B = GetPreDirection(Tuning,0);
%AutoRegressive Model
A = AutoRegression(Tuning,0,2);
%Residual Errors
[Q,R] = Residual(Tuning,A,B);
%save('Noise','Q','R');
%load('Noise')

%% Kalman Filter
%Prepare
% R = diag(diag(R));
T = [1 0 -1 0;0 1 0 -1];
H = B(:,2:3)*T;
zo = B(:,1);
%Kalman Filter
K = KalmanDecoder(A,Q,H,R);


%% Test
h1 = figure;
axis equal
for n = 1:9
    for d = 1:Tuning.Nd
    X = [];
    F = Tuning.Test{n,d}.Subspace;
    x = [0,0,0,0]';
        %Iterate
        for t = 1:Tuning.Nt
            z = F(:,t)-zo;
            [K,x] = KalmanUpdate(K,x,z); 
            X = [X,x];
        end
        figure();
        subplot(3,3,n)
        hold on
        plot(X(1,:),X(2,:))
    end
end







