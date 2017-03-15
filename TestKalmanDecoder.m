%% Test Kalman Decoder
clear
set(0,'DefaultFigureWindowStyle','docked')

%% Create State Space
%Simple Position and Velocity
A = eye(3);
Q = eye(3);
P = eye(3);
%% Regression on Average Trial Data
load('Data')
%Extract Data
Train = DataSet(trial,[1 50]);
Train = EliminateUnit(Train,[38,49,52,76],'Spikes');
%Convolve Firing Rate
v = -5:0.5:5;
w = gaussmf(-5:0.5:5,[2.5 0]);
Train = Convolution(Train,w,'Spikes','FiringRate');
%Get Preferred Direction
[~,B,E] = GetPreDirection(Train,[1,1]);
load('Noise');
%[Fe] = DataSet.TrialTest(B,w,trial,[52 100]);
% R = cov(Fe);
% save('Noise','R')

%% Generate Test Data
Test = DataSet(trial,[51 51]);
Test = EliminateUnit(Test,[38,49,52,76],'Spikes');
%Convolve Firing Rate
Test = Convolution(Test,w,'Spikes','FiringRate');

%% Create Filter Object
%Prepare Matrices
H = B(1:3,:);
H = H';
zo = -B(1,:)';
R = diag(diag(cov(E)));
%Filter
K = KalmanDecoder(A,Q,H,R);

%% Test
h1 = figure;
h2 = figure;
hold on
axis equal
for d = 1:Train.Nd
F = Test.Dir{d}.FiringRate;
v = [0,0,0]';
    %Iterate
    V = zeros(3,Test.Nt);
    X = zeros(3,Test.Nt+1);
    for t = 1:Test.Nt
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


   


 