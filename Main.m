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
Na = 2;
A = AutoRegression(Tuning,0,Na);
%Residual Errors
[Q,R] = Residual(Tuning,A,B);
save('Noise','Q','R');
%load('Noise')

%% Kalman Filter
%Prepare
T = [1 0 -1  0;
     0 1  0 -1];
T = [T,zeros(2,2*(Na-2))];
H = B(:,2:3)*T;
zo = B(:,1);
%Kalman Filter
K = KalmanDecoder(A,Q,H,R);


%% Test
h2 = figure;
h1 = figure;
axis equal
Out = [];
for n = 1:9
    for d = 1:Tuning.Nd
    X = [];
    F = Tuning.Test{n,d}.Subspace;
    x = zeros(2*Na,1);
        %Iterate
        for t = 1:Tuning.Nt
            z = F(:,t)-zo;
            [K,x] = KalmanUpdate(K,x,z); 
            X = [X,x];
        end
        figure(h1);
        subplot(3,3,n)
        hold on
        X = 2*X;
        plot(X(1,:),X(2,:),'r')
        Xreal = Tuning.Test{n,d}.Position;
        plot(Xreal(1,:),Xreal(2,:),'b');
        Out = [Out,[X;Xreal(:,1:end-1)]];
    end 
figure(h2);
subplot(3,3,n)
hold on
V = diff(X(1:2,:)')';
plot(V','r');
Vreal = Tuning.Test{n,d}.Velocity;
plot(Vreal','b');
end