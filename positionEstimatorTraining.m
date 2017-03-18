
function [model] = positionEstimatorTraining(training_data)
%% Create DataSet
disp('Creating Model')
disp('Preparing DataSet')
Tuning = DataSet(training_data,30);
disp('Tuning object added')

%% Preprocessing
disp('Starting preprocessing')
disp('Compute: offline convolution')
%Firing Rate
x = linspace(-5,5,20);
w = gaussmf(x,[2.5 0]);
Tuning = Convolution(Tuning,w);
disp('FiringRate field created')
%Averaging
dt = 20;
Tuning = MeanTimeStep(Tuning,dt);
%PCA
disp('Compute: PCA')
[Tuning,P] = PCA(Tuning);
disp('Subspace field created')

%% Regression
disp('Starting Regression')
%Preferred Direction
disp('Compute: preferred direction')
B = GetPreDirection(Tuning,0);
disp('B added')
%AutoRegressive Model
disp('Compute: autoregressive model')
Na = 2;
A = AutoRegression(Tuning,0,Na);
disp('A added')
%Residual Errors
disp('Compute: residuals')
[Q,R] = Residual(Tuning,A,B);
disp('Q,R added')

%% Kalman Filter Assembly 
disp('Starting Kalaman filter assembly')
%Observation
disp('Compute: observation and base firing rate')
[H,zo] = KalmanDecoder.Observation(B,Na,dt);
disp('H,zo added')
%Kalman Filter
disp('Creating Kalman filter')
K = KalmanDecoder(A,Q,H,R);
disp('K object added')
close 2 3 4 5 6 7
model = struct;
model.KFilter = K;
model.PCA = P;
model.BaseFiring = zo;
model.time = 320;
model.window = w;
end