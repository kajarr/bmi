%%% Team Members: WRITE YOUR TEAM MEMBERS' NAMES HERE
%%% BMI Spring 2015 (Update 17th March 2015)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         PLEASE READ BELOW            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function positionEstimator has to return the x and y coordinates of the
% monkey's hand position for each trial using only data up to that moment
% in time.
% You are free to use the whole trials for training the classifier.

% To evaluate performance we require from you two functions:

% A training function named "positionEstimatorTraining" which takes as
% input the entire (not subsampled) training data set and which returns a
% structure containing the parameters for the positionEstimator function:
% function modelParameters = positionEstimatorTraining(training_data)
% A predictor named "positionEstimator" which takes as input the data
% starting at 1ms and UP TO the timepoint at which you are asked to
% decode the hand position and the model parameters given by your training
% function:

% function [x y] = postitionEstimator(test_data, modelParameters)
% This function will be called iteratively starting with the neuronal data 
% going from 1 to 320 ms, then up to 340ms, 360ms, etc. until 100ms before 
% the end of trial.

% Place the positionEstimator.m and positionEstimatorTraining.m into a
% folder that is named with your official team name.

% Make sure that the output contains only the x and y coordinates of the
% monkey's hand.


function [x, y, model] = positionEstimator(test_data, model)
%% clock
[~,k] = size(test_data.spikes);
%% Model
dt = 20;
w = model.window;
P = model.PCA;
K = model.KFilter;
zo = model.BaseFiring;

%% Initialisation
if k == 320
%Initialise State
xo = test_data.startHandPos;
model.state = [[xo;0];[xo;0]];
end

% Reading 
S = test_data.spikes;
% Firing Rate
f = GetFiringRate(S,w,k,dt);
% PCA;
z = P*f-zo;
%KalmanFiltering
X = model.state;
[K,X] = KalmanUpdate(K,X,z);
%Appending model
model.state = X;
model.KFilter = K;
% Output
%r = [1;1;X(1);X(2)];
%p = M*r;
p = [X(1),X(2)];
x = p(1);
y = p(2);

end

%% Firing Rate
function [f] = GetFiringRate(S,w,k,dt)
%SubSelect
s = S(:,k-2*dt+1:end);
[Ns,~] = size(S);
C = zeros(Ns,2*dt);
    for t = 1:dt
        c = zeros(Ns,1);
        for l = 0:dt-1
            c = c + s(:,t+l)*w(:,l+1);
        end
        C(:,t+dt/2) = c; 
    end
    f = mean(C,2);
end
