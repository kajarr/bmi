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
[~,t] = size(test_data.spikes);
%% Model
dt = 20;
w = model.window;
P = model.PCA;
K = model.KFilter;
zo = model.BaseFiring;

%% Pre-Cue
if t == 320
%Initialise State
xo = test_data.startHandPos;
model.state = [xo;xo];
model.KFilter = K;
% Output
x = xo(1);
y = xo(2);
%% Post-Cue
else
% Reading 
s = test_data.spikes(:,t-dt+1:end);
% Firing Rate
f = GetFiringRate(s,w);
% PCA;
z = P*f-zo;
%KalmanFiltering
X = model.state;
[K,X] = KalmanUpdate(K,X,z);
%Appending model
model.state = X;
model.KFilter = K;
% Output
x = 10*X(1);
y = 10*X(2);
end
end

%% Firing Rate
function [f] = GetFiringRate(s,w)
Nw = length(w);
r = [];
for i = 1:98
    r = [r ; conv(w,s(i,:))];
end
f = r(:,Nw/2+1:end-Nw/2);
f = mean(f,2);
end
