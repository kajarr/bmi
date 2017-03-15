%Main File
close all
clear
set(0,'DefaultFigureWindowStyle','docked')
%% Mean Time
load('Data')

Train = DataSet(trial,[1 50]);
Train = EliminateUnit(Train,[38,49,52,76],'Spikes');

Test = DataSet(trial,[51 52]);
Test = EliminateUnit(Test,[38,49,52,76],'Spikes');

axis equal
hold on
for d = 1:Train.Nd
    X = Train.Dir{d}.Position;
    plot(X(1,:),X(2,:));
end
figure
for d = 1:Train.Nd
subplot(2,4,d)
plot(Train.Dir{d}.Velocity')
end

%% Firing Rate
x = -5:0.5:5;
w = gaussmf(-5:0.5:5,[5 0]);
Train = Convolution(Train,w,'Spikes','FiringRate');
Test = Convolution(Test,w,'Spikes','FiringRate');

figure
for d = 1:Train.Nd
    subplot(2,4,d)
    F = Train.Dir{d}.FiringRate;
    imagesc(F)
end
%% Preferred Direction
%Regression
[W,B,E] = GetPreDirection(Train,[1,1]);

%% Inverse the Weight
%Final Position Error 
Etrack = MeanTest(W,Test);

%% Statistic
[Fe] = DataSet.TrialTest(B,w,trial,[52 100]);
Q = cov(Fe);
figure
surf(V)



