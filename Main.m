%Main File
close all
clear all
set(0,'DefaultFigureWindowStyle','docked')
load('Data')
%% Mean Time
Train = DataSet(trial,[1 50]);
Train = EliminateUnit(Train,[38,49,52,76],'Spikes');

Test = DataSet(trial,[51 100]);
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
w = gaussmf(-5:0.5:5,[-5:0.5:5 0]);
Train = Convolution(Train,w,'Spikes','FiringRate');
%Train = BaseLineNormalisation(Train,'FiringRate',300);
Test = Convolution(Test,w,'Spikes','FiringRate');
%Test = BaseLineNormalisation(Test,'FiringRate',300);

figure
for d = 1:Train.Nd
    subplot(2,4,d)
    F = Train.Dir{d}.FiringRate;
    imagesc(F)
end
%% Preferred Direction
[W,E] = GetPreDirection(Train,[0,0]);

figure
bar(mean(abs(E)))
%% Inverse the Weight
%Final Position Error 
MeanTest(W,Test)



