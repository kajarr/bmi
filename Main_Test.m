%% Main_Test
close all
clear
set(0,'DefaultFigureWindowStyle','docked')

teamname = 'CDT';
[RMSE,Xreal,X] = testFunction_for_students_MTb(teamname);

[~,Nx] = size(X);
X = [ones(1,Nx);X];

M = Xreal*X'*(X*X')^-1;

            