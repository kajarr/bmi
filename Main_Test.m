%% Main_Test
close all
clear
set(0,'DefaultFigureWindowStyle','docked')

teamname = 'CDT';
N = eye(2);
save('Correction','N')

[RMSE,Xreal,X] = testFunction_for_students_MTb(teamname);

% Nnew = Xreal*X'*(X*X')^-1;
% Vr = diff(Xreal')';
% V = diff(X')';
% [~,Nv] = size(V);
% Vn = [];
% Vnr = [];
% for i = 1:Nv
%     Vn = [Vn,norm(V(:,i))];
%     Vnr = [Vnr,norm(Vr(:,i))];
% end
% Vr = [ones(1,Nv);Vr;Vnr];
% V = [V;Vn];
% M = V*Vr'*(Vr*Vr')^-1;






            