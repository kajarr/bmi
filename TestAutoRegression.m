%% Test Auto-Regression
close all
clear 
t = 1:100;
T = 100;
t = t/T;
P1 = [];
P2 = [];
figure
hold on
axis equal
for d = 1:8 
O = linspace(0,2*pi,8);
X = [0,cos(O(d))];
Y = [0,sin(O(d))];
x = X(1)+(X(2)-X(1));
y = Y(1)+(Y(2)-Y(1));
p = [x;y]*(10*t.^3 - 15*t.^4 + 6*t.^5);
p1 = awgn(p,30,'measured');
p2 = awgn(p,30,'measured');
plot(p1(1,:),p1(2,:));
P1 =[P1,p1];
P2 =[P2,p2];
end

n = 2;
[S,s] = DataSet.Shift(P1,n);
A = S*pinv(s);
[~,s] = DataSet.Shift(P2,n);
p_est = A*s;
plot(p_est(1,:),p_est(2,:));
