clear
close all
hold on
%Test Convolution
Nt = 361;
Nn = 100;
Test = round(rand(Nn,Nt));
stem(Test')
x = linspace(-5,5,20);
w = gaussmf(x,[1 0]);
Nw = 20;
dt = 20;
F = zeros(Nn,dt);
f = [];
for k = 40:dt:Nt
    F = [F,zeros(Nn,dt)];
    %Input
    S = Test(:,1:k);
    %SubSelect
    s = S(:,k-2*dt+1:end);
    C = zeros(Nn,2*dt);
    for t = 1:dt
        c = zeros(Nn,1);
        for l = 0:19
            c = c + s(:,t+l)*w(:,l+1);
        end
        C(:,t+dt/2) = c; 
    end
    f = [f,mean(C,2)];
    F(:,k-2*dt+1:end) = F(:,k-2*dt+1:end) + C; 
end
plot(F');