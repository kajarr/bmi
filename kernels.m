close all

x = -10:0.1:10;
gk = gaussk(x, 0, 1);
ek = expk(x, 0, 1);
ck2 = cauk2(x, -1, 2);

%Kernels
figure
plot(x,gk,'b')
hold on
plot(x,ek,'r')
hold on
plot(x,ck2,'m')
hold on
% plot(x, [gk ek ck2]);
title('Kernels Distribution')
xlabel('Time')
ylabel('Distribution')

%Gaussian kernel
function f = gaussk(x, mu, s)
ex = (-1/2)*((x - mu)/s).^2;
ar = (s * sqrt(2*pi));
f = (1/ar)*exp(ex);
end
%Exponential kernel
function f = expk(x, mu, s)
ex = (-sqrt(2))*abs((x-mu)/s);
ar = (s * sqrt(2));
f = (1/ar)*exp(ex);
end
%Causal kernel
function hw = cauk(x,mu,a)
q = ((a^2)*(x-mu));
w = exp((-a).*(x-mu));
    for i = 1:length(q)
        f(i)=q(i)*w(i);
    end
hw = zeros(size(f));
    for i = 1:length(f)
        if f(i) > 0.0
           hw(i) = f(i);
       else
           hw(i) = 0.0;
        end
    end
end
