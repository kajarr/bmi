%% Test Complex Filter
clear 
close all
load('TrainFilter');
x = Out(1:2,:)'/2;
d = Out(5:6,:)';

z = complex(x(:,1),x(:,2));
zd = complex(d(:,1),d(:,2));

hold on
plot(x)
plot(d)

%%
w = [0;0];
W = [];
mu = 1e-20;
y = zeros(3990,1);
for k = 2:3990
%Filter output
y(k) = w'*[z(k);y(k-1)];
e(k) = zd(k) - y(k);
w = w + mu*e(k)*conj([z(k);y(k-1)]);
W = [W;w];
end

figure
hold on
plot(real(W));
plot(imag(W));
figure
hold on
plot(real(e));
plot(imag(e));
figue
hold
on
plot(real(z));
plot(real(y));