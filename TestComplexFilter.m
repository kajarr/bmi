%% Test Complex Filter
clear 
close all
load('Complex');
x = x';
d = d';
subplot(1,2,1)
plot(x);
subplot(1,2,2)
plot(d);

figure
z = complex(x(:,1),x(:,2));
zd = complex(d(:,1),d(:,2));
hold on
plot(z)
plot(zd)

z = [z,ones(15,1)+1i*ones(15,1)];
z = z';
d = d';

%%
w = [0;0];
W = [];
mu = 0.005;
y = [];
e = [];
for k = 1:15
%Filter output
y(k) = w'*z(:,k);
e(k) = zd(k) - y(k);
w = w + mu*e(k)*conj(z(:,k));
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