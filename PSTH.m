set(0,'DefaultFigureWindowStyle','docked')
close all
clear
load('Data')

%Time Step (ms)
dt = 0.001;
Nt = 500;
%Trial Number
N = 100; 
%Gaussian Window
x = -5:0.5:5;
w = gaussmf(x,[0.5 0]);
Nw = length(x);
Nw_half = floor(Nw/2);

figure
%Direction 1 to 8
for d = 1:8
    %Create Direction Plot
    subplot(2,4,d);
    %For each trial super-pose activiy
    y = zeros(98,500);
    for tr = 1:N
    y = y + trial(tr,d).spikes(:,1:Nt);
    end
    %For each Neuron
    %Do Convolution
    for neu = 1:98
    hold on
    %Firing Rate
    r = conv(w,y(neu,:));
    %Trim
    R(neu,:) = r(Nw_half+1:end-Nw_half);
    end
    %Plot
    imagesc(R)
    xlim([1,500])
    %Average in Time for each Neuron
    F(:,d) = mean(R,2);
    %Save to folder
    Rall(:,:,d) = R;
end

%Normalize
F = F./max(max(F));
%Sort Using Maximum
[MaxRate,Direction] = max(F,[],2);

%Create Angular Position Space
O = linspace(0,2*pi,8);
count = zeros(8,1);

%Assemple prefered direction matrix
for neu = 1:98
d = Direction(neu);
count(d) = count(d)+1;
P(neu,:) = MaxRate(neu,:).*[cos(O(d)),sin(O(d))];
end

%Create Position Plot
figure
hold on
axis square
%For Each Direction
for d = 1:8
%Plot X and V
X = zeros(2,Nt+1);
V = zeros(2,Nt);
for t = 1:Nt
    %Regress
    V(:,t) = (P'*P)^-1*P'*(Rall(:,t,d)-F(:,d));
    X(:,t+1) = X(:,t) + V(:,t).*dt;
      
end
hold on
plot(X(1,:),X(2,:))
end

