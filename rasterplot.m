%Initiate
clear
close
load('Data')

%Direction
for k = 1:8
    
    subplot(2,4,k);
    y = zeros(98,500);
    
    %Trial
    for n = 1:100
    y = y + trial(n,k).spikes(:,1:500);
    end
    
    y = y/100;
    imagesc(y);
    ylabel('Neuron')
    xlabel('Time')
end





