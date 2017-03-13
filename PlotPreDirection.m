function [] = PlotPreDirection (P)
hold on
for i = 1:98
    plot([0 P(i,1)],[0 P(i,2)])
end
end