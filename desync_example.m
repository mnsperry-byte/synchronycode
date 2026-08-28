%params
r = 2.1;
x = 0:20;

% Find attractor point
N0 = 1.1;

for k = 1:10000
    N0 = N0*exp(r*(1-N0));
end

% Control trajectory
y = zeros(size(x));
y(1) = N0;

for i = 1:length(x)-1
    y(i+1) = y(i)*exp(r*(1-y(i)));
end


% Perturbed trajectory
h = zeros(size(x));
h(1) = N0;

for i = 1:length(x)-1

    if x(i) == 10
        h(i+1) = h(i)*exp(r*(1-h(i))+0.7);
    else
        h(i+1) = h(i)*exp(r*(1-h(i)));
    end

end

figure;

plot(x,y,'-o','LineWidth',1,'MarkerSize',5)
hold on
plot(x,h,'-o','LineWidth',1,'MarkerSize',5)

xline(10,'--k','Shock','HandleVisibility','off')

xlabel('t','FontSize',12)
ylabel('N','FontSize',12)

legend('Control','Perturbation','FontSize',12,'Location','northeast')

set(gca,'FontSize',12,'LineWidth',0.5)

ylim([0 2])