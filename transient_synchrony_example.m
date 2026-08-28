clear
close all
clc

%================= global parameters ======================%

L = 15000;                  % half domain length
N = 100000;                 % number of spatial points

x = linspace(-L,L-(2*L/N),N);
dx = x(2)-x(1);

numsteps = 400;            % number of time steps
num_sims = 10;              % simulations per r value

% R values to compare
rvals = [0.01 0.1 0.5 1 1.5 1.9 1.99];

% Dispersal kernel
alpha = 4;
k = exp(-x.^2/(2*alpha^2));
k = k/sum(k);

% Environmental driver correlation structures
alpha1 = 0.5;
alpha2 = 0.5;

sigma1 = 1e-3;
sigma2 = 1e-3;

lambda1 = 1;
phi = 1;

g1 = exp(-x.^2/(2*(lambda1/sqrt(2))^2));
g1 = g1/sqrt(sum(g1.^2));

g2 = exp(-phi*abs(x));
g2 = g2/sqrt(sum(g2.^2));


%================= storage ================================%

% Each row = one r value
% Each column = one time step

lEmp_mean = zeros(length(rvals),numsteps);

% Optional: store every individual simulation
lEmp_all = zeros(length(rvals),num_sims,numsteps);


%================= main loop ===============================%

for ir = 1:length(rvals)

    r = rvals(ir);

    fprintf('\nRunning r = %.4f\n',r)

    for isim = 1:num_sims

        fprintf('  Simulation %d/%d\n',isim,num_sims)

        % Make each simulation reproducible but different
        rng(isim)

        %---------------- simulation ----------------%

        n = zeros(numsteps+1,N);
        n(1,:) = 1.1;

        for t = 1:numsteps

            % Generate random environmental noise
            b1 = randn(1,N);
            b2 = randn(1,N);

            % Apply correlation structures
            z1 = sigma1*real(ifft(fft(g1).*fft(b1)));
            z2 = sigma2*real(ifft(fft(g2).*fft(b2)));

            % Environmental forcing + Ricker growth
            forcing = n(t,:).*exp(r*(1-n(t,:)) ...
                      + alpha1*z1 + alpha2*z2);

            % Dispersal
            n(t+1,:) = real(ifft(fft(k).*fft(forcing)));

        end


        %---------------- calculate synchrony scale ----------------%

        maxdist = 100;
        maxlag = round(maxdist/dx);
        s = (0:maxlag)*dx;

        lEmp = zeros(1,numsteps);

        for t = 1:numsteps

            data = n(t+1,:) - mean(n(t+1,:));

            % correlation function
            C = zeros(1,maxlag+1);

            for h = 0:maxlag

                C(h+1) = mean(data.*circshift(data,-h));

            end

            % truncation
            tail = C(round(0.8*length(C)):end);

            threshold = 3*std(tail);
            window = 5;

            idx_cut = length(C);

            for i = 1:length(C)-window+1

                if all(abs(C(i:i+window-1)) < threshold)

                    idx_cut = i;
                    break

                end

            end

            % calculate spatial scale
            l2 = sum((s(1:idx_cut).^2).*C(1:idx_cut)) ...
                 / sum(C(1:idx_cut));

            lEmp(t) = sqrt(abs(l2));

        end

        % Save this simulation
        lEmp_all(ir,isim,:) = lEmp;

    end

    % Average the 10 simulations at every time step
    lEmp_mean(ir,:) = mean(lEmp_all(ir,:,:),2);

end


%================= plot ================================%

figure
hold on

for ir = 1:length(rvals)

    plot(1:numsteps, ...
         lEmp_mean(ir,:), ...
         'LineWidth',2, ...
         'DisplayName',sprintf('r = %.2f',rvals(ir)))

end

xlabel('Time')
ylabel('Spatial synchrony scale')
legend('Location','best')
box on
set(gca,'FontSize',14)