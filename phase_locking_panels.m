%% parameters
L = 15000; %half domain length
N = 100000; %number of points for discretization

x = linspace(-L,L-(2*L/N),N); %domain
r = 2.6; %r value, corresponding to period 4 cycle in ricker map

%simulation length + shock time
numsteps = 100;
shock_time = 20;

% Columns: environmental shock variance
sigma_cases = [0.25 0.5 1];

% Rows: dispersal length scale
alpha_cases = [1 12.5 25 ];

%coefficient for environmental driver
alpha1 = 0.5;

% Spatial region to DISPLAY only
x_window = 1000;


% Spatially correlated noise for environmental shock

rng(1); %same random each time, so we can compare cases

lambda_shock = 20; %length scale of environmental shock correlation structure

g = exp(-x.^2/(2*lambda_shock^2));
g = g/sqrt(sum(g.^2));

shock = real(ifft(fft(g).*fft(randn(1,N))));

shock = shock - mean(shock);
shock = shock/std(shock);


% find determinstic attractor for later rounding use.
Ndet = 1.1;

for t = 1:10000
    Ndet = Ndet * exp(r*(1-Ndet));
end

orbit = zeros(1,200);

for t = 1:200
    Ndet = Ndet * exp(r*(1-Ndet));
    orbit(t) = Ndet;
end

states = unique(round(orbit(end-50:end),8));

fprintf('r = %.2f, number of phases = %d\n', ...
    r,length(states));


% simulations

results = cell(length(alpha_cases),length(sigma_cases));

for ia = 1:length(alpha_cases)

    alpha = alpha_cases(ia);

    k = exp(-x.^2/(2*alpha^2));
    k = k/sum(k);

    for is = 1:length(sigma_cases)

        sigma = sigma_cases(is);

        n = zeros(numsteps+1,N);

        n(1,:) = 1.1;

        for t = 1:numsteps

            if t == shock_time
                z = sigma * shock;
            else
                z = zeros(1,N);
            end

            forcing = n(t,:) .* ...
                exp(r*(1-n(t,:)) + alpha1*z);

            n(t+1,:) = real(ifft(fft(k).*fft(forcing)));

        end

        results{ia,is} = n;

    end
end


% phase rounding

phase_results = cell(size(results));

for ia = 1:length(alpha_cases)

    for is = 1:length(sigma_cases)

        n = results{ia,is};

        distances = zeros([size(n),length(states)]);

        for j = 1:length(states)
            distances(:,:,j) = abs(n - states(j));
        end

        [~,phase] = min(distances,[],3);

        phase_results{ia,is} = phase;

    end
end


%plotting
space_idx = abs(x) <= x_window;
time_idx = 70:100;
figure('Units','inches','Position',[3 3 20 14])
tl = tiledlayout(length(alpha_cases),length(sigma_cases), ...
    'TileSpacing','compact', ...
    'Padding','loose');
for ia = 1:length(alpha_cases)
    for is = 1:length(sigma_cases)
        nexttile
        imagesc( ...
            x(space_idx), ...
            time_idx, ...
            phase_results{ia,is}(time_idx,space_idx));
        set(gca,'YDir','normal')
        clim([0.5 length(states)+0.5])
        xlim([-x_window x_window])
        ylim([70 100])
        if is ~= 1
            set(gca,'YTickLabel',[])
        end
        if ia ~= length(alpha_cases)
            set(gca,'XTickLabel',[])
        end
        ax = gca;
        ax.FontSize = 24;
        ax.LineWidth = 1.5;
    end
end

%other plotting + labels can be used...