% r sweep configuration, 100 simulation per regime

bounds = [ ...
    0.001, ...   % Period 1
    2.0000, ...  % Period 2
    2.5265, ...  % Period 4
    2.6564, ...  % Period 8
    2.6846, ...  % Period 16
    2.6907, ...  % Period 32
    2.6920, ...  % Period 64
    2.6923, ...  % Period 128
    2.7 ...      % End
];


rvals = [];

for i = 1:length(bounds)-1

    regime_vals = linspace(bounds(i),bounds(i+1),100);

    rvals = [rvals regime_vals];

end


fprintf('Total R-values: %d\n',length(rvals))


%parameters

epsilon = 1e-6; %perturbation
tol = 1e-10; %alignment tolerance

burn = 5000; %burn in period
numsteps = 10000; %number of simulation time steps

relaxation_time = nan(size(rvals)); %storage




% simulations

for ir = 1:length(rvals)

    r = rvals(ir);

    fprintf('Running r = %.5f\n',r)

    % burn in period
    N0 = 1.2;
    for t = 1:burn
        N0 = N0*exp(r*(1-N0));
    end

    % steady trajectory
    steady = zeros(1,numsteps);
    N = N0;
    for t = 1:numsteps
        N = N*exp(r*(1-N));
        steady(t)=N;
    end

    % Perturbed trajectory

    N = N0 + epsilon;
    perturbed = zeros(1,numsteps);
    for t = 1:numsteps
        N = N*exp(r*(1-N));
        perturbed(t)=N;
    end

    % calculate the time it takes for alignment between steady and
    % perturbed
    err = abs(perturbed-steady);
    idx = find(err < tol,1,'first');
    if ~isempty(idx)
        relaxation_time(ir)=idx;
    end
end

%plotting relaxtion times against r valeus gives the relaxation times result in paper