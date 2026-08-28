%================= global parameters ======================%

L = 15000; % half domain length
N = 100000; % Number of points for discretization

x = linspace(-L,L-(2*L/N),N); % grid
dx = x(2)-x(1); %caldulate grid size

numsteps = 1000; %number of steps for propogation


% Dispersal kernel (Gaussian)
alpha = 4;
k = exp(-x.^2/(2*alpha^2));
k = k/sum(k);


% Environmental driver correlation structures
alpha1 = 0.5; 
alpha2 = 0.5;

sigma1 = 0.001;
sigma2 = 0.001;

lambda1 = 1;
phi = 1;

g1 = exp(-x.^2/(2*(lambda1/sqrt(2))^2));
g1 = g1/sqrt(sum(g1.^2));

g2 = exp(-phi*abs(x));
g2 = g2/sqrt(sum(g2.^2));


% ================ simulation infrastructure ==================%

%Regimes
bounds = [ ...
    0.001, ...   % Start of Period 1 (Transcritical bifurcation)
    2.0000, ...  % Start of Period 2
    2.5265, ...  % Start of Period 4
    2.6564, ...  % Start of Period 8
    2.6846, ...  % Start of Period 16
    2.6907, ...  % Start of Period 32
    2.6920, ...  % Start of Period 64
    2.6923, ...  % Start of Period 128
    2.7 ...      % End of Period 128 / Approximation of next flip
];

rvals = [];

for i = 1:length(bounds)-1

    regime_vals = linspace(bounds(i),bounds(i+1),5);

    rvals = [rvals regime_vals];

end   % 1 simulation for each regime

fprintf('Total R values = %d\n',length(rvals))


% main loop storage 
lTheory = nan(size(rvals)); %space for theoretical prediction
lEmp = nan(size(rvals)); %space for simulation result 
period_store = nan(size(rvals)); %space for assigning an empircally inferred period to each simulation


% ===================== main loop =========================

%shell 1: loop through r vals
for ir = 1:length(rvals)

    r = rvals(ir);

    fprintf('\nRunning r = %.4f\n',r)
    
    % shell 2: simulate rval dynamics
    %initialization
    n = zeros(numsteps+1,N);

    n(1,:) = 1.1; 


    for t = 1:numsteps

        %generate random noise each step
        b1 = randn(1,N); 
        b2 = randn(1,N);

        %apply correlation structure
        z1 = sigma1*real(ifft(fft(g1).*fft(b1)));
        z2 = sigma2*real(ifft(fft(g2).*fft(b2)));
        
        %define forcing
        forcing = n(t,:).*exp(r*(1-n(t,:)) ...
                  + alpha1*z1 + alpha2*z2);

        % convolve with dispersal to propogate
        n(t+1,:) = real(ifft(fft(k).*fft(forcing)));

    end


    % shell 3: simulation post processing, empirical part:
    % step 1: find and report determinsitic orbit for phase alignment

    Nval = 1.2;

    for i = 1:10000

        Nval = Nval*exp(r*(1-Nval));

    end


    maxPeriod = 128;

    orbit = zeros(maxPeriod+1,1);

    orbit(1)=Nval;
    
    
    for i=1:maxPeriod

        orbit(i+1)=orbit(i)*exp(r*(1-orbit(i)));

    end


    period=[];


    for p=1:maxPeriod

        if abs(orbit(end)-orbit(end-p)) < 1e-12

            period=p;

            break

        end

    end


    if isempty(period)

        fprintf('No periodic orbit\n')

        continue

    end


    %assign empircally inferred deterministic orbit period
    period_store(ir)=period;

    fprintf('Period = %d\n',period)


    % phase alignment order
    Nstar = orbit(end-period+1:end);


    % step 2: extract empircal length scale
    spatial_means=mean(n,2); %assumption: spatial mean should be close to determinstic attractor 

    search_range = 400:(numsteps-period);


    % first phase alignment
    [~,idx]=min(abs(spatial_means(search_range)-Nstar(1)));

    t0=search_range(idx);


    % Check alignment is close
    phase_error=zeros(period,1);

    for phase=1:period

        phase_error(phase)=abs(spatial_means(t0+phase-1)-Nstar(phase));

    end


    if max(phase_error)>1e-3

        fprintf('Bad phase alignment\n')

        continue

    end


    % for the first reference phase calculate the empircal correlation
    % function
    n_dec=n(t0:period:end,:);

    maxdist=100;

    maxlag=round(maxdist/dx);

    s=(0:maxlag)*dx;


    %correlation function
    frames=size(n_dec,1);

    C=zeros(1,maxlag+1);


    for t=1:frames

        data=n_dec(t,:)-mean(n_dec(t,:)); %subtract from mean

        for h=0:maxlag

            C(h+1)=C(h+1)+mean(data.*circshift(data,-h)); %two point periodicmnkk correlation

        end

    end


    C=C/frames;


    %truncation
    tail=C(round(0.8*length(C)):end);

    threshold=3*std(tail);

    window=5;

    idx_cut=length(C);


    for i=1:length(C)-window+1

        if all(abs(C(i:i+window-1))<threshold)

            idx_cut=i;

            break

        end

    end


    %calculate length scale
    l2=sum((s(1:idx_cut).^2).*C(1:idx_cut)) ...
       /sum(C(1:idx_cut));


    lEmp(ir)=sqrt(abs(l2));


    % shell 3: simulation post processing, analytical part (prediction)
    
    % derivative definitions
    fN=@(Nin) exp(r*(1-Nin)).*(1-r*Nin);

    fi=@(a,Nin) a.*Nin.*exp(r*(1-Nin));


    %calculate relevant Floquet multiplier
    F=prod(fN(Nstar));


    A=zeros(period,1);


    for i=1:period

        if i==period

            A(i)=1;

        else

            A(i)=prod(fN(Nstar(i+1:end)));

        end

    end


    % integrated correlation structures
    P01=sum(sigma1^2*exp(-x.^2/(2*lambda1^2)))*dx;

    P02=sum(sigma2^2*exp(-phi*abs(x)))*dx;


    %Environmental sensitivities squared
    F1_sq=fi(alpha1,Nstar).^2;

    F2_sq=fi(alpha2,Nstar).^2;


    %Weights for each driver
    Num1=A.^2 .* F1_sq * P01;

    Num2=A.^2 .* F2_sq * P02;


    Num=[Num1 Num2];

    w=Num/sum(Num(:));


    % dispersal term
    mu_grid=(0:period-1)';

    coeff_matrix=2*(period-mu_grid).*w;

    dispersal_coeff=sum(coeff_matrix(:));


    dispersal_term=((2*period*F^2)/(1-F^2) ...
                   +dispersal_coeff)*alpha^2;


    % environmental term
    environmental_term=sum(w(:,1)*lambda1^2 ...
                         +w(:,2)*(2/phi)^2);


    % calculate length scale
    lu2=dispersal_term+environmental_term;

    lTheory(ir)=sqrt(abs(lu2));

end


%Plotting Lemp and Ltheory on the same map gives the validation plot.

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


figure('Position',[100,100,950,500])


plot(1:length(rvals),lTheory,'r-o', ...
    'LineWidth',2, ...
    'MarkerSize',4, ...
    'DisplayName','Analytical')

hold on


plot(1:length(rvals),lEmp,'b-x', ...
    'LineWidth',2, ...
    'MarkerSize',4, ...
    'DisplayName','Empirical')


% Draw vertical lines separating period regimes
regime_boundaries = 5:5:35;


for boundary = regime_boundaries

    xline(boundary+0.5,'--k', ...
        'Alpha',0.36, ...
        'HandleVisibility','off', ...
        'LineWidth',3);

end


xticks(3:5:38);

xticklabels({'1','2','4','8','16','32','64','128'});


xlabel('Period ($n$)')

ylabel('Length Scale ($l$)')


legend('Location','best')


ax = gca;

ax.FontSize = 16;

ax.FontName = 'Times';

ax.TickLength = [0 0];


grid on