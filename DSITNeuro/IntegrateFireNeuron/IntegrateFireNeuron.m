%% Parameters for the LIF neuron (exercise 1a)
E_L = -0.070;               % leak potential (also resting potential)
Vth = -0.050;               % threshold potential (to produce spike)
Vreset = -0.065;            % reset potential (post-spike)
Cm = 2e-9;               % total membrane capacitance
Rm = 5e+6;             % membrane resistance
G_L = 1/Rm;               % total membrane conductance (leak conductance)
tau = Cm/G_L;                % membrane time constant

%% Time vector generation (exercise 1b)
dt = 0.0001;                % time-step
t = 0:dt:2;               % vector of time-points
ton = 0.5;                 % time to begin applied current (onset)
toff = 1.5;                % time to end applied current (offset)
non = round(ton/dt);        % time-point index of current onset
noff = round(toff/dt);      % time-point index of current offset

%% Membrane potential (Vm) and applied current (Iapp) vectors (exercise 1c)
V = E_L*ones(size(t));  % initialization of membrane potential vector
Iapp = 4.001e-9; % value of applied current 

%% &  Forward Euler method to update the membrane potential(exercise 1d)
Ntrials = length(Iapp);             % number of different trials
for trial = 1:Ntrials              % loop through different trials
    I = zeros(size(t));             % vector for current at each time-point
    I(non:noff) = Iapp(trial);      % add the applied current for the trial
    V = E_L*ones(size(t));          % initialize the membrane potential vector
    spikes = zeros(size(t));        % initialize a vector to record spikes
    
    for i = 2:length(t)            % loop through all time points
        V(i) = V(i-1) + dt*(I(i) +(E_L-V(i-1))/Rm)/Cm;
        if V(i) > Vth              % if potential is above threshold
            spikes(i) = 1;          % record the spike at that time-point
            V(i) = Vreset;          % reset the potential
        end
    end                            % end the loop & go to next time-point
    
    figure(1)
    subplot(1,3,1)  % Plot current vs time
    plot(t,I*1e9,'k')
    xlabel('Time (sec)')  
    ylabel('I_{app} (nA)')
    
    subplot(1,3,2)  % Plot membrane potential vs time
    plot(t,1000*V,'k');
    xlabel('Time (sec)')                   
    ylabel('V_m (mV)')               
    set(gca,'XTick',[0 0.5 1 1.5 2])
    axis([0 2 1000*(Vreset-0.005) 1000*(Vth+0.005)])    
    
    subplot(1,3,3) % Plot Spikes vs time
    plot(t,spikes,'k')          
    xlabel('Time (sec)')      
    ylabel('Spikes')     

end  
%% EXERCISE 2 Calculate the minimum applied current needed for the neuron to produce spikes. 

I_th = G_L * (Vth - E_L); %Calculate the threshold current

display(I_th); 
 
% print('The threshold current is',int2str(I_th),'nA.')

%% EXERCISE 2a+b Validate the calculated threshold current. 

% Now simulate trials, each with a different applied current
Iapp = [I_th-(I_th/100) I_th I_th+(I_th/100)];   % values of applied current steps
Ntrials = length(Iapp);             % number of different trials
for trial = 1:Ntrials               % loop through different trials
    I = zeros(size(t));             % vector for current at each time-point
    I(non:noff) = Iapp(trial);      % add the applied current for the trial
    V = E_L*ones(size(t));          % initialize the membrane potential vector
    spikes = zeros(size(t));        % initialize a vector to record spikes
    
    for i = 2:length(t)             % loop through all time points
        V(i) = V(i-1) + dt*(I(i) +((E_L-V(i-1))/(Rm)))/Cm;
        if V(i) > Vth               % if potential is above threshold
            spikes(i) = 1;          % record the spike at that time-point
            V(i) = Vreset;          % reset the potential
        end 
    end                             % end the loop & go to next time-point
    
    % Now generate a series of subplots to position the subfigures
    figure(2)
    subplot('position',[0.31*trial-0.18 0.74 0.22 0.22])
    
    plot(t,I*1e9,'k');              % Plot current (in nA) against time
    if ( trial == 1 )
        ylabel('I_{app} (nA)')
    end
    % The following line adds a title, which uses the command "strcat" to
    % combine the value of current in that trial with "nA" to make a
    % complete title.
    % The command "num2str" converts the number representing the applied
    % current into a string, which is something that can be printed.
    title(strcat(num2str(1e9*Iapp(trial)),'nA'))
    
    set(gca,'XTick',[0 0.5 1 1.5 2])
    
    axis([0 2 0 5])         % Sets x-axis then y-axis ranges
    
    % Next subfigure
    subplot('position',[0.31*trial-0.18 0.43 0.22 0.22])
    plot(t,1000*V,'k');                      % Plot membrane potential vs time
    if ( trial == 1 )
        ylabel('V_m (mV)')               % Label y-axis
    end
    set(gca,'XTick',[0 0.5 1 1.5 2])
    axis([0 2 1000*(Vreset-0.005) 1000*(Vth+0.005)])    % Set the axes
    
    subplot('position',[0.31*trial-0.18 0.12 0.22 0.22])
    plot(t,spikes,'k')          % Plots a line from 0 to 1 at each spike
    xlabel('Time (sec)')        % Label the x-axis
    if ( trial == 1)
        ylabel('Spikes')        % Label the y-axis
    end
    set(gca,'XTick',[0 0.5 1 1.5 2])
    set(gca,'YTick',[0 1])
    axis([0 2 0 1])           % Set the ranges of x- and y-axes
    
end;                            % Loop to the next trial with new current

annotation('textbox',[0.00 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',12,'FontWeight','Bold','String','A')
annotation('textbox',[0.37 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',12,'FontWeight','Bold','String','B')
annotation('textbox',[0.68 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',12,'FontWeight','Bold','String','C')


%% EXERCISE 3
% Calculate the maximum steady state membrane potential (Vmss) and the Imax
% so that firing rate=100Hz
tau = Rm * Cm;
Vmss_max = (Vreset- Vth * exp(0.01/tau))/(1- exp(0.01/tau));
I_max = (Vmss_max - E_L)/Rm;
I_range = I_th:((I_max-I_th)/100):I_max;

% Generate values for I so that 0<firing rate<100.
spikes = zeros(Ntrials,length(t));        % initialize a vector to record spikes
%Iapp_fr = sort(randsample(I_range, 10)); % for extremely rnadom values
Iapp_fr = I_th:((I_max-I_th)/9):I_max;
Ntrials = length(Iapp_fr);          % number of different trials
for trial = 1:Ntrials               % loop through different trials
    I = zeros(size(t));             % vector for current at each time-point
    I(non:noff) = Iapp_fr(trial);   % add the applied current for the trial
    V = E_L*ones(size(t));          % initialize the membrane potential vector
    
    for i = 2:length(t)            % loop through all time points
        V(i) = V(i-1) + dt*(I(i) +((E_L-V(i-1))/(Rm)))/Cm;
        if V(i) > Vth              % if potential is above threshold
            spikes(trial,i) = 1;    % record the spike at that time-point
            V(i) = Vreset;          % reset the potential
        end
    end                            % end the loop & go to next time-point
end  
% Calculate Vmss and firing rates
Vmss_fr= E_L + Iapp_fr * Rm;
ISI = tau * log((Vmss_fr-Vreset)./(Vmss_fr-Vth));
fir_rates_eq = 1 ./ ISI;

figure(3)
plot(Iapp_fr, sort(fir_rates_eq), '-o')
hold on;
plot(Iapp_fr, sum(spikes, 2), '-x')
hold off;
title('Firing-rate Curve')
xlabel('I (A)')
ylabel('Frequency (Hz)')
legend('ISI', 'experimental', 'Location', 'NorthWest')
ylim([0,100])
xlim([I_th, I_max])
%% EXERCISE 4 a Add noise term to the simulation of LIF model

% Adding the noise term to the total membrane potential at each time step:

% Generate noise vector
sigma_I = [0 1e-12 1e-5 1e-4 1e-3 0.01];
%sigma_I = [1e-12 1e-10 1e-8 1e-6 1e-5 1e-4 1e-3 1e-2] ;

% Generate the min and max current that will be applied to the model
Imin = I_th - I_th/100;
Imax = I_th + I_th/100;
I_range_min = Imin+(I_th-Imin)/10:(I_th-Imin)/10:I_th-(I_th-Imin)/10;
I_prev = sort(randsample(I_range_min, 3));
I_range_max = I_th + (Imax-I_th)/10:(Imax-I_th)/10:Imax-(Imax-I_th)/10;
I_after = sort(randsample(I_range_max, 3));

Iapp_fr = [Imin I_prev I_th I_after Imax]; % Creating 10 random values between Imin and Imax including the threshold current
Ntrials = length(Iapp_fr); % number of different trials
firing_rates = [];
count=1;
for sigma = sigma_I
    spikes = zeros(Ntrials,length(t));        % initialize a vector to record spikes
    for trial = 1:Ntrials               % loop through different trials
        I = zeros(size(t));             % vector for current at each time-point
        I(non:noff) = Iapp_fr(trial);   % add the applied current for the trial
        V = E_L*ones(size(t));          % initialize the membrane potential vector

        for i = 2:length(t)            % loop through all time points
            V(i) = V(i-1) + dt*(I(i) +((E_L-V(i-1))/(Rm)))/Cm + randn(1) * sigma * sqrt(dt);
            if V(i) > Vth              % if potential is above threshold
                spikes(trial,i) = 1;    % record the spike at that time-point
                V(i) = Vreset;          % reset the potential
            end
        end                            % end the loop & go to next time-point  
    end
    freq = sum(spikes, 2);
    firing_rates= [firing_rates,freq];
    count = count+1;  
end

figure(4)
legendCell = strcat('sigma=',string(num2cell(sigma_I)));
for i=1:length(sigma_I)
    plot(Iapp_fr, firing_rates(:,i), '-*')
    hold on;
    title('Firing-rate Curve')
    xlabel('I (A)')
    ylabel('Frequency (Hz)')
    legend(legendCell, 'Location', 'NorthWest')
    xlim([Imin, Imax])
end
hold off;

%% EXERCISE 4 b Add noise term to the simulation of LIF model and alter the timestep

%altering the timstep, a factor of ten smaller
dt_altered = dt/10;
t = 0:dt_altered:2;

% Adding the noise term to the total membrane potential at each time step:
% Generate noise vector
sigma_I = [0 1e-12 1e-5 1e-4 1e-3 0.01];
%sigma_I = [1e-12 1e-10 1e-8 1e-6 1e-5 1e-4 1e-3 1e-2] ;

% Generate the min and max current that will be applied to the model
Imin = I_th - I_th/100;
Imax = I_th + I_th/100;
I_range_min = Imin+(I_th-Imin)/10:(I_th-Imin)/10:I_th-(I_th-Imin)/10;
I_prev = sort(randsample(I_range_min, 3));
I_range_max = I_th + (Imax-I_th)/10:(Imax-I_th)/10:Imax-(Imax-I_th)/10;
I_after = sort(randsample(I_range_max, 3));

Iapp_fr = [Imin I_prev I_th I_after Imax]; % Creating 10 random values between Imin and Imax including the threshold current
Ntrials = length(Iapp_fr); % number of different trials
firing_rates = [];
count=1;
for sigma = sigma_I
    spikes = zeros(Ntrials,length(t));        % initialize a vector to record spikes
    for trial = 1:Ntrials               % loop through different trials
        I = zeros(size(t));             % vector for current at each time-point
        I(non:noff) = Iapp_fr(trial);   % add the applied current for the trial
        V = E_L*ones(size(t));          % initialize the membrane potential vector

        for i = 2:length(t)            % loop through all time points
            V(i) = V(i-1) + dt_altered*(I(i) +((E_L-V(i-1))/(Rm)))/Cm + randn(1) * sigma * sqrt(dt_altered);
            if V(i) > Vth              % if potential is above threshold
                spikes(trial,i) = 1;    % record the spike at that time-point
                V(i) = Vreset;          % reset the potential
            end
        end                            % end the loop & go to next time-point  
    end
    freq = sum(spikes, 2);
    firing_rates= [firing_rates,freq];
    count = count+1;  
end

figure(5)
legendCell = strcat('sigma=',string(num2cell(sigma_I)));
for i=1:length(sigma_I)
    plot(Iapp_fr, firing_rates(:,i), '-*')
    hold on;
    title('Firing-rate Curve')
    xlabel('I (A)')
    ylabel('Frequency (Hz)')
    legend(legendCell, 'Location', 'NorthWest')
    xlim([Imin, Imax])
end
hold off;