%% Exercise 1
clear
dt = 1e-6;          % time-step for integration (sec)
tmax=0.5;           % maximum time of simulation (sec)

istart = 0.05;      % time applied current starts (sec)
ilength=0.35;        % length of applied current pulse (sec)
Ibase = 0;          % Baseline current before pulse

%% Neuron parameters
V_L = -0.060;        % leak reversal potential (V)
E_Na = 0.045;           % reversal for sodium channels (V)
E_K = -0.082;           % reversal for potassium channels (V)

G_L = 30e-9;            % leak conductance (S)
G_Na = 12e-6;           % sodium conductance (S)
G_K = 3.6e-6;           % potassium conductance (S)

Cm = 100e-12;           % total membrane capacitance (F)

t=0:dt:tmax;            % vector of time points

Ievec = [0]*1e-9;      % series of applied currents (A)
Ntrials = length(Ievec);    % Number of trials to loop through

%% Now loop through trials, each one with a different current step
for trial = 1:Ntrials
    
    Ie= Ievec(trial);           % New applied current each trial
    
    V=zeros(size(t));           % membrane potential vector
    
    %% current clamp initialization
    Iapp=Ibase*ones(size(t));   % Applied current has a baseline
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        Iapp(i) = Ie;
    end
    
    V(1) = V_L;             % set the inititial value of voltage
    
    n=zeros(size(t));       % n: potassium activation gating variable
    n(1) = 0.35;            % start off near steady state when V is V_L
    m=zeros(size(t));       % m: sodium activation gating variable
    m(1) = 0.05;            % start off near steady state when V is V_L
    h=zeros(size(t));       % h: sodim inactivation gating variable
    h(1) = 0.75;            % start off near steady state when V is V_L
    
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 1:length(t)-1; % now see how things change through time
        
        I_L(i) = G_L*(V_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
        %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(1) 
    if ( trial == 1 )
        clf
    end
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('Iapp =  ',num2str(Ie*1e9),' nA')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})

end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')


%% Exercise 2

clear
dt = 1e-6;          % time-step for integration (sec)
tmax=0.3;           % maximum time of simulation (sec)

istart = 0.1;      % time applied current starts (sec)
ilength=0.1;        % length of applied current pulse (sec)
Ibase = 0;          % Baseline current before pulse

%% Neuron parameters
V_L = -0.060;        % leak reversal potential (V)
E_Na = 0.045;           % reversal for sodium channels (V)
E_K = -0.082;           % reversal for potassium channels (V)

G_L = 30e-9;            % leak conductance (S)
G_Na = 12e-6;           % sodium conductance (S)
G_K = 3.6e-6;           % potassium conductance (S)

Cm = 100e-12;           % total membrane capacitance (F)

t=0:dt:tmax;            % vector of time points

Ievec = [0.05 0.10 0.15 0.20 0.22]*1e-9;      % series of applied currents (A)
Ntrials = length(Ievec);    % Number of trials to loop through

%% Now loop through trials, each one with a different current step
for trial = 1:Ntrials
    
    Ie= Ievec(trial);           % New applied current each trial
    
    V=zeros(size(t));           % membrane potential vector
    
    %% current clamp initialization
    Iapp=Ibase*ones(size(t));   % Applied current has a baseline
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        Iapp(i) = Ie;
    end
    
    V(1) = V_L;             % set the inititial value of voltage
    
    n=zeros(size(t));       % n: potassium activation gating variable
    n(1) = 0.35;            % start off near steady state when V is V_L
    m=zeros(size(t));       % m: sodium activation gating variable
    m(1) = 0.05;            % start off near steady state when V is V_L
    h=zeros(size(t));       % h: sodim inactivation gating variable
    h(1) = 0.75;            % start off near steady state when V is V_L
    
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 1:length(t)-1; % now see how things change through time
                                                            
        I_L(i) = G_L*(V_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
        %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(2) 
    if ( trial == 1 )
        clf
    end
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('Iapp =  ',num2str(Ie*1e9),' nA')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})
end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.63 0.05 0.05],'String','C','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.47 0.05 0.05],'String','D','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.31 0.05 0.05],'String','E','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.15 0.05 0.05],'String','F','LineStyle','none','FontSize',16,'FontWeight','Bold')

%% Exercise 3
clear
dt = 2e-8;          % time-step for integration (sec)
tmax=0.5;           % maximum time of simulation (sec)
t=0:dt:tmax;

istart = 0.1;      % time applied current starts (sec)
ilength = 0.005;        % length of applied current pulse (sec)
step = [0.005 0.010 0.015 0.020 0.025]; 
Ibase = 0;          % Baseline current before pulse

%% Neuron parameters
V_L = -0.060;        % leak reversal potential (V)
E_Na = 0.045;           % reversal for sodium channels (V)
E_K = -0.082;           % reversal for potassium channels (V)

G_L = 30e-9;            % leak conductance (S)
G_Na = 12e-6;           % sodium conductance (S)
G_K = 3.6e-6;           % potassium conductance (S)

Cm = 100e-12;           % total membrane capacitance (F)

t=0:dt:tmax;            % vector of time points

Ievec = 0.22*1e-9;      % series of applied currents (A)
Ntrials = length(step);    % Number of trials to loop through


%% Now loop through trials, each one with a different current step
for trial = 1:Ntrials
    

    Ie= Ievec;           % New applied current each trial
    
    V=zeros(size(t));           % membrane potential vector
    
    delay = step(trial);
    
    %% current clamp initialization
    Iapp=Ibase*ones(size(t));   % Applied current has a baseline
    
    for s=1:10
        if s ==1
            istart(s) = 0.1;
        else 
            istart(s) = istart(s-1) + delay
        end
        for i=round(istart(s)/dt)+1:round((istart(s)+ilength)/dt) % make non-zero for duration of current pulse
            Iapp(i) = Ie;
        end
    end

    V(1) = V_L;             % set the inititial value of voltage
    
    n=zeros(size(t));       % n: potassium activation gating variable
    n(1) = 0.35;            % start off near steady state when V is V_L
    m=zeros(size(t));       % m: sodium activation gating variable
    m(1) = 0.05;            % start off near steady state when V is V_L
    h=zeros(size(t));       % h: sodim inactivation gating variable
    h(1) = 0.75;            % start off near steady state when V is V_L
    
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 1:length(t)-1; % now see how things change through time
        
        I_L(i) = G_L*(V_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
        %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(3) 
    if ( trial == 1 )
        clf
    end
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('delay =  ',num2str(delay),' s')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})
end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.63 0.05 0.05],'String','C','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.47 0.05 0.05],'String','D','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.31 0.05 0.05],'String','E','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.15 0.05 0.05],'String','F','LineStyle','none','FontSize',16,'FontWeight','Bold')

%% Exercise 4
clear
dt = 1e-6;          % time-step for integration (sec)
tmax=0.5;           % maximum time of simulation (sec)
t=0:dt:tmax;

istart = 0.1;      % time applied current starts (sec)
ilength = 0.005;        % length of applied current pulse (sec)
step = [0.020]; 
Ibase = 0.6*1e-9;          % Baseline current before pulse

%% Neuron parameters
V_L = -0.065;        % leak reversal potential (V)
E_Na = 0.045;           % reversal for sodium channels (V)
E_K = -0.082;           % reversal for potassium channels (V)

G_L = 30e-9;            % leak conductance (S)
G_Na = 12e-6;           % sodium conductance (S)
G_K = 3.6e-6;           % potassium conductance (S)

Cm = 100e-12;           % total membrane capacitance (F)

t=0:dt:tmax;            % vector of time points

Ievec = 0;      % series of applied currents (A)
Ntrials = length(step);    % Number of trials to loop through


%% Now loop through trials, each one with a different current step
for trial = 1:Ntrials
    

    Ie= Ievec;           % New applied current each trial
    
    V=zeros(size(t));           % membrane potential vector
    
    delay = step(trial);
    
    %% current clamp initialization
    Iapp=Ibase*ones(size(t));   % Applied current has a baseline
    
    for s=1:10
        if s ==1
            istart(s) = 0.1;
        else 
            istart(s) = istart(s-1) + delay
        end
        for i=round(istart(s)/dt)+1:round((istart(s)+ilength)/dt) % make non-zero for duration of current pulse
            Iapp(i) = Ie;
        end
    end

    V(1) = V_L;             % set the inititial value of voltage
    
    n=zeros(size(t));       % n: potassium activation gating variable
    n(1) = 0.35;            % start off near steady state when V is V_L
    m=zeros(size(t));       % m: sodium activation gating variable
    m(1) = 0.05;            % start off near steady state when V is V_L
    h=zeros(size(t));       % h: sodim inactivation gating variable
    h(1) = 0.5;            % start off near steady state when V is V_L
    
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 1:length(t)-1; % now see how things change through time
        
        I_L(i) = G_L*(V_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
        %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(4) 
    if ( trial == 1 )
        clf
    end
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('delay =  ',num2str(delay),' s')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})
end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')

%% Exercise 5

clear
dt = 1e-6;          % time-step for integration (sec)
tmax = 0.35;           % maximum time of simulation (sec)

istart = 0.1;      % time applied current starts (sec)
ilength = 0.005;        % length of applied current pulse (sec)
Ibase = 0.65*1e-9;          % Baseline current before pulse

%% Neuron parameters
V_L = -0.060;        % leak reversal potential (V)
E_Na = 0.045;           % reversal for sodium channels (V)
E_K = -0.082;           % reversal for potassium channels (V)

G_L = 30e-9;            % leak conductance (S)
G_Na = 12e-6;           % sodium conductance (S)
G_K = 3.6e-6;           % potassium conductance (S)

Cm = 100e-12;           % total membrane capacitance (F)

t=0:dt:tmax;            % vector of time points

Ievec = 1e-9;      % series of applied currents (A)
Ntrials = length(Ievec);    % Number of trials to loop through

%% Now loop through trials, each one with a different current step
for trial = 1:Ntrials
    
    Ie= Ievec(trial);           % New applied current each trial
    
    V=zeros(size(t));           % membrane potential vector
    
    %% current clamp initialization
    Iapp=Ibase*ones(size(t));   % Applied current has a baseline
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        Iapp(i) = Ie;
    end
    
    V(1) = V_L;             % set the inititial value of voltage
    
    n=zeros(size(t));       % n: potassium activation gating variable
    n(1) = 0.35;            % start off near steady state when V is V_L
    m=zeros(size(t));       % m: sodium activation gating variable
    m(1) = 0.05;            % start off near steady state when V is V_L
    h=zeros(size(t));       % h: sodim inactivation gating variable
    h(1) = 0.75;            % start off near steady state when V is V_L
    
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 1:length(t)-1; % now see how things change through time
                                                            
        I_L(i) = G_L*(V_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
        %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(5) 
    if ( trial == 1 )
        clf
    end
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('Iapp =  ',num2str(Ie*1e9),' nA')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})
end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')

%% Exercise 6

clear
dt = 1e-6;          % time-step for integration (sec)
tmax = 0.35;           % maximum time of simulation (sec)

istart = 0.1;      % time applied current starts (sec)
ilength = 0.005;        % length of applied current pulse (sec)
Ibase = 0.7*1e-9;          % Baseline current before pulse

%% Neuron parameters
V_L = -0.065;        % leak reversal potential (V)
E_Na = 0.045;           % reversal for sodium channels (V)
E_K = -0.082;           % reversal for potassium channels (V)

G_L = 30e-9;            % leak conductance (S)
G_Na = 12e-6;           % sodium conductance (S)
G_K = 3.6e-6;           % potassium conductance (S)

Cm = 100e-12;           % total membrane capacitance (F)

t=0:dt:tmax;            % vector of time points

Ievec = 1e-9;      % series of applied currents (A)
Ntrials = length(Ievec);    % Number of trials to loop through

%% Now loop through trials, each one with a different current step
for trial = 1:Ntrials
    
    Ie= Ievec(trial);           % New applied current each trial
    
    V=zeros(size(t));           % membrane potential vector
    
    %% current clamp initialization
    Iapp=Ibase*ones(size(t));   % Applied current has a baseline
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        Iapp(i) = Ie;
    end
    
    V(1) = V_L;             % set the inititial value of voltage
    
    n=zeros(size(t));       % n: potassium activation gating variable
    n(1) = 0;            % start off near steady state when V is V_L
    m=zeros(size(t));       % m: sodium activation gating variable
    m(1) = 0;            % start off near steady state when V is V_L
    h=zeros(size(t));       % h: sodim inactivation gating variable
    h(1) = 0;            % start off near steady state when V is V_L
    
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 1:length(t)-1; % now see how things change through time
                                                            
        I_L(i) = G_L*(V_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
        %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(6) 
    if ( trial == 1 )
        clf
    end
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('Iapp =  ',num2str(Ie*1e9),' nA')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})
end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')