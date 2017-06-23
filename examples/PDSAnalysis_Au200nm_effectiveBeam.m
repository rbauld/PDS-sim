%% Load data

Data_f(1).Filename = '.\10.00Hz.txt' ;
Data_f(2).Filename = '.\20.00Hz.txt' ;
Data_f(3).Filename = '.\30.00Hz.txt' ;

Data_f(1).freq = 10 ;
Data_f(2).freq = 20 ;
Data_f(3).freq = 30 ;

for ii=1:3
Data_f(ii).raw = importdata(Data_f(ii).Filename,'\t',1) ;    
end

DataRange = 1:length(Data_f(1).raw.data(:,2)) ;
FreqIndex = [1 2 3];

for ii=1:3
    % Scale mm to m
    Data_f(ii).data(:,1)=Data_f(ii).raw.data(DataRange,2).*10^-3 ;
    
    % Freq in column 2
    Data_f(ii).data(:,2)=ones(length(Data_f(ii).data(:,1)),1).*Data_f(ii).freq ;
    
    % Amplitude in column 3
    Data_f(ii).data(:,3)=sqrt(Data_f(ii).raw.data(DataRange,3).^2 +Data_f(ii).raw.data(DataRange,4).^2) ;
    
    % Complex representation in column 4
    Data_f(ii).data(:,4)=Data_f(ii).raw.data(DataRange,3)-1i.*Data_f(ii).raw.data(DataRange,4) ;
    
    % Normalized amplitude in column 5
    Data_f(ii).data(:,5)=Data_f(ii).data(:,3)./max(Data_f(ii).data(:,3)) ;
    
    % Phase in column 6
    Data_f(ii).data(:,6)=angle(Data_f(ii).data(:,4)) ;
end

PDS_sim = PDSExper()

%% Set experimental parameters

% Region 0, Fluid

PDS_sim.Region0_k0   = 0.057 ;
PDS_sim.Region0_rho0 = 1680 ;
PDS_sim.Region0_c0   = 1100 ;
PDS_sim.Region0_R1   = 0 ; % Reflection coeff between 0 & 1

% Region 1, Thin Film, Au
PDS_sim.Region1_k1   = 314 ;
PDS_sim.Region1_rho1 = 19300 ;
PDS_sim.Region1_c1   = 129 ;
PDS_sim.Region1_a1   = 5*10^7 ;
PDS_sim.Region1_L1   = 200*10^-9 ; % Thin film thickness
PDS_sim.Region1_R2   = 0 ; % Reflection coeff between 1 & 2

% Region 2, Substrate (Glass)

PDS_sim.Region2_k2   = 1.2 ;
PDS_sim.Region2_rho2 = 2400 ;
PDS_sim.Region2_c2   = 670 ;
PDS_sim.Region2_a2   = 0 ;
PDS_sim.Region2_L2   = 1e-3 ;
PDS_sim.Region2_Rth  = 0 ; % Interfacial Thermal resistance

% Region 3, Fluid

PDS_sim.Region3_k3 = 0.057 ;
PDS_sim.Region3_rho3 = 1680 ;
PDS_sim.Region3_c3 = 1100 ;
PDS_sim.n = 1.251 ;
PDS_sim.dNdT = -3.4e-4 ;

% Define other parameters

PDS_sim.PumpRadius = 50e-6 ;
PDS_sim.yshift = 0 ;
PDS_sim.ProbeRadius =100e-6 ;
PDS_sim.PhaseShift = 0 ;

PDS_sim.Power = 1
PDS_sim.xshift = 6.06e-3
PDS_sim.z = -100e-6 ;
PDS_sim.f = 10 ;



%% Plot frontside amplitude    

PDS_sim.xshift = 7.794e-3 ;
PDS_sim.Region1_k1 = 300 ;
PDS_sim.Region2_k2 = 1.0 ;

PDS_sim.P = 0.3 ;
ProbeSigma = 43e-6 ;
PDS_sim.ProbeRadius = 2*sqrt(2*log(2))*ProbeSigma ;
PDS_sim.PumpRadius = 40e-06 ;
PDS_sim.z = -1.08e-04 ;
PDS_sim.Region2_Rth = 2e-8 ; 

PDS_sim.plot_amp(Data_f)




%% Simulated annealing, amplitude
tic
clear Annealing
Annealing.Resid = [] ;
Annealing.Params = [] ;
NumberOfCycles = 200 ;
Temperature= [1 0.5 0.1].*1 ;


params(1) = 300 ;  % Thin film k
params(2) = 0.3 ; % Power
params(3) = -1.1006e-04 ; % z
params(4) = 1.0 ; % substrate k
params(5) = 2e-8 ; %Interfacial thermal resistance
params(6) = 2*sqrt(2*log(2))*ProbeSigma ;

%%
for jj=1:3
    %Define parameter perturbations
    param_perturb(1) = 0.*Temperature(jj) ;  % Thin film k
    param_perturb(2) = 1e-2.*Temperature(jj) ; % Power
    param_perturb(3) = 10e-6.*Temperature(jj) ; % z front
    param_perturb(4) = 0.0.*Temperature(jj) ; % substrate k
    param_perturb(5) = 0e-8.*Temperature(jj) ; % Rth
    param_perturb(6) = 10e-6.*Temperature(jj) ; % Probe beam
    
    for kk=1:NumberOfCycles            
        params_guess = params+param_perturb.*(rand(size(params)).*2-1) ;
        params_guess
        
        Residual = PDS_sim.resid_calc(Data_f, params_guess) ;

        Annealing.Resid = [Annealing.Resid ,Residual] ;
        Annealing.Params = [Annealing.Params; params_guess] ;

        % Choose best parameters at each iteration
        [value, Index] = sort(Annealing.Resid(:)) ;
        params = Annealing.Params(Index(1),:) ;
        Resid_final = Annealing.Resid(Index(1))
    end


end
toc

h = msgbox('Operation Completed','Success');

%% Load params and make test plot

%params = my_fit./my_scale ;

PDS_sim.Region1_k1 = params(1) ;
PDS_sim.P          = params(2) ;
PDS_sim.z          = params(3) ;
PDS_sim.Region2_k2 = params(4) ;
PDS_sim.Region2_Rth = params(5) ;
PDS_sim.ProbeRadius = params(6) ;

PDS_sim.plot_amp(Data_f)

%set(gca,'YScale','log')


