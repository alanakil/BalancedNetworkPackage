%% Simulation of a balanced network - Main script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

% Number of neurons in each population
Ne=4000;
Ni=1000;
N=Ne+Ni;
% Number of neurons in ffwd layer
Nx=1000;
% Fractions of E, I and X neurons.
qe=Ne/N;
qi=Ni/N;
qf=Nx/N;

% Time (in ms) for sim
T=10000;
% Time discretization
dt=.1;
% Number of time bins
Nt=round(T/dt);
time=dt:dt:T;

% Change connectivity.m to explore the effects of the weight matrix in the
% dynamics.
Connectivity

% Change the external layer at ExternalLayer.m to explore its effect in
% dynamics.
ExternalLayer

%% Extra stimulus: Istim is a time-dependent stimulus
% it is delivered to all neurons with weights given by Jstim.
% Specifically, the stimulus to neuron j at time index i is:
% Istim(i)*Jstim(j).
Istim=zeros(size(time)); 
Istim(time>T/2)=1; 
jestim=0; 
jistim=0;
Jstim=sqrt(N)*[jestim*ones(Ne,1); jistim*ones(Ni,1)]; 

% Synaptic timescales.
taux=8;
taue=6;
taui=4;

% Neuron parameters.
Cm=1;
gL=1/15;
EL=-72;
Vth=0;
Vre=-75;
DeltaT=2;
VT=-55;

% Maximum number of spikes for all neurons
% in simulation. Make it 50Hz across all neurons
% If there are more spikes, the simulation will
% terminate.
maxns=ceil(.05*N*T);

%% Simulation.m runs the actual simulation of the balanced network.
% Simulation runs the simulation with just EIF cells.
% Simulation_AdEX uses the AdEx model instead.
% Comment out whichever model the user doesn't want to use.

%Simulation
Simulation_AdEx

%% Analysis of the network (rasterplots, rates, and correlations). 

% Make a raster plot of first 500 neurons 
% s(1,:) are the spike times
% s(2,:) are the associated neuron indices
figure
plot(s(1,s(2,:)<500),s(2,s(2,:)<500),'k.')
xlim([ T-500 T]) % Only look at the last 500 ms of the simulation.
xlabel('time (ms)')
ylabel('Neuron index')

%% Plot the E and I rates over time.
% Mean rate of each neuron (excluding burn-in period)
Tburn=500; % burn-in period.
reSim=hist(s(2,s(1,:)>Tburn & s(2,:)<=Ne),1:Ne)/(T-Tburn); % Count E spikes
% between Tburn and T, then divide by the time difference.
riSim=hist(s(2,s(1,:)>Tburn & s(2,:)>Ne)-Ne,1:Ni)/(T-Tburn); % Same but for
% I spikes.

% Mean rate over E and I pops
reMean=mean(reSim); % Mean E rate of the simulation excluding burn-in period.
riMean=mean(riSim);% Mean I rate of the simulation excluding burn-in period.
disp(sprintf('\nMean E and I rates from sims: %.2fHz %.2fHz',1000*reMean,1000*riMean))

% Time-dependent mean rates
dtRate=100;
eRateT=hist(s(1,s(2,:)<=Ne),dtRate:dtRate:T)/(dtRate*Ne);
iRateT=hist(s(1,s(2,:)>Ne),dtRate:dtRate:T)/(dtRate*Ni);

% Plot time-dependent rates
figure
plot((dtRate:dtRate:T)/1000,1000*eRateT)
hold on
plot((dtRate:dtRate:T)/1000,1000*iRateT)
legend('r_e','r_i')
ylabel('rate (Hz)')
xlabel('time (s)')

%% Compute spike count covariances and correlations.

% Compute spike count covariances over windows of size
% winsize starting at time T1 and ending at time T2.
winsize=250; 
T1=250; % Burn-in period of 250 ms
T2=T;   % Compute covariances until end of simulation
C=SpikeCountCov(s,N,T1,T2,winsize);

% Get mean spike count covariances over each sub-pop
[II,JJ]=meshgrid(1:N,1:N);
mCee=mean(C(II<=Ne & JJ<=II));
mCei=mean(C(II<=Ne & JJ>Ne));
mCii=mean(C(II>Ne & JJ>II));

% Mean-field spike count cov matrix
% Compare this to the theoretical prediction
mC=[mCee mCei; mCei mCii]

% Compute spike count correlations from covariances.
% This takes a while, so make it optional
ComputeSpikeCountCorrs=0;
if(ComputeSpikeCountCorrs)
    
    % Get correlation matrix from cov matrix
    tic
    R=corrcov(C);
    toc
    
    mRee=mean(R(II<=Ne & JJ<=II & isfinite(R)));
    mRei=mean(R(II<=Ne & JJ>Ne & isfinite(R)));
    mRii=mean(R(II>Ne & JJ>II & isfinite(R)));

    % Mean-field spike count correlation matrix
    mR=[mRee mRei; mRei mRii]


end




