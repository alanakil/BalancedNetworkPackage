%% Simulation of the balanced network using AdEx model %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Record from a number of neurons.
% Indices of neurons to record currents, voltages
nrecord0=10; % Number to record from each population
Irecord=[randperm(Ne,nrecord0) randperm(Ni,nrecord0)+Ne];
numrecord=numel(Irecord); % total number to record

% Number of time bins to average over when recording
nBinsRecord=1;%20;
dtRecord=nBinsRecord*dt;
timeRecord=dtRecord:dtRecord:T;
Ntrec=numel(timeRecord);

% Integer division function
IntDivide=@(n,k)(floor((n-1)/k)+1);

% Random initial voltages.
V0=rand(N,1)*(VT-Vre)+Vre;

% Preallocate memory.
V=V0; % Vector of voltages of neurons. 
w = zeros(N,1); % Vecotr of adaptation currents.
Ie=zeros(N,1);
Ii=zeros(N,1);
Ix=zeros(N,1);
x=zeros(N,1);
IeRec=zeros(numrecord,Ntrec);
IiRec=zeros(numrecord,Ntrec);
IxRec=zeros(numrecord,Ntrec);
VRec=zeros(numrecord,Ntrec);
wRec=zeros(numrecord,Ntrec);
iFspike=1;
s=zeros(2,maxns);
nspike=0;
TooManySpikes=0;
tic

% Adaptive current parameters.
tauw = 150;
a = 0;
b = 2;
Ispike = [];

%% Start of actual simulation.
for i=1:numel(time)

    % Propogate ffwd spikes
    while(sx(1,iFspike)<=time(i) && iFspike<nspikeX)
        jpre=sx(2,iFspike);
        Ix=Ix+Jx(:,jpre)/taux;
        iFspike=iFspike+1;
    end
    
    % Euler update to V (AdEx)
    V = V+(dt/Cm)*(Istim(i)*Jstim+Ie+Ii+Ix+gL*(EL-V)+gL*DeltaT*exp((V-VT)/DeltaT));
    neuron_spiked = zeros(N,1);
    neuron_spiked(Ispike) = 1;
    w = w + dt/tauw * (a * (EL-V) - w + b * tauw * neuron_spiked);
    
    % Find which neurons spiked
    Ispike=find(V>=Vth);    
    
    % If there are spikes
    if(~isempty(Ispike))

        % Store spike times and neuron indices
        if(nspike+numel(Ispike)<=maxns)
            s(1,nspike+1:nspike+numel(Ispike))=time(i);
            s(2,nspike+1:nspike+numel(Ispike))=Ispike;
        else
            TooManySpikes=1;
            break;
        end
         
        % Update synaptic currents
        Ie=Ie+sum(J(:,Ispike(Ispike<=Ne)),2)/taue;    
        Ii=Ii+sum(J(:,Ispike(Ispike>Ne)),2)/taui;            
        
        % Update cumulative number of spikes
        nspike=nspike+numel(Ispike);
    end            
    
    % Euler update to synaptic currents
    Ie=Ie-dt*Ie/taue;
    Ii=Ii-dt*Ii/taui;
    Ix=Ix-dt*Ix/taux;
    
    % This makes plots of V(t) look better.
    % All action potentials reach Vth exactly. 
    % This has no real effect on the network sims
    V(Ispike)=Vth;
    
    % Store recorded variables
    ii=IntDivide(i,nBinsRecord); 
    IeRec(:,ii)=IeRec(:,ii)+Ie(Irecord);
    IiRec(:,ii)=IiRec(:,ii)+Ii(Irecord);
    IxRec(:,ii)=IxRec(:,ii)+Ix(Irecord);
    VRec(:,ii)=VRec(:,ii)+V(Irecord);
    wRec(:,ii)=wRec(:,ii)+w(Irecord);

    % Reset mem pot.
    V(Ispike)=Vre;
    
end
% Normalize recorded variables by # bins
IeRec=IeRec/nBinsRecord; 
IiRec=IiRec/nBinsRecord;
IxRec=IxRec/nBinsRecord;
VRec=VRec/nBinsRecord;
s=s(:,1:nspike); % Get rid of padding in s


tSim=toc; % Determine time of simulation.
disp(sprintf('\nTime for simulation: %.2f min',tSim/60))

figure; hold on
plot( time, wRec(1,:))
title('Adaptation current of one E cell')
xlabel('time (ms)')
ylabel('w (mV)')


