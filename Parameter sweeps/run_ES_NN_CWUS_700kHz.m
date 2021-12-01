USiRange = (10:10:200);
for iii = 1:20
clearvars -except USiRange iii 
clear functions %#ok<CLFUNC>
%clear all; clear functions; %#ok<CLFUNC,CLALL>
coder.extrinsic('nakeinterp1');
rng(3); 
% Code implements the cortex-basal ganglia-thalamus model used in  
% Tarnaud et al. (2021) [1]  

% The neuronal network is based on Kumaravelu et al. (2016) [2] and on the cortical network 
% from Plaksin et al. (2016) [3], Hayut et al. (2011) [4] and 
% Vierling-Claassen et al. (2010) [5] 

% The result is a completely biophysical HH-based network of the
% basal-ganglia-thalamus-cortex loop with ultrasonic coupling
% 
% Adaptations: 
% 1) The STN-Ca dynamics in Kumaravelu et al. (2016) [2] is modified according
% to Otsuka et al. (2004) [6].  
% 2) Also a connection from Thalamus to FS neurons is included (not present
% in Kumaravelu et al., but present in Hayut/Plaksin/Vierling-Claassen:
% strength is RthFSRS*(gThRS) -> Can be omitted by setting RthFSRS = 0

% Coupling with ultrasound is included via the SONIC-model (Lemaire et al. 2019 [7]).

% References
% [1] Tarnaud et al. (2021), Improved alpha-beta power reduction via 
% combined Electrical and Ultrasonic stimulation in a Parkinsonian Cortex-Basal 
% Ganglia-Thalamus Computational Model, Journal of Neural Engineering.
% [2] Kumaravelu et al. (2016) A biophysical model of the cortex-basal 
% ganglia-thalamus network in the 6-OHDA lesioned rat model of Parkinson's
% disease. Journal of computational neuroscience, 40(2), 207-229
% [3] Plaksin et aL. (2016) Cell-type-selective effects of intramembrane 
% cavitation as a unifying theoretical framework for ultrasonic 
% neuromodulation. Eneuro, 3(3).
% [4] Hayut et al. (2011) LTS and FS inhibitory interneurons, short-term synaptic
% plasticity, and cortical circuit dynamics. PLoS computational biology, 7(10), e1002248.
% [5] Vierling-Claassen (2010) Computational modeling of distinct neocortical oscillations 
% driven by cell-type selective optogenetic drive: separable resonant circuits controlled 
% by low-threshold spiking and fast-spiking interneurons. Frontiers in human neuroscience, 4, 198.
% [6] Otsuka, T. (2004). Conductance-based model of the voltage-dependent generation 
% of a plateau potential in subthalamic neurons. Journal of neurophysiology, 
% 92(1), 255-264.
% [7] Lemaire, T. et al. (2019). Understanding ultrasound neuromodulation using 
% a computationally efficient and interpretable model of intramembrane cavitation. 
% Journal of neural engineering, 16(4), 046007.

% --------------------------MATLAB PATH------------------------------------
Chronuxpath = '';                                            % Chronux path
DynasimPath = '';                                            % Dynasim path   

% Spectral analysis (as in Kumaravelu et al. (2016) modelDB code)
addpath(genpath(Chronuxpath));
% Dynasim toolbox (Sherfey et al. (2018))
addpath(genpath(DynasimPath));
% --------------------------MODIFICATION OF GAINS--------------------------
% Gains modified w.r.t. the original papers in Kumaravelu, Plaksin, Hayut
% and Vierling-Claassen, as ratios between the new and original gains.
% A. CORTEX
gCORRCTX = [1,1,1,1,1,1,1,1]; 
% 1. RS->RS 2. RS->FS 3. RS->LTS 4. FS->RS 5. FS->FS 6. FS->LTS 7. LTS->RS 8. LTS->FS
% B. Basal ganglia
gCORRBG = [1,1,1,5,1.5,1,1.5,1,0.75,1]; 
% 1. GPi->ThRT 2. STN->GPe (AMPA) 3. STN->GPe (NMDA) 4. GPe->GPe 5. StrMSN->GPe 6. STN->GPi 7. GPe->GPi 8. StrMSN->GPi 9. GPe->STN 10. StrMSN->StrMSN
% C. Cortex to and from basal ganglia 
gCORRCTXBG = [1,1,1.2,1.5,1.5];
% 1. RS->STN (AMPA) 2. RS->STN (NMDA) 3. RS->StrMSN 4. ThRT->RS 5. ThRT->FS
% -------------------------CORRECTION OF DELAYS----------------------------
CTXdelay = 1e-3;        % [ms] Cortical delays -  0 ms in Vierling claassen
% --------------------------DOPAMINE DEPLETIONS----------------------------
PDon = 1;           % 1 if PD is on  -> Also manually adjust the IM_StrMSN current gain and rndized input current to zero
% Gains affected by dopamine deplection, implemented as ratios between new
% and original gains
% A. CORTEX
if PDon, DACTX = [1,1,1,1,1,1,1,1]; else, DACTX = [1,1,1,1,1,1,1,1]; end 
% 1. RS->RS 2. RS->FS 3. RS->LTS 4. FS->RS 5. FS->FS 6. FS->LTS 7. LTS->RS 8. LTS->FS
% B. Basal ganglia
if PDon, DABG = [1,10,10,0.2,1,5,2,1,10,1]; else, DABG = [1,1,1,1,1,1,1,1,1,1]; end 
% 1. GPi->ThRT 2. STN->GPe (AMPA) 3. STN->GPe (NMDA) 4. GPe->GPe 5. StrMSN->GPe 6. STN->GPi 7. GPe->GPi 8. StrMSN->GPi 9. GPe->STN 10. StrMSN->StrMSN
% C. Cortex to and from basal ganglia 
if PDon, DACTXBG = [5,5,0.5,2,2]; else, DACTXBG = [1,1,1,1,1]; end 
% 1. RS->STN (AMPA) 2. RS->STN (NMDA) 3. RS->StrMSN 4. ThRT->RS 5. ThRT->FS

% Applied current
CorrIappPD = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
if PDon, CorrIappPD = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.02, -0.02, -0.02, 0]; end  % [A] Correction to applied currents due to PD 
% ---------------------------SOLVER OPTIONS--------------------------------
study_dir = mfilename;
addpath(genpath(pwd));
time_end = 11;   
%dt = min(25e-6,0.1*min(USdc/USprf,ESdc/ESprf));            % Discretisation time (s)
dt = 10e-6;
atol = 1e-6; rtol = 1e-6; % absolute and relative VSVO-tolerances
solver = 'euler';
randness = 0;           % Randomness factor for synaptic connections: 1 is full randomness, 0 is minimal randomness
% In more detail: randness=1 will take all connections uniformly at random,
% potentially resulting in some neurons having no synaptic afferents. 0
% will ensure always the same number of afferents chosen at random
% ---------------------------GENERAL---------------------------------------
Rg = 8.314;             % Universal gas constant (J/(K*mol))
Far = 96485.3329;       % Faraday constant (C/mol) 
Temp = 273.15+36;       % Surrounding medium temperature (K)
cCae = 2;               % Extracellular Ca2+ concentration (mol/m^3)
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Nbuffer = 10;            % Number of spikes monitored in buffer
T_trans = 1;            % [s] Transient time excluded from spectral analysis and FR calcs
T_max = 11;              % [s] Maximal time included in the spectral analysis and FR calcs
% ---------------------------ULTRASONIC DBS--------------------------------
modelName = 'STN';      % Insonicated nucleus name 
USps=3.5;                 % Start of insonication [s]
USpd = 5;               % Insonication duration [s]
USfreq = 700e3;         % Fundamental frequency [Hz] 
USdc = 1;    % Duty cycle [-]
USprf = 0; % Pulse repetition frequency [Hz]
USisppa = USiRange(iii);            % Spatial peak pulse averaged intensity [W/m^2] 

aBLS = 32e-9;           % Sonophore radius [m]  
fBLS = 1;               % Sonophore coverage fraction [-]
corrPec = 0;            % Correct electrostatic pressure [Boolean] 
proteinMode = 0;        % Protein mode [0,1, or 2]  - TODO - now only default 0
gateMultip = 1;         % Gate multiplier [-] - TODO - now only default 1
maxRate = 1e6;          % Max allowable rate [1/s]
% ---------------------------ELECTRICAL DBS--------------------------------
dbsA = [0,0,0*3,0,0,0*3,0,0,0,0,0,0*3,0,0,0]; % Intensity [A/m^2]
dbsPW = [0,0,0*300e-6,0,0,0*300e-6,0,0,0,0,0,100e-6,0,0,0]; % Pulse width [s]  
dbsF = [0,0,0*1,0,0,0*1,0,0,0,0,0,100,0,0,0]; % Frequency [Hz]
dbsPS = [0,0,0,0,0,0,0,0,0,0,0,3.5,0,0,0]; % Pulse start [s]
dbsPD = [0,0,0,0,0,0,0,0,0,0,0,5,0,0,0]; % Pulse duration [s]
%----------------------------Randomized input current----------------------
RndOn = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % Boolean: 1 if random input
rA = 1; rPS = 0; rPD = time_end;           % [A, s, s] Amplitude, pulse start, pulse duration 
rF = 20; rPW = 1/40; rsdF = 3; rsdPW = 1/80; % [Hz, s, Hz, s] Frequency, pulse width, standard deviation (sd) on frequency, sd on pulse width
randESi(0,rA,rPS,rPD,rF,rPW,rsdF,rsdPW,time_end);       % Initialise randESi
% ---------------------------NEURON POPULATIONS----------------------------
NNames = {'RS','RSFVX','RSIzh','FS','FSFVX','FSIzh',...
    'LTS','LTSCAX','ThTC','ThRE','ThRT','STN','GPi',...
    'GPe','StrMSN'}; 
NNumbers = [20,0,0,5,0,0,5,0,0,0,20,20,20,20,20];
Cm0 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];	% Rest capacitance (F/m^2)	
% Cell membrane areas (m^2) (If unknown: use 1 and express currents/gains per unit area)
cma=[11.88*10^(-9),11.88*10^(-9),1,10.17*10^(-9),10.17*10^(-9),1,25*10^(-9),...
    25*10^(-9),29*10^(-9),14*10^(-9),1,1,1,1,1];
RthFSRS = 1.4; % FS to RS thalamic input current ratio RTH
IthRS = 0.17*10^(-9); % Thalamic DC current input to the RS neuron (A) 
% Injected current (A/m^2) ---- Current implementation: DC-current injection
Iapp = [0*IthRS, 0, 0, 0*RthFSRS*IthRS, 0, 0, 0, 0, 0, 0, 0.012, 0, 0.03, 0.03, 0]./cma;
Iapp = Iapp+CorrIappPD./cma;
Vthresh = [-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10]; % Threshold [mV]
V0 = [-62,-62,-62,-62,-62,-62,-62,-62,-62,-62,-62,-62,-62,-62,-63.8];
Nmechanisms = {{'INa_RS','IK_RS','IM_RS','Ileak_RS'};
               {'INa_RSFVX','IK_RSFVX','IM_RSFVX','Ileak_RSFVX'};
               {};
               {'INa_FS','IK_FS','IM_FS','Ileak_FS'};
               {'INa_FSFVX','IK_FSFVX','Ileak_FSFVX'}; 
               {};
               {'INa_LTS','IK_LTS','IM_LTS','IT_LTS','Ileak_LTS'};
               {'INa_LTSCAX','IK_LTSCAX','IM_LTSCAX','IT_LTSCAX','Ileak_LTSCAX'};
               {'INa_ThTC','IK_ThTC','IT_ThTC','IKL_ThTC','Ih_ThTC','Ileak_ThTC','CA_ThTC'};
               {'INa_ThRE','IK_ThRE','IT_ThRE','Ileak_ThRE','CA_ThRE'};
               {'INa_ThRT','IK_ThRT','IT_ThRT','Ileak_ThRT'};
               {'INa_STN','IK_STN','IT_STN','IKCa_STN','IA_STN','IL_STN','Ileak_STN','CA_STN'};
               {'INa_GPi','IK_GPi','IT_GPi','ICa_GPi','Iahp_GPi','Ileak_GPi','CA_GPi'};
               {'INa_GPe','IK_GPe','IT_GPe','ICa_GPe','Iahp_GPe','Ileak_GPe','CA_GPe'};
               {'INa_StrMSN','IK_StrMSN','IM_StrMSN','Ileak_StrMSN'};}; 
eqnsHH = {'dV/dt = (10^3*Iapp+@current+@isyn+10^3*dbsA*(t>=dbsPS&(t<=(dbsPS+dbsPD))).*(mod(t-dbsPS,1/dbsF)<=dbsPW)+10^3*RndOn*randESi(t))/Cm; V(0) = V0+5*randn(N_pre,1)';
        'monitor V.spikes(Vthresh,Nbuffer)'
        };
eqnsHH_BLS = {'dQ/dt = (Iapp+0.001*@current+0.001*@isyn+dbsA*(t>=dbsPS&(t<=(dbsPS+dbsPD))).*(mod(t-dbsPS,1/dbsF)<=dbsPW)+RndOn*randESi(t)); Q(0) = (0.001*V0)*Cm+(0.005*Cm)*randn(N_pre,1)';
        'monitor Q.spikes(10^(-3)*Vthresh*Cm,Nbuffer)' 
        };
eqnsIF = {'dV/dt = (0.4*V.^2+50*V+1400-10*u+10^3*Iapp+@isyn+10^3*dbsA*(t>=dbsPS&(t<=(dbsPS+dbsPD))).*(mod(t-dbsPS,1/dbsF)<=dbsPW)+10^3*RndOn*randESi(t))/Cm; V(0) = c';
  'du/dt = a*(b*V-u); u(0) = b*c;';
  'if(V>VtrIF)(V=c; u=u+d)'; 
  'monitor V.spikes(Vthresh,Nbuffer)'
  };
eqnsHH = replace(eqnsHH,'Nbuffer',num2str(Nbuffer));
eqnsHH_BLS = replace(eqnsHH_BLS,'Nbuffer',num2str(Nbuffer));
eqnsIF = replace(eqnsIF,'Nbuffer',num2str(Nbuffer));

aIF= [0,0,20,0,0,100,0,0,0,0,0,0,0,0,0];
bIF = [0,0,0.2,0,0,0.2,0,0,0,0,0,0,0,0,0];
cIF = [0,0,-65,0,0,-65,0,0,0,0,0,0,0,0,0];
dIF = [0,0,8,0,0,2,0,0,0,0,0,0,0,0,0];
VtrIF = [0,0,30,0,0,30,0,0,0,0,0,0,0,0,0];
eqns = vertcat(repmat({eqnsHH},2,1),{eqnsIF},repmat({eqnsHH},2,1),{eqnsIF},repmat({eqnsHH},9,1));

% SONIC tables   --- Now: only STN-ultrasonic DBS
if (USisppa~=0)&&(USpd>0)&&((USps+USpd)>0)&&(USps<time_end)
indInsonic = find(contains(NNames,modelName));
eqns{indInsonic} = eqnsHH_BLS;                          % Use modified HH-BLS equation
Nmechanisms{indInsonic} = cellfun(@(X) [X '_BLS'],Nmechanisms{indInsonic},'UniformOutput',0); % Updated mechanisms
if fBLS < 1
if (corrPec)
SONIC = load(['SONIC-' modelName '-xfs-corrPec.mat']); %#ok<*UNRCH>
elseif (corrPec == 0)
SONIC = load(['SONIC-' modelName '-xfs.mat']);
end
else
SONIC = load(['SONIC-' modelName '.mat']);
end
SONICtable = SONIC.SONICtable;
% 2.1 SONIC functions (rate, Veff, Zeff, Cmeff, ngend)
QmRange = SONICtable.QmRange; USPaRange = SONICtable.USPaRange; 
USfreqRange = SONICtable.USfreqRange; aBLSRange = SONICtable.aBLSRange;
fBLSRange = SONICtable.fBLSRange;

Veff5D = SONICtable.Veff; Zeff5D = SONICtable.Zeff; Cmeff5D = SONICtable.Cmeff; ngend5D = SONICtable.ngend;

% rate 4D sonic tables
SONICfields = fieldnames(SONICtable);
SONICrates = sort(SONICfields(cellfun(@(X) contains(X,'a_')|contains(X,'apb_'),SONICfields)));
SONICgates = cellfun(@(X) X(3:end),SONICrates(cellfun(@(X) contains(X,'a_'),SONICrates)),'UniformOutput',0); 
rt = struct;
for i = 1:length(SONICrates)
rt.(SONICrates{i}) = min(SONICtable.(SONICrates{i}),maxRate);
end
f5Veff = @(Qm,USPa,USfreq,aBLS,fBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,fBLSRange,Veff5D,Qm,USPa,USfreq,aBLS,fBLS,'linear');
f5Zeff = @(Qm,USPa,USfreq,aBLS,fBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,fBLSRange,Zeff5D,Qm,USPa,USfreq,aBLS,fBLS,'linear');
f5Cmeff = @(Qm,USPa,USfreq,aBLS,fBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,fBLSRange,Cmeff5D,Qm,USPa,USfreq,aBLS,fBLS,'linear');
f5ngend = @(Qm,USPa,USfreq,aBLS,fBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,fBLSRange,ngend5D,Qm,USPa,USfreq,aBLS,fBLS,'linear');

f3Veff = @(Qm,USPa,USfreq) f5Veff(Qm,USPa,USfreq,aBLS,fBLS); 
f3Zeff = @(Qm,USPa,USfreq) f5Zeff(Qm,USPa,USfreq,aBLS,fBLS); %#ok<*NASGU>
f3Cmeff = @(Qm,USPa,USfreq) f5Cmeff(Qm,USPa,USfreq,aBLS,fBLS);
f3ngend = @(Qm,USPa,USfreq) f5ngend(Qm,USPa,USfreq,aBLS,fBLS);

f5rt = struct; f3rt = struct;
VecVeff0 = zeros(1,length(QmRange)); VecVeffPa = zeros(1,length(QmRange));
Vecrt0 = struct; VecrtPa = struct; f1rt0 = struct; f1rtPa = struct;
for i = 1:length(SONICrates)
f5rt.(SONICrates{i}) =  @(Qm,USPa,USfreq,aBLS,fBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,fBLSRange,rt.(SONICrates{i}),Qm,USPa,USfreq,aBLS,fBLS,'linear');   
f3rt.(SONICrates{i}) = @(Qm,USPa,USfreq) f5rt.(SONICrates{i})(Qm,USPa,USfreq,aBLS,fBLS); 
end

USPa = sqrt(2*rhol*c*USisppa);
USstep(0,USps,USpd,USdc,USprf); % Initialise USstep;

for i = 1:length(QmRange)
VecVeff0(i) = f3Veff(QmRange(i),0,USfreq);
VecVeffPa(i) = f3Veff(QmRange(i),USPa,USfreq);
end
for j = 1:length(SONICrates)
Vecrt0.(SONICrates{j}) = zeros(1,length(QmRange));
VecrtPa.(SONICrates{j}) = zeros(1,length(QmRange));
for i = 1:length(QmRange)
Vecrt0.(SONICrates{j})(i) = f3rt.(SONICrates{j})(QmRange(i),0,USfreq);
VecrtPa.(SONICrates{j})(i) = f3rt.(SONICrates{j})(QmRange(i),USPa,USfreq);
end
end
f1rt(0,'',0,QmRange,Vecrt0,VecrtPa);        % Initialise f1rt
f1Veff(0,0,QmRange,VecVeff0,VecVeffPa);     % Initialise f1Veff
end
%--------------------------Connection scheme-------------------------------
ibiexp={
  'tauD=0; tauR = 0; gSYN = 0; gPROB = 0; tau_delay = 0; ESYN = 0; netcon = zeros(N_pre,N_post)'  
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(V_post) = (gSYN*gPROB).*(sum(f(t-tspike_pre-tau_delay),1)*netcon).*(V_post-ESYN)'
  %'monitor Isyn'
  '@isyn += -Isyn(V_post)'
};
ibiexpP={
  % Plasticity (Varela et al. 1997)
  'tausf = inf; tausd1 = inf; tausd2 = inf; tausd3 = inf;'
  'sf = 0; sd1 = 1; sd2 = 1; sd3 = 1;'
  'dF/dt = (1-F)/tausf+(sf/dt)*sum(((t-tspike_pre)<=1.1*dt),1);'
  'dD1/dt = (1-D1)/tausd1+(D1.*sd1.^sum(((t-tspike_pre)<=1.1*dt),1)-D1)/dt;'
  'dD2/dt = (1-D2)/tausd2+(D2.*sd2.^sum(((t-tspike_pre)<=1.1*dt),1)-D2)/dt;'
  'dD3/dt = (1-D3)/tausd3+(D3.*sd3.^sum(((t-tspike_pre)<=1.1*dt),1)-D3)/dt;'
  'F(0)=ones(1,N_pre); D1(0) = ones(1,N_pre); D2(0) = ones(1,N_pre); D3(0) = ones(1,N_pre);'
  % Synapse
  'tauD=0; tauR = 0; gSYN = 0; gPROB = 0; tau_delay = 0; ESYN = 0; netcon = zeros(N_pre,N_post)'  
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(V_post) = (gSYN*gPROB).*(((F.*D1.*D2.*D3).*sum(f(t-tspike_pre-tau_delay),1))*netcon).*(V_post-ESYN)'
  %'monitor Isyn'
  '@isyn += -Isyn(V_post)'
};
ialpha={
  'tauS = 0; gSYN = 0; gPROB = 0; tau_delay = 0; ESYN = 0; netcon = zeros(N_pre,N_post)'
  'f(x) = (x/tauS).*(exp(-x/tauS)).*(x>0)'
  'Isyn(V_post) = exp(1)*(gSYN*gPROB).*(sum(f(t-tspike_pre-tau_delay),1)*netcon).*(V_post-ESYN)'
  %'monitor Isyn'
  '@isyn += -Isyn(V_post)'
};
ida = {
  'alphaDA = 2e3; betaDA = 1/(13e-3); sDAIC = 0; sDANoiseIC = 0;'
  'gSYN = 0; gPROB = 0; ESYN = -80; netcon = zeros(N_pre,N_post)'
  'ds/dt = alphaDA.*(1 + tanh(V_pre./4)).*(1-s) - betaDA.*s'
  's(0) = sDAIC+sDANoiseIC.*rand(1,N_pre)'
  'Isyn(V_post) = (gSYN*gPROB).*(s*netcon).*(V_post-ESYN)'
  %'monitor Isyn'
  '@isyn += -Isyn(V_post)'
};

% Structured connections - SC : synaptic connections, (s) -> reflective, pt -> partial target, (d): double matrix
% CUHtA/CLHtA : connection UPPER/LOWER half to all
dex = @(M,i,j) M(i,j);
SC = @(x,y,z) SynC(NNumbers(x),NNumbers(y),z);
SCs = @(x,y,z) circshift(SynC(NNumbers(x),NNumbers(y),z),1,1);
SCpt = @(x,y,z,pt) bsxfun(@times,SC(x,y,z),dex([ones(1,pt),zeros(1,NNumbers(y)-pt)],1,randperm(NNumbers(y))));
SCd = @(x,y,zU,zL) [SynC(NNumbers(x)/2,NNumbers(y)/2,zU),zeros(NNumbers(x)/2,NNumbers(y)/2);
                    zeros(NNumbers(x)/2,NNumbers(y)/2),SynC(NNumbers(x)/2,NNumbers(y)/2,zL)];
CUHtA = @(i,j) [ones(NNumbers(i)/2,NNumbers(j));zeros(NNumbers(i)/2,NNumbers(j))];
CLHtA = @(i,j) [zeros(NNumbers(i)/2,NNumbers(j));ones(NNumbers(i)/2,NNumbers(j))];
% Random connections - SCrand, (s) -> reflective
SCrand = @(x,y,P) double((randness).*(rand(NNumbers(x),NNumbers(y))<P)+...
    (~randness).*rndize(SC(x,y,ceil(P*NNumbers(x)))));
refl = @(A) A-diag(diag(A)); 
SCrands = @(x,y,P) double((randness).*refl(SCrand(x,y,P))+...
    (~randness).*rndize(SCs(x,y,ceil(P*(NNumbers(x)-1))),1));

CMechCTX = {'AMPAbiRSRS';'AMPAbiRSFS';'AMPAbiRSLTS';'GABAAbiFSRS';
    'GABAAbiFSFS';'GABAAbiFSLTS';'GABAAbiLTSRS';'GABAAbiLTSFS'};
CMech2CTX = {'AMPAbiRSRS';'AMPAbiRSFS';'AMPAbiRSLTS';'GABAAbiFSRS';
    'GABAAbiFSFS';'GABAAbiFSLTS';'GABAAbiLTSRS';'GABAAbiLTSFS'};
CSchemeCTX = {'RS->RS';'RS->FS';'RS->LTS';'FS->RS';
    'FS->FS';'FS->LTS';'LTS->RS';'LTS->FS'}; % in S-matrices: first idstr, second dstr ->netcon: upper half: idstr
CeqnsCTX = {ibiexp;ibiexpP;ibiexpP;ibiexp;ibiexp;ibiexp;ibiexp;ibiexp};

CMechBG = {'GABAAalGPiThRT';{'AMPAbiSTNGPe','NMDAbiSTNGPe'};'GABAAalGPeGPe';
    'GABAAalStrMSNGPe';'AMPAalSTNGPi';'GABAAalGPeGPi';'GABAAalStrMSNGPi';
    'GABAAbiGPeSTN';'GABAAdaStrMSNStrMSN'};
CMech2BG = {'GABAAalGPiThRT';'AMPAbiSTNGPe';'NMDAbiSTNGPe';'GABAAalGPeGPe';
    'GABAAalStrMSNGPe';'AMPAalSTNGPi';'GABAAalGPeGPi';'GABAAalStrMSNGPi';
    'GABAAbiGPeSTN';'GABAAdaStrMSNStrMSN';};
CSchemeBG = {'GPi->ThRT';'STN->GPe';'GPe->GPe';'StrMSN->GPe';
    'STN->GPi';'GPe->GPi';'StrMSN->GPi';'GPe->STN';'StrMSN->StrMSN'}; % in S-matrices: first idstr, second dstr ->netcon: upper half: idstr
CeqnsBG = {ialpha;ibiexp;ibiexp;ialpha;ialpha;ialpha;ialpha;ialpha;ibiexp;ida};
%% Cortex synaptical parameters
% Synaptic gain calculations from Plaksin et al. (2016) and convergence
% rates from Vierling-Claassen et al. (2010)
gRSRSp = 0.002*10^(-6)/cma(1); % Total maximal synaptic conductance used for RS to RS connection (S/m^2)
gFSRSp = 0.04*10^(-6)/cma(4); % RS to FS connection (S/m^2)
gLTSRSp = 0.09*10^(-6)/cma(7); % RS to LTS connection  (S/m^2)
gRSFSp = 0.015*10^(-6)/cma(1); % FS to RS connection (S/m^2)
gFSFSp = 0.135*10^(-6)/cma(4); % FS to FS connection (S/m^2)
gLTSFSp = 0.86*10^(-6)/cma(7); % FS to LTS connection (S/m^2)
gRSLTSp = 0.135*10^(-6)/cma(1); % LTS to RS connection (S/m^2)
gFSLTSp = 0.02*10^(-6)/cma(4); % LTS to FS connection (S/m^2)

cRSRS_VC = 76/36; cFSRS_VC = 93/6; cLTSRS_VC = 123/6;
cRSFS_VC = 95/36; cFSFS_VC = 15/6; cLTSFS_VC = 13/6;
cRSLTS_VC = 76/36; cFSLTS_VC = 22/6; cLTSLTS_VC = 0; % Vierling-Claassen convergence rates

% Vierling-Claassen individual rates (Valid for N_RS = 36, N_FS = 6, N_LTS = 6)
gRSRS = gRSRSp/(cRSRS_VC); % Maximal synaptic conductance used for RS to RS connection (S/m^2)
gFSRS = gFSRSp/(cFSRS_VC); % RS to FS connection (S/m^2)
gLTSRS = gLTSRSp/(cLTSRS_VC); % RS to LTS connection  (S/m^2)
gRSFS = gRSFSp/(cRSFS_VC); % FS to RS connection (S/m^2)
gFSFS = gFSFSp/(cFSFS_VC); % FS to FS connection (S/m^2)
gLTSFS = gLTSFSp/(cLTSFS_VC); % FS to LTS connection (S/m^2)
gRSLTS = gRSLTSp/(cRSLTS_VC); % LTS to RS connection (S/m^2)
gFSLTS = gFSLTSp/(cFSLTS_VC); % LTS to FS connection (S/m^2)
gLTSLTS = 0; % LTS to LTS connection (S/m^2)

% Rescaling of Vierling_Claassen gains with number of neurons 
gRSRS = gRSRS*(36-1)/(NNumbers(1)-1);
gFSRS = gFSRS*(36/NNumbers(1));
gLTSRS = gLTSRS*(36/NNumbers(1));
gRSFS = gRSFS*(6/NNumbers(4));
gFSFS = gFSFS*(6-1)/(NNumbers(4)-1);
gLTSFS = gLTSFS*(6/NNumbers(4));
gRSLTS = gRSLTS*(6/NNumbers(7));
gFSLTS = gFSLTS*(6/NNumbers(7));
gLTSLTS = gLTSLTS*(6-1)/(NNumbers(7)-1);


taurAMPArs = 0.1*10^(-3); % AMPA rise constant (s)
taudAMPArs = 3*10^(-3); % AMPA decay constant (s)
taurGABAAfs = 0.5*10^(-3); % GABA-A rise constant (s) for FS
taudGABAAfs = 8*10^(-3); % GABA-A decay constant (s) for FS
taurGABAAlts = 0.5*10^(-3); % GABA-A rise constant (s) for LTS
taudGABAAlts = 50*10^(-3); % GABA-A decay constant (s) for LTS
sfLTSRS = 0.2; % Short term synaptic  facilitation factor (-) (RS->LTS)
tausfLTSRS = 200*10^(-3); % Short term synaptic facilitation time constant (s) (RS->LTS)
sfFSRS = 0.5; % Short term synaptic facilitation factor (-) (RS->FS)
tausfFSRS = 94*10^(-3);  % Short term synaptic facilitation factor time constant (s) (RS->FS)
sd1FSRS = 0.46;  % Short term synaptic depression factor (-) (RS->FS)
sd2FSRS = 0.975; % Long term synaptic depression factor (-) (RS->FS)
tausd1FSRS = 380*10^(-3); tausd2FSRS = 9200*10^(-3); % Short and long term synaptic time constant (s) (RS-FS)

Anorm = @(tr,td) ((tr/td)^(-tr/(tr-td))-(tr/td)^(-td/(tr-td)))^(-1);

%% Connection parameters
NetConSTNGPe = SCpt(12,14,2,NNumbers(14)/2);

CparamsCTX = {{'netcon',SCrands(1,1,0.06),'tauR',taurAMPArs,'tauD',taudAMPArs,'tau_delay',CTXdelay,'gSYN',(DACTX(1)*gCORRCTX(1))*gRSRS,'gPROB',Anorm(taurAMPArs,taudAMPArs),'ESYN',0};
{'netcon',SCrand(1,4,0.43),'tauR',taurAMPArs,'tauD',taudAMPArs,'tau_delay',CTXdelay,'gSYN',(DACTX(2)*gCORRCTX(2))*gFSRS,'gPROB',Anorm(taurAMPArs,taudAMPArs),'ESYN',0,'sf',sfFSRS,'tausf',tausfFSRS,'sd1',sd1FSRS,'sd2',sd2FSRS,'tausd1',tausd1FSRS,'tausd2',tausd2FSRS}
{'netcon',SCrand(1,7,0.57),'tauR',taurAMPArs,'tauD',taudAMPArs,'tau_delay',CTXdelay,'gSYN',(DACTX(3)*gCORRCTX(3))*gLTSRS,'gPROB',Anorm(taurAMPArs,taudAMPArs),'ESYN',0,'sf',sfLTSRS,'tausf',tausfLTSRS};
{'netcon',SCrand(4,1,0.44),'tauR',taurGABAAfs,'tauD',taudGABAAfs,'tau_delay',CTXdelay,'gSYN',(DACTX(4)*gCORRCTX(4))*gRSFS,'gPROB',Anorm(taurGABAAfs,taudGABAAfs),'ESYN',-85};
{'netcon',SCrands(4,4,0.51),'tauR',taurGABAAfs,'tauD',taudGABAAfs,'tau_delay',CTXdelay,'gSYN',(DACTX(5)*gCORRCTX(5))*gFSFS,'gPROB',Anorm(taurGABAAfs,taudGABAAfs),'ESYN',-85};
{'netcon',SCrand(4,7,0.36),'tauR',taurGABAAfs,'tauD',taudGABAAfs,'tau_delay',CTXdelay,'gSYN',(DACTX(6)*gCORRCTX(6))*gLTSFS,'gPROB',Anorm(taurGABAAfs,taudGABAAfs),'ESYN',-85};
{'netcon',SCrand(7,1,0.35),'tauR',taurGABAAlts,'tauD',taudGABAAlts,'tau_delay',CTXdelay,'gSYN',(DACTX(7)*gCORRCTX(7))*gRSLTS,'gPROB',Anorm(taurGABAAlts,taudGABAAlts),'ESYN',-85};
{'netcon',SCrand(7,4,0.61),'tauR',taurGABAAlts,'tauD',taudGABAAlts,'tau_delay',CTXdelay,'gSYN',(DACTX(8)*gCORRCTX(8))*gFSLTS,'gPROB',Anorm(taurGABAAlts,taudGABAAlts),'ESYN',-85}};

CparamsBG = {{'netcon',SC(13,11,1),'tauS',5e-3,'tau_delay',5e-3,'gSYN',(DABG(1)*gCORRBG(1))*0.3360,'gPROB',1,'ESYN',-85};
           {'AMPAbiSTNGPe.netcon',NetConSTNGPe,'AMPAbiSTNGPe.tauR',0.4e-3,'AMPAbiSTNGPe.tauD',2.5e-3,'AMPAbiSTNGPe.gSYN',(DABG(2)*gCORRBG(2))*2.1772,'AMPAbiSTNGPe.gPROB',rand(1,NNumbers(12)),'AMPAbiSTNGPe.tau_delay',2e-3,'AMPAbiSTNGPe.ESYN',0,...
           'NMDAbiSTNGPe.netcon',NetConSTNGPe,'NMDAbiSTNGPe.tauR',2e-3,'NMDAbiSTNGPe.tauD',67e-3,'NMDAbiSTNGPe.gSYN',(DABG(3)*gCORRBG(3))*0.0099,'NMDAbiSTNGPe.gPROB',rand(1,NNumbers(12)),'NMDAbiSTNGPe.tau_delay',2e-3,'NMDAbiSTNGPe.ESYN',0};
           {'netcon',SCs(14,14,2),'tauS',5e-3,'tau_delay',1e-3,'gSYN',(DABG(4)*gCORRBG(4))*0.375,'gPROB',1,'ESYN',-85};
           {'netcon',CUHtA(15,14),'tauS',5e-3,'tau_delay',5e-3,'gSYN',(DABG(5)*gCORRBG(5))*1.5,'gPROB',1,'ESYN',-85};
           {'netcon',SCpt(12,13,2,NNumbers(13)/2),'tauS',5e-3,'tau_delay',1.5e-3,'gSYN',(DABG(6)*gCORRBG(6))*1.29,'gPROB',rand(1,NNumbers(12)),'ESYN',0};
           {'netcon',SC(14,13,2),'tauS',5e-3,'tau_delay',3e-3,'gSYN',(DABG(7)*gCORRBG(7))*1.5,'gPROB',1,'ESYN',-85};
           {'netcon',CLHtA(15,13),'tauS',5e-3,'tau_delay',4e-3,'gSYN',(DABG(8)*gCORRBG(8))*1.5,'gPROB',1,'ESYN',-85};
           {'netcon',SC(14,12,2),'tauR',0.4e-3,'tauD',7.7e-3,'gSYN',(DABG(9)*gCORRBG(9))*1.8605,'gPROB',1,'tau_delay',4e-3,'ESYN',-85};
           {'netcon',SCd(15,15,4,3),'alphaDA',2e3,'betaDA',76.9231,'gSYN',(DABG(10)*gCORRBG(10))*1,'gPROB',[(1/4).*ones(1,NNumbers(15)/2),(1/3).*ones(1,NNumbers(15)/2)],'ESYN',-80};
           }; 
       % PD : gSYN_gege : 0.375 -> 1.5, gSYN_RSIzhStrMSN :  0.301->0.1118,
       % I_M current of Striatum gain (IM_StrMSN.mech)
%% Synaptic coupling between CTX and BG
CparamsCTXBG={{'AMPAbiRSSTN.netcon',SC(1,12,2),'AMPAbiRSSTN.tauR',0.5e-3,'AMPAbiRSSTN.tauD',2.49e-3,'AMPAbiRSSTN.gSYN',(DACTXBG(1)*gCORRCTXBG(1))*1.2081,'AMPAbiRSSTN.gPROB',1,'AMPAbiRSSTN.tau_delay',5.9e-3,'AMPAbiRSSTN.ESYN',0,...
               'NMDAbiRSSTN.netcon',SC(1,12,2),'NMDAbiRSSTN.tauR',2e-3,'NMDAbiRSSTN.tauD',90e-3,'NMDAbiRSSTN.gSYN',(DACTXBG(2)*gCORRCTXBG(2))*0.0144,'NMDAbiRSSTN.gPROB',1,'NMDAbiRSSTN.tau_delay',5.9e-3,'NMDAbiRSSTN.ESYN',0};           
              {'netcon',SC(1,15,1),'tauS',5e-3,'tau_delay',5.1e-3,'gSYN',(DACTXBG(3)*gCORRCTXBG(3))*0.301,'gPROB',1,'ESYN',0};
              {'netcon',SC(11,1,1),'tauS',5e-3,'tau_delay',5.6e-3,'gSYN',(DACTXBG(4)*gCORRCTXBG(4))*0.645,'gPROB',1,'ESYN',0};
              {'netcon',SC(11,4,1),'tauS',5e-3,'tau_delay',5.6e-3,'gSYN',(DACTXBG(5)*gCORRCTXBG(5))*(RthFSRS*0.645),'gPROB',1,'ESYN',0};
};
CMechCTXBG = {{'AMPAbiRSSTN','NMDAbiRSSTN'};'AMPAalRSStrMSN';'AMPAalThRTRS';'AMPAalThRTFS'};
CMech2CTXBG = {'AMPAbiRSSTN';'NMDAbiRSSTN';'AMPAalRSStrMSN';'AMPAalThRTRS';'AMPAalThRTFS'};
CSchemeCTXBG = {'RS->STN';'RS->StrMSN';'ThRT->RS';'ThRT->FS'};
CeqnsCTXBG = {ibiexp;ibiexp;ialpha;ialpha;ialpha};

Cparams = vertcat(CparamsCTX,CparamsBG,CparamsCTXBG);
CMech = vertcat(CMechCTX,CMechBG,CMechCTXBG);
CMech2 = vertcat(CMech2CTX,CMech2BG,CMech2CTXBG);
CScheme = vertcat(CSchemeCTX,CSchemeBG,CSchemeCTXBG);
Ceqns = vertcat(CeqnsCTX,CeqnsBG,CeqnsCTXBG);
%% --------------------------CORE DYNASIM SOLVER----------------------------
% DynaSim options
save_data_flag = 0; 
verbose_flag =   1; 
overwrite_flag = 1; 
parfor_flag =    0; 
% If you want to run a simulation using the batch system of your cluster, make 
%   sure to set the following options to something like this: 
%    'cluster_flag', 1,... 
%    'memory_limit', '254G',... 
%    'num_cores', 16,... 
cluster_flag =   0; 
memory_limit =   '8G'; 
num_cores =      2; 
RndSeed = 'shuffle'; RndSeed = 0;
vary = {};

spec = []; j = 0;
for i = 1:length(NNames)
if NNumbers(i)~=0
j=j+1;
spec.populations(j).name = NNames{i};
spec.populations(j).size = NNumbers(i);
spec.populations(j).equations = eqns{i};
spec.populations(j).mechanism_list = Nmechanisms{i};
spec.populations(j).parameters = {'Iapp',Iapp(i),'Cm',Cm0(i),'Vthresh',Vthresh(i),'VtrIF',VtrIF(i),...
    'a',aIF(i),'b',bIF(i),'c',cIF(i),'d',dIF(i),'cCae',cCae,'Rg',Rg,'Temp',Temp,...
    'Far',Far,'V0',V0(i),'dbsA',dbsA(i),'dbsPW',dbsPW(i),'dbsF',dbsF(i),'dbsPS',dbsPS(i),'dbsPD',dbsPD(i),'RndOn',RndOn(i)};
end
end
for i = 1:length(CScheme)
spec.connections(i).mechanism_list = CMech{i};
spec.connections(i).direction = CScheme{i};
spec.connections(i).parameters = Cparams{i};
end
for i = 1:length(Ceqns)
spec.mechanisms(i).name = CMech2{i};
spec.mechanisms(i).equations = Ceqns{i};
end

% Ultrasound effects
% Synapses:
if (USisppa~=0)&&(USpd>0)&&((USps+USpd)>0)&&(USps<time_end)
USpostSyn = find(endsWith(CMech2,'STN')); % In Ceqns, the numbers where the insonicated nucleus is postsynaptic
for j = 1:length(USpostSyn)
spec.mechanisms(USpostSyn(j)).equations = replace(spec.mechanisms(USpostSyn(j)).equations,{'@isyn += -Isyn(V_post)','V_post'},{'@isyn += -Isyn(f1Veff({Q_post USstep(t)}))','X'});
end
end

[data] = dsSimulate(spec,... 
     'save_data_flag',save_data_flag,'study_dir',study_dir,... 
     'vary',vary,'cluster_flag',cluster_flag,... 
     'verbose_flag',verbose_flag,'overwrite_flag',overwrite_flag,... 
     'tspan',[0 time_end],'solver',solver,'dt',dt,... 
     'parfor_flag',parfor_flag,'random_seed',RndSeed,... 
     'memory_limit',memory_limit,'num_cores',num_cores);
 
%     'plot_functions',{@dsPlot,@dsPlot,@dsPlot},... 
%      'plot_options',{{'plot_type','waveform','format','png'},... 
%                      {'plot_type','rastergram','format','png'},... 
%                      {'plot_type','power','format','png','xlim',[0 40]}}); 

if (USisppa~=0)&&(USpd>0)&&((USps+USpd)>0)&&(USps<time_end)
data.labels = horzcat(data.labels,'STN_V');
data.([modelName '_V_spikes']) = data.([modelName '_Q_spikes']); 
for j = 1:size(data.([modelName '_Q']),2)
for i = 1:length(data.time)
data.([modelName '_V'])(i,j) = f1Veff(double(data.([modelName '_Q'])(i,j)),USstep(data.time(i)));
end
end
end 

% Spectral power
dt1=0.01*10^-3;
params.Fs = 1/dt1; %Hz
params.fpass = [1 100];
params.tapers = [3 5];
params.trialave = 1;

GPi_APs = cell2struct(cellfun(@(X) X(X>=T_trans&X<=T_max),cellfun(@(X) data.time(~~X),mat2cell(data.GPi_V_spikes,size(data.GPi_V_spikes,1),ones(size(data.GPi_V_spikes,2),1)),'UniformOutput',false),'uniformoutput',false),'times');
STN_APs = cell2struct(cellfun(@(X) X(X>=T_trans&X<=T_max),cellfun(@(X) data.time(~~X),mat2cell(data.STN_V_spikes,size(data.STN_V_spikes,1),ones(size(data.STN_V_spikes,2),1)),'UniformOutput',false),'uniformoutput',false),'times');
GPe_APs = cell2struct(cellfun(@(X) X(X>=T_trans&X<=T_max),cellfun(@(X) data.time(~~X),mat2cell(data.GPe_V_spikes,size(data.GPe_V_spikes,1),ones(size(data.GPe_V_spikes,2),1)),'UniformOutput',false),'uniformoutput',false),'times');

[gpi_alpha_beta_area gpi_S gpi_f]=spectAnalysis.make_Spectrum(GPi_APs,params);
[stn_alpha_beta_area stn_S stn_f]=spectAnalysis.make_Spectrum(STN_APs,params);
[gpe_alpha_beta_area gpe_S gpe_f]=spectAnalysis.make_Spectrum(GPe_APs,params);

% Spectrogram
movingwin = [1,0.1];            % [s] - Moving window (first: width of movingwin, second: step size of window)
[gpi_S_gram,gpi_t_gram,gpi_f_gram]=mtspecgrampt(GPi_APs,movingwin,params);
[stn_S_gram,stn_t_gram,stn_f_gram]=mtspecgrampt(STN_APs,movingwin,params);
[gpe_S_gram,gpe_t_gram,gpe_f_gram]=mtspecgrampt(GPe_APs,movingwin,params);

% FR calculations
FiringRate = struct; % [Hz]
j = 0;
for i = 1:length(NNames)
if NNumbers(i)~=0
j = j+1;
FiringRate(j).('Population') = NNames(i);
FiringRate(j).('FR') = double(sum(data.([NNames{i} '_V_spikes'])(data.time>=T_trans&data.time<=T_max,:))./(T_max-T_trans));
FiringRate(j).('MeanFR') = mean(double(sum(data.([NNames{i} '_V_spikes'])(data.time>=T_trans&data.time<=T_max,:))./(T_max-T_trans)));
FiringRate(j).('stdFR') = std(double(sum(data.([NNames{i} '_V_spikes'])(data.time>=T_trans&data.time<=T_max,:))./(T_max-T_trans)));
end
end

%save(['data_PDon_dbs' num2str(dbsFreqRange(iii)) '.mat'],'data','-v7.3');
%save('data_PDon.mat','data','-v7.3');
%save('data_healthy.mat','data','-v7.3');


save(['data_PDon_USi' num2str(USiRange(iii)) '.mat'],'data','-v7.3');
% fieldNM = fieldnames(data);
% FN_spikes = fieldNM(contains(fieldNM,'_spikes')|contains(fieldNM,'time'));
% for i = 1:length(FN_spikes)
% data_spikes.(FN_spikes{i}) = data.(FN_spikes{i});
% end
% save(['data_PDon_USi' num2str(InPhaseUSiRange(iii)) '_ESi' num2str(InPhaseESiRange(jjj)) 'inPhase_spikes.mat'],'data_spikes','-v7.3');

rmdir run_ES_NN_CWUS_700kHz s
end
