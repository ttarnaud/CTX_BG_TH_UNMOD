% --------------------------MATLAB PATH------------------------------------
Chronuxpath = '';                                            % Chronux path
DynasimPath = '';                                            % Dynasim path   

% Spectral analysis (as in Kumaravelu et al. (2016) modelDB code)
addpath(genpath(Chronuxpath));
% Dynasim toolbox (Sherfey et al. (2018))
addpath(genpath(DynasimPath));

%% 1. Analysis of electrical DBS
% Params
T_trans = 3.5;
T_max = 8.5; 
tauES = 100;        % Electrostimulation pulsewidth [us]

dbs_specP = cell(16,1); dbs_specGram = cell(16,1); dbs_FiringRate = cell(16,1);
for i = 1:16
fileName = ['DBS ES/tauES_' num2str(tauES) 'us/data_PDon_dbs' num2str(i*10) '.mat'];
[specP,specGram,FiringRate]=funPostProc(fileName,T_trans,T_max);
dbs_specP{i} = specP; dbs_specGram{i} = specGram; dbs_FiringRate{i} = FiringRate;
end

dbs_data = struct; dbs_data.specP = dbs_specP; dbs_data.specGram = dbs_specGram;
dbs_data.FiringRate = dbs_FiringRate;
save(['dbs_data_tauES' num2str(tauES) 'us.mat'],'dbs_data');

%% 2. Parkinson's disease 
T_trans=3.5;
T_max = 8.5;
[PDon_specP,PDon_specGram,PDon_FiringRate]=funPostProc('PD/data_PDon.mat',T_trans,T_max);

PDon_data = struct; PDon_data.specP = PDon_specP; PDon_data.specGram = PDon_specGram;
PDon_data.FiringRate = PDon_FiringRate;
save('PDon_data.mat','PDon_data');

%% 3. Continuous wave ultrasound
T_trans = 3.5;
T_max = 8.5;
fUS = 700;          % Ultrasound frequency [kHz]

US_specP = cell(20,1); US_specGram = cell(20,1); US_FiringRate = cell(20,1);
for i = 1:20
fileName = ['CW US/fUS' num2str(fUS) 'kHz/data_PDon_USi' num2str(i*10) '.mat'];
[specP,specGram,FiringRate]=funPostProc(fileName,T_trans,T_max);
US_specP{i} = specP; US_specGram{i} = specGram; US_FiringRate{i} = FiringRate;
end

US_data = struct; US_data.specP = US_specP; US_data.specGram = US_specGram;
US_data.FiringRate = US_FiringRate;
save(['US_data_fUS' num2str(fUS) 'kHz.mat'],'US_data');

%% 4. Ultrasonic DBS
% Params
T_trans = 3.5;
T_max = 8.5; 
fUS = 700;          % Ultrasound frequency [kHz]
Ius = 400;          % Ultrasonic intensity [W/m^2]

USdbs_specP = cell(16,1); USdbs_specGram = cell(16,1); USdbs_FiringRate = cell(16,1);
for i = 1:16
fileName = ['DBS US/fUS' num2str(fUS) 'kHz/data_PDon_USdbs' num2str(i*10) '_spikes_Ius' num2str(Ius) '.mat'];
[specP,specGram,FiringRate]=funPostProc(fileName,T_trans,T_max,'spikes');
USdbs_specP{i} = specP; USdbs_specGram{i} = specGram; USdbs_FiringRate{i} = FiringRate;
end

USdbs_data = struct; USdbs_data.specP = USdbs_specP; USdbs_data.specGram = USdbs_specGram;
USdbs_data.FiringRate = USdbs_FiringRate;
save(['USdbs_data_Ius' num2str(Ius) '_fUS' num2str(fUS) 'kHz.mat'],'USdbs_data');

%% 5. Analysis of EUS DBS  --- Temporally outphase
% Params
T_trans = 3.5;
T_max = 8.5; 
fUS = 700;          % Ultrasound frequency [kHz]
tauES = 300;        % Electrostimulation pulsewidth [us]

EUSdbs_TempOutphase_specP = cell(16,1); EUSdbs_TempOutphase_specGram = cell(16,1); EUSdbs_TempOutphase_FiringRate = cell(16,1);
for i = 1:16
fileName = ['Outphase/fUS' num2str(fUS) 'kHz_tauES' num2str(tauES) 'us/data_PDon_EUSdbs' num2str(i*5) 'maxOutPhase_spikes.mat'];
[specP,specGram,FiringRate]=funPostProc(fileName,T_trans,T_max,'spikes');
EUSdbs_TempOutphase_specP{i} = specP; EUSdbs_TempOutphase_specGram{i} = specGram; EUSdbs_TempOutphase_FiringRate{i} = FiringRate;
end

EUSdbs_TempOutphase_data = struct; EUSdbs_TempOutphase_data.specP = EUSdbs_TempOutphase_specP; EUSdbs_TempOutphase_data.specGram = EUSdbs_TempOutphase_specGram;
EUSdbs_TempOutphase_data.FiringRate = EUSdbs_TempOutphase_FiringRate;
save(['EUSdbs_TempOutphase_data_fUS' num2str(fUS) 'kHz_tauES' num2str(tauES) 'us.mat'],'EUSdbs_TempOutphase_data');

%% 6. Analysis of EUS DBS  --- inphase
% Params
T_trans = 3.5;
T_max = 8.5; 
fUS = 700;          % Ultrasound frequency [kHz]
tauES = 100;        % Electrostimulation pulsewidth [us]

EUSdbs_inphase_specP = cell(11,11); EUSdbs_inphase_specGram = cell(11,11); EUSdbs_inphase_FiringRate = cell(11,11);
for i = 1:11
    for j = 1:11
fileName = ['Inphase/fUS' num2str(fUS) 'kHz_tauES' num2str(tauES) 'us/data_PDon_USi' num2str((i-1)*50) '_ESi' num2str((j-1)*0.2) 'inPhase_spikes.mat'];
[specP,specGram,FiringRate]=funPostProc(fileName,T_trans,T_max,'spikes');
EUSdbs_inphase_specP{i,j} = specP; EUSdbs_inphase_specGram{i,j} = specGram; EUSdbs_inphase_FiringRate{i,j} = FiringRate;
    end
end

EUSdbs_inphase_data = struct; EUSdbs_inphase_data.specP = EUSdbs_inphase_specP; EUSdbs_inphase_data.specGram = EUSdbs_inphase_specGram;
EUSdbs_inphase_data.FiringRate = EUSdbs_inphase_FiringRate;
save(['EUSdbs_inphase_data_fUS' num2str(fUS) 'kHz_tauES' num2str(tauES) 'us.mat'],'EUSdbs_inphase_data');

%%
function varargout = funPostProc(fileName,T_trans,T_max,varargin) 
if (nargin == 4) && strcmp(varargin{1},'spikes')
load (fileName,'data_spikes');
data = data_spikes;
elseif (nargin == 3)
load(fileName,'data');
else 
error('Incorrect number of input or input values')
end
% Specifications
NNames = {'RS','RSFVX','RSIzh','FS','FSFVX','FSIzh',...
    'LTS','LTSCAX','ThTC','ThRE','ThRT','STN','GPi',...
    'GPe','StrMSN'}; 
NNumbers = [20,0,0,5,0,0,5,0,0,0,20,20,20,20,20];

% Spectral power
dt1=0.01*10^-3;
params.Fs = 1/dt1; %Hz
params.fpass = [1 100];
params.tapers = [3 5];
params.trialave = 1;

GPi_APs = cell2struct(cellfun(@(X) X(X>=T_trans&X<=T_max),cellfun(@(X) data.time(~~X),mat2cell(data.GPi_V_spikes,size(data.GPi_V_spikes,1),ones(size(data.GPi_V_spikes,2),1)),'UniformOutput',false),'uniformoutput',false),'times');
STN_APs = cell2struct(cellfun(@(X) X(X>=T_trans&X<=T_max),cellfun(@(X) data.time(~~X),mat2cell(data.STN_V_spikes,size(data.STN_V_spikes,1),ones(size(data.STN_V_spikes,2),1)),'UniformOutput',false),'uniformoutput',false),'times');
GPe_APs = cell2struct(cellfun(@(X) X(X>=T_trans&X<=T_max),cellfun(@(X) data.time(~~X),mat2cell(data.GPe_V_spikes,size(data.GPe_V_spikes,1),ones(size(data.GPe_V_spikes,2),1)),'UniformOutput',false),'uniformoutput',false),'times');

[gpi_alpha_beta_area, gpi_S, gpi_f]=spectAnalysis.make_Spectrum(GPi_APs,params);
[stn_alpha_beta_area, stn_S, stn_f]=spectAnalysis.make_Spectrum(STN_APs,params);
[gpe_alpha_beta_area, gpe_S, gpe_f]=spectAnalysis.make_Spectrum(GPe_APs,params);

specP = {gpi_alpha_beta_area gpi_S gpi_f;
         stn_alpha_beta_area stn_S stn_f;
         gpe_alpha_beta_area gpe_S gpe_f};
% Spectrogram
movingwin = [1,0.1];            % [s] - Moving window (first: width of movingwin, second: step size of window)
[gpi_S_gram,gpi_t_gram,gpi_f_gram]=mtspecgrampt(GPi_APs,movingwin,params);
[stn_S_gram,stn_t_gram,stn_f_gram]=mtspecgrampt(STN_APs,movingwin,params);
[gpe_S_gram,gpe_t_gram,gpe_f_gram]=mtspecgrampt(GPe_APs,movingwin,params);

specGram = {gpi_S_gram,gpi_t_gram,gpi_f_gram;
            stn_S_gram,stn_t_gram,stn_f_gram;
            gpe_S_gram,gpe_t_gram,gpe_f_gram};

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

varargout = {specP,specGram,FiringRate};
end


