clear all; clc;
% --------------------------MATLAB PATH------------------------------------
Chronuxpath = '';                                            % Chronux path
DynasimPath = '';                                            % Dynasim path   

% Spectral analysis (as in Kumaravelu et al. (2016) modelDB code)
addpath(genpath(Chronuxpath));
% Dynasim toolbox (Sherfey et al. (2018))
addpath(genpath(DynasimPath));

% plot 1: Spectograms, Firing rates and power spectral density
% plot 2: Rastergrams
% plot 3: EUS Spectograms, Firing rates and PSD
% plot 4: EUS in phase PSD, Firing rate
plotN = 1;

if (plotN == 1) || isnan(plotN)
tauES = 300;                % Electrical pulse width [us]
Ius = 500;                  % Ultrasonic intensity [W/m^2]
fUS = 700;                  % Ultrasonic frequency [kHz]
% Electrical DBS
load(['cdata/dbs_data_tauES' num2str(tauES) 'us.mat']);
load(['cdata/USdbs_data_Ius' num2str(Ius) '_fUS' num2str(fUS) 'kHz.mat']);
load('cdata/PDon_data.mat');
load(['cdata/US_data_fUS' num2str(fUS) 'kHz.mat']);
USintens = (10:10:200);
f = (10:10:160);         % Hz

% 1. Spectral power by ES/US pulsed DBS
ESdbs_gpi_alpha_beta_P = cellfun(@(X) X{1,1},dbs_data.specP);
ESdbs_stn_alpha_beta_P = cellfun(@(X) X{2,1},dbs_data.specP);
ESdbs_gpe_alpha_beta_P = cellfun(@(X) X{3,1},dbs_data.specP);
USdbs_gpi_alpha_beta_P = cellfun(@(X) X{1,1},USdbs_data.specP);
USdbs_stn_alpha_beta_P = cellfun(@(X) X{2,1},USdbs_data.specP);
USdbs_gpe_alpha_beta_P = cellfun(@(X) X{3,1},USdbs_data.specP);
% 2. Continuous wave US dbs
US_gpi_alpha_beta_P = cellfun(@(X) X{1,1},US_data.specP);
US_stn_alpha_beta_P = cellfun(@(X) X{2,1},US_data.specP);
US_gpe_alpha_beta_P = cellfun(@(X) X{3,1},US_data.specP);
% 3. No stim PD
alpha_beta_P_PDon = cell2mat(PDon_data.specP(:,1));

% 4. Firing rate calcs
ESdbs_FR_STN = cellfun(@(X) X(5).MeanFR,dbs_data.FiringRate);
ESdbs_FR_GPi = cellfun(@(X) X(6).MeanFR,dbs_data.FiringRate);
ESdbs_FR_GPe = cellfun(@(X) X(7).MeanFR,dbs_data.FiringRate);
ESdbs_stdFR_STN = cellfun(@(X) X(5).stdFR,dbs_data.FiringRate);
ESdbs_stdFR_GPi = cellfun(@(X) X(6).stdFR, dbs_data.FiringRate);
ESdbs_stdFR_GPe = cellfun(@(X) X(7).stdFR,dbs_data.FiringRate);

USdbs_FR_STN = cellfun(@(X) X(5).MeanFR,USdbs_data.FiringRate);
USdbs_FR_GPi = cellfun(@(X) X(6).MeanFR,USdbs_data.FiringRate);
USdbs_FR_GPe = cellfun(@(X) X(7).MeanFR,USdbs_data.FiringRate);
USdbs_stdFR_STN = cellfun(@(X) X(5).stdFR,USdbs_data.FiringRate);
USdbs_stdFR_GPi = cellfun(@(X) X(6).stdFR,USdbs_data.FiringRate);
USdbs_stdFR_GPe = cellfun(@(X) X(7).stdFR,USdbs_data.FiringRate);

PDon_FR_STN = PDon_data.FiringRate(5).MeanFR;
PDon_FR_GPi = PDon_data.FiringRate(6).MeanFR;
PDon_FR_GPe = PDon_data.FiringRate(7).MeanFR;
PDon_stdFR_STN = PDon_data.FiringRate(5).stdFR;
PDon_stdFR_GPi = PDon_data.FiringRate(6).stdFR;
PDon_stdFR_GPe = PDon_data.FiringRate(7).stdFR;

US_FR_STN = cellfun(@(X) X(5).MeanFR,US_data.FiringRate);
US_FR_GPi = cellfun(@(X) X(6).MeanFR,US_data.FiringRate);
US_FR_GPe = cellfun(@(X) X(7).MeanFR,US_data.FiringRate);
US_stdFR_STN = cellfun(@(X) X(5).stdFR,US_data.FiringRate);
US_stdFR_GPi = cellfun(@(X) X(6).stdFR, US_data.FiringRate);
US_stdFR_GPe = cellfun(@(X) X(7).stdFR,US_data.FiringRate);

col = linspecer(3);
figure('units','normalized','position',[0 0 1 1]); set(gcf,'color','w');
subplot(4,3,1);
hold on;
p1 = plot(f,ESdbs_stn_alpha_beta_P,'linewidth',2,'linestyle','-','color',col(1,:)); 
q1 = plot(f,USdbs_stn_alpha_beta_P,'linewidth',2,'linestyle','--','color',col(1,:)); 
m1 = plot(f,alpha_beta_P_PDon(2).*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(1,:));
p2 = plot(f,ESdbs_gpi_alpha_beta_P,'linewidth',2,'linestyle','-','color',col(2,:)); 
q2 = plot(f,USdbs_gpi_alpha_beta_P,'linewidth',2,'linestyle','--','color',col(2,:)); 
m2 = plot(f,alpha_beta_P_PDon(1).*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(2,:));
p3 = plot(f,ESdbs_gpe_alpha_beta_P,'linewidth',2,'linestyle','-','color',col(3,:)); box off;
q3 = plot(f,USdbs_gpe_alpha_beta_P,'linewidth',2,'linestyle','--','color',col(3,:)); box off;
m3 = plot(f,alpha_beta_P_PDon(3).*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(3,:));
ylabel('\alpha-\beta SE [-]');
g1 = plot(nan,nan,':','color','k','linewidth',2); % No stim
g2 = plot(nan,nan,'-','color','k','linewidth',2);   % ES
g3 = plot(nan,nan,'--','color','k','linewidth',2); % US
l1 = legend([p1,p2,p3],{'STN','GPi','GPe'},'orientation','horizontal'); legend boxoff;
ah1=axes('position',get(gca,'position'),'visible','off');
l2 = legend(ah1,[g1,g2,g3],{'PD','ES','US'},'orientation','horizontal','position',[0.7912 0.8127 0.1055 0.0417]); legend boxoff;
Nmb1 = text(0.02,1,'(a)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,2);
hold on;
o1 = plot(USintens,US_stn_alpha_beta_P,'linewidth',2,'linestyle','--','color',col(1,:)); 
p1 = plot(USintens,alpha_beta_P_PDon(2).*ones(length((USintens)),1),'linewidth',2,'linestyle',':','color',col(1,:));
o2 = plot(USintens,US_gpi_alpha_beta_P,'linewidth',2,'linestyle','--','color',col(2,:)); 
p2 = plot(USintens,alpha_beta_P_PDon(1).*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(2,:));
o3 = plot(USintens,US_gpe_alpha_beta_P,'linewidth',2,'linestyle','--','color',col(3,:)); box off;
p3 = plot(USintens,alpha_beta_P_PDon(3).*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(3,:));
Nmb3 = text(0.02,1,'(b)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,4);
hold on;
u1 = plot(f,ESdbs_FR_STN,'linewidth',2,'linestyle','-','color',col(1,:)); 
w1 = plot(f,USdbs_FR_STN,'linewidth',2,'linestyle','--','color',col(1,:)); 
v1 = plot(f,PDon_FR_STN.*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(1,:));
u2 = plot(f,ESdbs_FR_GPi,'linewidth',2,'linestyle','-','color',col(2,:));
w2 = plot(f,USdbs_FR_GPi,'linewidth',2,'linestyle','--','color',col(2,:));
v2 = plot(f,PDon_FR_GPi.*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(2,:));
u3 = plot(f,ESdbs_FR_GPe,'linewidth',2,'linestyle','-','color',col(3,:)); 
w3 = plot(f,USdbs_FR_GPe,'linewidth',2,'linestyle','--','color',col(3,:)); 
v3 = plot(f,PDon_FR_GPe.*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(3,:));
a = plot(f,f,'linewidth',1,'linestyle',':','color','k');
ylabel('\muFR [Hz]');
Nmb4 = text(0.02,1,'(d)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,5);
hold on;
x1 = plot(USintens,US_FR_STN,'linewidth',2,'linestyle','--','color',col(1,:)); 
y1 = plot(USintens,PDon_FR_STN.*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(1,:));
x2 = plot(USintens,US_FR_GPi,'linewidth',2,'linestyle','--','color',col(2,:));
y2 = plot(USintens,PDon_FR_GPi.*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(2,:));
x3 = plot(USintens,US_FR_GPe,'linewidth',2,'linestyle','--','color',col(3,:)); 
y3 = plot(USintens,PDon_FR_GPe.*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(3,:));
Nmb5 = text(0.02,1,'(e)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,7);
hold on;
b1 = plot(f,ESdbs_stdFR_STN,'linewidth',2,'linestyle','-','color',col(1,:)); 
d1 = plot(f,USdbs_stdFR_STN,'linewidth',2,'linestyle','--','color',col(1,:)); 
c1 = plot(f,PDon_stdFR_STN.*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(1,:));
b2 = plot(f,ESdbs_stdFR_GPi,'linewidth',2,'linestyle','-','color',col(2,:));
d2 = plot(f,USdbs_stdFR_GPi,'linewidth',2,'linestyle','--','color',col(2,:));
c2 = plot(f,PDon_stdFR_GPi.*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(2,:));
b3 = plot(f,ESdbs_stdFR_GPe,'linewidth',2,'linestyle','-','color',col(3,:)); 
d3 = plot(f,USdbs_stdFR_GPe,'linewidth',2,'linestyle','--','color',col(3,:)); 
c3 = plot(f,PDon_stdFR_GPe.*ones(length(f),1),'linewidth',2,'linestyle',':','color',col(3,:));
ylabel('\sigmaFR [Hz]'); xlabel('f_{DBS} [Hz]');
Nmb6 = text(0.02,1,'(g)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,8);
hold on;
r1 = plot(USintens,US_stdFR_STN,'linewidth',2,'linestyle','--','color',col(1,:)); 
t1 = plot(USintens,PDon_stdFR_STN.*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(1,:));
r2 = plot(USintens,US_stdFR_GPi,'linewidth',2,'linestyle','--','color',col(2,:));
t2 = plot(USintens,PDon_stdFR_GPi.*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(2,:));
r3 = plot(USintens,US_stdFR_GPe,'linewidth',2,'linestyle','--','color',col(3,:)); 
t3 = plot(USintens,PDon_stdFR_GPe.*ones(length(USintens),1),'linewidth',2,'linestyle',':','color',col(3,:));
xlabel('I_{US} [W/m^2]');
Nmb7 = text(0.02,1,'(h)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,3);
hold on;
plot(PDon_data.specP{2,3},PDon_data.specP{2,2},'linewidth',1.5,'linestyle',':','color',col(1,:)); box off;
plot(PDon_data.specP{1,3},PDon_data.specP{1,2},'linewidth',1.5,'linestyle',':','color',col(2,:)); box off;
plot(PDon_data.specP{3,3},PDon_data.specP{3,2},'linewidth',1.5,'linestyle',':','color',col(3,:)); box off;

plot(dbs_data.specP{16}{2,3},dbs_data.specP{16}{2,2},'linewidth',1.5,'linestyle','-','color',col(1,:)); box off;
plot(dbs_data.specP{16}{1,3},dbs_data.specP{16}{1,2},'linewidth',1.5,'linestyle','-','color',col(2,:)); box off;
plot(dbs_data.specP{16}{3,3},dbs_data.specP{16}{3,2},'linewidth',1.5,'linestyle','-','color',col(3,:)); box off;
ylabel('PSD [s]');
title('160 Hz DBS versus PD');
Nmb8 = text(0.02,1,'(c)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(4,3,6);
hold on;
plot(PDon_data.specP{2,3},PDon_data.specP{2,2},'linewidth',1.5,'linestyle',':','color',col(1,:)); box off;
plot(PDon_data.specP{1,3},PDon_data.specP{1,2},'linewidth',1.5,'linestyle',':','color',col(2,:)); box off;
plot(PDon_data.specP{3,3},PDon_data.specP{3,2},'linewidth',1.5,'linestyle',':','color',col(3,:)); box off;

plot(US_data.specP{20}{2,3},US_data.specP{20}{2,2},'linewidth',1.5,'linestyle','--','color',col(1,:)); box off;
plot(US_data.specP{20}{1,3},US_data.specP{20}{1,2},'linewidth',1.5,'linestyle','--','color',col(2,:)); box off;
plot(US_data.specP{20}{3,3},US_data.specP{20}{3,2},'linewidth',1.5,'linestyle','--','color',col(3,:)); box off;
ylim([0 150]);
ylabel('PSD [s]'); xlabel('f [Hz]');
title('200 W/m^2 US versus PD');
Nmb9 = text(0.02,1,'(f)','fontsize',18,'fontweight','bold','units','normalized');
hold off;


% Spectograms
PSDdbs160 = dbs_data.specGram{16};
PSDus200 = US_data.specGram{20};
PSDpdon = PDon_data.specGram;
caxisLIMS = [floor(min(min(PSDpdon{1,1}(:,PSDpdon{1,3}>=8&PSDpdon{1,3}<=50)))),ceil(max(max(PSDpdon{1,1}(:,PSDpdon{1,3}>=8&PSDpdon{1,3}<=50))))];

subplot(4,3,10);
pcolor(PSDdbs160{1,2},PSDdbs160{1,3}(PSDdbs160{1,3}>=8&PSDdbs160{1,3}<=50),PSDdbs160{1,1}(:,PSDdbs160{1,3}>=8&PSDdbs160{1,3}<=50)'); shading interp;
xlabel('t [s]'); ylabel('f [Hz]'); title('GPi 160 Hz DBS');
caxis(caxisLIMS);
Nmb10 = text(0.02,1,'(j)','fontsize',18,'fontweight','bold','units','normalized');


subplot(4,3,11);
pcolor(PSDus200{1,2},PSDus200{1,3}(PSDus200{1,3}>=8&PSDus200{1,3}<=50),PSDus200{1,1}(:,PSDus200{1,3}>=8&PSDus200{1,3}<=50)'); shading interp;
xlabel('t [s]'); title('GPi 200 W/m^2 US');
caxis(caxisLIMS);
Nmb10 = text(0.02,1,'(k)','fontsize',18,'fontweight','bold','units','normalized');


subplot(4,3,9);
pcolor(PSDpdon{2,2},PSDpdon{2,3}(PSDpdon{2,3}>=8&PSDpdon{2,3}<=50),PSDpdon{2,1}(:,PSDpdon{2,3}>=8&PSDpdon{2,3}<=50)'); shading interp;
xlabel('t [s]'); ylabel('f [Hz]'); title('STN PD');
colorbar;
Nmb11 = text(0.02,1,'(i)','fontsize',18,'fontweight','bold','units','normalized');

subplot(4,3,12);
pcolor(PSDpdon{1,2},PSDpdon{1,3}(PSDpdon{1,3}>=8&PSDpdon{1,3}<=50),PSDpdon{1,1}(:,PSDpdon{1,3}>=8&PSDpdon{1,3}<=50)'); shading interp;
xlabel('t [s]'); title('GPi PD');
c = colorbar; ylabel(c,'PSD [s]');  caxis(caxisLIMS);
Nmb12 = text(0.02,1,'(l)','fontsize',18,'fontweight','bold','units','normalized');

set(findobj('type','axes'),'fontsize',18);


 
% % Spectral density
% figure; set(gcf,'color','w'); 
% subplot(1,3,1); plot(stn_f,stn_S); xlabel('f [Hz]'); ylabel('Power spectral density'); title('STN'); box off;
% subplot(1,3,2); plot(gpe_f,gpe_S); xlabel('f [Hz]'); title('GPe'); box off;
% subplot(1,3,3); plot(gpi_f,gpi_S); xlabel('f [Hz]'); title('GPi'); box off;
% 
% % Spectogram
% figure; set(gcf,'color','w');
% subplot(1,3,1); pcolor(stn_t_gram,stn_f_gram(stn_f_gram<=50&stn_f_gram>=8),stn_S_gram(:,stn_f_gram<=50&stn_f_gram>=8)'); xlabel('t [s]'); ylabel('Power spectral density'); shading interp
% subplot(1,3,2); pcolor(gpe_t_gram,gpe_f_gram(gpe_f_gram<=50&gpe_f_gram>=8),gpe_S_gram(:,gpe_f_gram<=50&gpe_f_gram>=8)'); xlabel('t [s]'); shading interp
% subplot(1,3,3); pcolor(gpi_t_gram,gpi_f_gram(gpi_f_gram<=50&gpi_f_gram>=8),gpi_S_gram(:,gpi_f_gram<=50&gpi_f_gram>=8)'); xlabel('t [s]'); shading interp
% 
end

if (plotN == 2) || isnan(plotN)
tauES = 300;                % Electrical pulse width [us]
fUS = 700;                  % Ultrasonic frequency [kHz]

% Rastergrams
 
PDondata = load('PD/data_PDon.mat');

figure; set(gcf,'color','w');
subplot(1,3,1);
hold on;
dsPlot(PDondata.data,'plot_type','rastergram','xlim',[1 11],'lock_gca',1);
title('PD');
xlabel('Time [s]');
hold off;

clear PDondata; % 3 GB in mem, load one by one
DBS160data = load(['DBS ES/tauES_' num2str(tauES) 'us/data_PDon_dbs160.mat']);
subplot(1,3,2);
hold on;
dsPlot(DBS160data.data,'plot_type','rastergram','xlim',[1 11],'lock_gca',1);
xlabel('Time [s]');
title('DBS 160 Hz');
hold off;

clear DBS160data;
US200data = load(['CW US/fUS' num2str(fUS) 'kHz/data_PDon_USi200.mat']);
US200data.data.labels = horzcat({US200data.data.labels{1:21},US200data.data.labels{end},US200data.data.labels{22:62},'STN_V_spikes',US200data.data.labels{63:end-1}}); % Reshuffle plot order
US200data.data = orderfields(US200data.data,[(1:23) 72 (24:64) 71 (65:70)]);
subplot(1,3,3);
hold on;
dsPlot(US200data.data,'plot_type','rastergram','xlim',[1 11],'lock_gca',1);
title('US 200 W/m^2');
xlabel('Time [s]');
hold off;
clear US200data;
set(findobj('type','axes'),'fontsize',18);
end


if (plotN == 3) || isnan(plotN)
tauES = 100;                % Electrical pulse width [us]
Ius = 500;                  % Ultrasonic intensity [W/m^2]
fUS = 700;                  % Ultrasonic frequency [kHz]

% Electro-ultrasonic stimulation
load(['cdata/EUSdbs_TempOutphase_data_fUS' num2str(fUS) 'kHz_tauES' num2str(tauES) 'us.mat']);
load(['cdata/dbs_data_tauES' num2str(tauES) 'us.mat']);
load(['cdata/USdbs_data_Ius' num2str(Ius) '_fUS' num2str(fUS) 'kHz.mat']);
fEUS = 0.5*(10:10:160);         % Hz
fDBS = (10:10:160);             % Hz


% 1. Spectral power by ES/US/EUS pulsed DBS
ESdbs_gpi_alpha_beta_P = cellfun(@(X) X{1,1},dbs_data.specP);
ESdbs_stn_alpha_beta_P = cellfun(@(X) X{2,1},dbs_data.specP);
ESdbs_gpe_alpha_beta_P = cellfun(@(X) X{3,1},dbs_data.specP);
USdbs_gpi_alpha_beta_P = cellfun(@(X) X{1,1},USdbs_data.specP);
USdbs_stn_alpha_beta_P = cellfun(@(X) X{2,1},USdbs_data.specP);
USdbs_gpe_alpha_beta_P = cellfun(@(X) X{3,1},USdbs_data.specP);

EUSdbs_TempOutphase_gpi_alpha_beta_P = cellfun(@(X) X{1,1},EUSdbs_TempOutphase_data.specP);
EUSdbs_TempOutphase_stn_alpha_beta_P = cellfun(@(X) X{2,1},EUSdbs_TempOutphase_data.specP);
EUSdbs_TempOutphase_gpe_alpha_beta_P = cellfun(@(X) X{3,1},EUSdbs_TempOutphase_data.specP);


% 4. Firing rate calcs
ESdbs_FR_STN = cellfun(@(X) X(5).MeanFR,dbs_data.FiringRate);
ESdbs_FR_GPi = cellfun(@(X) X(6).MeanFR,dbs_data.FiringRate);
ESdbs_FR_GPe = cellfun(@(X) X(7).MeanFR,dbs_data.FiringRate);
ESdbs_stdFR_STN = cellfun(@(X) X(5).stdFR,dbs_data.FiringRate);
ESdbs_stdFR_GPi = cellfun(@(X) X(6).stdFR, dbs_data.FiringRate);
ESdbs_stdFR_GPe = cellfun(@(X) X(7).stdFR,dbs_data.FiringRate);

USdbs_FR_STN = cellfun(@(X) X(5).MeanFR,USdbs_data.FiringRate);
USdbs_FR_GPi = cellfun(@(X) X(6).MeanFR,USdbs_data.FiringRate);
USdbs_FR_GPe = cellfun(@(X) X(7).MeanFR,USdbs_data.FiringRate);
USdbs_stdFR_STN = cellfun(@(X) X(5).stdFR,USdbs_data.FiringRate);
USdbs_stdFR_GPi = cellfun(@(X) X(6).stdFR,USdbs_data.FiringRate);
USdbs_stdFR_GPe = cellfun(@(X) X(7).stdFR,USdbs_data.FiringRate);

EUSdbs_TempOutphase_FR_STN = cellfun(@(X) X(5).MeanFR,EUSdbs_TempOutphase_data.FiringRate);
EUSdbs_TempOutphase_FR_GPi = cellfun(@(X) X(6).MeanFR,EUSdbs_TempOutphase_data.FiringRate);
EUSdbs_TempOutphase_FR_GPe = cellfun(@(X) X(7).MeanFR,EUSdbs_TempOutphase_data.FiringRate);
EUSdbs_TempOutphase_stdFR_STN = cellfun(@(X) X(5).stdFR,EUSdbs_TempOutphase_data.FiringRate);
EUSdbs_TempOutphase_stdFR_GPi = cellfun(@(X) X(6).stdFR, EUSdbs_TempOutphase_data.FiringRate);
EUSdbs_TempOutphase_stdFR_GPe = cellfun(@(X) X(7).stdFR,EUSdbs_TempOutphase_data.FiringRate);

col = linspecer(3);
figure('units','normalized','position',[0 0 1 1]); set(gcf,'color','w');
subplot(1,3,1);
hold on;
p1 = plot(fEUS,EUSdbs_TempOutphase_stn_alpha_beta_P,'linewidth',2,'linestyle',':','marker','*','color',col(1,:)); 
p2 = plot(fEUS,EUSdbs_TempOutphase_gpi_alpha_beta_P,'linewidth',2,'linestyle',':','marker','*','color',col(2,:)); 
p3 = plot(fEUS,EUSdbs_TempOutphase_gpe_alpha_beta_P,'linewidth',2,'linestyle',':','marker','*','color',col(3,:)); box off;
axes1p = gca; ylimsp = ylim;
ylabel('\alpha-\beta SE [-]');
l1 = legend([p1,p2,p3],{'STN','GPi','GPe'},'orientation','horizontal'); legend boxoff;
ah1=axes('position',get(gca,'position'),'visible','off');
Nmb1 = text(0.02,1,'(a)','fontsize',18,'fontweight','bold','units','normalized');
hold off;
axes2p = axes('Position',axes1p.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','ylim',ylimsp);
hold on;
p1x = plot(fDBS,ESdbs_stn_alpha_beta_P,'Parent',axes2p,'linewidth',2,'linestyle','-','color',col(1,:));
p2x = plot(fDBS,ESdbs_gpi_alpha_beta_P,'Parent',axes2p,'linewidth',2,'linestyle','-','color',col(2,:));
p3x = plot(fDBS,ESdbs_gpe_alpha_beta_P,'Parent',axes2p,'linewidth',2,'linestyle','-','color',col(3,:));

p1y = plot(fDBS,USdbs_stn_alpha_beta_P,'Parent',axes2p,'linewidth',2,'linestyle','--','color',col(1,:));
p2y = plot(fDBS,USdbs_gpi_alpha_beta_P,'Parent',axes2p,'linewidth',2,'linestyle','--','color',col(2,:));
p3y = plot(fDBS,USdbs_gpe_alpha_beta_P,'Parent',axes2p,'linewidth',2,'linestyle','--','color',col(3,:));
ylim(ylimsp);
hold off;

subplot(1,3,2);
hold on;
u1 = plot(fEUS,EUSdbs_TempOutphase_FR_STN,'linewidth',2,'linestyle',':','marker','*','color',col(1,:)); 
u2 = plot(fEUS,EUSdbs_TempOutphase_FR_GPi,'linewidth',2,'linestyle',':','marker','*','color',col(2,:));
u3 = plot(fEUS,EUSdbs_TempOutphase_FR_GPe,'linewidth',2,'linestyle',':','marker','*','color',col(3,:)); 
a = plot(fEUS,2*fEUS,'linewidth',1,'linestyle',':','color','k');
ylabel('\muFR [Hz]'); xlabel('f_{EUS} [Hz]');
Nmb4 = text(0.02,1,'(b)','fontsize',18,'fontweight','bold','units','normalized');
hold off;
axes1u = gca; ylimsu = ylim;
axes2u = axes('Position',axes1u.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','ylim',ylimsu);
hold on;
u1x = plot(fDBS,ESdbs_FR_STN,'Parent',axes2u,'linewidth',2,'linestyle','-','color',col(1,:));
u2x = plot(fDBS,ESdbs_FR_GPi,'Parent',axes2u,'linewidth',2,'linestyle','-','color',col(2,:));
u3x = plot(fDBS,ESdbs_FR_GPe,'Parent',axes2u,'linewidth',2,'linestyle','-','color',col(3,:));

u1y = plot(fDBS,USdbs_FR_STN,'Parent',axes2u,'linewidth',2,'linestyle','--','color',col(1,:));
u2y = plot(fDBS,USdbs_FR_GPi,'Parent',axes2u,'linewidth',2,'linestyle','--','color',col(2,:));
u3y = plot(fDBS,USdbs_FR_GPe,'Parent',axes2u,'linewidth',2,'linestyle','--','color',col(3,:));
ylim(ylimsu);
xlabel('f_{DBS} [Hz]');
hold off;


subplot(1,3,3);
hold on;
b1 = plot(fEUS,EUSdbs_TempOutphase_stdFR_STN,'linewidth',2,'linestyle',':','marker','*','color',col(1,:)); 
b2 = plot(fEUS,EUSdbs_TempOutphase_stdFR_GPi,'linewidth',2,'linestyle',':','marker','*','color',col(2,:));
b3 = plot(fEUS,EUSdbs_TempOutphase_stdFR_GPe,'linewidth',2,'linestyle',':','marker','*','color',col(3,:)); 
ylabel('\sigmaFR [Hz]');
Nmb6 = text(0.02,1,'(c)','fontsize',18,'fontweight','bold','units','normalized');
ylimsb = ylim;
hold off;
axes1b = gca;
axes2b = axes('Position',axes1b.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','ylim',ylimsb);
hold on;
b1x = plot(fDBS,ESdbs_stdFR_STN,'Parent',axes2b,'linewidth',2,'linestyle','-','color',col(1,:));
b2x = plot(fDBS,ESdbs_stdFR_GPi,'Parent',axes2b,'linewidth',2,'linestyle','-','color',col(2,:));
b3x = plot(fDBS,ESdbs_stdFR_GPe,'Parent',axes2b,'linewidth',2,'linestyle','-','color',col(3,:));

b1y = plot(fDBS,USdbs_stdFR_STN,'Parent',axes2b,'linewidth',2,'linestyle','--','color',col(1,:));
b2y = plot(fDBS,USdbs_stdFR_GPi,'Parent',axes2b,'linewidth',2,'linestyle','--','color',col(2,:));
b3y = plot(fDBS,USdbs_stdFR_GPe,'Parent',axes2b,'linewidth',2,'linestyle','--','color',col(3,:));

n1 = plot(nan,nan,'color','k','linestyle',':','linewidth',2,'marker','*');
n2 = plot(nan,nan,'color','k','linestyle','-','linewidth',2);
n3 = plot(nan,nan,'color','k','linestyle','--','linewidth',2);
ylim(ylimsb);
l2 = legend([n1,n2,n3],{'cEUS','ES','US'},'orientation','horizontal'); legend boxoff;
hold off;
axes2b.Position = axes1b.Position;
axes2u.Position = axes1u.Position;
axes2p.Position = axes1p.Position;
axes2b.YTickLabel = {};
axes2u.YTickLabel = {};
axes2p.YTickLabel = {};
set(findobj('type','axes'),'fontsize',17);

l1.Position = [0.2826 0.9641 0.1757 0.0395];


end


if (plotN == 4) || isnan(plotN)
tauES = 100;                % Electrical pulse width [us]
fUS = 700;                  % Ultrasonic frequency [kHz]

% Electro-ultrasonic stimulation - in phase
load(['cdata/EUSdbs_inphase_data_fUS' num2str(fUS) 'kHz_tauES' num2str(tauES) 'us.mat']);
InPhaseUSiRange = linspace(0,500,11);   % W/m^2
InPhaseESiRange = linspace(0,2,11);     % A/m^2


% 1. Spectral power by EUS pulsed DBS in phase
EUSdbs_inphase_gpi_alpha_beta_P = cellfun(@(X) X{1,1},EUSdbs_inphase_data.specP);
EUSdbs_inphase_stn_alpha_beta_P = cellfun(@(X) X{2,1},EUSdbs_inphase_data.specP);
EUSdbs_inphase_gpe_alpha_beta_P = cellfun(@(X) X{3,1},EUSdbs_inphase_data.specP);


% 4. Firing rate calcs
EUSdbs_inphase_FR_STN = cellfun(@(X) X(5).MeanFR,EUSdbs_inphase_data.FiringRate);
EUSdbs_inphase_FR_GPi = cellfun(@(X) X(6).MeanFR,EUSdbs_inphase_data.FiringRate);
EUSdbs_inphase_FR_GPe = cellfun(@(X) X(7).MeanFR,EUSdbs_inphase_data.FiringRate);
EUSdbs_inphase_stdFR_STN = cellfun(@(X) X(5).stdFR,EUSdbs_inphase_data.FiringRate);
EUSdbs_inphase_stdFR_GPi = cellfun(@(X) X(6).stdFR, EUSdbs_inphase_data.FiringRate);
EUSdbs_inphase_stdFR_GPe = cellfun(@(X) X(7).stdFR,EUSdbs_inphase_data.FiringRate);

figure('units','normalized','position',[0 0 1 1]); set(gcf,'color','w');
subplot(3,3,1);
hold on;
p1 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stn_alpha_beta_P); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stn_alpha_beta_P,'showtext','on','linecolor','k');
c1 = colorbar('orientation','horizontal','location','northoutside'); ylabel(c1,'\alpha-\beta SE [-]');
Nmb1 = text(4.4010e-04,1.1015,'(a)','fontsize',18,'fontweight','bold','units','normalized');
ylabel('I_{US} [W/m^2]');
hold off;
subplot(3,3,4);
hold on;
p2 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_gpi_alpha_beta_P); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_gpi_alpha_beta_P,'showtext','on','linecolor','k');
c2 = colorbar('orientation','horizontal','location','northoutside'); %ylabel(c2,'GPi \alpha-\beta SE [-]');
Nmb4 = text(4.4010e-04,1.1015,'(d)','fontsize',18,'fontweight','bold','units','normalized');
ylabel('I_{US} [W/m^2]'); 
hold off;
subplot(3,3,7);
hold on;
p3 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_gpe_alpha_beta_P); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_gpe_alpha_beta_P,'showtext','on','linecolor','k');
c3 = colorbar('orientation','horizontal','location','northoutside'); %ylabel(c3,'GPe \alpha-\beta SE [-]');
ylabel('I_{US} [W/m^2]'); 
xlabel('I_{ES} [\muA/cm^2]'); 
Nmb7 = text(4.4010e-04,1.1015,'(g)','fontsize',18,'fontweight','bold','units','normalized');
hold off;

subplot(3,3,2);
hold on;
u1 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_FR_STN); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_FR_STN,'showtext','on','linecolor','k');
d1 = colorbar('orientation','horizontal','location','northoutside'); ylabel(d1,'\muFR [Hz]');
Nmb2 = text(4.4010e-04,1.1015,'(b)','fontsize',18,'fontweight','bold','units','normalized');
hold off;
subplot(3,3,5);
hold on;
u2 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_FR_GPi); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_FR_GPi,'showtext','on','linecolor','k');
d2 = colorbar('orientation','horizontal','location','northoutside'); %ylabel(d2,'GPi \muFR [Hz]');
Nmb5 = text(4.4010e-04,1.1015,'(e)','fontsize',18,'fontweight','bold','units','normalized');
hold off;
subplot(3,3,8);
hold on;
u3 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_FR_GPe); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_FR_GPe,'showtext','on','linecolor','k');
d3 = colorbar('orientation','horizontal','location','northoutside'); %ylabel(d3,'GPe \muFR [Hz]');
Nmb8 = text(4.4010e-04,1.1015,'(h)','fontsize',18,'fontweight','bold','units','normalized');
xlabel('I_{ES} [\muA/cm^2]');
hold off;


subplot(3,3,3);
hold on;
b1 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stdFR_STN); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stdFR_STN,'showtext','on','linecolor','k');
Nmb3 = text(4.4010e-04,1.1015,'(c)','fontsize',18,'fontweight','bold','units','normalized');
e1 = colorbar('orientation','horizontal','location','northoutside'); ylabel(e1,'\sigmaFR [Hz]');
hold off;
subplot(3,3,6);
hold on;
b2 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stdFR_GPi); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stdFR_GPi,'showtext','on','linecolor','k');
Nmb6 = text(4.4010e-04,1.1015,'(f)','fontsize',18,'fontweight','bold','units','normalized');
e2 = colorbar('orientation','horizontal','location','northoutside'); %ylabel(e2,'GPi \sigmaFR [Hz]');
hold off;
subplot(3,3,9);
hold on;
b3 = pcolor(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stdFR_GPe); shading interp;
contour(100*InPhaseESiRange,InPhaseUSiRange,EUSdbs_inphase_stdFR_GPe,'showtext','on','linecolor','k');
Nmb9 = text(4.4010e-04,1.1015,'(i)','fontsize',18,'fontweight','bold','units','normalized');
e3 = colorbar('orientation','horizontal','location','northoutside'); %ylabel(e3,'GPe \sigmaFR [Hz]');
xlabel('I_{ES} [\muA/cm^2]');
hold off;
STNtxt = text(-3.0244,4.5101,'STN','fontsize',18,'fontweight','bold','units','normalized');
GPitxt = text(-3.0244, 2.5790,'GPi','fontsize',18,'fontweight','bold','units','normalized');
GPetxt = text(-3.0244, 0.6273,'GPe','fontsize',18,'fontweight','bold','units','normalized');
set(findobj('type','axes'),'fontsize',18);
end
