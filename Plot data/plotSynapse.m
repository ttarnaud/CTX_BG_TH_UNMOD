clear all; clc; close all;
load('PD/data_PDon.mat'); 
ii = 1; % First cell
xlims = [5.9 6.2];      % [s]
tspike_pre_RS = cell(1,size(data.RS_V_spikes,2));
for i = 1:size(data.RS_V_spikes,2)
tspike_pre_RS{:,i} = data.time(~~data.RS_V_spikes(:,i))';
end
tau_r = 0.1e-3; tau_d = 3e-3; tau_del = 1e-3;  % [s]
Anorm = @(tr,td) ((tr/td)^(-tr/(tr-td))-(tr/td)^(-td/(tr-td)))^(-1);
biexp = @(x) Anorm(tau_r,tau_d)*(exp(-double(x>0).*x/tau_d)-exp(-double(x>0).*x/tau_r)).*double(x>0);
RS_s = sum(biexp(data.time-tspike_pre_RS{:,ii}-tau_del),2);

figure; set(gcf,'color','w');
subplot(2,3,1); plot(data.time,data.LTS_RS_AMPAbiRSLTS_F(:,ii).*data.LTS_RS_AMPAbiRSLTS_D1(:,ii).*data.LTS_RS_AMPAbiRSLTS_D2(:,ii).*...
data.LTS_RS_AMPAbiRSLTS_D3(:,ii).*RS_s,'linewidth',1.5); box off; ylabel('Relative synaptic conductance');
subplot(2,3,2); plot(data.time,data.FS_RS_AMPAbiRSFS_F(:,ii).*data.FS_RS_AMPAbiRSFS_D1(:,ii).*data.FS_RS_AMPAbiRSFS_D2(:,ii).*...
data.FS_RS_AMPAbiRSFS_D3(:,ii).*RS_s,'linewidth',1.5); box off;
subplot(2,3,3); plot(data.time,data.StrMSN_StrMSN_GABAAdaStrMSNStrMSN_s(:,ii),'linewidth',1.5); box off;
set(findobj('type','axes'),'xlim',xlims)

data2.RS_V_spikes = data.RS_V_spikes(:,ii); data2.time = data.time; data2.model = data.model;
data3.StrMSN_V_spikes = data.StrMSN_V_spikes(:,ii); data3.time = data.time; data3.model = data.model;
subplot(2,3,4); dsPlot(data2,'plot_type','rastergram','xlim',xlims,'lock_gca',1);
subplot(2,3,5); dsPlot(data2,'plot_type','rastergram','xlim',xlims,'lock_gca',1);
subplot(2,3,6); dsPlot(data3,'plot_type','rastergram','xlim',xlims,'lock_gca',1);
set(findobj('type','axes'),'xlim',xlims)
set(findobj('type','axes'),'fontsize',18)