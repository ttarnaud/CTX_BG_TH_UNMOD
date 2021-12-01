This code implements the cortex-basal ganglia-thalamus model used in Tarnaud et al. (2021) [1].  
The neuronal network is based on Kumaravelu et al. (2016) [2] and on the cortical network 
from Plaksin et al. (2016) [3], Hayut et al. (2011) [4] and Vierling-Claassen et al. (2010) [5] 

The result is a completely biophysical HH-based network of the basal-ganglia-thalamus-cortex loop with ultrasonic coupling
Coupling with ultrasound is included via the SONIC-model (Lemaire et al. 2019 [6]).

1 Matlab scripts
--------------------------------------------------
Install Dynasim: https://github.com/DynaSim/DynaSim [7].
Install chronux: http://chronux.org/  [8]

Dynasim and chronux should be added to the MATLAB path
Simulations were performed on MATLAB/R2017a, with the attached Dynasim version (Dynasim-master)

Compile nakeinterp1 to .mex 

1.1 Main script 
---------------
- run_ES_NN_PDon.m - the main Parkinsonian network
Parameter sweeps performed for [1] can be found in the subfolder /Parameter sweeps

1.2 Membrane dynamics 
---------------
models/-folder contains the membrane dynamics .mech-files

SONIC-STN.mat contains the SONIC table [6] for the subthalamic nucleus for ultrasonic neuromodulation

1.3 Subfunctions
---------------
- linspecer (by Jonathan Lansey [9]: line colours) 
- spectAnalysis (spectral analysis class by Kumaravelu et al. (2016) [2] from modelDB 206232)
- nakeinterp1 (fast interpolation, Bruno Luong). 

Helper functions: 
USStep, SynC, rndize, randESI, fVCa, funtauu, f1Veff and f1rt

2 Plot data folder
-----------------
This folder contains the scripts used to gerenate the plots from the simulation data

- PostProc: will extract firing rate/spectrogram/... data from data_....mat into ..._data.mat (cdata) files. (parameters tauES, fUS and Ius can be altered)
  
2.1 Plot functions
------------------
- plotResults - Plot the neuronal network figures in the paper
- plotSynapse (Fig. 1 - plots typical synaptic currents)

2.2 Subfunctions
------------------
- linspecer (by Jonathan Lansey [9]: line colours) 
- spectAnalysis (spectral analysis class by Kumaravelu et al. (2016) [2] from modelDB 206232)


2.3 Data folders
------------------------------------------------
The following folder-structure is used in the plot-function scripts:

CW - continuous wave 
DBS - pulsed 

US - Ultrasound
ES - Electrostimulation

Inphase/Outphase : combined electroultrasonic stimulation (cEUS) in phase and alternating, respectively
PD : Parkinson's disease (no stimulation)

fUS500kHz or fUS700kHz - 500 kHz or 700 kHz fundamental ultrasound frequency 
tauES100us or tauES300us - 100 us or 300 us electrical pulsewidth

USi : ultrasound intensity (followed by number in W/m^2)
ESi : electrostimulation current (followed by A/m^2)

_spikes : reduced data (containing spike timing) is saved

-> Files starting with data_.... contain the full data, while ..._data-files combine the data of the former. 

For example: save the data_... files obtained by running the run_ES_NN_CWUS_700kHz.m parameter sweep, under the folder: CW US/fUS700kHz/ and place the ..._data file obtained by postproc.m under cdata/, before running plotResults.m

--------------------------------------------------


References
-------------
[1] Tarnaud et al. (2021), Improved alpha-beta power reduction via combined Electrical and Ultrasonic stimulation in a Parkinsonian Cortex-Basal 
Ganglia-Thalamus Computational Model, Journal of Neural Engineering.

[2] Kumaravelu et al. (2016) A biophysical model of the cortex-basal ganglia-thalamus network in the 6-OHDA lesioned rat model of Parkinson's
disease. Journal of computational neuroscience, 40(2), 207-229

[3] Plaksin et aL. (2016) Cell-type-selective effects of intramembrane cavitation as a unifying theoretical framework for ultrasonic 
neuromodulation. Eneuro, 3(3).

[4] Hayut et al. (2011) LTS and FS inhibitory interneurons, short-term synaptic plasticity, and cortical circuit dynamics. PLoS computational biology, 7(10), e1002248.

[5] Vierling-Claassen (2010) Computational modeling of distinct neocortical oscillations driven by cell-type selective optogenetic drive: separable resonant circuits controlled 
by low-threshold spiking and fast-spiking interneurons. Frontiers in human neuroscience, 4, 198.

[6] Lemaire, T. et al. (2019). Understanding ultrasound neuromodulation using a computationally efficient and interpretable model of intramembrane cavitation. 
Journal of neural engineering, 16(4), 046007.

[7] Sherfey J. et al. (2018) DynaSim: A MATLAB Toolbox for Neural Modeling and Simulation. Front. Neuroinform. 12:10. doi: 10.3389/fninf.2018.00010

[8] Observed Brain Dynamics, Partha Mitra and Hemant Bokil, Oxford University Press, New York, 2008.

[9] Jonathan C. Lansey (2021). Beautiful and distinguishable line colors + colormap 
(https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap), MATLAB Central File Exchange. Retrieved October 18, 2021.

Contact
------------------
Email for bugs/questions: thomas.tarnaud@ugent.be
