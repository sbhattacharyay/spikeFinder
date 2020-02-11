%% Author: Shubhayu Bhattacharyay 
% Department of Biomedical Engineering 
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% January 2020;
%% Loading calcium traces and spike trains from training and testing data
load_calc_trains = {};
load_spike_trains = {};
for i = 1:10
    % load dataset
    dataset = num2str(i);
    calcium_train = csvread([dataset '.train.calcium.csv']);
    spike_train = csvread([dataset '.train.spikes.csv']);
    load_calc_trains{i,1} = calcium_train(2:end,:);
    load_spike_trains{i,1} = spike_train(2:end,:);
end
load_calc_tests = {};
cd 'C:\Users\Shubhayu\Documents\Shubhayu_Projects\Johns Hopkins University\Spike Inference Study\spikefinder.test'
for i = 1:5
    % load test dataset
    dataset = num2str(i);
    load_calcium_test = csvread([dataset '.test.calcium.csv']);
    load_calc_tests{i,1} = load_calcium_test(2:end,:);
end
cd 'C:\Users\Shubhayu\Documents\Shubhayu_Projects\Johns Hopkins University\Spike Inference Study\spikefinder.train'
%% Calculating number of recordings in the dataset and splitting neuronal recordings:
neuron_counts = cellfun(@(x) length(x(1,:)),load_calc_trains);
no_of_neurons = sum(neuron_counts);
dataset_idx_splits = cumsum(neuron_counts);
datasetNo = ["I","II","III","IV","V","VI","VII","VIII","IX","X"]';
scanMethod = ["3D AOD","galvo","resonant","galvo","resonant"]';
indicator = ["OGB-1","OGB-1","GCamp6s","OGB-1","GCamp6s"]';
idxLists = 1:174;

calc_trains = {};
spike_trains = {};
counter = 1;

ogb1_IdxRange1 = 1:dataset_idx_splits(1);
ogb1_IdxRange2 = dataset_idx_splits(1)+1:dataset_idx_splits(2);
ogb1_IdxRange3 = dataset_idx_splits(3)+1:dataset_idx_splits(4);
ogb1_Indices = [ogb1_IdxRange1 ogb1_IdxRange2 ogb1_IdxRange3]';
gcamp_Indices = setdiff(idxLists,ogb1_Indices);

for i = 1:length(load_calc_trains)
    curr_dataset_calc = cell2mat(load_calc_trains(i));
    curr_dataset_spike = cell2mat(load_spike_trains(i));
    [~,col_len] = size(curr_dataset_calc);
    for j = 1:col_len
        calc_trains{counter,1} = curr_dataset_calc(:,j);
        spike_trains{counter,1} = curr_dataset_spike(:,j);
        counter = counter + 1;
    end
end
%% Selecting sample spike train and calcium time-series for algorithm testing
indizes = cellfun(@(x) find(~isnan(x) & ((x~=0 | circshift(x,1)~=0))),calc_trains,'UniformOutput',false);
traces = cellfun(@(x,y) x(y),calc_trains,indizes,'UniformOutput',false);
traces = cellfun(@(x) (x-median(x))/std(x),traces,'UniformOutput',false);
spikes = cellfun(@(x,y) x(y),spike_trains,indizes,'UniformOutput',false);

ffx =cellfun(@(x,y) xcorr(x,y,'unbiased'),traces,spikes,'UniformOutput',false);
filterX =cellfun(@(x,y) y((numel(x)-20):(numel(x)+250)),traces,ffx,'UniformOutput',false); % window of 2 sec

for i = 1:length(filterX)
    curr=cell2mat(filterX(i));
    curr(21)=(curr(20)+curr(22))/2;
    filterX{i}=curr;
end

predictions=cellfun(@(x,y) conv(x,y,'same'),spikes,filterX,'UniformOutput',false);
predictions=cellfun(@(x) circshift(x,115),predictions,'UniformOutput',false);
predictions=cellfun(@(x)(x-median(x))/std(x),predictions,'UniformOutput',false);

corr_coeff=cellfun(@(x,y) corr(x,y),predictions,traces);

[~,selectionIdx]=max(corr_coeff);
[~,g_model]=max(corr_coeff(gcamp_Indices));
[~,o_model]=max(corr_coeff(ogb1_Indices));

g_model = gcamp_Indices(g_model);
o_model = ogb1_Indices(o_model);

%% Removing low-fidelity recordings and separate by indicator type
highFi_Idx = corr_coeff>=.75;
highFi_Idx = (idxLists(highFi_Idx))';

highFi_ogb1 = ogb1_Indices(corr_coeff(ogb1_Indices) >= 0.75);
highFi_gcamp = setdiff(highFi_Idx,highFi_ogb1);

%% mean OGB-1 and GCaMP filters

ogb1_Filters = filterX(highFi_ogb1);
gcamp_Filters = filterX(highFi_gcamp);

combined_ogb1_Filter = [];
combined_gcamp_Filter = [];

for i = 1:length(ogb1_Filters)
    temp = cell2mat(ogb1_Filters(i));
    combined_ogb1_Filter = [combined_ogb1_Filter temp];
end

for i = 1:length(gcamp_Filters)
    temp = cell2mat(gcamp_Filters(i));
    combined_gcamp_Filter = [combined_gcamp_Filter temp];
end

norm_gcampFilters = normalize(combined_gcamp_Filter);
norm_ogb1Filters = normalize(combined_ogb1_Filter);

t_decay = (-20:250)/100;% 100Hz sampling rate,

mean_gcampFilter = mean(norm_gcampFilters,2);
mean_ogb1Filter = mean(norm_ogb1Filters,2);

mean3 = mean(combined_gcamp_Filter,2);
mean4 = mean(combined_ogb1_Filter,2);

std_gcampFilter = std(norm_gcampFilters,0,2);
std_ogb1Filter = std(norm_ogb1Filters,0,2);

std3 = std(combined_gcamp_Filter,0,2);
std4 = std(combined_ogb1_Filter,0,2);

f = figure;
ax = gca;
h=plot(t_decay,mean_gcampFilter,'k','LineWidth',4);
hold on
plot(t_decay,mean_gcampFilter-std_gcampFilter,'-.k','LineWidth',4);
plot(t_decay,mean_gcampFilter+std_gcampFilter,'-.k','LineWidth',4);
xline(0,'--r','LineWidth',3);
xlim([t_decay(1) t_decay(end)]);
set(gca,'YTickLabel',[]);
xlabel('$\Delta t$ (sec)','Interpreter','latex');
ylabel('$\tilde{C}^{(g)}$','Interpreter','latex');
t = text;
lgd = legend;
IEEEFigureSettings(f,h,ax,t,lgd);
set(lgd,'visible','off');
set(t,'visible','off');

f = figure;
ax = gca;
h=plot(t_decay,mean_ogb1Filter,'k','LineWidth',4);
hold on
plot(t_decay,mean_ogb1Filter-std_ogb1Filter,'-.k','LineWidth',4);
plot(t_decay,mean_ogb1Filter+std_ogb1Filter,'-.k','LineWidth',4);
xline(0,'--r','LineWidth',3);
xlim([t_decay(1) t_decay(end)]);
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$\tilde{C}^{(o)}$','Interpreter','latex');
set(gca,'YTickLabel',[]);
t = text;
lgd = legend;
IEEEFigureSettings(f,h,ax,t,lgd);
set(lgd,'visible','off');
set(t,'visible','off');
%% Plot correlations
f = figure;
vlines_idx = dataset_idx_splits + 0.5;
scatter(1:length(corr_coeff),corr_coeff,'filled');
hold on
ylim([0 1]);
xlim([0 175]);
set(gca, 'XTick', []);
ylabel('Correlation')
xlabel('Recording Index')
lgd = legend;
ax = gca;
t = text;
h=plot([0 174],[0.75 0.75],'r--');
for j=1:length(vlines_idx)
    xline(vlines_idx(j),':','LineWidth',4);
end
IEEEFigureSettings(f,h,ax,t,lgd)
set(lgd,'visible','off');
set(t,'visible','off');
%set(h,'visible','off');

%% Plot reconstructions of model recordings

f = figure;
ax = gca;
t = (0:length(predictions{g_model})-1)/100;% 100Hz sampling rate,
plt1=plot(t,traces{g_model},'k');
hold on
plt2=plot(t,predictions{g_model},'r-.');
spikeIdx = find(spikes{g_model}~=0);
h = stem(t(spikeIdx),-1.2*ones(length(t(spikeIdx)),1),'Marker','none','LineWidth',2);
hbase = h.BaseLine;
hbase.Visible = 'off';
set(h,'BaseValue',-2);
xlabel('$t$ (sec)','Interpreter','latex')
yol = t(end);
xlim([0 t(end)])
ylim([-2.5 8])
lgd = legend('$F^{(g^*)}$','$\hat{F}^{(g^*)}$','Location','northeast','Interpreter','latex');
t = text;
IEEEFigureSettings(f,plt1,ax,t,lgd)
IEEEFigureSettings(f,plt2,ax,t,lgd)
lgd.Location='northeast';
%set(lgd,'visible','off');
set(t,'visible','off');

f = figure;
ax = gca;
t = (0:length(predictions{o_model})-1)/100;% 100Hz sampling rate,
plt1=plot(t,traces{o_model},'k');
hold on
plt2=plot(t,predictions{o_model},'r-.');
spikeIdx = find(spikes{o_model}~=0);
h = stem(t(spikeIdx),-1.2*ones(length(t(spikeIdx)),1),'Marker','none','LineWidth',2);
hbase = h.BaseLine;
hbase.Visible = 'off';
set(h,'BaseValue',-2);
xlabel('$t$ (sec)','Interpreter','latex')
xlim([0 yol])
ylim([-2.5 8])
lgd = legend('$F^{(o^*)}$','$\hat{F}^{(o^*)}$','Location','northeast','Interpreter','latex');
t = text;
IEEEFigureSettings(f,plt1,ax,t,lgd)
IEEEFigureSettings(f,plt2,ax,t,lgd)
lgd.Location='northeast';
%set(lgd,'visible','off');
set(t,'visible','off');

% f = figure;
% ax = gca;
% t_decay = (-20:250)/100;% 100Hz sampling rate,
% h=plot(t_decay,filterX{selectionIdx},'k');
% xline(0,'--r','LineWidth',3);
% xlabel('Time after spike (seconds)')
% ylabel('Norm. Ca^{2+} Fluorescence')
% xlim([t_decay(1) t_decay(end)]);
% set(gca, 'YTick', []);
% t = text;
% lgd = legend;
% IEEEFigureSettings(f,h,ax,t,lgd);
% set(lgd,'visible','off');
% set(t,'visible','off');
%% Finding optimal L-tap value for filter
hiFi_tracesGCAMP = traces{highFi_gcamp};
hiFi_tracesOGB1 = traces{highFi_ogb1};
hiFi_spikesGCAMP = spikes{highFi_gcamp};
hiFi_spikesOGB1 = spikes{highFi_ogb1};

hiFi_tracesGCAMP_ffx = ffx{highFi_gcamp};
hiFi_tracesOGB1_ffx = ffx{highFi_ogb1};
L = 2:2:400; %keep L even
corr_values = zeros(length(L),1);

for j = 1:length(L)
    curr_filterX =select_ffx((numel(select_trace)-20):(numel(select_trace)+L(j)));
    curr_filterX(21) = (curr_filterX(20)+curr_filterX(22))/2;
    curr_prediction=conv(select_spikes,curr_filterX,'same');
    curr_prediction=circshift(curr_prediction,(L(j)-20)/2);
    curr_prediction=(curr_prediction-median(curr_prediction))/std(curr_prediction);
    curr_corr = corr(curr_prediction,select_trace);
    corr_values(j)=curr_corr;
end

f = figure;
ax = gca;
h=plot(L/100,corr_values,'k');
xlabel('Length of filter (seconds)');
ylabel('Correlation');
hold on
[maxCorr,optIdx]=max(corr_values);
optL = L(optIdx);
scatter(optL/100,maxCorr,400,'r','filled','h');
t = text;
lgd = legend;
IEEEFigureSettings(f,h,ax,t,lgd);
set(lgd,'visible','off');
set(t,'visible','off');
%% Plotting of optimal filter and reconstruction

new_filterX =select_ffx((numel(select_trace)-20):(numel(select_trace)+optL));
new_filterX(21) = (new_filterX(20)+new_filterX(22))/2;
new_prediction=conv(select_spikes,new_filterX,'same');
new_prediction=circshift(new_prediction,(optL-20)/2);
new_prediction=(new_prediction-median(new_prediction))/std(new_prediction);

f = figure;
ax = gca;
t = (0:length(new_prediction)-1)/100;% 100Hz sampling rate,
plt1=plot(t,new_prediction,'r-.');
hold on
plt2=plot(t,select_trace);
spikeIdx = find(select_spikes~=0);
h = stem(t(spikeIdx),-1.4*ones(length(t(spikeIdx)),1),'Marker','none','LineWidth',2);
hbase = h.BaseLine;
hbase.Visible = 'off';
set(h,'BaseValue',-.6);
xlabel('Time (seconds)')
ylabel('Norm. Ca^{2+} Fluorescence')
xlim([0 t(end)])
lgd=legend('recreated','original','Location','northeast');
t = text;
IEEEFigureSettings(f,plt1,ax,t,lgd);
IEEEFigureSettings(f,plt2,ax,t,lgd);
lgd.Location='northeast';
%set(lgd,'visible','off');
set(t,'visible','off');

f = figure;
ax = gca;
t_decay = (-20:optL)/100;% 100Hz sampling rate,
h=plot(t_decay,new_filterX);
xline(0,'--r','LineWidth',3);
xlabel('Time after spike (seconds)')
ylabel('Norm. Ca^{2+} Fluorescence')
xlim([t_decay(1) t_decay(end)]);
set(gca, 'YTick', []);
t = text;
lgd = legend;
IEEEFigureSettings(f,h,ax,t,lgd);
set(lgd,'visible','off');
set(t,'visible','off');

%% Transient and Frequency Response Properties
model_ogb1_spikes=spikes{o_model}';
model_ogb1_trace=traces{o_model}';

model_gcamp_spikes=spikes{g_model}';
model_gcamp_trace=traces{g_model}';
%% Cross-validation to prepare for system identification
% 3-fold cross validation for OGB1
K1 = 3;
ogb1_testCV = crossvalind('Kfold',highFi_ogb1,K1);
% 5-fold cross validation for GCaMP variants
K2 = 5;
gcamp_testCV = crossvalind('Kfold',highFi_gcamp,K2);

gcamp_testsIdx = {};
ogb1_testsIdx = {};

gcamp_trainsIdx = {};
ogb1_trainsIdx = {};

for i = 1:K1
    ogb1_testsIdx{i}=highFi_ogb1(ogb1_testCV == i);
    ogb1_trainsIdx{i}=highFi_ogb1(ogb1_testCV ~= i);
end

for j = 1:K2
    gcamp_testsIdx{j}=highFi_gcamp(gcamp_testCV == j);
    gcamp_trainsIdx{j}=highFi_gcamp(gcamp_testCV ~= j);
end

Ts = 0.01;

gcampTrain1=iddata(traces(gcamp_trainsIdx{1}),spikes(gcamp_trainsIdx{1}),Ts);
gcampTrain2=iddata(traces(gcamp_trainsIdx{2}),spikes(gcamp_trainsIdx{2}),Ts);
gcampTrain3=iddata(traces(gcamp_trainsIdx{3}),spikes(gcamp_trainsIdx{3}),Ts);
gcampTrain4=iddata(traces(gcamp_trainsIdx{4}),spikes(gcamp_trainsIdx{4}),Ts);
gcampTrain5=iddata(traces(gcamp_trainsIdx{5}),spikes(gcamp_trainsIdx{5}),Ts);

gcampTest1=iddata(traces(gcamp_testsIdx{1}),spikes(gcamp_testsIdx{1}),Ts);
gcampTest2=iddata(traces(gcamp_testsIdx{2}),spikes(gcamp_testsIdx{2}),Ts);
gcampTest3=iddata(traces(gcamp_testsIdx{3}),spikes(gcamp_testsIdx{3}),Ts);
gcampTest4=iddata(traces(gcamp_testsIdx{4}),spikes(gcamp_testsIdx{4}),Ts);
gcampTest5=iddata(traces(gcamp_testsIdx{5}),spikes(gcamp_testsIdx{5}),Ts);

ogb1Train1=iddata(traces(ogb1_trainsIdx{1}),spikes(ogb1_trainsIdx{1}),Ts);
ogb1Train2=iddata(traces(ogb1_trainsIdx{2}),spikes(ogb1_trainsIdx{2}),Ts);
ogb1Train3=iddata(traces(ogb1_trainsIdx{3}),spikes(ogb1_trainsIdx{3}),Ts);

ogb1Test1=iddata(traces(ogb1_testsIdx{1}),spikes(ogb1_testsIdx{1}),Ts);
ogb1Test2=iddata(traces(ogb1_testsIdx{2}),spikes(ogb1_testsIdx{2}),Ts);
ogb1Test3=iddata(traces(ogb1_testsIdx{3}),spikes(ogb1_testsIdx{3}),Ts);
%% System Identification

%% Stochastic Gradient Descent for System Identification
% System identification approach to estimate the waveform of the pulses
% Learn a basis for the pulse shape
lms = dsp.LMSFilter;

% Generate own signal given hypothesis shape:
% 1. Assume decaying exponential --> parameter estimation

% Use spike times to convolve with hypothesis shape
% Solve for optimal parameter
% Deisgnate a length for shape filter
% implement gradient descent to solve for parameters

% Search sys id methods/examples in MATLAB

% Generate my own signal, try out to solve for opt parameters
% Start with exact gradient descent
% System identification
% No assumptions about the shape
%% Gram-Schmidt for recursive component estimation

% Test both different shapes and differnet decay factors


% Subtract converged estimations from y to find additional components
% Error vs components
% https://arxiv.org/pdf/0912.1637.pdf

n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for i = 2:k
    U(:,i) = V(:,i);
    for j = 1:i-1
        U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
    end
    U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
end


%% Spike Detection with block-sparse estimation
% Solve a (weighted) block sparse coding problem (there must be libraries
% to do it)

%% Encoder-decoder architecture
% Use an auto-encoder (e.g., U-net) and a loss that balances reconstruction
% accuracy VS classification accuracy. Could initialize first layer to
% filters learned in a).