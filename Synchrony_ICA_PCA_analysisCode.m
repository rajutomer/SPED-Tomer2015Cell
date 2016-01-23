%% Code to perform synchrony, PCA analysis, plotting and labelling of 3D fish volumes
% Author: Raju Tomer (rajutomer@gmail.com), with inputs from Sethuraman
% Sankaran on dF/F traces filtering (as described in Rajasethupathy et al Nature 2015).
%%

%% Input Parameters %%
data_set = 2; % 1: whole brain, 4x; 2: whole nervous system; 3: Brain + spinal cord with 10x
data_dir = 'D:\Segmentation\Data\';
out_dir = 'D:\Segmentation\Data\PCA_ICA\';
use_global_cutoff = 0; % if to use global cutoff for filtering traces

if (data_set == 1)
    fname = [data_dir 'ImO4x_12Hz_StD_TL0-467_DC10_VW0_LOC000D_CM0_CHN00_T0000_.stack.tif.fnuclei.txt']; %O4x12Hz, whole brain
    base_dir = 'D:\Segmentation\O4x_12Hz';
    fname_fnuc = [base_dir '\ImO4x_12Hz_StD_TL0-467_DC10_VW0_LOC000D_CM0_CHN00_T0000_.stack.tif.fnuclei.tif'];
    base_out_name = 'O4x12Hz_';
    tstart = 101; % Brain, 12 Hz, remove the first 100 time points to remove initial "light" illumination induced behavior
    vol_rate = 12; % volumes per second
elseif (data_set==2)
    fname = [data_dir 'ImO4x_StD_TL0-2500_DC10_VW0_LOC000D_CM0_CHN00_T0000_.stack.fnuclei_AMat.txt']; %O4x, whole nervous system
    base_dir = 'D:\Segmentation\O4x_39s';
    fname_fnuc = [base_dir '\ImO4x_StD_TL0-2500_DC10_VW0_LOC000D_CM0_CHN00_T0000_.stack.fnuclei.tif'];
    base_out_name = 'O4x39s_PC';
    tstart = 401; % NS, 6.23 Hz, remove the first 400 time points to remove initial "light" illumination induced behavior
    vol_rate = 6.23; % volumes per second
else
    fname = [data_dir 'StD_TL0-1239_DC10_VW0_LOC000D_CM0_CHN00_T0000_.stack_fnuclei_AMat.txt']; data_set=3; %O10x, whole brain
    base_dir = 'D:\Segmentation\O10x';
    fname_fnuc = [base_dir '\StD_TL0-1239_DC10_VW0_LOC000D_CM0_CHN00_T0000_.stack_fnuclei.tif'];
    base_out_name = 'O10x_PC';
    tstart = 251; % Brain, 4.14 Hz, remove the first 300 time points to remove initial "light" illumination induced behavior
    vol_rate = 4.14; % volumes per second
end

%% load data
full_dat = load(fname);
time_series = full_dat(:,2:end);

%% normalize by subtracting mean and dividing by mean of a trace
time_series_norm = time_series;
for i = 1:size(time_series,1)
    time_series_norm(i,:) = (time_series(i,:) - mean(time_series(i,:)))./mean(time_series(i,:));
end

%% filter to identify active cells
time_series_norm_filt = [];
cell_ids_filt = [];
cutoff = [];
cutoff_v = [];
if (use_global_cutoff == 0)
    for i = 1:size(time_series_norm,1)
        cutoff = calculate_iterative_noise(time_series_norm(i,:));
        disp(cutoff);
        if (max(time_series_norm(i,tstart:end))>cutoff)
            time_series_norm_filt = [time_series_norm_filt; time_series_norm(i,tstart:end)];
            cell_ids_filt = [cell_ids_filt; full_dat(i,1)];
            cutoff_v = [cutoff_v; cutoff]; 
        end
    end
else
    cutoff = calculate_iterative_noise(time_series_norm(:));
	disp(cutoff);
    for i = 1:size(time_series_norm,1)
        if (max(time_series_norm(i,tstart:end))>cutoff)
            time_series_norm_filt = [time_series_norm_filt; time_series_norm(i,tstart:end)];
            cell_ids_filt = [cell_ids_filt; full_dat(i,1)];
        end
    end
end

time_series_norm_filt_counts = [];
if (use_global_cutoff > 0)
    time_series_norm_filt_counts = time_series_norm_filt > cutoff;
else
    sz = size(time_series_norm_filt);
    time_series_norm_filt_counts = time_series_norm_filt > repmat(cutoff_v, [1 sz(2)]) ;
end    
fract_active_cells = sum(time_series_norm_filt_counts);
t=1:size(fract_active_cells,2);
t=t/vol_rate;
plot(t,fract_active_cells); xlim([0,max(t)]);    

%% ICA calculations
time_series_norm_filt_Z = time_series_norm_filt;
for i = 1:size(time_series_norm_filt,1)
    time_series_norm_filt_Z(i,:) = (time_series_norm_filt(i,:) - mean(time_series_norm_filt(i,:)))./std(time_series_norm_filt(i,:));
end

no_of_eig = 20; % for 10x 25 eig, 20 ICs
[icasig, A, W] = fastica (time_series_norm_filt_Z, 'numOfIC', 10, 'lastEig', no_of_eig); % Note that, the initialy seeding is choosing randomly, therefore the results may vary in every run even with same set of parameters.
icasig_pos = [];
pos = [];
for i = 1:size(icasig,1)
    if ( ( min(icasig(i,:)) + max(icasig(i,:)) ) > 0)
        pos = [pos; i];
        icasig_pos = [icasig_pos; icasig(i,:)] ;
    end
end

icasig_pos = icasig_pos./max(icasig_pos(:));

t=1:size(icasig_pos,2);
t=t/vol_rate;

sz1 = size(icasig_pos,1);
hold on
for i = 1:sz1
    plot(t, (icasig_pos(i,:) + sz1 - i));
end
xlim([0 max(t)]);
hold off

%% Code to label Fish imaging volume
% requires DIPlib and DIPimage (http://www.diplib.org/)
info.Data = imfinfo(fname_fnuc);
[N_slices, ~]  = size(info.Data);
tmp = imread(fname_fnuc, 1);
sz = size(tmp); clear tmp
fnuclei = zeros(sz(1), sz(2), N_slices, 'uint16');
for i = 1:N_slices
    fnuclei(:,:,i) = imread(fname_fnuc,i);
end
fnuclei = dip_image(fnuclei, 'uint16');
msr = measure(fnuclei,[],({'center'}));

id = msr.id ;
maxid = max(id);

for i = 1:size(pos,1)
    pos(i)
    labs = zeros(1,maxid);
    labs(cell_ids_filt) = abs(W(pos(i),:));
    labs = [0 labs];
    fnuclei_lab1 = lut(fnuclei,labs);
    fname_stack1 = [out_dir base_out_name '_ICA_' num2str(pos(i)) '_Eig' num2str(no_of_eig) '_ICs.ics'];
    writeim(fnuclei_lab1,fname_stack1,'ICSv2',0,[]);
end

%% PCA calculations
time_series_norm_filt_Z_tr = time_series_norm_filt_Z';
[coeff_weightage, tr_score, lat] = princomp(time_series_norm_filt_Z_tr);

%Plot graphs
s=size(tr_score);
t=1:s(1);
t=t/vol_rate;
plot(t,tr_score(:,1),'red',t,tr_score(:,2),'green',t,tr_score(:,3),'blue'); xlim([0 max(t)]);

%Code to label Fish imaging volume
%requires DIPlib and DIPimage (http://www.diplib.org/)
info.Data = imfinfo(fname_fnuc);
[N_slices, ~]  = size(info.Data);
tmp = imread(fname_fnuc, 1);
sz = size(tmp); clear tmp
fnuclei = zeros(sz(1), sz(2), N_slices, 'uint16');
for i = 1:N_slices
    fnuclei(:,:,i) = imread(fname_fnuc,i);
end
fnuclei = dip_image(fnuclei, 'uint16');
msr = measure(fnuclei,[],({'center'}));

id = msr.id ;
maxid = max(id);

labs = zeros(1,maxid);
labs(cell_ids_filt) = abs(coeff_weightage(:,1));
labs = [0 labs];
fnuclei_lab1 = lut(fnuclei,labs);
fname_stack1 = [out_dir base_out_name 'Stack_absCoeff_PC1.ics'];
writeim(fnuclei_lab1,fname_stack1,'ICSv2',0,[])

labs = zeros(1,maxid);
labs(cell_ids_filt) = abs(coeff_weightage(:,2));
labs = [0 labs];
fnuclei_lab2 = lut(fnuclei,labs);
fname_stack2 = [out_dir base_out_name '_Stack_absCoeff_PC2.ics'];
writeim(fnuclei_lab2,fname_stack2,'ICSv2',0,[])

labs = zeros(1,maxid);
labs(cell_ids_filt) = abs(coeff_weightage(:,3));
labs = [0 labs];
fnuclei_lab3 = lut(fnuclei,labs);
fname_stack3 = [out_dir base_out_name 'Stack_absCoeff_PC3.ics'];
writeim(fnuclei_lab3,fname_stack3,'ICSv2',0,[])
