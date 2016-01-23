function [] = proj_In_Time()
%% max_proj_in_time - generates maximum intensity, avg intensity and standard deviation projections along time points
% Author: Raju Tomer (rajutomer@gmail.com)
%%

%%
%%%% Parameters
no_of_threads = 3; % no. of workers for parallel processing
base_dir = 'D:\SPED_data\Exp1'; % data base directory
data_dir = [base_dir '\LOC000_dc']; % data dir
out_dir = [base_dir '\Combos']; % output dir
N_IT_Arr = [10 20]; %which iterations deconvolved data to process, parse filenames
start_TP = 0; % start time point
stop_TP = 2500; % stop time point
%%%%%%%%%%%%%%%
%%

%%
parpool(no_of_threads);
for N_IT = N_IT_Arr
    list = dir([data_dir '\DC' num2str(N_IT) '*.tif']);
    finf = imfinfo([data_dir '\' list(1).name]);
    stack_size = [finf(1).Height finf(1).Width numel(finf)];
    max_proj_stack = zeros(stack_size, 'uint16');
    avg_proj_stack = zeros(stack_size, 'uint16');
    std_stack = zeros(stack_size, 'uint16');
    ofname_max = [out_dir '\MAX_TL' num2str(start_TP) '-' num2str(stop_TP) '_' list(1).name '.tif'];
    ofname_avg = [out_dir '\Avg_TL' num2str(start_TP) '-' num2str(stop_TP) '_' list(1).name '.tif'];
    ofname_std = [out_dir '\StD_TL' num2str(start_TP) '-' num2str(stop_TP) '_' list(1).name '.tif'];
    parfor j = 1:stack_size(3)
        j
        slice_across_time = zeros([stack_size(1) stack_size(2) (stop_TP - start_TP + 1)], 'single');
        for i = (start_TP + 1):(stop_TP + 1)
            im = imread([data_dir '\' list(i).name],j);
            slice_across_time(:,:,i) = single(im);
        end
        max_proj_stack(:,:,j) = uint16(max(slice_across_time, [], 3));
        avg_proj_stack(:,:,j) = uint16(mean(slice_across_time, 3));
        std_stack(:,:,j) = uint16(std(slice_across_time, 0, 3));
    end
    for k = 1:stack_size(3)
        imwrite(uint16(max_proj_stack(:,:,k)), ofname_max, 'tiff', 'Compression', 'none', 'WriteMode', 'append');
        imwrite(uint16(avg_proj_stack(:,:,k)), ofname_avg, 'tiff', 'Compression', 'none', 'WriteMode', 'append');
        imwrite(uint16(std_stack(:,:,k)), ofname_std, 'tiff', 'Compression', 'none', 'WriteMode', 'append');
    end
end
delete(gcp('nocreate'));
