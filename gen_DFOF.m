function [] = gen_DFOF()
%% gen_DFOF generate DFOF : Delta F over F is calculated by subtracting and 
% dividing by baseline signal image calculated using gen_Baseline_For_DFOF.
% Author: Raju Tomer (rajutomer@gmail.com)
%%

%% Parameters
no_of_threads = 12;
stack_size = [2048 632 39]; % image stack size
base_dir = 'D:\SPED_data\Exp1';
dir_data = [base_dir '\LOC000_dc'];
dir_out = [base_dir '\DFOF'];
avg_img_fname = [base_dir '\ref_im\Avg_medF_DC10_Range1-2500.tif']; % base line F, generated using gen_Baseline_For_DFOF
DC_IT = 10; % used for parsing file names
start_image_ind = 1; % start ind for first time point
end_image_ind = 2500; % end ind for last time point
rOffset = 1; % Ratio offset added for increasing visualization range
denom_offset = 10; % offset added to baseline F signal to avoid divisions by zero.
scale_range = 5000; %for scaling ratios to image intensities
sig_im_threshold = [1 23]; % a threshold used to remove noise outside the sample.
if_med_filt = 1; %whether to do median filtering?
med_filt_kernel = 3; %median filter kernel
write_max_proj = 1; % whether to write maximum intensity projection image
if_tiff_stack = 1; % image images are stored as tiff stack?
if_reverse_stack = 0; % whether to reverse the slice order in the stack
if_raw = 0; % if the the data is raw or deconvolved, needed for appropriate file parsing and saving
if_sub_bg = 0; % whether to subtract camera background? 
camera_bg = 100; % camera background
%%

parpool(no_of_threads);
if (if_raw == 1)
    list = dir([dir_data '\*.stack']);
else
    list = dir([dir_data '\DC' num2str(DC_IT) '*.tif']);
end

%Read Avg
avg_img = zeros(stack_size(2), stack_size(1), stack_size(3), 'uint16');
for k = 1:stack_size(3)
    avg_img(:,:,k) = imread(avg_img_fname, k);
end
if (if_sub_bg == 1)
	avg_img = avg_img - camera_bg;
    avg_img(avg_img < 0) = 0;
    %avg_img = avg_img + 1;
end
avg_img = double(avg_img);

parfor i = start_image_ind:end_image_ind
    fname = [dir_data '\' list(i).name];
    im =[];
    if (if_tiff_stack > 0)
        im = zeros([stack_size(2) stack_size(1) stack_size(3)], 'uint16');
        for r = 1:stack_size(3)
            im(:,:,r) = imread(fname,r);
        end
    else
        fid = fopen(fname, 'r', 'l' );
        im = fread( fid, prod(stack_size), 'uint16' );
        fclose( fid );
        im = reshape(im,stack_size);
        im = permute(im, [2 1 3]);    
    end
	if (if_reverse_stack == 1)
        im = flip(im,3);
	end
    if (if_sub_bg == 1)
        im = im - camera_bg;
        im(im < 0 ) = 0;
    end
    im = uint16(((double(im) - 0.6*avg_img)./(avg_img + denom_offset) + rOffset).*scale_range);
    im (im < 0) = 0;
    if (sig_im_threshold(1) == 1)
        sig_im_threshold(2)
        im(avg_img <= sig_im_threshold(2)) = 0;
    end
    if (if_med_filt == 1)
        for rj = 1:stack_size(3)
            im(:,:,rj) = medfilt2(im(:,:,rj), [med_filt_kernel med_filt_kernel]);
        end
    end
    ofname = [dir_out '\F_' list(i).name '.tif'];
    ofname_mip = [dir_out '\MAX_F_' list(i).name '.tif'];
    i
    for j = 1:stack_size(3)
        imwrite(im(:,:,j), ofname, 'tiff', 'Compression', 'none', 'WriteMode', 'append');
    end
    if (write_max_proj == 1)
        imwrite(max(im,[],3), ofname_mip, 'tiff', 'Compression', 'None');
    end
end

delete(gcp('nocreate'));
