function [] = deconv_Time_Series_Data()
%% deconv_Time_Series_Data - FUNCTION to perform deconvolution of SPED time series stacks. To be used after align_PSF_DataStack.
% Aligned 3D PSF (output of align_PSF_DataStack) is used to deconvolve time
% series image stacks. Deconvolution of image slices is performed using coressponding 2D
% PSF images from 3D PSF.
% Author: Raju Tomer (rajutomer@gmail.com)
%%
%%%%%%%%% Parameters %%%%%%%%%%%
no_of_threads = 12; % no. of threads to use for deconvolution calculations. Each thread process a different time point.
base_dir = 'D:\SPED_data\Exp1';
out_dir = [base_dir '\Data_dc'];
data_dir = [base_dir '\Data'];
psf_fname = [base_dir '\O4xPSF_AvgBeadV2_AvgOf10imagesLo40_Hi100b65k_Cent_114Rev_073m1.tif'];
data_ind_arr = 1:39; % this and the parameter below, psf_ind_arr, defines the z-axis mapping of data and PSF indices
psf_ind_arr = 170:5:360; % same as above
stack_size = [2048 316 40]; % individual data stack size informations, needed to read images stored in raw binary format
start_TP = 0; % range of time points to process, starting time point.
stop_TP = 4000; % range of time points to process, ending time point.
N_IT_arr = [10]; % vector containing number of iterations to use for Richardson Lucy deconvolution
camera_bg = 100; % empty image background
doMedFiltering = 0; % whether to do median filtering
doCropping = 0; % wether to crop the data, if 1 then cropping_rect should be defined
cropping_rect = []; %a 4-element vector with the form [XMIN YMIN XMAX YMAX] from imageJ
upsample_im = 0; % if to upsample the data before deconvolution
upsample_scale = 2; % if above is 1, the umsampling fold factor
write_mip=1; % whether to write maximum intensity projection images of deconvolved time points
tiff_sw = 0; % if the data is in tiff stack (1) or in raw (binary format (0)
reverse_stack = 1; % 1 for reversing image stacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if (numel(data_ind_arr) ~= numel(psf_ind_arr))
    disp('Number of Indices for PSF and data stack do not match, check that');
    return;
end

if (doCropping > 0)
    cropping_rect = cropping_rect + 1;
end

parpool(no_of_threads);

if (tiff_sw > 0)
    list_TPs = dir([data_dir '\*.tif*']);
else
    list_TPs = dir([data_dir '\*.stack*']);
end

kernel_for_PSF_taper = fspecial('gaussian',[20,20],4);

parfor tp = (start_TP + 1):(stop_TP + 1)
    fname_TP = [data_dir '\' list_TPs(tp).name];
    im = [];
    if (tiff_sw == 1)
        im = zeros(stack_size, 'uint16');
        for r = 1:stack_size(3)
            im(:,:,r) = imread(fname_TP,r);
        end
    else    
        fid = fopen(fname_TP, 'r', 'l' );
        im = fread( fid, prod(stack_size), 'uint16' );
        fclose( fid );
        im = reshape(im,stack_size);
        im = permute(im, [2 1 3]);
    end
    if (reverse_stack == 1)
        im = flip(im,3);
    end
    if (doCropping > 0)
        im = im(cropping_rect(2):cropping_rect(4), cropping_rect(1):cropping_rect(3),:)
    end
    for N_IT = N_IT_arr
        ofname = [out_dir '\DC' num2str(N_IT) '_' list_TPs(tp).name '.tif'];
        if (exist(ofname, 'file') == 2)
            disp('skipping');
            ofname
            continue
        end
        ofname_mip = [out_dir '\MAX_DC' num2str(N_IT) '_' list_TPs(tp).name '.tif'];
        ofname_mip_tmp = [out_dir '\Tmp_MAX_DC' num2str(N_IT) '_' list_TPs(tp).name '.tif'];
        first_slice=1;
        tic
        for i = 1:numel(data_ind_arr)
            data_ind = data_ind_arr(i);
            psf_ind = psf_ind_arr(i);    
            im_raw = double(im(:,:,data_ind));
            im_raw = im_raw - camera_bg;
            im_raw(im_raw < 0) = 0;
            if (doMedFiltering > 0)
                im_raw = medfilt2(im_raw);
            end
            if (upsample_im == 1)
                disp('UpSampling Data');
                im_raw = imresize(im_raw, upsample_scale, 'bicubic');
                disp('UpSampling Data Done');
            end
            psf2D = imread(psf_fname, psf_ind);
            psf2D = psf2D - min(psf2D(:));
            psf2D = double( psf2D) ./ sum(sum( double(psf2D)));
            psf2D = edgetaper( psf2D, kernel_for_PSF_taper);
            psf2D = padarray( psf2D, [20,20] );
            data_ind
            N_IT
        	im_dc = deconvlucy( im_raw, psf2D, N_IT );
            im_dc = im_dc - min(min(im_dc));
            imwrite(uint16(im_dc), ofname, 'tiff', 'Compression', 'none', 'WriteMode', 'append');
            if (first_slice == 1)
                max_im = im_dc;
                first_slice = 0;
            else
                max_im = max(im_dc, max_im);
            end            
            if (write_mip == 1)
                imwrite(uint16(max_im), ofname_mip_tmp, 'tiff', 'Compression', 'None');
            end
        end
        toc
        if (write_mip == 1)
            imwrite(uint16(max_im), ofname_mip, 'tiff', 'Compression', 'None');
            delete(ofname_mip_tmp);
        end
    end
end 
delete(gcp('nocreate'));
