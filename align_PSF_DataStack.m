function [] = align_PSF_DataStack()
% align_PSF_DataStack - FUNCTION to perform alignment (along z-axis) of
% the SPED empirical PSF with the raw SPED image stack.
% A defined set of slices in the data stack are deconvolved by a defined
% set of 2D PSF images sampled along the z-axis.
% Output is a series of deconvolved images which are then manually
% inspected to determine the global z-axis alignment of the PSF with the raw data stack.
% Author: Raju Tomer (rajutomer@gmail.com), with initial inputs from Sean Quirin.
%%

%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
no_of_threads = 12; % no of threads to use for computation.
if_tif_ser = 0; % if the stack is a series of tiff files (0) otherwise a 3D tiff stack assumed
base_dir = 'D:\SPED_data\Exp1'; % base directory containing data sets
out_dir = [base_dir '\match']; % output dir
data_fname = [base_dir '\VW0_LOC000D_CM0_CHN00_T0000_.stack_rev.tif']; % file name for stack file, if tiff series this should be full dir name
sped_psf_fname = [base_dir '\O4xPSF.tif']; %SPED psf fine 3d tiff stack
dataInd = 1:10:39; % indices for 2d images in raw data stack, deconvolution is performed for these slices only
N_IT_arr = [10 15]; % vector containing number of different iterations to use for Richardson Lucy deconvolution
psf_step_mult = 10; % this and the parameter below define the indices in PSF stack to use for deconvolution
psf_ind_arr = 2:80; % this vector is multiplied by psf_step_mult for indices in PSF stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
parpool(no_of_threads); % start the thread pool
kernel_for_PSF_taper = fspecial('gaussian',[20,20],4); % for minimizing edge artefacts
if (if_tif_ser == 1)
    list = dir([data_fname '\*.tif']);
    fprintf(['Number of images, ' num2str(numel(list)) '\n']);
else
    info.Data = imfinfo(data_fname);
    [N_im, ~]  = size(info.Data);
    fprintf(['Number of images, ' num2str(N_im) '\n']);
end

for N_IT = N_IT_arr
    for ind = dataInd
        if (if_tif_ser == 1)
            im_raw = double( imread( [data_fname '\' list(ind).name] ) );
        else
            im_raw = double( imread( data_fname, ind ) );
        end
        parfor j=psf_ind_arr
            k = psf_step_mult*j;
            psf2D = imread(sped_psf_fname, k);
            psf2D = psf2D - min(psf2D(:));
            psf2D = double( psf2D) ./ sum(sum( double(psf2D)));
            psf2D = edgetaper( psf2D, kernel_for_PSF_taper);
            psf2D = padarray( psf2D, [30,30] );
            tic
            im_dc = deconvlucy( im_raw, psf2D, N_IT );
            toc
            im_dc = im_dc - min(min(im_dc));
            if k < 10
                indnum = ['00' num2str(k)];
            elseif k < 100
                indnum = ['0' num2str(k)];
            else
                indnum = ['' num2str(k)];
            end
            imwrite(uint16(im_dc),[out_dir '\SliceNo' num2str(ind) '_iter' num2str(N_IT) 'TestedWith_' indnum '.tif']);
        end
    end
end
delete(gcp('nocreate'));
