%% Author: Raju Tomer (rajutomer@gmail.com)
% Uses DIPimage and DIPlib libraries. Code partly based on PointCloudXplore
% toolbox (http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=pcx)
%%

%% Parameters
aspect = [0.585 0.585 5]; % for 10X
%aspect = [0.73 0.73 5]; % for 4x, data was up-sampled 2-fold before deconvolution
bg = 3; % background level determined by inspetcing images in ImageJ
fname = 'D:\Segmentation\O10x\StD_TL0-1239_DC10_VW0_LOC000D_CM0_CHN00.stack.PIC';
%%

%% read data
ar = aspect(1)/aspect(3);
image_data = readim(fname,'PIC') ; % read image
nuc_image = dip_image(squeeze(image_data),'uint16');
clear image_data
nuc_image = nuc_image - bg ;
nuc_image(nuc_image<1) = 0;
%%

%% Generate signal mask
nuc_image_sm = nuc_image;
nuc_image_sm = gaussf_iir(nuc_image_sm,1*[1,1,ar]);
image_shell = nuc_image_sm > gaussf_iir(nuc_image_sm,10*[1,1,ar]);
image_shell(nuc_image<1) = 0;
%%

%% Normalized convolution for local intensity correction
local_thresh = nuc_image_sm;
local_thresh(~image_shell) = 0;
local_thresh = gaussf_iir(local_thresh,5*[1,1,ar]);
shell_sm = gaussf_iir(+image_shell,5*[1,1,ar]);
local_thresh = local_thresh(image_shell)/shell_sm(image_shell);
clear shell_sm
nuc_signal = newim(image_shell,'bin');
nuc_signal(image_shell) = nuc_image_sm(image_shell)>local_thresh;

ac = nuc_image_sm;
ac(~nuc_signal) = 0;
ac = gaussf_iir(ac,5*[1,1,ar]);
shell_sm = gaussf_iir(+nuc_signal,5*[1,1,ar]);
mask = shell_sm>1e-2;
ac = ac(mask)/shell_sm(mask);
nuc_image = dip_image(nuc_image,'sfloat');
nuc_image(mask) = nuc_image(mask) / ac;
nuc_image(~mask) = 0;   % Final dnaimg
nuc_image_sm(mask) = nuc_image_sm(mask) / ac;
nuc_image_sm(~mask) = 0;
[tmp,th] = threshold(nuc_image_sm(image_shell),'isodata');
clear tmp
nuc_signal = nuc_image_sm > 0.9;   % Final DNA mask
%%

%% Find local maxima
nuc_image_sm = -nuc_image_sm;
nuclei = dip_localminima(nuc_image_sm,nuc_signal,3,10,3,0);
nuclei = dip_image(nuclei,'uint16');
msrO = measure(nuclei,[],'center');
no_nuc = size(msrO,1);
id = msrO.id;
no_nuc
maxid = max(id);
lab = zeros(1,maxid);
lab(id) = id;
clear id maxid msrO
lab = uint16([0,lab]);
nuclei = lut(nuclei,lab);
%%

%% nuclei modelling
tmp = ~nuc_signal;
dnaimg_model = nuc_image;
dnaimg_model(tmp) = 4095;
offset = erosion(dnaimg_model,parabolic_size*[1,1,ar],'parabolic');
offset = offset(nuc_signal);
dnaimg_model(tmp) = 0;
range = dilation(dnaimg_model,parabolic_size*[1,1,ar],'parabolic');
range = max(range(nuc_signal)-offset,1); % Avoid division by zero.
tmp = (dnaimg_model(nuc_signal)-offset)/range;
tmp = 1-clip(tmp,1,0);
dnaimg_model(nuc_signal) = exp(3*tmp);
fnuclei = nuclei;
fnuclei = dip_growregionsweighted(fnuclei,dnaimg_model,nuc_signal,aspect,5,[]);
%%

%% filter by size
msr = measure(fnuclei,[],({'size','center'}));
lower_sz_thresh = 10;
to_rm = msr.size < lower_sz_thresh;
id = msr.id ;
sum(to_rm)
maxid = max(id);
lab = zeros(1,maxid);
lab(id) = id;
lab(id(to_rm)) = 0;
lab = uint16([0,lab]);
nuclei = lut(nuclei,lab);
%%

%% Final segmented volume
fnuclei = nuclei;
fnuclei = dip_growregionsweighted(fnuclei,dnaimg_model,nuc_signal,aspect,5,[]);
%%