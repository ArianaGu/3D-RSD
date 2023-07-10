clear;
close all;
%% Load scene data
scene_index = 0;
load(['./data/scene_' int2str(scene_index) '_raw.mat']);

%% Set up parameters
% Virtual wavelength
sc = 3; 
sample_spacing = 0.01;
lambda = 2 * sample_spacing * sc;    
ill_cycle = 5;
peak_ratio = 1e-1;

% Gate
z_gate = 0;
rect_data(:, 1:z_gate) = 0;

% Downsample
K = 0;
ts = ts*2^K;
for k = 1:K
    rect_data = rect_data(:,1:2:end) + rect_data(:,2:2:end);
end

% Truncation
tag_maxdepth = 2;

%% Generate illumination
t_ill = Illumination_kernel(lambda, ill_cycle*lambda, ts, 'gaussian'); % * cycle respect to the virtual wavelength
t_ill = Illumination_pad(t_ill, tag_maxdepth, ts);
[fre_mask, lambda_loop, omega_space, weight] = Illumination_decomp(t_ill, ts, peak_ratio);

%% Generate Fourier components
total_rect_data = total_rect_data(:,1:length(t_ill));
total_u_total = generate_Fcomponents(total_rect_data, fre_mask);

%% Build 3D grid
% Transverse
v_apt_Sz = 2; 
M = 2 * round(round(v_apt_Sz/(sample_spacing))/2);
N = M;
% Axial
depth_spacing = 0.01;
pad_size = 100; % decided by depth range
total_pad_size = 2*pad_size+1;

%% Recon
tic

% 1. 3D RSD
padding_volume = zeros(M, N, total_pad_size, length(lambda_loop));
padding_volume = fit_space(padding_volume, camera_pos, total_u_total, M, N, sample_spacing, depth_spacing, pad_size);
u_volume = f_3D(padding_volume, lambda_loop, weight,pad_size, sample_spacing, depth_spacing, laser_pos_base);
u_volume = permute(u_volume,[2,1,3]);

% 2. BP RSD
% u_volume = f_3D_precise(total_u_total, M, N, lambda_loop, weight, sample_spacing,camera_pos,laser_pos_base);

% 3. BP mask-R RSD
% u_volume = f_3D_R(total_u_total, M, N, lambda_loop, weight, sample_spacing,camera_pos,laser_pos_base, 0.3);

mgn_volume = abs(u_volume);

toc

%% Plot
[img, depth_img] = max(mgn_volume,[],3);
depth_img = depth_img*0.02 + 0.8;

% Denoise
index = find(img < 0.3*max(img,[], 'all'));
depth_img(index) = 0;

% Visualization
flip_and_show(img, 'RSD');
flip_and_show(depth_img, 'depth map');

% load chirp
% sound(y,Fs)
% save(['./data/scene_' int2str(scene_index) '_recon_noise_05.mat'], 'img');

%% print SSIM with respect to BP recon
rsd_img = img / max(img(:));

load(['./BP/volume/scene_' int2str(scene_index) '.mat'])
W = abs(W_c);
W = permute(W, [2 3 1]);

% Denoise
for i = 1:size(W, 3)
    W(:,:,i) = W(:,:,i)*(i*sampling_grid_spacing + minimalpos(3))^2;
end

% Compute SSIM
[bp_img, ~] = max(W,[],3);
bp_img = bp_img / max(bp_img(:));
ssim(rsd_img, bp_img)


%% functions
function [] = flip_and_show(img, title_str)
    img = imresize(img,2);
    img = flip(img,1);
    img = flip(img,2);

    figure;
    imagesc(img);
    colormap 'hot';
%     title(title_str);
    axis image;
    axis off;
%     colorbar;
end

function  new_volume = fit_space(padding_volume, camera_pos, u_total, M, N, sample_spacing, depth_spacing, pad_size)
    n = 0;
    new_volume = zeros(size(padding_volume), 'single');
    clear padding_volume
    for index = 1: size(camera_pos,1)
        x = round(camera_pos(index, 1)/sample_spacing) + M/2;
        y = round(camera_pos(index, 2)/sample_spacing) + N/2;
        z = round(camera_pos(index, 3)/depth_spacing) + pad_size;
        if x>0 && x<M+1 && y>0 && y<N+1 && z>0 && z<2*pad_size+1
            temp = squeeze(u_total(index, :));
            new_volume(x, y, z, :) =  squeeze(new_volume(x, y, z, :))' + temp;
        else
            n = n+1;
        end
    end
    display(n)
end

%% Functions from generate_fourierDhist
function[frequency_index, L, O, W] = Illumination_decomp(t, ts, ratio)
    %{
    Illumination Frequency Decomposition
    Parameters input: 
        1. t: temporal illumination kernel function
        2. ts: temporal sampling rate, unit second
        3. ratio: ratio of the energy coefficient for the fitering step

    Parameters output: 
        1. frequency_index: frequency mask
        2. l: lambda sequence 
        3. o: omega sequence
        4. w: linear weight coefficient
%}
    c_light = 299792458; % speed of light
    
    % Calculate frequency sampling
    fs = 1/ts;
    
    % Get length of the signal
    N = length(t); % N refers to the total discrete sample 
   
    % Energy for the illumination kernel in Fourier domain
    P_wave = abs(fft(t)/N);
    P_wave = P_wave(1:round(N/2+1)); 
    
    % Calculate the ratio based on the Energy peak value
    coeff_ratio = P_wave./max(P_wave(:));
    
    % Find  the filtered frequency index (only index, no unit associate)
    frequency_index = find(coeff_ratio>=ratio); 
    
    % Calculate the final synethsis quantity
    F = fs.*frequency_index./N; % Note that F is converted to analog frequency
    O = 2 * pi * F;
    L = c_light * 2 * pi ./ O;
    W = P_wave(frequency_index); % weight from the gaussian spectrum
    
%     P_wave_c = fft(t)/N;
%     W = P_wave_c(frequency_index);
end

function[kernel] = Illumination_kernel(lambda, length, ts, type)
%{
    Parameters input: 
        1. lambda: virtual wavelength, unit meter
        2. length: length of the signal, unit meter
        3. ts: temporal samplingm unit second
        4. type: type of the illumination kernel

    Parameters output: 
        1. kernel: output 1d illumination kernel
%}
    c_light = 299792458; % speed of light
    bin_res = c_light * ts; % map temporal sampling to spatial meter, unit m
    
    switch type
        case 'gaussian'
            v_s = (lambda/bin_res);
            length = round(length/bin_res);
            sin_pattern = length/v_s; 
            
            cmp_sin = exp(1j * 2 * pi * (sin_pattern * linspace(1,length,length)')/length);
            gauss_wave = gausswin(length, 1/0.3);
            
            kernel = cmp_sin.*gauss_wave;
    end
end

function[t_out] = Illumination_pad(t_in, Maxdepth, ts)
%{
    Parameters input: 
        1. t_in: temproal signal input
        2. Maxdepth: unit meter, maximum depth for signal of interest in reconstruction
        3. ts: temporal sampling rate, unit second

    Parameters output: 
        1. kernel: output 1d illumination kernel
%}
    c_light = 299792458; % speed of light
    bin_res = c_light * ts; % map temporal sampling to spatial meter, unit m
    
    nMax = round(2 * Maxdepth / bin_res); % synethsis illumination kernel length
    
    if (length(t_in)<nMax)
        % Pad the virtual illumination kernel to required size
        t_out = padarray(t_in, round(nMax - length(t_in)), 'post');
    end
    
end


function u_total = generate_Fcomponents(rect_data, fre_mask)
    num_component = length(fre_mask);
    u_total = zeros(size(rect_data,1), num_component);

    for ii = 1 : size(rect_data,1)
        t_tmp = squeeze(rect_data(ii, :));
        f_tmp = fft(t_tmp);
        f_slice = f_tmp(fre_mask);
        u_total(ii,:) = f_slice;
    end
end

