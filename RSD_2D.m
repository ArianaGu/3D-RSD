clear;
close all;
%% Load scene data
scene_index = 10;
load(['./data/scene_' int2str(scene_index) '_total.mat']);

% enable phase mapping compensation
tic
% total_u_total = phase_mapping(total_u_total, camera_pos);

toc

u_total = reshape(total_u_total, 190, 190, []);
sampling_spacing = 0.01;
aperturefullsize = [1.9, 1.9];

%% Set up parameters
% Virtual wavelength
sc = 3; 
sample_spacing = 0.01;
lambda = 2 * sample_spacing * sc;    
ill_cycle = 5;
peak_ratio = 1e-1;

% Truncation
tag_maxdepth = 2;

%% Generate illumination
t_ill = Illumination_kernel(lambda, ill_cycle*lambda, ts, 'gaussian'); % * cycle respect to the virtual wavelength
t_ill = Illumination_pad(t_ill, tag_maxdepth, ts);
[fre_mask, lambda_loop, omega_space, weight] = Illumination_decomp(t_ill, ts, peak_ratio);

%% Build 3D grid
depth_min = 0.8; % depth minimal volumn
depth_max = 1.8; % depth maximum volumn
v_apt_Sz = 2.4; % virtual aperture maximum size, automatical pad to symmetry (square), unit meter, standard size 2
d_offset = 0; % electrical delay as an constant offset 


%% Basic parameters for wavefront cube
c_light = 299792458; % define speed of light

% Pad Virtual Aperture Wavefront
% additional if physical dimension is odd number
if(mod(size(u_total,1),2)==1)
    [tmp_x, tmp_y, tmp_z] = size(u_total);
    tmp3D = zeros(round(tmp_x/2)*2, round(tmp_y/2)*2, tmp_z);
    aperturefullsize = [size(tmp3D,1)-1, size(tmp3D,2)-1] * sample_spacing;
    tmp3D(1:size(u_total,1), 1:size(u_total,2),:) = u_total;
    u_total = tmp3D;
end

% perallocated memeory for padding wavefront
u_tmp = zeros(2 * round(round(v_apt_Sz/(sample_spacing))/2), 2 * round(round(v_apt_Sz/(sample_spacing))/2), size(u_total,3));

% create Virtual Aperture by zero padding
for index = 1 : size(u_total,3)
    tmp = squeeze(u_total(:,:,index));
    [u_tmp(:,:,index), apt_tmp, ~]= Create_VirtualAperture(tmp, aperturefullsize, v_apt_Sz, 0);
end

aperturefullsize = apt_tmp; % update virtual aperture size
u_total = u_tmp; % update wavefront cube

% create depth slice for the volume
depth_loop = depth_min : 2 * sample_spacing : depth_max;


%% Reconstruction using fast RSD
nZ = length(depth_loop);
u_volume = zeros([size(u_total,1), size(u_total,2), nZ]);

tic
for i = 1 : length(depth_loop)
    depth = depth_loop(i);
    u_tmp = zeros(size(u_total,1), size(u_total,2));    
        
    for spectrum_index = 1 : length(lambda_loop)
        u_field = squeeze(u_total(:,:,spectrum_index));
        lambda = lambda_loop(spectrum_index);
        omega = omega_space(spectrum_index);
        
        % Define time of arrival and apply fast RSD reconstruction formula
        u1 = u_field.*exp(1j * omega * (depth + d_offset)/c_light * ones(size(u_field))); 
        u2_RSD_conv = Camera_Focusing(u1, aperturefullsize, lambda, depth, 'RSD convolution', 0);     
        u_tmp = u_tmp + weight(spectrum_index) * u2_RSD_conv;

    end

    u_volume(:,:,i) = u_tmp;
end

toc

%% Visualization part
mgn_volume = abs(u_volume);

% Maximum projection along depth direction
[img, ~] = max(mgn_volume,[],3);

% Upsampling for visualization
% img = imresize(img,2);

% Adjust alignment
rsd_img = img(13:212,9:208);
rsd_img = rsd_img / max(rsd_img(:));

figure;
imagesc(imresize(rsd_img,2));
colormap 'hot';
axis image;
axis off;

%% print SSIM with respect to BP recon
load(['./BP/volume/scene_' int2str(scene_index) '.mat'])
W = abs(W_c);
W = permute(W, [2 3 1]);
% Denoise
for i = 1:size(W, 3)
    W(:,:,i) = W(:,:,i)*(i*sampling_grid_spacing + minimalpos(3))^2;
end

[bp_img, ~] = max(W,[],3);
bp_img = bp_img / max(bp_img(:));
ssim(rsd_img, bp_img)

% cross correlation analysis to decide alignment adjustment
crossCorr_fft = fft2(bp_img) .* conj(fft2(rsd_img));
crossCorr = ifftshift(ifft2(crossCorr_fft));
figure()
imagesc(crossCorr)
[maxVal, maxIndex] = max(abs(crossCorr(:)));
[row, col] = ind2sub(size(crossCorr), maxIndex);
fprintf('The maximum value is %f at row %d and column %d.\n', maxVal, row, col);

%% Alternative functions
function mapped_u_total = phase_mapping(u_total, camera_pos)
    mapped_u_total = zeros(size(u_total));
    z = camera_pos(:,3);
    for lambda = 1:size(u_total, 2)
        mapped_u_total(:, lambda) = u_total(:, lambda).*exp(1j*2*pi/lambda*z);
    end
end

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

%{
    Parameters input: 
        1. uin: input virtual wavefront
        2. L: aperture size, unit meter
        3. lambda: virtual wavelength, unit meter
        4. depth: propagation (focusing) depth, unit meter
        5. method: string for choosing variety wave propagation kernel 
        6. alpha: coefficient for accuracness, ideally >= 1, 0.1 is
        acceptable for small amplitude error (large phase error)

    Parameters output: 
        1. uout: output wavefront at the destination plane
%}
function[uout] = Camera_Focusing(uin, L, lambda, depth, method, alpha)
    if (alpha~= 0)
        [u1_prime, aperturefullsize_prime, pad_size, Nd] = tool_fieldzeropadding(uin, L, lambda, depth, alpha);
    else
        u1_prime = uin;
        aperturefullsize_prime = L;
        pad_size = 0;
        Nd = size(u1_prime,1);
    end
    
    switch method        
        case 'RSD convolution'
            uout = fliplr(propRSD_conv(u1_prime, aperturefullsize_prime, lambda, depth));
        
    end
    uout = uout(pad_size+1: pad_size+Nd, pad_size+1: pad_size+Nd);

% Padding the Wavefront
function[u2, phyoutapt, pad_size, sM] = tool_fieldzeropadding(u1, phyinapt, lambda, depth, alpha)
    % % Input: 
    % u1: input complex field
    % phyinapt: Input physical aperture dimension
    % lambda : wavelength, unit m
    % depth: unit m
    % alpha: uncertainty parameter
    % % Output: 
    % u2: output complex field
    % phyoutapt: Output physical aperture dimension
    % pad_size: one sided padding size
    % sM: ola square symmetry, parameter for redoing the padding

    % This program input x-y sampling interval is the same, this program
    % also check if input field is symmetry (square)
    
    [M, N] = size(u1);
    delta = phyinapt(1)/(M-1);
    
    if(M == N)

    else
        [u1, ~] = tool_symmetry(u1, phyinapt);
    end
    
    [sM, ~] = size(u1); % update new discret index, sM output for redo the padding
    
    N_uncertaint = (lambda * abs(depth))/(delta^2);
    pad_size = (N_uncertaint - sM)/2;

    if(pad_size > 0)
        pad_size = round(alpha * pad_size);
    else
        pad_size = 0;
    end
    
    u2 = padarray(u1, [pad_size,pad_size], 0);
    phyoutapt = delta * (size(u2)-1);
    
    function[u2, phyoutapt] = tool_symmetry(u1, phyinapt)
        % % Input: 
        % u1: input complex field
        % phyinapt: Input physical aperture dimension
        % % Output: 
        % u2: output complex field
        % phyoutapt: Output physical aperture dimension

        % This program input x-y sampling interval is the same

        [M, N] = size(u1);
        delta = phyinapt(1)/(M-1);

        if(M > N)    
            square_pad = ceil(0.5 * (M - N));
            u2 = padarray(u1, [0,square_pad],0);
            phyoutapt = delta.*size(u2);
        elseif (M < N)
            square_pad = ceil(0.5 * (N - M));
            u2 = padarray(u1, [square_pad,0],0);
            phyoutapt = delta.*size(u2);
        else
            u2 = u1;
            phyoutapt = phyinapt;
        end

    end
end


% RSD convolution Focusing
function[u2] = propRSD_conv(u1, L,lambda, z)
	% RSD in the convolution foramt
    
	% assumes same x and y side lengths and uniform sampling
	% u1 - source plane field
	% L - source and observation plane side lengths
	% lambda - wavelength
	% z - propagation distance
	% u2 - observation plane field
    
    [M,N] = size(u1);			% get input field array size, physical size

    % Extract each physical dimension
    l_1 = L(1);
    l_2 = L(2);
     
    % Spatial sampling interval
    dx = l_1/(M-1);						% sample interval x direction
    dy = l_2/(N-1);						% sample interval y direction

    % spatial sampling resolution needs to be equal in both dimension
    z_hat = z/dx;                   
    mul_square = lambda * z_hat/(M * dx);

    % center the grid coordinate
    m = linspace(0, M-1, M); m = m - M/2;
    n = linspace(0, N-1, N); n = n - N/2;
    
    [g_m, g_n] = meshgrid(n,m);  % coordinate mesh			

    %  Convolution Kernel Equation including the firld drop off term
    h = exp(1j * 2 * pi * (z_hat.^2 * sqrt(1 + (g_m.^2/z_hat.^2 + g_n.^2./z_hat.^2)) ./(mul_square * M)))./sqrt(1 + g_m.^2/z_hat^2 + g_n.^2/z_hat^2);  
    
    % Convolution or multiplication in Fourier domain
    H = fft2(h);
    U1 = fft2(u1);
    U2 = U1.*H;
    u2 = ifftshift(ifft2(U2));
end

end

%{
    Parameters input: 
        1. uin: captured virtual wavefront, unitless
        2. apt_in: captured physical aperture sizem, uniter meter (array x-y)
        3. maxSz: virtual aperture physical maximum size, unit meter (real valed scalar)
        4. posL: virtual point source location on the aperture (origional
        physical captured) index [x y] unitless

    Parameters output: 
        1. uout: output virtual wavefront, square size, unitless
        2. apt_out: output virtual aperture physical maximum size, unit meter (array x-y)
        3. posO: output new virtual point source location on the apeture,
        index [x y] unitless

    Description:
        This function works for creating virtual aperture, independent from
        RSD propagation. Input can be symmetry or not, output dimension
        follow by the input sampling resolution and the maxSz requirement.
        
        Padding process refers to the center of the captured wavefront
        center.
%}
function [uout, apt_out, posO] = Create_VirtualAperture(uin, apt_in, maxSz, posL)
    % Get input wavefront dimension
    [M, N] = size(uin);
    
    % Calculate the spatial sampling density, unit meter
    delta = apt_in(1)/(M-1);
    
    % Calculate the ouptut size
    M_prime = round(maxSz/delta);
    
    % symmetry square padding
    if(M == N)

    else
        [uin, ~] = tool_symmetry(uin, apt_in);
    end
    
    % Padding size difference
    diff = size(uin) - [M, N] ;
    posL = posL + diff/2; % update the virtual point source location after symmetry
    
    % Update the symmetry aperture size
    [M, ~] = size(uin);
    
    % Virtual aperture extented padding on the boundary
    % difference round to nearest even number easy symmetry padding
    dM = 2 * ceil((M_prime - M)/2);
    
    % symmetry padding on the boundary
    if(dM > 0)
        % Using 0 boundary condition
        uout = padarray(uin, [dM/2 dM/2], 0, 'both');
        posO = posL + dM/2; % update the virtual point source location after padding
    else
        uout = uin;
        posO = posL;
    end
    
    % update the virtual aperture size
    apt_out = delta * size(uout);
    
    
function[u2, phyoutapt] = tool_symmetry(u1, phyinapt)
    % % Input: 
    % u1: input complex field
    % phyinapt: Input physical aperture dimension
    % % Output: 
    % u2: output complex field
    % phyoutapt: Output physical aperture dimension

    % This program input x-y sampling interval is the same

    [M, N] = size(u1);
    delta = phyinapt(1)/(M-1);

    if(M > N)    
        square_pad = ceil(0.5 * (M - N));
        u2 = padarray(u1, [0,square_pad],0);
        phyoutapt = delta.*size(u2);
    elseif (M < N)
        square_pad = ceil(0.5 * (N - M));
        u2 = padarray(u1, [square_pad,0],0);
        phyoutapt = delta.*size(u2);
    else
        u2 = u1;
        phyoutapt = phyinapt;
    end

end
    
end


