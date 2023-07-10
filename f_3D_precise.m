function u_volume = f_3D_precise(u_total, M, N, lambda_loop, weight, sample_spacing,camera_pos,laser_pos)

depth_min = 0.8; % depth minimal volumn
depth_max = 1.8; % depth maximum volumn
depth_resolution = 0.02;
depth_loop = depth_min: depth_resolution: depth_max;
u_volume = zeros([M, N,length(depth_loop)]);

m = linspace(0, M-1, M); m = (m - (M-1)/2)*sample_spacing;
n = linspace(0, N-1, N); n = (n - (N-1)/2)*sample_spacing;
[p_x, p_y] = meshgrid(n,m);

%% Pad the original plane
depth1 = 0.8;
for spectrum_index = 1 : length(lambda_loop)
    lambda = lambda_loop(spectrum_index);
    u1 = u_total(:, spectrum_index)*weight(spectrum_index);
    u_mid = propRSD_3D(u1, M, N, camera_pos,sample_spacing,lambda, depth1);
    
    % 2nd propogation
    for i = 1 : length(depth_loop)
        depth = depth_loop(i)-depth1;
        u2_RSD_conv = RSD(u_mid, sample_spacing, lambda, depth);
        r1 = sqrt((p_x - laser_pos(1)).^2 + (p_y - laser_pos(2)).^2 + (depth_loop(i) - laser_pos(3))^2);
        p1 = exp(1j*2*pi/lambda.*r1);
        u2_RSD_conv = u2_RSD_conv.* p1;
        u_volume(:,:,i) = u_volume(:,:,i) + u2_RSD_conv;
    end
    
end
end


%% Alternative functions
function[u2] = propRSD_3D(u1,M, N, camera_pos,sample_spacing,lambda, z)
    % center the grid coordinate
    m = linspace(0, M-1, M); m = (m - M/2)*sample_spacing;
    n = linspace(0, N-1, N); n = (n - N/2)*sample_spacing;
    
    [g_m, g_n] = meshgrid(n,m);  % coordinate mesh
    
    xc = camera_pos(:,1);
    yc = camera_pos(:,2);
    zc = camera_pos(:,3);
    
    u2 = zeros(M, N, size(u1,2));
    for x = 1:M
        for y = 1:N
            xv = g_m(x,y);
            yv = g_n(x,y);
            zv = z;
            r = sqrt((xv-xc).^2 + (yv-yc).^2 + (zv-zc).^2);
            RSD = exp(1j * 2*pi/lambda * r)./r;
            u2(x,y) = sum(u1.*RSD, 'all');
        end
    end
end

function uout = RSD(uin, sample_spacing, lambda, depth)
    [M,N] = size(uin);
    m = linspace(0, M-1, M); m = (m - M/2)*sample_spacing;
    n = linspace(0, N-1, N); n = (n - N/2)*sample_spacing;
    [g_m, g_n] = meshgrid(n,m);
    r = sqrt(g_m.^2+g_n.^2+depth^2);
    h = exp(1j*2*pi/lambda*r);
    H = fft2(h);
    U1 = fft2(uin);
    U2 = U1.*H;
    uout = ifftshift(ifft2(U2));
end

function[uout] = Camera_Focusing(uin, sample_spacing, lambda, depth)
u1_prime = uin;
uout = propRSD_conv(u1_prime, sample_spacing, lambda, depth);

% RSD convolution Focusing
    function[u2] = propRSD_conv(u1, dx,lambda, z)
        
        
        [M,N] = size(u1);			% get input field array size, physical size
        
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

function [] = flip_and_show(img, title_str)
img = squeeze(img);
img = imresize(img,2);
img = flip(img,1);
img = flip(img,2);

figure;
imagesc(img);
colormap 'hot';
title(title_str);
axis image;
axis off;
colorbar;
end