function u_volume = f_3D(padding_volume, lambda_loop, weight,pad_size, sample_spacing, depth_spacing, laser_pos)

M = size(padding_volume, 1);
N = size(padding_volume, 2);
depth_min = 0.8; % depth minimal volumn
depth_max = 1.8; % depth maximum volumn
depth_resolution = 0.02;
depth_loop = depth_min: depth_resolution: depth_max;
u_volume = zeros([M, N,length(depth_loop)]);

m = linspace(0, M-1, M); m = (m - (M-1)/2)*sample_spacing;
n = linspace(0, N-1, N); n = (n - (N-1)/2)*sample_spacing;
d = -depth_spacing*pad_size: depth_spacing: depth_spacing*pad_size;
[g_x, g_y, g_z] = meshgrid(n,m,d);  % coordinate mesh
[p_x, p_y] = meshgrid(n,m);

%% Pad the original plane
depth1 = 0.8;
for spectrum_index = 1 : length(lambda_loop)
    lambda = lambda_loop(spectrum_index);
    padding_volume_f = squeeze(padding_volume(:,:,:,spectrum_index));
    
    r = sqrt(g_x.^2 + g_y.^2 + (g_z + depth1).^2);
    G = exp(1j*2*pi/lambda .* r)./r;
    
    fft_pad1 = round(M/2);
    fft_pad2 = round(N/2);
%     padding_volume_f = padarray(padding_volume_f, [fft_pad1, fft_pad2], 0, 'both');
%     G = padarray(G, [fft_pad1, fft_pad2], 0, 'both');

    Pf = fftn(padding_volume_f);
    Gf = fftn(G);
    uout = ifftshift(ifftn(Pf.*Gf));
    u_mid = uout(:, :, pad_size + 1);
%     u_mid = u_mid(fft_pad1+1:M+fft_pad1, fft_pad2+1:N+fft_pad2);
    
    % 2nd propogation
    for i = 1 : length(depth_loop)
        depth = depth_loop(i)-depth1;
        u2_RSD_conv = RSD(u_mid, sample_spacing, lambda, depth);
        r1 = sqrt((p_x - laser_pos(1)).^2 + (p_y - laser_pos(2)).^2 + (depth_loop(i) - laser_pos(3))^2);
        p1 = exp(1j*2*pi/lambda.*r1)./r1;
        u2_RSD_conv = u2_RSD_conv.* p1;
        u_volume(:,:,i) = u_volume(:,:,i) + u2_RSD_conv*weight(spectrum_index);
    end
    
end
end


%% Alternative functions
function uout = RSD(uin, sample_spacing, lambda, depth)
    [M,N] = size(uin);
    m = linspace(0, M-1, M); m = (m - M/2)*sample_spacing;
    n = linspace(0, N-1, N); n = (n - N/2)*sample_spacing;
    [g_m, g_n] = meshgrid(n,m);
    r = sqrt(g_m.^2+g_n.^2+depth^2);
    h = exp(1j*2*pi/lambda*r)./r;
    H = fft2(h);
    U1 = fft2(uin);
    U2 = U1.*H;
    uout = ifftshift(ifft2(U2));
end

