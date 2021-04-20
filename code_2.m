
% 2020-09-15 --- Atlanta
% 2021-03-03 --- Shanghai
% Optimal Gaussian Illumination

% settings of optical cofiguration
rho_o = 1;  
rho_c = 1;

% settings of coordinate size 
Ma = 500;
Mb = 500;
M = Ma*Mb;
N = 50;

% non-uniform coordinate
a = 2/Ma : 2/Ma : 2;
b = 1/Mb : 1/Mb : 1;
[etam, rhom] = meshgrid( b.^2, a );    % non-uniform coordinate 

% settings of Gaussian width
sigma = 3 : 0.1 : 4;
Leng = length(sigma);
rms = ones(1,Leng);

% calculate and evaluate Gaussian POTF for different sigma values
for ii = 1 : Leng
    
    r = -rho_c : 0.01 : rho_c;
    y = exp( sigma(ii) * r.^2 );    % Gaussian curve
    p = polyfit( r, y, 8 );     % Gaussian curve using 8-th polynomial approximation
    
    % calculate Gaussian POTF
    G = getPOTFnP_Circular_Eight( rhom, etam, rho_o, rho_c, p(9), p(7), p(5), p(3), p(1) );
    G = imag( G );
       
    G = abs( G );   % add absolute operation to the POTF
       
    % remove zero points that are outside the Spatial frequency coverage
    G(G == 0) = [];
    
    G = G * sum(sum(G)) / sum(sum(G.^2));    % self-adaptive range
    
    rms(ii) = sqrt( sum(sum(abs(G - 1).^2)) / length(G(:)) );    % calculate RMS value
        
end

% show the rms versus sigma
figure(1);
enorm = min(rms);
best_sigma = sigma(rms == enorm);
plot( sigma, rms, '.-', 'linewidth', 0.5, 'markersize', 16 );
title( num2str([best_sigma enorm]), 'fontsize', 12 );
pause(0.1);

% show the optimized Gaussian POTF in uniform coordinate
[eta0, rho0] = meshgrid(b(1:2:end), a);          % uniform coordinate  
y = exp( best_sigma * r.^2 );    % Gaussian curve
p = polyfit( r, y, 8 );              % Gaussian curve using 8-th polynomial approximation
T = getPOTFnP_Circular_Eight(rho0, eta0, rho_o, rho_c, p(9), p(7), p(5), p(3), p(1));
T = flipud(imag(T));
G = abs(T);
G(G == 0) = [];
T = T * sum(sum(G)) / sum(sum(G.^2));
figure(2);
imshow( abs(T), [0, 2.5] );
title( num2str([best_sigma enorm]), 'fontsize', 12 );
axis equal; colormap hot;
pause(0.1);

% generate the relaxed optimal illumination
[X, Y] = meshgrid( -2 : 4/256 : 2, -2 : 4/256 : 2 );
Source = zeros(257);
for ii = 1 : 257
    for jj = 1 : 257
        R = sqrt( X(ii,jj)^2 + Y(ii,jj)^2 );
        if  R <= rho_c
            Source( ii, jj ) = exp( best_sigma * R.^2 );
        end
    end
end
figure(3);
imshow( Source, [] );
axis equal; colormap gray;
pause(0.1);
