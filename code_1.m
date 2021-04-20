
% 2020-07-08 --- Atlanta
% 2021-03-03 --- Shanghai
% POTF shape Given & Solve, based on numerical Matrix | least square method

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

% spanned matrix of a large number of basic POTFs
Omega = zeros(M, N);   
Tu{N} = zeros( Ma, Mb );  % basic uniform POTF

% generate matrix Omega
for ii = 1 : N
    
    k = ii/N * rho_c;

    Tu{ii} = imag( getPOTFnP_Circular_Zero( rhom, etam, rho_o, k, 1 ));    % basic uniform POTF

    Omega(:, ii) = reshape( Tu{ii}, M, 1);     % obtain the Omega matrix
    
end


% Mark the all-zero column or all-zero row
Find_zero_column = sum(abs(Omega),1) == 0;
Find_zero_row = sum(abs(Omega),2) == 0;

% ensure matrix Omega full-rank
Omega( Find_zero_row, :) = [];
Omega( :, Find_zero_column ) = [];

% POTF target
t0 = ones(length(Omega), 1);

% solve the optimal illumination profile
S = [];
s = (Omega'*Omega) \ (Omega'*t0);
t = Omega * s;
for ii = length(s) : -1 : 1
    S(ii) = sum(s(ii:end));
end
S = S/max(max(S));

% add back the zero-point and calculate the optimized POTF
temp = zeros(M, 1);
temp( ~Find_zero_row ) = t;
T = reshape( temp, Ma, Mb );

% calculate the RMS between the target POTF and the optimized POTF
enorm = sqrt( sum(abs(t - 1).^2) / length(t) );
disp([N enorm]);

% show the optimized POTF in uniform coordinate
[eta0, rho0] = meshgrid(b(1:2:end), a);          % uniform coordinate    
img = interp2( etam, rhom, T, eta0, rho0);
figure(1);
imshow( flipud(img), [0, 2.5] );
title( num2str([N enorm]), 'fontsize', 12 );
axis equal; colormap hot;
pause(0.1);

% show the optimized illumination profile
figure(2);
K = rho_c/N : rho_c/N : rho_c;
K( Find_zero_column ) = [];
plot( K, S, '.-', 'linewidth', 0.5, 'markersize', 16 );
xlim([0, rho_o]);
ylim([0, 1]);
title( num2str([N enorm]), 'fontsize', 12 );
pause(0.1);

% generate the relaxed optimal illumination
[X, Y] = meshgrid( -2 : 4/256 : 2, -2 : 4/256 : 2 );
Source = zeros(257);
for ii = 1 : 257
    for jj = 1 : 257
        R = sqrt( X(ii,jj)^2 + Y(ii,jj)^2 );
        if  R <= rho_c
            Source( ii, jj ) = interp1( [0 K], [1 S], R, 'linear' );
        end
    end
end
figure(3);
imshow( Source, [] );
axis equal; colormap gray;
pause(0.1);
