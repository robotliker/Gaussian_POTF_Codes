
% 2020-08-20 | Tomography
% modified at 2020-12-07, Atlanta
% 2021-03-03 Shanghai

% 2D Version with rotation Tomography

% function [PSF_2D, OTF_2D] = POTF_2D_in_TDPM_SSBPM (OTF_type, source_type, Leng, ...
%     object_2D,noil,useOF,NAo,NAc,SAMPLING_RATE,lambda,NAci)

% performs a simulation on a central line scatterer to generate the cross-sectional intensity (idatapsf)
% Based on SS-BPM described by Eq. (23) and (24) in Jenkins_2015b.

% PSF_2D is the intensity point spread function
% OTF_2D is the optical transfer function.
% (if OTF_type=='A', this will be AOTF; if OTF_type=='P', this will be POTF;
% if OTF_type=='O', this will be the simulation of an object (phantom), and only PSF_2D is meaningful);

% OTF_type={'P','A','O'} indicated what type of object is used.
% 'P' uses a point RI object, so the output OTF_2D is POTF.
% 'A' uses a point absorptive object, so the output OTF_2D is AOTF.
% 'O' uses object_2D as the object, so the output PSF_2D is 3D intensity scattered by the ojbect.

% source_type={'Gaus','disk','annular'} select the type of source function.
% Leng is the size of POTF and PSF needed.
% Lengs=[Lengs1,Lengs2] is the size of source. It can be different from the size of object.
% SAMPLING_RATE is the sample distance of the camera.
% SAMPLING_RATEs is the sample distance of the source. Usually SAMPLING_RATEs==SAMPLING_RATE

% Leng is the size of POTF and PSF needed.
% useof={'1','OF','OF2'} indicate whether use (inverse) OF or not


% ---------------------------------------------------------------------------------------- %
% Source_Optimal_FigX means relaxed optimal illumination
% ---------------------------------------------------------------------------------------- %


clear;

load( 'USAF_S.mat' );

sample = USAF_S;     % select the circular object
RIrange = 0.002;     % RI range of the object
TomoAngle = 0 : 60 : 120;   % angles for tomography

OTF_type = 'P';
useOF = '1';

Leng = 257;   % To make a good alignment, this value must be odd, such as 65 or 255
sample = imresize( sample, [Leng Leng], 'nearest' );

noil = 1;
NAo = 0.7;
NAc = 0.5;
NAci = 0;
lambda = 550 * 1e-9;
SAMPLING_RATE = lambda / 4;
GaussianFactor = 'Source_Optimal_Fig9';      % This value can be 'Source_Optimal_Fig6'
% GaussianFactor = 3.6;     % This value can be set arbitrary

k0 = 2*pi / lambda;                % wavenumber in the air
df = 1 / ( Leng * SAMPLING_RATE );

if mod( Leng, 2 ) == 0
    t = Leng / 2;
    [ fX , fY ] = meshgrid( -t : t-1 , -t : t-1 );
else
    t = (Leng-1) / 2;
    [ fX , fY ] = meshgrid( -t : t , -t : t );
end

fX = fX * df;
fY = fY * df;
rho = sqrt( fX.^2 + fY.^2 );         % spatial frequency of the light
dz = SAMPLING_RATE * 1;     % defocus distance

mask_obj = double( rho <= ( NAo / lambda ) ) ;       % objective lens


% Select the source distribution function
r = rho ./ ( 1/lambda );


%---% Gaussian
if isnumeric( GaussianFactor )
    source = double( rho <= NAc / lambda  &  rho >= NAci / lambda );
    source = source .* exp( GaussianFactor * r.^2 );
end

%---% User Defined
if ~isnumeric( GaussianFactor )
    load( [ GaussianFactor  '.mat' ] );
    source( 129, : ) = [];
    source( :, 129 ) = [];
end
% size matching
source = imresize( source, [Leng Leng], 'nearest' );

source = source ./ max(max( source));
source_curve = source( fix( Leng/2 ) + 1, : );

%--------------------------------------------------------------------------------------------------%

delta_n = 1e-6;
object_delta = noil * ones( Leng, Leng );        % RI distribution in x-z plane.
object_delta( t + 1, t + 1 ) = noil + delta_n;

[ A, B ] = find( abs( source ) > eps );
% A is the row number and B is the column number of a source point

%--------------------------------------------------------------------------------------------------%

idata_Source = zeros( Leng, Leng );

% for processing the source: row by row
for aa = 1 : length( A )
    % Modified Split-Step Beam Propagation Method, described in Jenkins_2015b
    fprintf( '\b\b\b\b\b\b\b\b%8d' , length( A ) - aa );
    idata2 = zeros( Leng , Leng );
    
    % incident light
    % don't understand
    templ = ifft( circshift( fft(ones(Leng, 1)) , A(aa) - 1 ));  
    

    if strcmp( useOF , 'OF' ) % use receprocal obliquity factor modification
        OF = 1 / cos( asin( rho( A( aa ) , B( aa ) ) / ( noil / lambda ) ) );
        % use receprocal obliquity factor modification and crop at sqrt(2)
    elseif strcmp( useOF , 'OF2' )
        OF = 1 / cos( asin( rho( A( aa ) , B( aa ) ) / ( noil / lambda ) ) );
        OF( OF > sqrt( 2 ) ) = sqrt( 2 );
        % do not use receprocal obliquity factor modification
    elseif strcmp( useOF , '1' )
        OF=1;
    end
    
    %---% forward propagation in the object space
    for z = 1 : Leng
        %---%  propagation by deltaz
        % the spatial frequency along z-axis for a single 
        fz = real( sqrt( (1/lambda)^2 - rho( : , B( aa ) ).^2 ) );
        templ = ifft( fft( templ ) .* exp( 1i * noil * 2*pi * fz * dz  ) );
        %---% applied phase delay
        templ = templ .* exp( 1i * ( k0 * OF * ( object_delta( : , z ) - noil ) * dz ) );
    end

    templ = ifft( fft( templ ).* mask_obj( : , B( aa ) ) .* exp( 1i * noil * 2*pi * fz * dz  ) );
    templ = ifft( fft( templ ) .* mask_obj( : , B( aa ) ) );
    
    % backwrad propagation without object, which represnets propagation in the image space
    for a = Leng : -1 : 1
        templr = ifft( fft( templ ) .* exp( -1i * noil * 2*pi * fz * ( Leng - a ) * dz  ) );
        idata2( : , a ) = abs( templr ) .^ 2;
    end
    
    idata_Source = idata_Source + source( A( aa ) , B( aa ) ) * idata2;
    
end




%---% object figure --- Fig 9
idata_Object{length(TomoAngle)} = zeros(Leng);

% scan the tomography angles
for ii  = 1  :  length( TomoAngle )

    tep = imrotate( sample( 1 : Leng, 1 : Leng ) , TomoAngle( ii ), 'crop', 'bilinear' );
    
    object_2D = tep;
    object = ones( Leng ) * noil + RIrange * object_2D;
    
    idata_Object{ ii } = zeros( Leng, Leng );

    
    for aa = 1 : length( A )
        % Modified Split-Step Beam Propagation Method, described in Jenkins_2015b
        fprintf( '\b\b\b\b\b\b\b\b%8d' , length( A ) - aa );
        idata2 = zeros( Leng , Leng );
        
        % incident light
        templ = ifft( circshift ( fft( ones( Leng, 1 ) ) , ( A( aa ) - 1 )));   % incident light
        
        
        if strcmp( useOF , 'OF' ) % use receprocal obliquity factor modification
            OF = 1 / cos( asin( rho( A( aa ) , B( aa ) ) / ( noil / lambda ) ) );
            % use receprocal obliquity factor modification and crop at sqrt(2)
        elseif strcmp( useOF , 'OF2' )
            OF = 1 / cos( asin( rho( A( aa ) , B( aa ) ) / ( noil / lambda ) ) );
            OF( OF > sqrt( 2 ) ) = sqrt( 2 );
            % do not use receprocal obliquity factor modification
        elseif strcmp( useOF , '1' )
            OF=1;
        end
        
        %---% forward propagation in the object space
        for z = 1 : Leng
            %---%  propagation by deltaz
            fz = real( sqrt( (1/lambda)^2 - rho( : , B( aa ) ).^2 ) );
            templ = ifft( fft( templ ) .* exp( 1i * noil * 2*pi * fz * dz  ) );
            %---% applied phase delay
            tep = object;
            templ = templ .* exp( 1i * ( k0 * OF * ( tep( : , z ) - noil ) * dz ) );
        end
        
        templ = ifft( fft( templ ).* mask_obj( : , B( aa ) ) .* exp( 1i * noil * 2*pi * fz * dz  ) );
        templ = ifft( fft( templ ) .* mask_obj( : , B( aa ) ) );
        
        % backwrad propagation without object, which represnets propagation in the image space
        for a = Leng : -1 : 1
            templr = ifft( fft( templ ) .* exp( -1i * noil * 2*pi * fz * ( Leng - a ) * dz  ) );
            idata2( : , a ) = abs( templr ) .^ 2;
        end
        
        idata_Object{ ii } = idata_Object{ ii } + source( A( aa ) , B( aa ) ) * idata2;
        
    end
    
end

disp( 'end' );



% object RI recovery using POTF method
%---%
tau = 0.02;
number = length( TomoAngle );
fI{number} = zeros(Leng);
fT{number} = zeros(Leng);
RO{number} = zeros(Leng);
Rshow = [];

for ii = 1 : number
    
    Ib = idata_Object{ii};
    Ib_offset = (Ib - mean(mean(Ib)));
    fIb = fftshift(fft2(ifftshift( Ib_offset )));
    
    tep = imrotate( fIb , -TomoAngle( ii ) );
    m = length( tep );
    index = ceil( m/2 ) - ceil( Leng/2 ) + 1 : ceil( m/2 ) - ceil( Leng/2 ) + Leng;
    fI{ii} = tep( index, index );
    
    Is = idata_Source;
    Is_offset = (Is - mean(mean(Is)));
    fIs = fftshift(fft2(ifftshift( Is_offset )));
    tep = imrotate( fIs , -TomoAngle( ii ) );
    fT{ii} = imag( tep( index, index ) );
    
    RO{ii} = delta_n * fI{ii}.*fT{ii} ./ ( fT{ii}.^2 + tau * max(max(fT{ii}.^2)));
    RO{ii} = fftshift( ifft2( ifftshift( RO{ii} )));
    RO{ii} = real(RO{ii}) + imag(RO{ii});
    
    Rshow = [Rshow RO{ii}];
    
end


% RI recovery using tomography
a = 0;
b = 0;
for ii = 1 : number
    a = a + delta_n * fI{ii} .* fT{ii};
    b = b + fT{ii}.^2;
end
r = a ./ (b + tau * max(max(b)));
R = fftshift( ifft2( ifftshift( r )));
R = imag(R) + real(R);
Rshow = [Rshow R];

% show result
%---%
figure(1);
O = sample * RIrange;
set(gca,'position',[0 0 1 1]);
imshow([O R], []);
pause(0.01);
err = nmse(O, R, 'offset');
Info.sigma = GaussianFactor;
Info.nmse = err;
disp( Info );

        
