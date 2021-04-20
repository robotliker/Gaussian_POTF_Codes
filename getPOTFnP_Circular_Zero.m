

% For Example
% [eta, rho] = meshgrid( -2:0.02:2, -2:0.01:2 );
% a = getPOTFnP_Circular( rho, eta, 0.9, 0.3 );

% 
function OTF = getPOTFnP_Circular_Zero( rho, eta, cc, cu, p0 )


%---% normalized illumination 
% Esum = pi/2 * ( 2 * cu^2  -  c * cu^4 );
%---%

%---%
a = 1 + rho.^2 ./ eta.^2;
b = 1 - ( rho.^2 + eta.^2 ) / 4;
%---%
ec = sqrt( 1 - cc^2 );
eu = sqrt( 1 - cu^2 );
xx = -sqrt( b ./ a .* ( b > 0 ) );
%---%

%---%
flag = eta > 0;
eta = abs( eta ); 
rho = abs( rho );
%---%

%---%
m1 = 2 * (a - 1) ./ eta ./ a ./ sqrt( a );
m1( rho == 0 ) = 0;
m2 = -eta ./ sqrt( a ) / 2;
m2( isnan( m2 ) ) = 0;
%---%


%---%
h0 = p0;
%---%


X1a = eta ./ rho .* ( -eta/2 - ec );
X1b = eta ./ rho .* ( eta/2 - eu );
X1 = min( X1a, X1b );

X2a = eta ./ rho .* ( -eta/2 - ec );
X2b = eta ./ rho .* ( -eta/2 - eu );
X2 = min( X2a, X2b );

Y1 = real( sqrt( b - a .* X1.^2 ) );
Y2 = real( sqrt( b - a .* X2.^2 ) );
    
      
%---% calcuate the POTF
theta1 = asin( Y1 ./ sqrt( b ) );
theta2 = asin( Y2 ./ sqrt( b ) );

%---% calulate the POTF contributed by the left circle

%---%
s2t1 = sin( 2*theta1 );
s2t2 = sin( 2*theta2 );
%---%
u01 = b/4 .* ( 2*theta1 + s2t1 );
u02 = b/4 .* ( 2*theta2 + s2t2 );
v01 = theta1;
v02 = theta2;
%---%

%---%
S1 = m1 .* ( h0.* u01 ) ...
     + m2 .* ( h0.* v01 );
 %---%
S2 = m1 .* ( h0.* u02 ) ...
     + m2 .* ( h0.* v02 );
%---%


AOTF =  S1 + S2;
          
POTF =  S1 - S2;

%---% calculate the POTF contributed by the right circle

%---%
POTF = POTF .* flag - ( 1 - flag ) .* POTF;
POTF( isnan( POTF ) ) = 0;
AOTF( isnan( AOTF ) ) = 0;

%---%
OTF = AOTF + 1i * POTF;
OTF = OTF / pi / 4;
% OTF = OTF ./ Esum;
    
end


