

% For Example
% [eta, rho] = meshgrid( -2:0.02:2, -2:0.01:2 );
% a = getPOTFnP_Circular( rho, eta, 0.9, 0.3 );

% 
function OTF = getPOTFnP_Circular_Eight( rho, eta, cc, cu, p0, p1, p2, p3, p4 )


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
a2 = a.^2;
a3 = a.^3;
a4 = a.^4;
%---%
v1 = b ./ a;
v2 = v1.^2;
v3 = v1.^3;
v4 = v1.^4;
%---%
r1 = rho;
r2 = rho.^2;
r3 = rho.^3;
r4 = rho.^4;
r5 = rho.^5;
r6 = rho.^6;
r7 = rho.^7;
r8 = rho.^8;
%---%
c1 = sqrt( a );
c2 = a;
c3 = a.^( 3/2 );
c4 = a2;
c5 = a.^( 5/2 );
c6 = a3;
c7 = a.^( 7/2 );
c8 = a4;
%---%
d1 = a - 1;
d2 = d1.^2;
d3 = d1.^3;
d4 = d1.^4;
%---%
h0 = p0 +...
         p1 * ( v1 + r2 / 4 ) +...
         p2 * ( v2 + 3/2 * v1 .* r2 + r4 / 16 ) +...
         p3 * ( v3 + 15/4 * v2 .* r2 + 15/16 * v1 .* r4 + r6 / 64 ) +...
         p4 * ( v4 + 7 * v3 .* r2 + 35/8 * v2 .* r4 + 7/16 * v1 .* r6 + r8 / 256 );
%---%
h1 = p1 * ( r1 ) +...
         p2 * ( 2 * v1 .* r1 + r3 / 2 ) +...
         p3 * ( 3 * v2 .* r1 + 5/2 * v1 .* r3 + 3/16 * r5 ) +...
         p4 * ( 4 * v3 .* r1 + 7 * v2 .* r3 + 7/4 * v1 .* r5 + r7 / 16 );
h1 = h1 ./ -c1;
%---%
h2 = p1 * d1 +...
         p2 * ( 2 * v1 .* d1 + 1/2 * r2 .* ( a - 3 ) ) +...
         p3 * ( 3 * v2 .* d1 + 3/2 * v1 .* r2 .* ( 3*a - 5 ) + 3/16 * r4 .* ( a - 5 ) ) +...
         p4 * ( 4 * v3 .* d1 + 3 * v2 .* r2 .* ( 5*a - 7 ) + 5/4 * v1 .* r4 .* ( 3*a - 7 ) + 1/16 * r6 .* ( a - 7 ) );
h2 = h2 ./ c2;
%---%
h3 = p2 * ( 2 * r1 .* d1 ) +...
         p3 * ( 6 * v1 .* r1 .* d1 + 1/2 * r3 .* ( 3*a - 5 ) ) +...
         p4 * ( 12 * v2 .* r1 .* d1 + 2 * v1 .* r3 .* ( 5*a - 7 ) + 1/4 * r5 .* ( 3*a - 7 ) );
h3 = h3 ./ -c3;
%---%
h4 = p2 * d2 +...
         p3 * ( 3 * v1 .* d2 + ( 3/4 * ( a - 3 ).^2 - 3 ) .* r2 ) +...
         p4 * ( 6 * v2 .* d2 + ( ( 3*a - 5 ).^2 - 4 ) .* v1 .* r2 + ( 3/8 * ( a - 5 ).^2 - 5 ) .* r4 );
h4 = h4 ./ c4;
%---%
h5 = p3 * ( 3 * r1 .* d2 ) +...
         p4 * ( 12 * v1 .* r1 .* d2 + ( 3*a2 - 10*a + 7 ) .* r3 );
h5 = h5 ./ -c5;
%---%
h6 = p3 * d3 +...
         p4 * ( 4 * v1 .* d3 + ( a.^3 - 9*a2 + 15*a - 7 ) .* r2 );
h6 = h6 ./ c6;
%---%
h7 = p4 * ( 4 * r1 .* d3 );
h7 = h7 ./ -c7;
%---%
h8 = p4 * d4;
h8 = h8 ./ c8;
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
s4t1 = sin( 4*theta1 );
s4t2 = sin( 4*theta2 );
s6t1 = sin( 6*theta1 );
s6t2 = sin( 6*theta2 );
s8t1 = sin( 8*theta1 );
s8t2 = sin( 8*theta2 );
s10t1 = sin( 10*theta1 );
s10t2 = sin( 10*theta2 );
%---%
u01 = b/4 .* ( 2*theta1 + s2t1 );
u02 = b/4 .* ( 2*theta2 + s2t2 );
v01 = theta1;
v02 = theta2;
%---%
u11 = b .* Y1 - Y1.^3 / 3;
u12 = b .* Y2 - Y2.^3 / 3;
v11 = Y1;
v12 = Y2;
%---%
u21 = b.^2 / 32 .* ( 4*theta1 - s4t1 );
u22 = b.^2 / 32 .* ( 4*theta2 - s4t2 );
v21 = b .* ( 2*theta1 - s2t1 ) / 4;
v22 = b .* ( 2*theta2 - s2t2 ) / 4;
%---%
u31 = b .* Y1.^3 / 3  -  Y1.^5 / 5;
u32 = b .* Y2.^3 / 3  -  Y2.^5 / 5;
v31 = Y1.^3 / 3;
v32 = Y2.^3 / 3;
%---% 
u41 = b.^3 / 192 .* ( 12 * theta1 - 3 * s2t1 - 3 * s4t1 + s6t1 );
u42 = b.^3 / 192 .* ( 12 * theta2 - 3 * s2t2 - 3 * s4t2 + s6t2 );
v41 = b.^2 / 32   .* ( 12 * theta1 - 8 * s2t1 + s4t1 );
v42 = b.^2 / 32   .* ( 12 * theta2 - 8 * s2t2 + s4t2 );
%---%
u51 = b/5 .* Y1.^5  -  Y1.^7 / 7;
u52 = b/5 .* Y2.^5  -  Y2.^7 / 7;
v51 = Y1.^5 / 5;
v52 = Y2.^5 / 5;
%---%
u61 = b.^4 / 1024 .* ( 40 * theta1 - 16 * s2t1 - 8 * s4t1 + 16/3 * s6t1 - s8t1 );
u62 = b.^4 / 1024 .* ( 40 * theta2 - 16 * s2t2 - 8 * s4t2 + 16/3 * s6t2 - s8t2 );
v61 = b.^3 / 192   .* ( 60 * theta1 - 45 * s2t1 + 9 * s4t1 - s6t1 );
v62 = b.^3 / 192   .* ( 60 * theta2 - 45 * s2t2 + 9 * s4t2 - s6t2 );
%---%
u71 = b/7 .* Y1.^7  -  Y1.^9 / 9;
u72 = b/7 .* Y2.^7  -  Y2.^9 / 9;
v71 = Y1.^7 / 7;
v72 = Y2.^7 / 7;
%---%
u81 = b.^5 / 10240 .* ( 280 * theta1 - 140 * s2t1 - 40 * s4t1 + 130/3 * s6t1 - 15 * s8t1 + 2 * s10t1 );
u82 = b.^5 / 10240 .* ( 280 * theta2 - 140 * s2t2 - 40 * s4t2 + 130/3 * s6t2 - 15 * s8t2 + 2 * s10t2 );
v81 = b.^4 / 1024   .* ( 280 * theta1 - 224 * s2t1 + 56 * s4t1 - 32/3 * s6t1 + s8t1 );
v82 = b.^4 / 1024   .* ( 280 * theta2 - 224 * s2t2 + 56 * s4t2 - 32/3 * s6t2 + s8t2 );
%---%


%---%
S1 = m1 .* ( h0.* u01 + h1 .* u11 + h2.* u21 + h3.* u31 + h4.*u41 + h5.* u51 + h6.* u61 + h7.* u71 + h8.* u81 ) ...
     + m2 .* ( h0.* v01 + h1 .* v11 + h2.* v21 + h3.* v31  + h4.* v41+ h5.* v51 + h6.* v61  + h7.* v71 + h8.* v81 );
 %---%
S2 = m1 .* ( h0.* u02 - h1 .* u12 + h2.* u22 - h3.* u32 + h4.*u42  - h5.* u52 + h6.* u62 - h7.* u72 + h8.* u82 ) ...
     + m2 .* ( h0.* v02 - h1 .* v12 + h2.* v22 - h3.* v32  + h4.* v42 - h5.* v52 + h6.* v62  - h7.* v72 + h8.* v82 );
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


