
%==%
function Z = getCircle( N, x0, y0, r0 )

Z = zeros( N );
x0 = fix( N * x0 );
y0 = fix( N * y0 );
r0 = N*r0/2;
for i = 1:N;
    for j = 1:N;
        if (i-x0).^2 + (j-y0).^2 < r0.^2;
            Z(i,j) = 1;
        end
    end
end
end
