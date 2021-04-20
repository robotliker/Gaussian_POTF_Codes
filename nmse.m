

function nmse = nmse(A,B,C)

flag = isnan(A) | isnan(B);
A( flag ) = [];
B( flag ) = [];

if strcmp(C,'offset')
    E = (B-A) - mean(mean(B-A));
else
    E = B-A;
end

D = A - mean(mean(A));

nmse = sum(sum(E.^2)) / sum(sum(D.^2));

end

