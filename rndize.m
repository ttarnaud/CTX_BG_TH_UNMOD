function Out = rndize(A,refl)
% refl = 1: don't introduce self-connections during randomization
if (nargin == 1), refl = 0; end
% Function permutes columns of A
Out = zeros(size(A));
if (~refl)
for j=1:size(A,2)
    temp = A(:,j);
    Out(:,j) = temp(randperm(size(A,1)));  
end
else
for j=1:size(A,2)
    temp = [A(1:j-1,j);A(j+1:end,j)];
    temp = temp(randperm(size(temp,1)));
    Out(:,j) = [temp(1:j-1); 0; temp(j:end)];  
end    
end
end