function d = divisor(n)
%   Check whether input is positive integer and scalar.
if ~isscalar(n)
    error('divisor:NonScalarInput','Input must be a scalar.');
end
if (n < 1) || (floor(n) ~= n)
  error('divisor:PositiveIntegerOnly', 'Input must be a positive integer.'); 
end

% Find prime factors of number :
pf = factor(n);         % Prime factors of n
upf = unique(pf);       % Unique

% Calculate the divisors
d = upf(1).^(0:1:sum(pf == upf(1)))';
for f = upf(2:end)
    d = d*(f.^(0:1:sum(pf == f)));
    d = d(:);
end
d = sort(d)'; 

end