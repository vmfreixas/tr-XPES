function doubfac = doubfac( n )
% calculate n!!

% Initialize
%doubfac = 1;                         

%while n>0
%    doubfac = doubfac*n;
%    n = n-2;
%end

doubfac = prod(2*n-1:-2:1);

end

