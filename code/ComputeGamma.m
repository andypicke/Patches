function gam=ComputeGamma(n2,dtdz,chi,eps)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
% 10/27/16 - APickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
gam =  n2 .* chi ./2 ./ eps ./ (dtdz.^2);

return
%%
