function hf=ShadePatchDepths(p1,p2,xl)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Shade depth regions of patches (overturns) on a plot
%
%
%
%----------------
% 12/10/16 - A.P.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

Np=length(p1); % number of patches

for ip=1:Np
    
    hf=fill([xl(1) xl(2) xl(2) xl(1)],[p1(ip) p1(ip) p2(ip) p2(ip)],0.75*[1 1 1],'edgecolor','none');
    hold on
end




%%