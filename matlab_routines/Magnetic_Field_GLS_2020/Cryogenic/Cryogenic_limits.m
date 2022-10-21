% Cryogenic_limits.m
% Defines the safe ranges of currents that can be used with Cryogenic
% magnet.
%
    cryogenic_limits = zeros(4,2);
    cryogenic_limits(:,1) = 0;
    
    cryogenic_limits(1,2) =  81.30;
    cryogenic_limits(2,2) =  89.40;
    cryogenic_limits(3,2) = 117.55;
    cryogenic_limits(4,2) = 114.80;