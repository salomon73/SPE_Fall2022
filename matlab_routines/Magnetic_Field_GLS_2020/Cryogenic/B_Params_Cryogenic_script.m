function void = B_Params_script(magnet,Ips)
%
% Small script to run the more general B_Params_170 which returns beam
% parametersog the 170GHz gyrotron tube.
% The only difference is that the currents delivered by the power supplies
% are input here. They are rearranged to run B_Params_170.
% Ips = [Ig1, Ig2, Imain, Iadd1, Iadd2]
%
% Example: 
% >> magnet = 'asg_modified_201206'
% >> Ips    = [9.7, 7.7, 8.6, 3.54,87.47]
% >> B_Params_script(magnet,Ips)
%
%   J.-P. Hogge, 12.003.2008
%   Version: 1.0

    I = [Ips(1) Ips(2) (Ips(3)+Ips(5))  Ips(4)+Ips(5)  Ips(5) ]
    B_Params_170(magnet,I)