function Aeff = ModeArea(x,y,Field)

%% Description : 
%% This function computes the effective area of the considered optical field, with : $$ A_{eff} =  {\Big({\int{\int{\vert E\vert^{2}}} dxdy\Big)^{2}}\over{\int{\int{\vert E\vert^{4}dxdy}}}}\ \ \ \ \ \ \ \ \ \  \textbf{(3)}$$
%% Inputs : 
%% * x, y : Transverse coordinates of the simulation (m)
%% * Field : Matrix of the optical field of the considered mode (must be complex !)
%%
%% Outputs : 
%% * Aeff : Effective mode area (m²)

Aeff = (trapz(x,trapz(y,abs(Field).^2)).^2)./((trapz(x,trapz(y,abs(Field).^4))));