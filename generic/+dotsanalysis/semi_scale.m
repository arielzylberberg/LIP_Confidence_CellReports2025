function [um,ym] = semi_scale(uscoh,y,ischoice_flag)

if nargin<3 || isempty(ischoice_flag)
    ischoice_flag = false;
end

u = uscoh;
I = u<0;
u(I) = -1 * u(I);
um = unique(u);

if (ischoice_flag)
    y(I) = 1 - y(I);
end
[um,ym] = curva_media(y,u,[],0);