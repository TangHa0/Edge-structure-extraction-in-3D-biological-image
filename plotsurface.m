function plotsurface(surf,coloring,alpha)  %alphaΪ͸����
% function plotsurface(surf,coloring,alpha)  ����plotsurface(surf����ɫ����)
% Plot surface (face & vertex list) with given coloring ����(��Ͷ����б�)���л���

% default: uniform coloring 
if ~exist('coloring','var')
    coloring=ones(length(surf.vertices),1);  %ones(m,n)    ����һ��m��n��������Ԫ�ض���1�ľ���  ������������
end

% default: facealpha = 0.5 ͸����
if ~exist('alpha','var');
    alpha=0.5;
end

patch('Faces',surf.faces,'Vertices',surf.vertices,'FaceVertexCData',coloring,'FaceColor','interp','edgecolor', 'interp','FaceAlpha',alpha); 
%patch�Ǹ��ײ��ͼ�κ���������������Ƭͼ�ζ���һ����Ƭ���������䶥������ȷ����һ����������Ρ��û�����ָ����Ƭ�������ɫ�͵ƹ�


xlabel('x'); ylabel('y'); zlabel('z'); view(3); daspect([1 1 1]);

shading flat;

end



% despect ÿ����ĵ�λ����Ϊ1��1��1��Ҳ������ͬ��
%%%����ǻ���ͼ��ĺ���



