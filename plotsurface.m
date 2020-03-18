function plotsurface(surf,coloring,alpha)  %alpha为透明度
% function plotsurface(surf,coloring,alpha)  函数plotsurface(surf、着色、α)
% Plot surface (face & vertex list) with given coloring 利用(面和顶点列表)进行绘制

% default: uniform coloring 
if ~exist('coloring','var')
    coloring=ones(length(surf.vertices),1);  %ones(m,n)    生成一个m行n列且所有元素都是1的矩阵  这里是列向量
end

% default: facealpha = 0.5 透明度
if ~exist('alpha','var');
    alpha=0.5;
end

patch('Faces',surf.faces,'Vertices',surf.vertices,'FaceVertexCData',coloring,'FaceColor','interp','edgecolor', 'interp','FaceAlpha',alpha); 
%patch是个底层的图形函数，用来创建补片图形对象。一个补片对象是由其顶点坐标确定的一个或多个多边形。用户可以指定补片对象的颜色和灯光


xlabel('x'); ylabel('y'); zlabel('z'); view(3); daspect([1 1 1]);

shading flat;

end



% despect 每个轴的单位长度为1比1比1，也就是相同的
%%%这个是绘制图像的函数



