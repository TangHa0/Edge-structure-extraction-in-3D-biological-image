function labels = segment(edge,sz,edgesize)
% function labels = segment(edge,sz,edgesize)  函数标签
%
% Segment domain of given size into regions enclosed by given edge set 将给定大小的域分割成由给定边集合包围的区域
% using bwlabeln (image processing toolbox)  使用bwlabeln(图像处理工具箱)
%
% Output: 
%   labels: Segmentation label for every voxel 标签:每个体素的分割标签。
%
% Input:
%   edge: Edge set as obtained from edgedetect  边缘:从 edge detect 得到的边缘集。
%   sz: Size of image domain, (3x1)-vector  sz:图像域的大小，(3x1)-矢量。
%   edgesize: Thickness of edge set in voxels  在voxels中设置的边缘厚度
%            (large values help to close holes in the edge set)
%            较大的值有助于关闭边缘中的孔

% default parameter: edges with thickness 1 厚度默认值为1

if ~exist('edgesize','var')
    edgesize=1;
end

% create label space 创建标签空间
labels=ones(sz);

% return if no edge is given  如果没有给出边，则返回
if isempty(edge.vertices)
    return;
end

% pixels corresponding to edge itself 像素匹配边缘本身
edgepix=unique(round(edge.vertices),'rows');
thin_e=zeros(sz); thin_e(sub2ind(sz,edgepix(:,2),edgepix(:,1),edgepix(:,3)))=1;

% fatten edge set 边缘厚度加厚（有助于关闭孔）
fat_e=thin_e;

for i=1:edgesize
    fat_e=grow(fat_e);
end

% set label to 0 on fat edge set  将标签设置为0。
labels(fat_e>0)=0;

% find connected components using bwlabeln (image processing toolbox)  使用bwlabeln(图像处理工具箱)找到连接的组件
labels=bwlabeln(labels,6); 

end

function x_new = grow(x)
% Function to grow level set of 1 by one pixel (in binary image)   函数以1 / 1像素(二进制图像)的形式增长

x_new=x; 

for i=-1:1
    
    for j=-1:1
        
        for k=-1:1
            
            x_new = x_new + circshift(x,[i j k]) + circshift(x,[-i -j -k]);
        end
    end
end

end