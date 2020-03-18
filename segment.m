function labels = segment(edge,sz,edgesize)
% function labels = segment(edge,sz,edgesize)  ������ǩ
%
% Segment domain of given size into regions enclosed by given edge set ��������С����ָ���ɸ����߼��ϰ�Χ������
% using bwlabeln (image processing toolbox)  ʹ��bwlabeln(ͼ��������)
%
% Output: 
%   labels: Segmentation label for every voxel ��ǩ:ÿ�����صķָ��ǩ��
%
% Input:
%   edge: Edge set as obtained from edgedetect  ��Ե:�� edge detect �õ��ı�Ե����
%   sz: Size of image domain, (3x1)-vector  sz:ͼ����Ĵ�С��(3x1)-ʸ����
%   edgesize: Thickness of edge set in voxels  ��voxels�����õı�Ե���
%            (large values help to close holes in the edge set)
%            �ϴ��ֵ�����ڹرձ�Ե�еĿ�

% default parameter: edges with thickness 1 ���Ĭ��ֵΪ1

if ~exist('edgesize','var')
    edgesize=1;
end

% create label space ������ǩ�ռ�
labels=ones(sz);

% return if no edge is given  ���û�и����ߣ��򷵻�
if isempty(edge.vertices)
    return;
end

% pixels corresponding to edge itself ����ƥ���Ե����
edgepix=unique(round(edge.vertices),'rows');
thin_e=zeros(sz); thin_e(sub2ind(sz,edgepix(:,2),edgepix(:,1),edgepix(:,3)))=1;

% fatten edge set ��Ե��ȼӺ������ڹرտף�
fat_e=thin_e;

for i=1:edgesize
    fat_e=grow(fat_e);
end

% set label to 0 on fat edge set  ����ǩ����Ϊ0��
labels(fat_e>0)=0;

% find connected components using bwlabeln (image processing toolbox)  ʹ��bwlabeln(ͼ��������)�ҵ����ӵ����
labels=bwlabeln(labels,6); 

end

function x_new = grow(x)
% Function to grow level set of 1 by one pixel (in binary image)   ������1 / 1����(������ͼ��)����ʽ����

x_new=x; 

for i=-1:1
    
    for j=-1:1
        
        for k=-1:1
            
            x_new = x_new + circshift(x,[i j k]) + circshift(x,[-i -j -k]);
        end
    end
end

end