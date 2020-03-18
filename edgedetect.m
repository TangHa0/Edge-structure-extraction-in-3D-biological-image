
function [edge,conn_comp,edge_strength] = edgedetect(f,filter,thresholds,min_edge_functional)

%{
 function [edge,conn_comp,edge_strength] = edgedetect(f,filter,thresholds,min_edge_functional)
 ����min_edge_functional��maxedge��һ��Сֵ�ض�(�ڵ���isosurface֮ǰ��������min_edge_functional��ģ��С��min_edge_function��ֵ������Ϊ0)��
���������Ĵ�ͬ�����������ࡰ�ٱߡ������ݶȴ�С�ľֲ�����ֵ��������Щ��������maxedge���㼯(��˱���ֵ���⵽)��
����ֻ����ֵ������ɾ�����⽫������������ʱ�������������һ��ѡ�����Ƴ�С���ȵı仯�����һ�ж��㹻�죬�ҾͰ�����Ϊ0��
     
%}


% ���: 
%   edge:  ��Ե:��Ե/��Ծ����(��Ͷ����б�)
%   conn_comp: Connected components of edge set ��Ե������ͨ���������
%   edge_strength: Gradient magnitude |nabla f| (of smoothed data) �ݶȴ�С


%{
 
����:
   f: ����ͼ��ǰ����ģ�����ά���񡣺���ʹ����ʵ�ĺ˴Ź���ͼ��
   filter: ��˹�˲���      ˫���� filter(1)������׼ƫ����˲�(2)��˹����˲����Ĵ��ڴ�С��
   thresholds:  ˫��ֵ
   min_edge_functional: ��ֵ�ж� ���ٶ���׼ȷ��֮���ѡһ


%}



%���������ȷ�Խ�����֤
% Ĭ������ĸ���
if nargin <= 3    % �����������
  min_edge_functional=0;
end

if nargin <= 2
  thresholds=0;  %��ֵ
end

if nargin <= 1
  filter=0;   %������
end

if nargin == 0
  error('Missing input arguments.');
end


% ��֤�������  ��������
% ndims:����ά����С
% isvector:ȷ�������Ƿ�Ϊ����
% isscalar:ȷ�������Ƿ�Ϊ����
if ndims(f)~=3 || ~isvector(thresholds) || ~isvector(filter) ...
               || ~isscalar(min_edge_functional)
    error('Wrong input argument structure.');
end

% standard window size: 2*sigma   ��˹�˲�ƽ������Ϊ������׼��
if isscalar(filter)    % isscalar��ȷ�������Ƿ�Ϊ����
    window=ceil(2*filter);  %ceil�������������������
else
    sigma=filter(1); window=filter(2); 
end

% ����ͼ�ε�ƽ������˹�˲�Ƕ�뵽��ά�˲�����
if(sigma > 0)    
    filt=gaussfilter([window;window;window],(diag(sigma).^2));
    f_smooth=imfilter(f,filt,NaN);  %��ά�˲�  NaNֵ�þ���
else
    f_smooth=f;
end

% calculate gradient magnitude  �����ݶȷ�ֵ
[gx,gy,gz]=gradient(f_smooth);   % һά�����ϵĲ�֣� x y z ����������������
f_gradmag=sqrt(gx.^2+gy.^2+gz.^2);  %�ݶȷ�ֵ�ļ���
    
%   ����һ�ײ��ģ��
filt_dx=[-1/2, 0, 1/2];
filt_dy=shiftdim(filt_dx,1);    %shiftdimά��ת�� ����һ��ά��
filt_dz=shiftdim(filt_dx,-1);

%   ������ײ��ģ��
filt_dxx=[1, -2, 1];

filt_dyy=shiftdim(filt_dxx,1);

filt_dzz=shiftdim(filt_dxx,-1);

filt_dxy=[1/4, 0, -1/4; 0 0 0; -1/4, 0, 1/4];

filt_dyz=shiftdim(filt_dxy,-1);

filt_dxz=shiftdim(filt_dyz,1);


% define differential filter masks - third order  
filt_dxxx=[-1/2, 1, 0, -1, 1/2];

filt_dyyy=shiftdim(filt_dxxx,1);

filt_dzzz=shiftdim(filt_dxxx,-1);

filt_dxxy=[-1/2, 1, -1/2; 0 0 0; 1/2, -1, 1/2];

filt_dxxz=permute(filt_dxxy,[3 2 1]);

filt_dyyx=permute(filt_dxxy,[2 1 3]);

filt_dyyz=permute(filt_dxxy,[2 3 1]);

filt_dzzx=permute(filt_dxxy,[3 1 2]);

filt_dzzy=permute(filt_dxxy,[1 3 2]);

filt_dxyz=(1/2)*cat(3,-filt_dxy,zeros(3,3),filt_dxy);


% ������˹�˲���ͼ����������

% ����һ�ײ��
f_dx=imfilter(f_smooth,filt_dx);

f_dy=imfilter(f_smooth,filt_dy);

f_dz=imfilter(f_smooth,filt_dz);



%������ײ��
f_dxx=imfilter(f_smooth,filt_dxx);

f_dyy=imfilter(f_smooth,filt_dyy);

f_dzz=imfilter(f_smooth,filt_dzz);

f_dxy=imfilter(f_smooth,filt_dxy);

f_dyz=imfilter(f_smooth,filt_dyz);

f_dxz=imfilter(f_smooth,filt_dxz);


%���ײ��
f_dxxx=imfilter(f_smooth,filt_dxxx);

f_dyyy=imfilter(f_smooth,filt_dyyy);

f_dzzz=imfilter(f_smooth,filt_dzzz);

f_dxxy=imfilter(f_smooth,filt_dxxy);

f_dxxz=imfilter(f_smooth,filt_dxxz);

f_dyyx=imfilter(f_smooth,filt_dyyx);

f_dyyz=imfilter(f_smooth,filt_dyyz);

f_dzzx=imfilter(f_smooth,filt_dzzx);

f_dzzy=imfilter(f_smooth,filt_dzzy);

f_dxyz=imfilter(f_smooth,filt_dxyz);


% �������ڷǼ���ֵ���ƵĲ��
maxedge = f_dx.^2.*f_dxx + f_dy.^2.*f_dyy  ...
    + f_dz.^2.*f_dzz + 2*f_dx.*f_dy.*f_dxy ...
    + 2*f_dy.*f_dz.*f_dyz + 2*f_dx.*f_dz.*f_dxz;


% ģ����С�ڵ�ֵ�жϵĶ���Ϊ0���������Ļ������õ�ֵ�ж�
maxedge(abs(maxedge)<min_edge_functional)=0;
%abs Ϊ����ֵ

% find potential edges by taking zero pass of above function  ͨ��ȡ��������������ҵ�Ǳ�ڵı�Ե
edge=isosurface(maxedge,0);  
%fv = isosurface(X,Y,Z,V,isovalue) ���� isovalue ��ָ���ĵ�ֵ��ֵ���������� V �����ֵ�����ݡ�����ֵ�����Ӿ���ָ��ֵ�ĵ㣬��ȸ�������������ͬ�ĵ�ķ�ʽ������ͬ��
%���� X��Y �� Z ��ʾ Cartesian ���������V ������Щ����㴦�Ķ�Ӧֵ���������飨X��Y �� Z�������ǵ����Ĳ��ҷ��� meshgrid ���ɵĸ�ʽ��V �������� X��Y �� Z ��С��ͬ����ά�����顣
%struct fv ������ֵ�����Ͷ��㣬�����Խ�����ֱ�Ӵ��ݵ� patch ����



% return if threshold is zero or edge set is empty  �����ֵΪ���߼�Ϊ�գ��򷵻�
if ~thresholds(1) || isempty(edge.vertices)
    return; 
end

% third order sign condition to filter out extrema which are not maxima  ���׷����������˳������ֵ�ļ�ֵ��
signcond = f_dx.^3.*f_dxxx + f_dy.^3.*f_dyyy + f_dz.^3.*f_dzzz ...
         + 3*f_dx.^2.*f_dy.*f_dxxy + 3*f_dx.^2.*f_dz.*f_dxxz ...
         + 3*f_dy.^2.*f_dx.*f_dyyx + 3*f_dy.^2.*f_dz.*f_dyyz ...
         + 3*f_dz.^2.*f_dx.*f_dzzx + 3*f_dz.^2.*f_dy.*f_dzzy ...
         + 6*f_dx.*f_dy.*f_dz.*f_dxyz;

     
% Ӧ����������ɾ���Ǽ�ֵ�ļ�ֵ 
edge_sign=interp3(signcond,edge.vertices(:,1),edge.vertices(:,2),edge.vertices(:,3),'linear');
ind=(edge_sign<0); edge=reducesurface(edge,ind);



% ����ֵ�ı�Ե�����ϲ����ݶȴ�С�� �ݶȷ�ֵ��ֵ
edge_strength=interp3(f_gradmag,edge.vertices(:,1),edge.vertices(:,2),edge.vertices(:,3),'linear');

%  ���ݶȷ��ȵ�����ֵ���޵ı�Ե����
ind=(edge_strength>thresholds(2)); edge=reducesurface(edge,ind); edge_strength=edge_strength(ind);

% remove connected components of edge set where gradient magnitude  �Ƴ���Ե�������������ͨ���������������ͨ�����ݶȴ�С��Զ���ᳬ���������ֵ��
% is never above the upper threshold
[ind,conn_comp]=hysteresis(edge,edge_strength,thresholds);
edge=reducesurface(edge,ind); edge_strength=edge_strength(ind); 
[~,~,conn_comp]=unique(conn_comp(ind));

     end
  
     
     

function [ind,conn_comp] = hysteresis(edge,edge_strength,thresholds)
% Apply upper hysteresis threshold on face & vertex list ����Ͷ����б���Ӧ�� ���ͺ���ֵ��

% initialize accepted vertex list  ��ʼ�������б�
ind=true(size(edge.vertices,1),1); 
%szdim = size(A,dim) ����ά�� dim �ĳ���

% find connected components of edge set  ���ұ�Ե������ͨ����
conn_comp=conncomp(edge);

%conncomp(G) �� bin ��ʽ����ͼ G ����ͨ������bin ���ָʾͼ�е�ÿ���ڵ������ķ�����
% G ������ͼ������·�������������ڵ�����ͬһ������


% loop over connected components  ������ͨ����
for i=1:max(conn_comp)
    
    % remove component if maximal edge strength is below upper threshold
    % �������Ե��ֵ����������ֵ����ɾ������ͨ����
    if max(edge_strength(conn_comp==i)) < thresholds(1)
        ind(conn_comp==i)=false;
    end
    
end

end


%ɾ��һЩ��
function surf_new = reducesurface(surf,ind)
% Reduce list of vertices and faces (triangles) to given indices  �����������б������Σ����ٵ�����������

% reduce vertex list ���ٶ����б�
surf_new.vertices=surf.vertices(ind,:);

% corresponding faces ��Ӧ����
num=(1:sum(ind))'; indlist=NaN(size(surf.vertices,1),1);
indlist(ind)=num; new_elem=indlist(surf.faces);
surf_new.faces=new_elem(~isnan(sum(new_elem,2)),:);

end

function cc = conncomp(surf,label)
% Find connected components of given surface (list of faces and vertices) ���Ҹ�������������������Ͷ����б�

% initialize label list if not given a prior ��ʼ����ǩ�б�
if ~exist('label','var')
   label=(1:size(surf.vertices,1))';
end

% loop until all nodes of every face have the same label ѭ����ֱ��ÿ��������нڵ㶼����ͬ�ı�ǩ��
while max(std(label(surf.faces),0,2))~=0
    
    % smallest label of every face (vectorized) ÿ�������С��ǩ��ʸ������
    m=repmat(min(label(surf.faces),[],2),[3 1]);

    % set label of node to smallest of face (randomize to break cycles) ���ڵ�ı�ǩ����Ϊ�����Сֵ(�����������)
    j=surf.faces(:); p=randperm(length(m));
    label(j(p))=m(p);    
end

% relabel  �����ñ�ǩ����
[~,~,cc]=unique(label);

end

function filter = gaussfilter(sz,sigma)
% Calculate Gaussian filter tensor with given dimensions sz
% ��������ߴ�Ϊsz�ĸ�˹�˲�������Э�������
% and covariance matrix sigma

% initialize  ��ʼ��
dims=length(sz); ind_cell=cell(dims,1); filter=zeros(sz');

% calculate midpoint   �����е�
midpoint=(sz+1)/2;

% define (unnormalized) gaussian function ����(�ǹ淶)��˹����
gauss=@(x) exp(-x'*(sigma\x)/2);

% loop over all elements of tensor  ������������Ԫ�ؽ���ѭ��
for i=1:numel(filter)  % numel ������������м���Ԫ��
    
    % get indices of element with linear index i �õ���������i��Ԫ��������
    [ind_cell{:}]=ind2sub(sz',i); ind=cell2mat(ind_cell);   %ind2sub  �����������±�
    
    % calculate gaussian value at size i  ����ߴ�Ϊi�ĸ�˹ֵ��
    filter(i)=gauss(ind-midpoint);
    
end

% get l1 norm of filter �õ��˲�����l1������
% ������Ϊŷ����þ���
n=filter;

for d=1:dims
    n=sum(n,dims-d+1);
end

% normalize  ��һ��
filter=filter/n;

end