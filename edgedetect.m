
function [edge,conn_comp,edge_strength] = edgedetect(f,filter,thresholds,min_edge_functional)

%{
 function [edge,conn_comp,edge_strength] = edgedetect(f,filter,thresholds,min_edge_functional)
 参数min_edge_functional是maxedge的一个小值截断(在调用isosurface之前，所有与min_edge_functional的模量小于min_edge_function的值都设置为0)。
带有噪声的大同质区域包含许多“假边”，即梯度大小的局部极大值。所有这些都包含在maxedge的零集(因此被等值面检测到)，
并且只在阈值步骤中删除，这将显著增加运行时。所以我添加了一个选项来移除小幅度的变化。如果一切都足够快，我就把它设为0。
     
%}


% 输出: 
%   edge:  边缘:边缘/跳跃表面(面和顶点列表)
%   conn_comp: Connected components of edge set 边缘集的连通组件（检测后）
%   edge_strength: Gradient magnitude |nabla f| (of smoothed data) 梯度大小


%{
 
输入:
   f: 输入图像，前期用模拟的三维网格。后期使用真实的核磁共振图像。
   filter: 高斯滤波器      双参数 filter(1)包含标准偏差和滤波(2)高斯卷积滤波器的窗口大小。
   thresholds:  双阈值
   min_edge_functional: 低值切断 在速度与准确率之间二选一


%}



%对输入的正确性进行验证
% 默认输入的个数
if nargin <= 3    % 输入变量个数
  min_edge_functional=0;
end

if nargin <= 2
  thresholds=0;  %阈值
end

if nargin <= 1
  filter=0;   %过滤器
end

if nargin == 0
  error('Missing input arguments.');
end


% 验证输入参数  （错误处理）
% ndims:返回维数大小
% isvector:确定输入是否为向量
% isscalar:确定输入是否为标量
if ndims(f)~=3 || ~isvector(thresholds) || ~isvector(filter) ...
               || ~isscalar(min_edge_functional)
    error('Wrong input argument structure.');
end

% standard window size: 2*sigma   高斯滤波平滑窗口为两倍标准差
if isscalar(filter)    % isscalar：确定输入是否为标量
    window=ceil(2*filter);  %ceil：朝正无穷大四舍五入
else
    sigma=filter(1); window=filter(2); 
end

% 进行图形的平滑，高斯滤波嵌入到多维滤波器中
if(sigma > 0)    
    filt=gaussfilter([window;window;window],(diag(sigma).^2));
    f_smooth=imfilter(f,filt,NaN);  %多维滤波  NaN值得矩阵
else
    f_smooth=f;
end

% calculate gradient magnitude  计算梯度幅值
[gx,gy,gz]=gradient(f_smooth);   % 一维方向上的差分， x y z 三个方向三个分量
f_gradmag=sqrt(gx.^2+gy.^2+gz.^2);  %梯度幅值的计算
    
%   定义一阶差分模板
filt_dx=[-1/2, 0, 1/2];
filt_dy=shiftdim(filt_dx,1);    %shiftdim维数转换 左移一个维度
filt_dz=shiftdim(filt_dx,-1);

%   定义二阶差分模板
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


% 经过高斯滤波的图像用来求差分

% 计算一阶差分
f_dx=imfilter(f_smooth,filt_dx);

f_dy=imfilter(f_smooth,filt_dy);

f_dz=imfilter(f_smooth,filt_dz);



%计算二阶差分
f_dxx=imfilter(f_smooth,filt_dxx);

f_dyy=imfilter(f_smooth,filt_dyy);

f_dzz=imfilter(f_smooth,filt_dzz);

f_dxy=imfilter(f_smooth,filt_dxy);

f_dyz=imfilter(f_smooth,filt_dyz);

f_dxz=imfilter(f_smooth,filt_dxz);


%三阶差分
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


% 计算用于非极大值抑制的差分
maxedge = f_dx.^2.*f_dxx + f_dy.^2.*f_dyy  ...
    + f_dz.^2.*f_dzz + 2*f_dx.*f_dy.*f_dxy ...
    + 2*f_dy.*f_dz.*f_dyz + 2*f_dx.*f_dz.*f_dxz;


% 模量级小于低值切断的都设为0，如果够快的话都不用低值切断
maxedge(abs(maxedge)<min_edge_functional)=0;
%abs 为绝对值

% find potential edges by taking zero pass of above function  通过取上述函数的零点找到潜在的边缘
edge=isosurface(maxedge,0);  
%fv = isosurface(X,Y,Z,V,isovalue) 基于 isovalue 中指定的等值面值处的体数据 V 计算等值面数据。即等值面连接具有指定值的点，与等高线连接仰角相同的点的方式大致相同。
%数组 X、Y 和 Z 表示 Cartesian 轴对齐网格。V 包含这些网格点处的对应值。坐标数组（X、Y 和 Z）必须是单调的并且符合 meshgrid 生成的格式。V 必须是与 X、Y 和 Z 大小相同的三维体数组。
%struct fv 包含等值面的面和顶点，您可以将它们直接传递到 patch 命令



% return if threshold is zero or edge set is empty  如果阈值为零或边集为空，则返回
if ~thresholds(1) || isempty(edge.vertices)
    return; 
end

% third order sign condition to filter out extrema which are not maxima  三阶符号条件过滤出非最大值的极值。
signcond = f_dx.^3.*f_dxxx + f_dy.^3.*f_dyyy + f_dz.^3.*f_dzzz ...
         + 3*f_dx.^2.*f_dy.*f_dxxy + 3*f_dx.^2.*f_dz.*f_dxxz ...
         + 3*f_dy.^2.*f_dx.*f_dyyx + 3*f_dy.^2.*f_dz.*f_dyyz ...
         + 3*f_dz.^2.*f_dx.*f_dzzx + 3*f_dz.^2.*f_dy.*f_dzzy ...
         + 6*f_dx.*f_dy.*f_dz.*f_dxyz;

     
% 应用三阶条件删掉非极值的极值 
edge_sign=interp3(signcond,edge.vertices(:,1),edge.vertices(:,2),edge.vertices(:,3),'linear');
ind=(edge_sign<0); edge=reducesurface(edge,ind);



% 在阈值的边缘顶点上插入梯度大小。 梯度幅值插值
edge_strength=interp3(f_gradmag,edge.vertices(:,1),edge.vertices(:,2),edge.vertices(:,3),'linear');

%  除梯度幅度低于阈值下限的边缘顶点
ind=(edge_strength>thresholds(2)); edge=reducesurface(edge,ind); edge_strength=edge_strength(ind);

% remove connected components of edge set where gradient magnitude  移除边缘的连接组件（连通分量），在这个连通分量梯度大小永远不会超过上面的阈值。
% is never above the upper threshold
[ind,conn_comp]=hysteresis(edge,edge_strength,thresholds);
edge=reducesurface(edge,ind); edge_strength=edge_strength(ind); 
[~,~,conn_comp]=unique(conn_comp(ind));

     end
  
     
     

function [ind,conn_comp] = hysteresis(edge,edge_strength,thresholds)
% Apply upper hysteresis threshold on face & vertex list 在面和顶点列表上应用 上滞后阈值。

% initialize accepted vertex list  初始化顶点列表
ind=true(size(edge.vertices,1),1); 
%szdim = size(A,dim) 返回维度 dim 的长度

% find connected components of edge set  查找边缘集的连通分量
conn_comp=conncomp(edge);

%conncomp(G) 以 bin 形式返回图 G 的连通分量。bin 编号指示图中的每个节点所属的分量。
% G 是无向图，则有路径相连的两个节点属于同一分量。


% loop over connected components  遍历连通分量
for i=1:max(conn_comp)
    
    % remove component if maximal edge strength is below upper threshold
    % 如果最大边缘幅值低于上限阈值，则删除此连通分量
    if max(edge_strength(conn_comp==i)) < thresholds(1)
        ind(conn_comp==i)=false;
    end
    
end

end


%删除一些边
function surf_new = reducesurface(surf,ind)
% Reduce list of vertices and faces (triangles) to given indices  将顶点和面的列表（三角形）减少到给定的索引

% reduce vertex list 减少顶点列表
surf_new.vertices=surf.vertices(ind,:);

% corresponding faces 相应的面
num=(1:sum(ind))'; indlist=NaN(size(surf.vertices,1),1);
indlist(ind)=num; new_elem=indlist(surf.faces);
surf_new.faces=new_elem(~isnan(sum(new_elem,2)),:);

end

function cc = conncomp(surf,label)
% Find connected components of given surface (list of faces and vertices) 查找给定曲面的连接组件（面和顶点列表）

% initialize label list if not given a prior 初始化标签列表。
if ~exist('label','var')
   label=(1:size(surf.vertices,1))';
end

% loop until all nodes of every face have the same label 循环，直到每个面的所有节点都有相同的标签。
while max(std(label(surf.faces),0,2))~=0
    
    % smallest label of every face (vectorized) 每个面的最小标签（矢量化）
    m=repmat(min(label(surf.faces),[],2),[3 1]);

    % set label of node to smallest of face (randomize to break cycles) 将节点的标签设置为面的最小值(随机打乱周期)
    j=surf.faces(:); p=randperm(length(m));
    label(j(p))=m(p);    
end

% relabel  重新用标签标明
[~,~,cc]=unique(label);

end

function filter = gaussfilter(sz,sigma)
% Calculate Gaussian filter tensor with given dimensions sz
% 计算给定尺寸为sz的高斯滤波张量和协方差矩阵
% and covariance matrix sigma

% initialize  初始化
dims=length(sz); ind_cell=cell(dims,1); filter=zeros(sz');

% calculate midpoint   计算中点
midpoint=(sz+1)/2;

% define (unnormalized) gaussian function 定义(非规范)高斯函数
gauss=@(x) exp(-x'*(sigma\x)/2);

% loop over all elements of tensor  对张量的所有元素进行循环
for i=1:numel(filter)  % numel 计算出数组中有几个元素
    
    % get indices of element with linear index i 得到线性索引i的元素索引。
    [ind_cell{:}]=ind2sub(sz',i); ind=cell2mat(ind_cell);   %ind2sub  线性索引的下标
    
    % calculate gaussian value at size i  计算尺寸为i的高斯值。
    filter(i)=gauss(ind-midpoint);
    
end

% get l1 norm of filter 得到滤波器的l1范数。
% 范数即为欧几里得距离
n=filter;

for d=1:dims
    n=sum(n,dims-d+1);
end

% normalize  归一化
filter=filter/n;

end