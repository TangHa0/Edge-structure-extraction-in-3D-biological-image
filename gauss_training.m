load mri;
%其中第三个维度为颜色维度，应该剔除

%使用squeeze()函数删除单一的维度
D2 = squeeze(D);

%数据类型强制转换为double型
D1 = double(D2);

filter = [2; 3];           
sigma=filter(1); window=filter(2);

filt=gaussfilter([window;window;window],(diag(sigma).^2));
    f_smooth=imfilter(D1,filt,NaN);  %多维滤波  NaN值得矩阵
    
sigma=filter(1); window=filter(2);

    



voxelSurf(f_smooth,true);





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

