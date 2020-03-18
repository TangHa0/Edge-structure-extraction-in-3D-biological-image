load mri;
%���е�����ά��Ϊ��ɫά�ȣ�Ӧ���޳�

%ʹ��squeeze()����ɾ����һ��ά��
D2 = squeeze(D);

%��������ǿ��ת��Ϊdouble��
D1 = double(D2);

filter = [2; 3];           
sigma=filter(1); window=filter(2);

filt=gaussfilter([window;window;window],(diag(sigma).^2));
    f_smooth=imfilter(D1,filt,NaN);  %��ά�˲�  NaNֵ�þ���
    
sigma=filter(1); window=filter(2);

    



voxelSurf(f_smooth,true);





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

