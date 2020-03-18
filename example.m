%ģ���ò��Ժ�����ģ��3Dͼ�񣬷���ʵ����ͼ��
%��ά100x180x100����ά����
[X,Y,Z] = meshgrid(1:180,1:100,1:100);

background = sin((X+Y+Z)*pi/100);

blobs = double(sqrt((X-60).^2+(Y-60).^2+(Z-50).^2)<=20) ...
      + double(sqrt((X-140).^2+(Y-60).^2+(Z-50).^2)<=10) ...
      - double(abs(X-100)<=60 & abs(Y-20)<=5 & abs(Z-50)<=20);
  
f = background + blobs + 0.2*randn(size(X));


%��ȡmatlabͼƬ�������е��Դ�mriͼƬ 128x128x1x7
load mri;
%���е�����ά��Ϊ��ɫά�ȣ�Ӧ���޳�

%ʹ��squeeze()����ɾ����һ��ά��
D2 = squeeze(D);

%��������ǿ��ת��Ϊdouble��
D1 = double(D2);

%Z���n��
image_num = 18;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%{
prefix = 'I0000';
fnum = 417:433;
ext = '_anon.dcm';

%ϵ���еĵ�һ���ļ����ļ��� (nobkpt)
fname = [prefix num2str(fnum(1)) ext];

%examine file header (nobkpt)
info = dicominfo(fname)

%��Ԫ��������ȡ��С��Ϣ�� (nobkpt)
voxel_size = [info.PixelSpacing; info.SliceThickness]'

%����Ƭ�ѵ���������ȶѵ��������γɾ�������ά�ȵľ���
hWaitBar = waitbar(0,'Reading DICOM files');
for i=length(fnum):-1:1
  fname = [prefix num2str(fnum(i)) ext];
  J(:,:,i) = uint8(dicomread(fname));
  waitbar((length(fnum)-i)/length(fnum))
end
delete(hWaitBar)
whos J

J1 = double (J);
%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ϊ�������ݼ�




%��˹�˲����Ĵ��ڴ�С����˹�˵ĳ��� 5x5

%�����Ե������
%��˹�˲����� ��׼ƫ��  �Լ� 5X5��ģ�崰�ڴ�С
%��˹�˲����Ĵ��ڴ�С��window������˹�˵ĳ��� 5x5x5

filter = [2; 7];           

%��ֵ һ��һ�ͣ�����ִ��˫��ֵɸѡ
thresholds = [0.2; 0.1];   

%��ֵ�жϣ�ֵԽ���ٶ�Խ�� ����׼ȷ��Խ��
min_edge_functional = 1e-5; 

% ���ͼ��ı�Ե���Ҷȷ����׼���Ծ�ĵط��������ĺ��� 3D Canny����
[edge,conncomp,edgestrength] = edgedetect(D1,filter,thresholds,min_edge_functional);

% �ָ���� �������ı�Ե���
edgesize = 1;   

% �ָ�ͼ����ʹ�õı�Ե�����Ѿ�����ͨ��ɸѡ��
seglabelsresult = segment(edge,size(D1),edgesize);

% �Ҷ�ͼ����ʾ
%figure; imagesc(D(:,:,image_num)); axis image; set(gca,'YDir','normal');
%colormap gray; colorbar; xlabel('x'); ylabel('y'); title('3D-data [z=50]');

%��Եǿ��ͼ����ʾ
%figure; title('��Եǿ��');  %��Ե�ݶ�
%plotsurface(edge,edgestrength,0.4); colorbar; view(130,30);

%��Ե3D͸�����ӻ�
figure; title('�����ı�Ե');   %��Ե�������Ӳ����������˳��ˣ�
plotsurface(edge,conncomp,0.4); camlight headlight; view(130,30);  %view��������ͼƬ�ĵ�һ����ʾ�ӽǣ��ӽ��ǿ������������


%2D�ָ�ͼ
%figure; imagesc(seglabelsresult(:,:,8)); axis image; set(gca,'YDir','normal');


%���ص���ά���ӻ���ͼ
%voxelSurf(J1,true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ����Ϊ�ָ��


% Edge detection parameters ��Ե���Ĳ���
filter = [1.4; 3];            % Options for Gaussian filter, standard deviation and window size ѡ���˹�˲�������׼ƫ��ʹ��ڴ�С
thresholds = [1; 0];    % Upper and lower hysteresis thresholds �ߵ�������ֵ
min_edge_functional = 1e-4; % Low value cutoff for edge detection functional:  ��ֵ�ж�
                            % larger value -> faster, less accurate Խ���ֵ->Խ�죬��Խ��׼ȷ

% Find jump surfaces/edges of f
[edge,conncomp,edgestrength] = edgedetect(seglabelsresult,filter,thresholds,min_edge_functional);

% Segmentation parameters �ָ����
edgesize = 1;   % Thickness of the voxelized edge surface �������ı�Ե���



%figure; title('��Եǿ��');  %��Ե�ݶ�
%plotsurface(edge,edgestrength,0.4); colorbar; view(130,30);

%figure; title('�ָ�');   %��Ե�������Ӳ���
%plotsurface(edge,conncomp,0.4); camlight headlight; view(130,30); 








%%%%%%%%%%%%%%%%%%%%%%%%�������Ƭͼ




 

%voxelSurf(J1,true);



filter = [2; 5];           

%��ֵ һ��һ�ͣ�����ִ��˫��ֵɸѡ
thresholds = [0.2; 0.1];   

%��ֵ�жϣ�ֵԽ���ٶ�Խ�� ����׼ȷ��Խ��
min_edge_functional = 1e-3; 

% ���ͼ��ı�Ե���Ҷȷ����׼���Ծ�ĵط��������ĺ��� 3D Canny����
%[edge,conncomp,edgestrength] = edgedetect(J1,filter,thresholds,min_edge_functional);

% �ָ���� �������ı�Ե���
edgesize = 1;   

% �ָ�ͼ����ʹ�õı�Ե�����Ѿ�����ͨ��ɸѡ��
%eglabelsresult = segment(edge,size(J1),edgesize);

% �Ҷ�ͼ����ʾ
%figure; imagesc(J1(:,:,50)); axis image; set(gca,'YDir','normal');
%colormap gray; colorbar; xlabel('x'); ylabel('y'); title('3D-data [z=50]');

%��Ե3D͸�����ӻ�
%figure; title('�����ı�Ե');   %��Ե�������Ӳ����������˳��ˣ�
%plotsurface(edge,conncomp,0.4); camlight headlight; view(130,30);  %view��������ͼƬ�ĵ�һ����ʾ�ӽǣ��ӽ��ǿ������������

