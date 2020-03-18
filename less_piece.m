prefix = 'I0000';
fnum = 417:422;
ext = '_anon.dcm';

%first filename in series (nobkpt)
fname = [prefix num2str(fnum(1)) ext];

%examine file header (nobkpt)
info = dicominfo(fname)

%extract size info from metadata (nobkpt)
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

 


%��˹�˲����Ĵ��ڴ�С����˹�˵ĳ��� 5x5

%�����Ե������
%��˹�˲����� ��׼ƫ��  �Լ� 5X5��ģ�崰�ڴ�С
%��˹�˲����Ĵ��ڴ�С��window������˹�˵ĳ��� 5x5x5

filter = [2; 7];           

%��ֵ һ��һ�ͣ�����ִ��˫��ֵɸѡ
thresholds = [1.1; 1];   

%��ֵ�жϣ�ֵԽ���ٶ�Խ�� ����׼ȷ��Խ��
min_edge_functional = 1e-4; 

% ���ͼ��ı�Ե���Ҷȷ����׼���Ծ�ĵط��������ĺ��� 3D Canny����
[edge,conncomp,edgestrength] = edgedetect(J1,filter,thresholds,min_edge_functional);

% �ָ���� �������ı�Ե���
edgesize = 1;   

% �ָ�ͼ����ʹ�õı�Ե�����Ѿ�����ͨ��ɸѡ��
%seglabelsresult = segment(edge,size(J1),edgesize);

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


