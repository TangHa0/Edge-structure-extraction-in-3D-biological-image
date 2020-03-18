prefix = 'I0000';
fnum = 417:422;
ext = '_anon.dcm';

%first filename in series (nobkpt)
fname = [prefix num2str(fnum(1)) ext];

%examine file header (nobkpt)
info = dicominfo(fname)

%extract size info from metadata (nobkpt)
voxel_size = [info.PixelSpacing; info.SliceThickness]'



%把切片堆叠按纵向深度堆叠起来，形成具有三个维度的矩阵
hWaitBar = waitbar(0,'Reading DICOM files');

for i=length(fnum):-1:1
    
  fname = [prefix num2str(fnum(i)) ext];
  J(:,:,i) = uint8(dicomread(fname));
  waitbar((length(fnum)-i)/length(fnum))
  
end
delete(hWaitBar)
whos J

J1 = double (J);

 


%高斯滤波器的窗口大小即高斯核的长宽 5x5

%定义边缘检测参数
%高斯滤波器的 标准偏差  以及 5X5的模板窗口大小
%高斯滤波器的窗口大小（window）即高斯核的长宽 5x5x5

filter = [2; 7];           

%阈值 一高一低，用来执行双阈值筛选
thresholds = [1.1; 1];   

%低值切断，值越大速度越快 ，但准确率越低
min_edge_functional = 1e-4; 

% 检测图像的边缘（灰度发生阶级跳跃的地方），核心函数 3D Canny算子
[edge,conncomp,edgestrength] = edgedetect(J1,filter,thresholds,min_edge_functional);

% 分割参数 立体表面的边缘厚度
edgesize = 1;   

% 分割图像所使用的边缘集（已经过连通性筛选）
%seglabelsresult = segment(edge,size(J1),edgesize);

% 灰度图像显示
%figure; imagesc(D(:,:,image_num)); axis image; set(gca,'YDir','normal');
%colormap gray; colorbar; xlabel('x'); ylabel('y'); title('3D-data [z=50]');

%边缘强度图像显示
%figure; title('边缘强度');  %边缘梯度
%plotsurface(edge,edgestrength,0.4); colorbar; view(130,30);

%边缘3D透明可视化
figure; title('检测出的边缘');   %边缘集的连接部件（经过滤除了）
plotsurface(edge,conncomp,0.4); camlight headlight; view(130,30);  %view仅仅定义图片的第一个显示视角，视角是可以随意调整的


%2D分割图
%figure; imagesc(seglabelsresult(:,:,8)); axis image; set(gca,'YDir','normal');


%体素的三维可视化绘图
%voxelSurf(J1,true);


