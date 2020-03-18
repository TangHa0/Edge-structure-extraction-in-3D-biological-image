%模拟用测试函数（模拟3D图像，非真实立体图）
%三维100x180x100的三维网格
[X,Y,Z] = meshgrid(1:180,1:100,1:100);

background = sin((X+Y+Z)*pi/100);

blobs = double(sqrt((X-60).^2+(Y-60).^2+(Z-50).^2)<=20) ...
      + double(sqrt((X-140).^2+(Y-60).^2+(Z-50).^2)<=10) ...
      - double(abs(X-100)<=60 & abs(Y-20)<=5 & abs(Z-50)<=20);
  
f = background + blobs + 0.2*randn(size(X));


%读取matlab图片工具箱中的自带mri图片 128x128x1x7
load mri;
%其中第三个维度为颜色维度，应该剔除

%使用squeeze()函数删除单一的维度
D2 = squeeze(D);

%数据类型强制转换为double型
D1 = double(D2);

%Z轴第n个
image_num = 18;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%{
prefix = 'I0000';
fnum = 417:433;
ext = '_anon.dcm';

%系列中的第一个文件的文件名 (nobkpt)
fname = [prefix num2str(fnum(1)) ext];

%examine file header (nobkpt)
info = dicominfo(fname)

%从元数据中提取大小信息。 (nobkpt)
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
%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上为输入数据集




%高斯滤波器的窗口大小即高斯核的长宽 5x5

%定义边缘检测参数
%高斯滤波器的 标准偏差  以及 5X5的模板窗口大小
%高斯滤波器的窗口大小（window）即高斯核的长宽 5x5x5

filter = [2; 7];           

%阈值 一高一低，用来执行双阈值筛选
thresholds = [0.2; 0.1];   

%低值切断，值越大速度越快 ，但准确率越低
min_edge_functional = 1e-5; 

% 检测图像的边缘（灰度发生阶级跳跃的地方），核心函数 3D Canny算子
[edge,conncomp,edgestrength] = edgedetect(D1,filter,thresholds,min_edge_functional);

% 分割参数 立体表面的边缘厚度
edgesize = 1;   

% 分割图像所使用的边缘集（已经过连通性筛选）
seglabelsresult = segment(edge,size(D1),edgesize);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  以下为分割部分


% Edge detection parameters 边缘检测的参数
filter = [1.4; 3];            % Options for Gaussian filter, standard deviation and window size 选择高斯滤波器，标准偏差和窗口大小
thresholds = [1; 0];    % Upper and lower hysteresis thresholds 高低两个阈值
min_edge_functional = 1e-4; % Low value cutoff for edge detection functional:  低值切断
                            % larger value -> faster, less accurate 越大的值->越快，但越不准确

% Find jump surfaces/edges of f
[edge,conncomp,edgestrength] = edgedetect(seglabelsresult,filter,thresholds,min_edge_functional);

% Segmentation parameters 分割参数
edgesize = 1;   % Thickness of the voxelized edge surface 立体表面的边缘厚度



%figure; title('边缘强度');  %边缘梯度
%plotsurface(edge,edgestrength,0.4); colorbar; view(130,30);

%figure; title('分割');   %边缘集的连接部件
%plotsurface(edge,conncomp,0.4); camlight headlight; view(130,30); 








%%%%%%%%%%%%%%%%%%%%%%%%另外的切片图




 

%voxelSurf(J1,true);



filter = [2; 5];           

%阈值 一高一低，用来执行双阈值筛选
thresholds = [0.2; 0.1];   

%低值切断，值越大速度越快 ，但准确率越低
min_edge_functional = 1e-3; 

% 检测图像的边缘（灰度发生阶级跳跃的地方），核心函数 3D Canny算子
%[edge,conncomp,edgestrength] = edgedetect(J1,filter,thresholds,min_edge_functional);

% 分割参数 立体表面的边缘厚度
edgesize = 1;   

% 分割图像所使用的边缘集（已经过连通性筛选）
%eglabelsresult = segment(edge,size(J1),edgesize);

% 灰度图像显示
%figure; imagesc(J1(:,:,50)); axis image; set(gca,'YDir','normal');
%colormap gray; colorbar; xlabel('x'); ylabel('y'); title('3D-data [z=50]');

%边缘3D透明可视化
%figure; title('检测出的边缘');   %边缘集的连接部件（经过滤除了）
%plotsurface(edge,conncomp,0.4); camlight headlight; view(130,30);  %view仅仅定义图片的第一个显示视角，视角是可以随意调整的

