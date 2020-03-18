function [hh,TT,X,Y,Z,CC,AA] = voxelSurf(data,internal,varargin)

%在3d数组数据的位置绘制体素立方体，其中的条目为>0。
%数据是一个正整数的索引数组，0被当作空气，索引1到n。
%的颜色取决于colormap(推荐的行)。
%第二个参数指定是否绘制内部的面。
%使用第三个参数来指定坐标轴。
%voxel matrix数据的最大推荐尺寸:512x512x512 voxels。
%警告:计算将需要大量的RAM用于大型体素矩阵。
%(~16GB用于512x512x512像素矩阵)。


%这是一种完全矢量化的方法来显示三维索引整数数据中的体素。 颜色是根据颜色表进行的，建议使用颜色表（线）。
%该方法比其他现有的体素绘图仪快得多，并且可以在具有相当好的图形卡的计算机上轻松处理高达512x512x512数据点的3D阵列。
%（它只需要大量的RAM阵列）。
%另一个很好的特点是，它通过将体系边界处的体素面切割成两个三角形，在两个体系的边界处显示两种颜色。

[sY,sX,sZ]=size(data);

colordata=single(real(data)); %color data
alphadata=single(imag(data)); %alpha data

usealpha=false;

if  isreal(data) == false
   usealpha=true;
end
    
if internal ==false
    colordata=single(real(data) > 0);
    alphadata=single(imag(data) > 0);
end

axes=[1 sX 1 sY 1 sZ];
if nargin > 2 
   axes=varargin{1}; 
end

a=1;%only used if there are no complex values in data
if nargin > 3 
   a=varargin{2}; 
end

if nargin > 4 %specify to use C as transparency and not as color index
    %C作为透明度 而不是颜色索引
    
end


colorpadded=zeros(sY+2,sX+2,sZ+2,'single');
colorpadded(2:end-1,2:end-1,2:end-1)=colordata;
if usealpha
    alphapadded=zeros(sY+2,sX+2,sZ+2,'single');
    alphapadded(2:end-1,2:end-1,2:end-1)=alphadata;
end

C=single(real(data(:))); %stores colors
if usealpha
    A=single(imag(data(:))); %stores alpha values
end


%this is the total ids of all triangle points:
%所有体素三角形的ID总和
Pointlabels=uint32(1:(sY+1)*(sX+1)*(sZ+1));
Pointlabels=reshape(Pointlabels,sY+1,sX+1,sZ+1);

%there are 6 sides of the cubes, hence the following section is done six
%times . there are 3 different kinds of triangles, A = triangle with the
%color of the cube on the "inside", B = triangle with color of the cube on
%the oposite side, C= special case for triangles bordering air so that the
%cubes bordering air are not just half drawn.
%{
立方体有6个边，因此下面的部分完成6次。 有三种不同的三角形，A =内部立方体颜色的三角形，B =反方向立方体颜色的三角形，C =与空气接界的三角形的特殊情况 
不只是一半绘制。
%}


h=[1   2   2;
   3   2   2;
   2   1   2;
   2   3   2;
   2   2   1;
   2   2   3];



h2=[ 1 2 1  1;
     1 2 2  1;
     1 2 1  2;
     1 2 2  2];
 
 
 
h3=[h2(:,1),h2(:,3),h2(:,4);
    h2(:,2),h2(:,3),h2(:,4);
    h2(:,3),h2(:,1),h2(:,4);
    h2(:,3),h2(:,2),h2(:,4);
    h2(:,3),h2(:,4),h2(:,1);
    h2(:,3),h2(:,4),h2(:,2)];



cellA=cell(2,2,2);
cellA{1,1,1}=Pointlabels(1:end-1,1:end-1,1:end-1);
cellA{2,1,1}=Pointlabels(2:end-0,1:end-1,1:end-1);
cellA{1,2,1}=Pointlabels(1:end-1,2:end-0,1:end-1);
cellA{2,2,1}=Pointlabels(2:end-0,2:end-0,1:end-1);
cellA{1,1,2}=Pointlabels(1:end-1,1:end-1,2:end-0);
cellA{2,1,2}=Pointlabels(2:end-0,1:end-1,2:end-0);
cellA{1,2,2}=Pointlabels(1:end-1,2:end-0,2:end-0);
cellA{2,2,2}=Pointlabels(2:end-0,2:end-0,2:end-0);


TT=[];%三角形矢量的索引
CC=[];%向量的颜色值
AA=[];%向量的透明度

for s=1:6
    D =colorpadded(h(s,1):end+h(s,1)-3,h(s,2):end+h(s,2)-3,h(s,3):end+h(s,3)-3);
    
    
    SA=(D < colordata ) ;%| ( D2 < alphadata);
    SB=(D > colordata & colordata > 0);% | (D2 > alphadata & alphadata >0);
    SC=(D < (colordata>0));% | ( D2 < (alphadata>0));

    ss=4*(s-1)+1;
    tmp1=cellA{h3(ss  ,1),h3(ss  ,2),h3(ss  ,3)}(:).*uint32(SA(:));
    tmp2=cellA{h3(ss+1,1),h3(ss+1,2),h3(ss+1,3)}(:).*uint32(SA(:));
    tmp3=cellA{h3(ss+2,1),h3(ss+2,2),h3(ss+2,3)}(:).*uint32(SA(:));
    
    tmp4=cellA{h3(ss+1,1),h3(ss+1,2),h3(ss+1,3)}(:).*uint32(SB(:));
    tmp5=cellA{h3(ss+3,1),h3(ss+3,2),h3(ss+3,3)}(:).*uint32(SB(:));
    tmp6=cellA{h3(ss+2,1),h3(ss+2,2),h3(ss+2,3)}(:).*uint32(SB(:));
    
    tmp7=cellA{h3(ss+1,1),h3(ss+1,2),h3(ss+1,3)}(:).*uint32(SC(:));
    tmp8=cellA{h3(ss+3,1),h3(ss+3,2),h3(ss+3,3)}(:).*uint32(SC(:));
    tmp9=cellA{h3(ss+2,1),h3(ss+2,2),h3(ss+2,3)}(:).*uint32(SC(:));
    
    if usealpha
        D2=alphapadded(h(s,1):end+h(s,1)-3,h(s,2):end+h(s,2)-3,h(s,3):end+h(s,3)-3);
        SA2=( D2 < alphadata & D == colordata );

        tmp10=cellA{h3(ss  ,1),h3(ss  ,2),h3(ss  ,3)}(:).*uint32(SA2(:));
        tmp11=cellA{h3(ss+1,1),h3(ss+1,2),h3(ss+1,3)}(:).*uint32(SA2(:));
        tmp12=cellA{h3(ss+2,1),h3(ss+2,2),h3(ss+2,3)}(:).*uint32(SA2(:));

        tmp13=cellA{h3(ss+1,1),h3(ss+1,2),h3(ss+1,3)}(:).*uint32(SA2(:));
        tmp14=cellA{h3(ss+3,1),h3(ss+3,2),h3(ss+3,3)}(:).*uint32(SA2(:));
        tmp15=cellA{h3(ss+2,1),h3(ss+2,2),h3(ss+2,3)}(:).*uint32(SA2(:));
        
        TT=[TT;[[tmp1(tmp1~=0),tmp2(tmp2~=0),tmp3(tmp3~=0)];
        [tmp4(tmp4~=0),tmp5(tmp4~=0),tmp6(tmp6~=0)];
        [tmp7(tmp7~=0),tmp8(tmp8~=0),tmp9(tmp9~=0)];
        [tmp10(tmp10~=0),tmp11(tmp11~=0),tmp12(tmp12~=0)];  
        [tmp13(tmp13~=0),tmp14(tmp14~=0),tmp15(tmp15~=0)];]];
        %配置颜色和透明度:
        CC=[CC;[C(SA(:)); C(SB(:));C(SC(:));C(SA2(:));C(SA2(:))]];
        AA=[AA;[A(SA(:)); A(SB(:));A(SC(:));A(SA2(:));A(SA2(:))]];
    else
        TT=[TT;[[tmp1(tmp1~=0),tmp2(tmp2~=0),tmp3(tmp3~=0)];
        [tmp4(tmp4~=0),tmp5(tmp4~=0),tmp6(tmp6~=0)];
        [tmp7(tmp7~=0),tmp8(tmp8~=0),tmp9(tmp9~=0)];]];
        %assign colors accordingly:
        CC=[CC;[C(SA(:)); C(SB(:));C(SC(:))]]; 
    end
    
end




%free some memory before creating vertex positions:
%在创建顶点位置之前释放一些内存
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9

[XX,YY,ZZ]=meshgrid(...
    linspace(axes(1)-(axes(2)-axes(1))/2/(sX-1),axes(2)+(axes(2)-axes(1))/2/(sX-1),sX+1),...
    linspace(axes(3)-(axes(4)-axes(3))/2/(sY-1),axes(4)+(axes(4)-axes(3))/2/(sY-1),sY+1),...
    linspace(axes(5)-(axes(6)-axes(5))/2/(sZ-1),axes(6)+(axes(6)-axes(5))/2/(sZ-1),sZ+1));

%X, Y and Z coordinates of those indices:
X=XX(:);
Y=YY(:);
Z=ZZ(:);
%draw it all:
if usealpha
    hh = trisurf(TT,X,Y,Z,CC,'EdgeColor','none','FaceVertexAlphaData',AA,'AlphaDataMapping','direct','FaceAlpha','flat'); 
else
    hh = trisurf(TT,X,Y,Z,CC,'EdgeColor','none','FaceAlpha',a); 
end
%create a lightsource in the middle of all six sides of the data cube,
%在数据立方体的所有六个面中间创建一个光源，
%outset by one cube diameter:
%从一个立方体直径开始


material([0.2 0.7 0.3]) 
daspect([1 1 1]);
pbaspect([1 1 1]);
camproj('perspective')


%重新绘制光源
delete(findall(gcf,'Type','light'))
light('Position',[   (axes(2)-axes(1))/2      (axes(4)-axes(3))/2      axes(5)-(axes(6)-axes(5))],'Style','local')
light('Position',[   (axes(2)-axes(1))/2      (axes(4)-axes(3))/2      axes(6)+(axes(6)-axes(5))],'Style','local')
light('Position',[   (axes(2)-axes(1))/2    axes(3)-(axes(4)-axes(3))     (axes(6)-axes(5))/2   ],'Style','local')
light('Position',[   (axes(2)-axes(1))/2    axes(4)+(axes(4)-axes(3))     (axes(6)-axes(5))/2   ],'Style','local')
light('Position',[axes(1)-(axes(2)-axes(1))   (axes(4)-axes(3))/2         (axes(6)-axes(5))/2   ],'Style','local')
light('Position',[axes(2)+(axes(2)-axes(1))   (axes(4)-axes(3))/2         (axes(6)-axes(5))/2   ],'Style','local')

