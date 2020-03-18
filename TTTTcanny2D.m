I = imread('QQ20180604-213843.png'); 
x=rgb2gray(I);
edge(x,'canny'); 
