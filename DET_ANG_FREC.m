function [sumaANG,A]=DET_ANG_FREC(stack,TamX,rad,extra)

sumANG=zeros(1,4);

    
    
%%PROCESAMIENTO
cubo=Recorte_Cubo(stack,TamX,1);
% figure(2);imshow(abs(cubo(:,:,1)))

for j=1:12
    im=cubo(:,:,j);
%     cubo(:,:,j)=(im-min(im(:)))/(max(im(:))-min(im(:)))
cubo(:,:,j)=mat2gray(im);
end



[tx,ty,tz]=size(cubo);%
tx=pow2(nextpow2(tx));%

for i=1:12
    
    cubo2(:,:,i)=padarray(cubo(:,:,i),[tx/2 tx/2]);
%             cubo2(:,:,i)=cubo(:,:,i);

    TF_Im_Todas(:,:,i)=(fftshift(fft2(cubo2(:,:,i))));
    
end

sTF=size(TF_Im_Todas(:,:,1));
cont=1;
A=zeros(4,3);
daf=floor(sTF(1)/2)+1;
mask(:,:,1)=ones(sTF(1),sTF(2))-circ(rad,size(TF_Im_Todas(:,:,1)),[daf daf])-circ(rad+extra,size(TF_Im_Todas(:,:,1)),[daf daf]);
mask(:,:,2)=ones(sTF(1),sTF(2))-circ(rad,size(TF_Im_Todas(:,:,1)),[daf daf])-circ(rad+extra,size(TF_Im_Todas(:,:,1)),[daf daf]);
mask(:,:,3)=ones(sTF(1),sTF(2))-circ(rad,size(TF_Im_Todas(:,:,1)),[daf daf])-circ(rad+extra,size(TF_Im_Todas(:,:,1)),[daf daf]);
mask(:,:,4)=ones(sTF(1),sTF(2))-circ(rad,size(TF_Im_Todas(:,:,1)),[daf daf])-circ(rad+extra,size(TF_Im_Todas(:,:,1)),[daf daf]);
mask (mask == 1)=8;
mask (mask==-1)=8;
mask(mask==0)=1;
mask(mask==8)=0;


mask(:,1:daf,1)=0;
mask(daf-4:end,:,2)=0;mask(:,1:daf,2)=0;
mask(daf+1:end,:,3)=0;
mask(daf+1:end,:,4)=0;

for i=1:3:12
    
    TF_Bloqueada=TF_Im_Todas(:,:,i).*(mask(:,:,cont));
   imtool(abs(TF_Bloqueada))
   
    [valor,point]=((max(abs(TF_Bloqueada(:)))));
    [R,C] = ind2sub(size(TF_Bloqueada),point);
    A(cont,:)=([C-daf daf-R 0]);
    cont=cont+1;
end

% imtool(abs(TF_Im_Todas(:,:,12)))
B=[1 0 0];


angle1_rad=atan2(norm(cross(A(1,:),B)),dot(A(1,:),B));
angle2_rad=atan2(norm(cross(A(2,:),B)),dot(A(2,:),B));
angle3_rad=atan2(norm(cross(A(3,:),B)),dot(A(3,:),B));
angle4_rad=atan2(norm(cross(A(4,:),B)),dot(A(4,:),B));

theta=zeros(1,4) ;

theta(1,1) = rad2deg(angle1_rad);
ang=theta(1,1);

if ang> 90
    theta(1,1)=ang-180;
end

theta(1,2) = rad2deg(angle2_rad)
theta(1,3) = rad2deg(angle3_rad)
theta(1,4) = rad2deg(angle4_rad)
sumaANG=theta(1,:)+sumANG(1,:) ;





end
