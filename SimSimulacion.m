% Imagenes Simuladas
ImOriginal=rgb2gray(imread('USAFTESIS.jpg')); % Lectura de la imagen que se utilizará en la simulación
Im=double(imresize(ImOriginal,[512 512]))/255; clear ImOriginal; %Cambia el tamaño de la imagen a 512x512

pxsize=0.05; %Tamaño de pixel


%Coordenadas Espacio Real
n=512;m=n;
Lx=n*pxsize;
Ly=m*pxsize;

[X1,Y1]=meshgrid(-Lx/2:pxsize:Lx/2-pxsize,-Ly/2:pxsize:Ly/2-pxsize);
%Coordenadas Espacio Recíproco
fx=-(1/(2*pxsize)):1/Lx:1/(2*pxsize)-1/Lx;
fy=-(1/(2*pxsize)):1/Ly:1/(2*pxsize)-1/Ly;
[Fx,Fy]=meshgrid(fx,fy);

theta=[0 135 90 45 ]; % ángulos del patrón de iluminación
nangles=4;
nphases=3;
phi=linspace(0,pi,nphases+1); %
phi=phi(1:nphases);           %  El vector phi contiene las fases de cada patrón de iluminación (PI)
H=OTF2D(n,m,0,0,pxsize);      %  Generar OTF, el tamaño de ésta se controla en la función OTF2D  
PSF=otf2psf(H);               %  Generar PSF del sistema
%plot(-Lx/2:pxsize:Lx/2-pxsize,abs(PSF(257,:)));axis([-3 3 0 inf]) 




ImF=fftshift(fft2(Im)).*H; % Acción del sistema óptico sobre la imagen en el espacio recíproco
ImOTF=ifft2(ImF); % Imagen formada por el sistema óptico
ki=0.65*2;        % Frecuecia espacial del patrón de iluminación

%Generación del patrón de iluminación
for q=1:nangles                          
    for i=1:nphases
        k = rotxy(deg2rad(theta(q)))*[ki; 0]; 
        kx = k(1); ky = k(2); 
        Ir(:,:)=(1+cos(2*pi*(kx*X1+ky*Y1)+phi(i)))/2; % Patrón de iluminación
        temporal=fftshift(fft2(Im.*Ir)).*H; % Acción del sistema óptico sobre el objeto iluminaco por el PI
        stack(:,:,(q-1)*3+i)=abs(ifft2(temporal)); %Stack es un cubo de datos donde se guardan las imágenes simuladas
        imshow(stack(:,:,(q-1)*3+i),[])
        pause(0.2)
    end
end

%%
% temp=stack;clear stack;
% 
% for j=1:12
%   stack(:,:,j)=temp(:,:,j)';  
% end
% clear temp;
%%
[theta,A]=DET_ANG_FREC(stack(:,:,:),512,60,20);% La función DET_ANG_FREC determina los ángulos
% y periodos de los patrones de iluminación utilizados. Mediante la
% localizacióon de los máximos en el espacio recíproco

%%


%%PROCESAMIENTO
pxsize=0.05;
cubo=(Recorte_Cubo(stack,512,1));

for j=1:12
    im=cubo(:,:,j);
%   cubo(:,:,j)=(im-min(im(:)))/(max(im(:))-min(im(:)))
    cubo(:,:,j)=mat2gray(im);
end



[tx,ty,tz]=size(cubo);%Tamaño del cubo de datos de las imágenes

tx=pow2(nextpow2(tx))%

%Zero-Padding de las imagenes a procesar,esto es, aumnentar el tamaño de la imagen
%mediante la incorporación de un borde de zeros, en este caso la imagen
%pasa de 512x512 a 1024x1024. El fin de esto es aumentar la densidad de
%puntos en el espacio de frecuencias.
for i=1:12
    
    cubo2(:,:,i)=padarray(cubo(:,:,i),[tx/2 tx/2]); %Agrega borde de ceros (zero-padding)
    
    TF_Im_Todas(:,:,i)=(fftshift(fft2(cubo2(:,:,i)))); % Cubo con las transformadas de Fourier de todas las imágenes a tratar.
    
end
% clear cubo 


[n,m]=size(TF_Im_Todas(:,:,1));

Lx=n*pxsize;
Ly=m*pxsize;


%Nuevas Coordenadas Espacio Real (Recordemos que el zero-padding cambió el tamaño de las imágenes)
[X1,Y1]=meshgrid(-Lx/2:pxsize:Lx/2-pxsize,-Ly/2:pxsize:Ly/2-pxsize);

%Nuevas Coordenadas Espacio Recíproco
fx=-(1/(2*pxsize)):1/Lx:1/(2*pxsize)-1/Lx;
fy=-(1/(2*pxsize)):1/Ly:1/(2*pxsize)-1/Ly;
[Fx,Fy]=meshgrid(fx,fy);

nangles=4;
nphases=3;


ki=[ norm(A(1,:)) norm(A(2,:)) norm(A(3,:)) norm(A(4,:)) ]*(1/Lx); %Vector con la  frecuencia del PI por cada ángulo.



%DECONVOLUCIÓN 
SUM1=0;
SUM2=0;
ind=[0 -1 1]; % Estos números representarán en el algoritmo laa componentea (-1, 0, -1) sobre la cual se aplican las funciones
sp = zeros(n,m,nphases*nangles); % Cubo de datos que guardará las componentes separadas

p=0;
qq=0;
FASES=zeros(1,nphases);
clear Ir
% FASES=(phi);
for q=1:nangles
 k = rotxy(deg2rad(theta(q)))*[ki(q); 0]; 
 kx = k(1); ky = k(2); 

     Ir(:,:)=exp(-1i*2*pi*(-kx*X1+ky*Y1)); % Factor de desplazamiento

for l=1:3
    %Determinación de la fase del patrón mediante las autocorrelaciones en el espacio recíproco
    A1=TF_Im_Todas(:,:,l+qq).*conj(OTF2D(n,m,0,0,pxsize));
    A2=fft2(ifft2((TF_Im_Todas(:,:,l+qq))).*Ir(:,:)).*conj(OTF2D(n,m,0,0,pxsize)); 
    AAA=(A1.*conj(A2)); % Autocorrelación 
    FASES(1,l)=-angle(sum(AAA(:)));%

    
end
 rad2deg(mod(FASES,2*pi))% Deja las fases expresadas entre 0 y 2pi
%  Fases(:,q)=wrapTo360(radtodeg(FASES))
qq=qq+3;



  M(:,:)=0.5*[2 exp(-1i*FASES(1,1)) exp(1i*FASES(1,1)) ;2 exp(-1i*(FASES(1,2)))  exp(1i*(FASES(1,2)));...
       2 exp(-1i*(FASES(1,3)))  exp(1i*(FASES(1,3)))];
    
  
  Inv_M=inv(M); % Inversa de la matriz de fases, esta se utiliza para separar las 3 componentes presentes en cada imagen.
  
    for mm=1:nphases
        
     temp_separated = zeros(n,m,nphases); 
     Ir(:,:)=exp(-1i*2*pi*(-kx*ind(mm)*X1+ky*ind(mm)*Y1)); 
     
             for k = 1:nphases %!!!!! 
                
              temp_separated(:,:,k) = Inv_M(mm,k).*TF_Im_Todas(:,:,k+p); 
              
           
              sp(:,:,(q-1)*3+mm) = (sp(:,:,(q-1)*3+mm)+temp_separated(:,:,k));  % Guarda las componetnes separadas
              
              
              end
%               gca= imshow(log(abs(sp(:,:,(q-1)*3+mm))),[]);
%               filename=['sp',num2str((q-1)*3+mm)];
%               fpath='D:\Desktop\Gimp SIM';              
%               saveas(gca, fullfile(fpath, filename), 'tiff');
              
             
   Cmov1(:,:)= gauss_atenua(n,m,ind(mm)*kx,ind(mm)*ky,pxsize).*sp(:,:,(q-1)*3+mm); 
    Cmov(:,:)=fft2((ifft2(Cmov1 ).*Ir(:,:))).*conj(OTF2D(n,m,ind(mm)*kx,ind(mm)*ky,pxsize)); 
     T=((q-1)*3+mm);
     
     if   ( T==2||T==5||T==8||T==11 )  
    %Determinación FASE GLOBAL entre componentes 0 y 1,-1  
%          xcore=xcorr2(sp(:,:,(q-1)*3+mm-1),Cmov(:,:));
          xcore=sp(:,:,(q-1)*3+mm-1).*conj(Cmov(:,:));
          dim_xcore=size(xcore);
%          offset=angle(xcore(floor(dim_xcore(1)/2)+1,floor(dim_xcore(2)/2)+1));
         
         [M2,I2] = max(xcore(:));

         [I_row, I_col] = ind2sub(size(xcore),I2);
         offset=angle(xcore(I_row,I_col));
     elseif (T==1||T==4||T==7||T==10)
         offset=0;
     end
     
     SUM1=(exp(-1i*ind(mm)*offset)*Cmov) + SUM1;
%       
     SUM2= (OTF2D(n,m,ind(mm)*kx,ind(mm)*ky,pxsize).*conj(OTF2D(n,m,ind(mm)*kx,ind(mm)*ky,pxsize))).*gauss_atenua(n,m,ind(mm)*kx,ind(mm)*ky,pxsize) +  SUM2;
            
        
   
    end
    
    p=p+3;
    
end

%Función de Apodización
% 
% Apo=(OTF2D(n,m,0,0,pxsize));
% 
% Apo (Apo ~= 0)    = 0.1;
% Apo=bwdist(~Apo);
% maxOTF=max(max(abs(Apo)));
% minOTF=min(min(abs(Apo)));
% Apo=(Apo-minOTF)/maxOTF;

k_r = sqrt(Fx.^2+Fy.^2);
k_max = (0.6*4)^2;
bhs = cos(pi*k_r/(2*k_max));
indi = find( k_r > k_max ); 
bhs(indi) = 0;


S=(  SUM1./(SUM2+0.3*(max(abs(SUM2(:)))))) ;%.*(Apo.^0.4).*bhs
RS=(ifft2(S));
RS=(abs(RS));
Wide=(cubo2(:,:,1)+cubo2(:,:,2)+cubo2(:,:,3))/3;
% Wide=padarray(abs(ImOTF),[tx/2 tx/2]);
FIG=figure(1)
subplot(1,3,1)
imshow(Wide,[]);title(['Original, ','NA=1.2, ','ResMáx=',num2str(0.61*0.61/1.2),'um.'])
subplot(1,3,2)
imshow(RS);title(['Reconstrucción , ','FrecPatrón=',num2str(ki(1)^-1),'um.'])
subplot(1,3,3)
imshow(log(abs(S)),[]);

