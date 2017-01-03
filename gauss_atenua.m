function G0=gauss_atenua(n,m,off_x,off_y,pxsize)

d=0.2;
a=0.8;
kx=off_x;
ky=off_y;
fx=-(1/(2*pxsize)):1/(n*pxsize):1/(2*pxsize)-1/(n*pxsize);
fy=-(1/(2*pxsize)):1/(m*pxsize):1/(2*pxsize)-1/(m*pxsize);

[Fx,Fy]=meshgrid(fx,fy);



G0=1-a*exp(-( (sqrt((Fx-kx).^2+(Fy+ky).^2)).^2)/(2*d^2));




end