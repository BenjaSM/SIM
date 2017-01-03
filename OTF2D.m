 
function [ otf ] = OTF2D(n,m,kx,ky,pxsize)

fx=-(1/(2*pxsize)):1/(n*pxsize):1/(2*pxsize)-1/(n*pxsize);
fy=-(1/(2*pxsize)):1/(m*pxsize):1/(2*pxsize)-1/(m*pxsize);
r0=0.91;
[X, Y] = meshgrid(fx,fy);

[angle,mag] = cart2pol(X-kx,Y+ky);

idx = mag <= 2*r0;

otf = zeros(n,m);

dum =  2/pi*(acos(mag/(2*r0)) - (mag/(2*r0)).*sqrt(1-(mag/(2*r0)).^2));
otf(idx) = dum(idx);


end