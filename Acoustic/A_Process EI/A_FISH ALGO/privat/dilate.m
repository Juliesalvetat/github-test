function matrix2=dilate(matrix1,mask)

%mask: masque de dilatation (de 0 et 1)
%matrix1: matrice à éroder (de 0 et 1)

smask=sum(mask);
n=size(mask,1);
m=size(mask,2);

cor=xcorr2(matrix1,mask);
matrix2=double(cor>=1);

%on retire l'extension de dimension due à la correlation
%matrix2=matrix2(1+ceil((n-1)/2):end-floor((n-1)/2),1+ceil((m-1)/2):end-floor((m-1)/2));
matrix2=matrix2(1+floor((n-1)/2):end-ceil((n-1)/2),1+floor((m-1)/2):end-ceil((m-1)/2));