function matrix2=erode(matrix1,mask)

%mask: masque d'erosion (de 0 et 1)
%matrix1: matrice à éroder (de 0 et 1)

smask=sum(sum(mask,1),2);
n=size(mask,1);
m=size(mask,2);

cor=xcorr2(matrix1,mask);
matrix2=double(cor>=smask);

%on retire l'extension de dimension due à la correlation
%matrix2=matrix2(1+ceil((n-1)/2):end-floor((n-1)/2),1+ceil((m-1)/2):end-floor((m-1)/2));
matrix2=matrix2(1+floor((n-1)/2):end-ceil((n-1)/2),1+floor((m-1)/2):end-ceil((m-1)/2));