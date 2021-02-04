function rgb = ek500()
%  ek500 - return a 13x3 array of RGB values representing the Simrad EK500 color table

    rgb = [[255 255 255];
           [159 159 159];
           [095 095 095];
           [000 000 255];
           [000 000 127];
           [000 191 000];
           [000 127 000];
           [255 255 000];
           [255 127 000];
           [255 000 191];
           [255 000 000];
           [166 083 060];
           [120 060 040]];

   rgb = rgb / 255.;
   
end