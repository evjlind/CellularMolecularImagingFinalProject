% Read the image and convert to gray-scale
close all
I = imread("blood1-1.png");
figure
imshow(I)
Igray = rgb2gray(I);
Igray = imadjust(Igray);
figure
imshow(I)
% Apply adaptive image threshold
Igray = medfilt2(Igray,[10 10]);
T = adaptthresh(Igray,0.85);
BW = imbinarize(Igray,T);
figure
imshow(Igray)
%%
% Calculate the average background RGB color
Imask = immultiply(I,repmat(BW,[1 1 3]));
Imask = reshape(Imask,[],3);
idx = any(Imask == 0,2);
bgColor = mean(Imask(~idx,:));
% Erase areas connected border
BW = ~BW;
BW = imclearborder(BW);
BW = imfill(BW,'holes');
% Erase noise (small regions)
se = strel('disk',3);
BW = imerode(BW,se);
% Erase some regions where its average RGB is similar to the backgroud RGB
s = regionprops('table',BW,{'Area','PixelIdxList'});
I2 = reshape(I,[],3);
s.AvgRGB = zeros(height(s),1);
for kk = 1:height(s)
regColor = mean(I2(s.PixelIdxList{kk},:));
s.AvgRGB(kk,:) = norm(regColor - bgColor);
end
for kk = 1:height(s)
if s.AvgRGB(kk) <= 10 || s.Area(kk) < 200
  BW(s.PixelIdxList{kk}) = false;
end
end
% Apply watershed to separate overlapped cells
D = bwdist(~BW);
D = -D;
D = imgaussfilt(D, 10);
D = medfilt2(D,[15 15]);
L = watershed(D);
L(~BW) = 0;
BW2 = L > 0;
% Visualize the result
label = label2rgb(L,'jet',[0.5 0.5 0.5]);
s = regionprops('table',BW2,{"BoundingBox",'Eccentricity',"MajorAxisLength","MinorAxisLength","Area","Circularity"});
figure
imshowpair(I,label,'blend')
hold on
sickle = 0;
normal = 0;
for kk = 1:height(s)
if 2500>s.Area(kk)
if s.Area(kk)>250
val = 0.3 - s.Eccentricity(kk)*0.5-s.Circularity(kk)*0.55+(s.MajorAxisLength(kk)/s.MinorAxisLength(kk))*0.4;
if val>0.5
sickle_certainty(kk) = (val-0.5);
rectangle(...
'Position',   s.BoundingBox(kk,:),...
'EdgeColor',  'r', ...
'LineWidth', 2);
text(s.BoundingBox(kk,1),s.BoundingBox(kk,2)+5,...
  strcat('sickle ',num2str(kk)),...
  'FontSize', 6)
sickle = sickle+1;
else
 rectangle(...
'Position',   s.BoundingBox(kk,:),...
'EdgeColor',  'b')
text(s.BoundingBox(kk,1),s.BoundingBox(kk,2)+5,...
  strcat('normal ',num2str(kk)),...
  'FontSize', 6)
normal = normal+1;
normal_certainty(kk) = abs((val+0.5));
end
end
end
end
pct = sickle/(normal+sickle)*100;
disp("number of normal cells")
disp(normal)
disp("number of sickle cells")
disp(sickle)
disp("percentage of sickle cells:")
disp(pct)
