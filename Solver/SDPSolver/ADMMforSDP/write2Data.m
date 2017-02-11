%this function writes data in the file
function [] = write2Data(fileName,d1,d2,d3,d4,d5,d6,d7)

fileID = fopen(fileName,'w');

for i = 1:length(d1)
  fprintf(fileID,'%10.9f %10.9f %f %f %f %f %f\n',d1(i),d2(i),d3(i),d4(i),d5(i),d6(i),d7(i));
end
fclose(fileID);
