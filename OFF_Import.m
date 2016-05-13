function [v,t] = OFF_Import(filename)
% OFF_Import is a tool designed to import into MATLAB both binary and ASCII OFF files.
% 
% 
% SYNOPSIS:
% 
% INPUT:
% 
% filename: string representing the name fo the file
% 

%open file
fid=fopen(filename, 'r'); %Open the file, assumes OFF ASCII format.
if fid == -1
    error('File could not be opened, check name or path.')
end


vnum=0;

if (feof(fid) ~= 0) 
    error ('File empty?')
end

tline = fgetl(fid);                 % reads a line of data from file.
fword = sscanf(tline, '%s ');       % make the line a character string
if ~strncmpi(fword, 'off',3)       % Checking if a "OFF" line
    exit('File does not start with "OFF" tag.')
end

tline = fgetl(fid);                 % reads a line of data from file.
while (tline(1) == '#' && feof(fid)==0)
    tline = fgetl(fid);
end
sizes =sscanf(tline,'%i');
vnum=sizes(1)
tnum=sizes(2)
v=zeros(vnum,3);
vdeg=zeros(vnum);
t=zeros(tnum,3);
vcount = 0;
tcount = 0;
while (feof(fid) == 0 && vcount < vnum)  % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    if (tline(1) == '#')
        continue;
    end
    vcount = vcount+1;
    v(vcount,:) = sscanf(tline, '%f')';       % make the line a vector
end

if (vcount ~= vnum)
    exit('Vertex count not matching');
end

while (feof(fid) == 0 && tcount < tnum)  % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    if (tline(1) == '#')
        continue;
    end
    tcount = tcount+1;
    temp = sscanf(tline, '%i')';       % make the line a vector
    if ( temp(1) ~= 3 )
        exit('Faces are not triangles?')
    end
    t(tcount,:) = temp(2:4)';
end

if (tcount ~= tnum)
    exit('Triangle count not matching');
end

if (min(t(:)) == 0 && max(t(:)) == vnum-1)
    t = t+1;
end



