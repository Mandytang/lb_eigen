function [tn, vn] = get_normals(v,t)

vnum = size(v,1);
tnum = size(t,1);

tn = zeros(tnum,3);
vn = zeros(vnum,3);

for ti = 1:tnum
    
    % get local matrices for this element:
    i1 = t(ti,1);
    i2 = t(ti,2);
    i3 = t(ti,3);
    v1 = v(i1,:);
    v2 = v(i2,:);
    v3 = v(i3,:);
    v2mv1 = v2-v1;
    v3mv1 = v3-v1;
    
    cr = cross(v2mv1,v3mv1);
    cr = (1.0/norm(cr,2)) * cr;
    tn(ti,:) = cr;
    
    vn(i1,:) = vn(i1,:) + cr;
    vn(i2,:) = vn(i2,:) + cr;
    vn(i3,:) = vn(i3,:) + cr;
  
end

for vi = 1:vnum
    vn(vi,:) = vn(vi,:) ./ norm(vn(vi,:),2);
end
