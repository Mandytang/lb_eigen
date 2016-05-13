function a = get_area(v,t)
vnum = size(v,1);
tnum = size(t,1);

a = 0;

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
    a = a + norm(cr,2)/2;
end

