function [A,B] = computeAB(v,t)

vnum = size(v,1);
tnum = size(t,1);
vdeg = accumarray(t(:),1);
maxdeg = max(vdeg);

A = spalloc(vnum,vnum,vnum*maxdeg);
B = spalloc(vnum,vnum,vnum*maxdeg);

% linear local matrices on unit triangle:
tB = (ones(3,3) + eye(3) )./24.0;

tA00 = [ 0.5 -0.5  0.0;
        -0.5  0.5  0.0;
         0.0  0.0  0.0];
 
tA11 = [ 0.5  0.0 -0.5;
         0.0  0.0  0.0;
        -0.5  0.0  0.5];
    
tA0110 = [ 1.0 -0.5 -0.5;
          -0.5  0.0  0.5;
          -0.5  0.5  0.0];

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
    a0 = v3mv1 * v3mv1';
    a1 = v2mv1 * v2mv1';
    a0110 = v2mv1 * v3mv1';
    cr = cross(v2mv1,v3mv1);
    vol = norm(cr,2);
    localB = vol * tB;
    localA = (1.0/vol) * (a0*tA00 + a1*tA11 - a0110 * tA0110);

    % add to global matrices:
    A([i1,i2,i3],[i1,i2,i3]) = A([i1,i2,i3],[i1,i2,i3]) + localA;
    B([i1,i2,i3],[i1,i2,i3]) = B([i1,i2,i3],[i1,i2,i3]) + localB;
    
end
