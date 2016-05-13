function closed = is_closed(v,t)

vnum = size(v,1);
tnum = size(t,1);
vdeg = accumarray(t(:),1);
maxdeg = max(vdeg);

adj = spalloc(vnum,vnum,vnum*maxdeg);

for ti = 1:tnum
    for ei = 1:3
        v1 = t(ti,ei);
        v2 = t(ti,mod(ei,3)+1);
        
        adj(v1,v2) = adj(v1,v2) +1;
        adj(v2,v1) = adj(v2,v1) +1;
    end
end

m = max(adj(:));
if (m > 2)
    error('Mesh not manifold');
end

if ~ isempty (find(adj == 1))
    closed = 0;
    return 
else
    closed = 1;
    return 
end
