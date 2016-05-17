% Martin Reuter
% Laplace-Beltrami Eigenstuff Exercises
% %% Theory 1D
% % Exercise 1
% fprintf('Exercise 1\n');
% % Plot the Eigenfunctions of Laplace Operator on the Line Segment [0,a]
% % neumann
% a = 1;
% x = 0:0.02:a;
% evn1 = ones(1,length(x)); % constant for n=0 
% evn2 = cos(1*x*pi/a); % n=1
% evn3 = cos(2*x*pi/a); % n=2
% figure; 
% plot(x,[evn1;evn2;evn3]);
% legend('evn1','evn2','evn3');
% figure;  
% plot(evn2,evn3,'k.-'); 
% % dirichlet
% evd1 = sin(1*x*pi/a);
% evd2 = sin(2*x*pi/a);
% figure; 
% plot(x,[evd1;evd2]);
% legend('evn1','evn2');
% figure; 
% plot(evd1,evd2,'k.-'); 
% % note that -1 *efunc is also an efunc
% fprintf('Press Enter to continue...\n');
% pause;

%% Linear FEM on line
% Exercise 2
fprintf('Exercise 2\n');
% create spiral
phi = 0:0.2:4*pi; 
vnum = length(phi); 
x = phi.*cos(phi);
y = phi.*sin(phi); 
figure;
plot(x,y,'b.-');
% integrals of local form functions
a = [1 -1; -1 1];
b = [1/3 1/6; 1/6 1/3];
fprintf('Press Enter to continue...\n');
pause;

% Exercise 3
fprintf('Exercise 3\n');
% allocate sparse (filled with zero) 
A = spalloc(vnum,vnum,2*vnum); 
B = spalloc(vnum,vnum,2*vnum); 
for i = 2:vnum
    % global index of first and second vertex of this edge:
    i1 = i-1; 
    i2 = i;
    % length of segment
    W = sqrt((x(i2)-x(i1))^2+(y(i2)-y(i1))^2);
    % fill local stuff into global A matrix
    A(i1,i1) = A(i1,i1)+a(1,1)/W;
    A(i1,i2) = A(i1,i2)+a(1,2)/W;
    % same as line before: symmetric
    A(i2,i1) = A(i2,i1)+a(2,1)/W;
    A(i2,i2) = A(i2,i2)+a(2,2)/W;
    % or simply
    % A([i1,i2],[i1,i2]) = A([i1,i2],[i1,i2]) + (1.0/W) * a;
    B([i1,i2],[i1,i2]) = B([i1,i2],[i1,i2])+W*b;
end
fprintf('Press Enter to continue...\n');
pause;

% Exercise 4
fprintf('Exercise 4\n');
% compute efunc evals:
opts.issym = 1;
N = 3; % need the first 2 non zero 
[evecsn, evalsn] = eigs(A,B,N,'SM',opts);
% reverse order (smallest first)
evalsn = diag(evalsn);
evalsn = evalsn(end:-1:1)
evecsn = evecsn(:,end:-1:1);
% plot first non-zero efunc as color
figure;
surface([x(:) x(:)],[y(:) y(:)],[zeros(length(x),2)],[evecsn(:,2) evecsn(:,2)],...
    'FaceColor','none','EdgeColor','flat','Marker','none','LineWidth',2);
fprintf('Press Enter to continue...\n');
pause;

% Exercise 5
fprintf('Exercise 5\n');
% embedding (onto 2 and 3, as first is constant)
figure
plot(evecsn(:,2),evecsn(:,3),'b.-');
fprintf('Press Enter to continue...\n');
pause;

% Exercise 6
fprintf('Exercise 6\n');
% dirichlet (remove first and last row and column) 
Ad = A(2:end-1,2:end-1);
Bd = B(2:end-1,2:end-1);
% compute efunc evals
opts.issym = 1;
N = 2;
[evecsd, evalsd] = eigs(Ad,Bd,N,'SM',opts);
% reverse order (smallest first) 
evalsd = diag(evalsd);
evalsd = evalsd(end:-1:1)
evecsd = evecsd(:,end:-1:1);
% zero for boundary nodes
evecsd = [zeros(1,size(evecsd,2));evecsd;zeros(1,size(evecsd,2))];
fprintf('Press Enter to continue...\n');
pause;

% Exercise 7
fprintf('Exercise 7\n');
% embedding (onto 1 and 2) 
figure;
plot(evecsd(:,1),evecsd(:,2),'b.-');
% embedding (onto 3neuman and 2dirichlet) 
% to create circle
figure;
plot(evecsn(:,3),evecsd(:,2),'b.-');
axis equal;
fprintf('Press Enter to continue...\n');
pause;

%% 2D Square Theory
% Exercise 8
fprintf('Exercise 8\n');
% square Dirichlet
a = 1; 
x = 0:0.02:a; 
[X,Y] = meshgrid(x,x); 
% one dim eigenspace
Z = sin(1*X*pi./a).*sin(1*Y*pi./a); 
figure;
surface(X,Y,zeros(size(X)),Z,'FaceColor','interp','EdgeColor','none','Marker','none'); 
axis square
% 2dim eigenspace space
Z = sin(2*X*pi./a).*sin(1*Y*pi./a); 
figure;
surface(X,Y,zeros(size(X)),Z,'FaceColor','interp','EdgeColor','none','Marker','none'); 
axis square;
Z = sin(1*X*pi./a).*sin(2*Y*pi./a); 
figure;
surface(X,Y,zeros(size(X)),Z,'FaceColor','interp','EdgeColor','none','Marker','none'); 
axis square;
% one dim eigenspace
Z = sin(2*X*pi./a).*sin(2*Y*pi./a); 
figure;
surface(X,Y,zeros(size(X)), Z,'FaceColor','interp','EdgeColor','none','Marker','none'); 
axis square;
fprintf('Press Enter to continue...\n');
pause;

%% Meshes
% Exercise 9
fprintf('Exercise 9\n');
[v,t] = OFF_Import('gorilla18.off'); 
% bbox
cmin = min(v);
cmax = max(v);
%plot mesh (color height function)
figure;
h = trisurf(t,v(:,1),v(:,2),v(:,3),'EdgeColor','none','FaceColor','interp');
light('Position',[10 0 0],'Style','infinite'); 
light('Position',[-10 0 0],'Style','infinite');
lighting phong
axis tight equal off
zoom(1)
% vector with vertex degrees
vdeg = accumarray(t(:),1); 
maxdeg = max(vdeg)
is_closed(v,t) % also checks if manifold, look at source
fprintf('Press Enter to continue...\n');
pause;

% Exercise 10
fprintf('Exercise 10\n');
% construct matrices (for closed meshes or Neumann Bcond)
[A,B] = computeAB(v,t);
% solve sparse symmetric generalized eigenvalue problem
opts.issym = 1;
N = 800; % will take a while
[evecs,evals] = eigs(A,B,N,'SM',opts);
% evecs are normalized with respect to B 
% evecs' * B * evecs = id
% why with respect to B?
% reverse order (smallest first)
evals = diag(evals);
evals = evals(end:-1:1);
evecs = evecs(:,end:-1:1);
fprintf('Press Enter to continue...\n');
pause;

% Exercise 11
fprintf('Exercise 11\n');
% plot eigenvalues
figure;
plot((1:N),evals,'k-');
% plot first non-constant eigenfunction
figure;
trisurf(t,v(:,1),v(:,2),v(:,3),evecs(:,2),'EdgeColor','none','FaceColor','interp'); 
light; 
lighting phong; 
axis tight equal off; 
zoom(1)
% same with
figure;
trisurf(t,v(:,1),v(:,2),v(:,3),evecs(:,3),'EdgeColor','none','FaceColor','interp');
figure;
trisurf(t,v(:,1),v(:,2),v(:,3),evecs(:,4),'EdgeColor','none','FaceColor','interp');
% project into 3d spectral domain
figure;
trisurf(t,evecs(:,2),evecs(:,3),evecs(:,4),evecs(:,2),'EdgeColor','none','FaceColor','interp'); 
light; 
lighting phong; 
axis tight equal off;
zoom(1)
fprintf('Press Enter to continue...\n');
pause;

% Exercise 12
fprintf('Exercise 12\n');
% add noise to right hemisphere
[tn, vn] = get_normals(v,t);
pos = find(v(:,1)>0);
sigma = 0.7;
vnoise = v;
vnoise(pos,:) = v(pos,:)+vn(pos,:).*(normrnd(0,sigma,size(pos,1),1)*ones(1,3));
figure;
trisurf(t,vnoise(:,1),vnoise(:,2),vnoise(:,3),0.5*ones(size(vnoise,1),1),... 
    'EdgeColor','none','FaceColor','flat'); 
axis equal;
light('Position',[10 0 0],'Style','infinite'); 
light('Position',[-10 0 0],'Style','infinite');
lighting flat
material dull
zoom(1)
axis tight equal off
fprintf('Press Enter to continue...\n');
pause;

% Exicise 13
fprintf('Exercise 13\n');
% filter geometry
Nmax = 800;
% represent vnoise as lin comb of eigenvectors
% by projecting xyz coords onto eigenbasis (inner product) 
c = (vnoise'*B)*evecs(:,1:Nmax);
% compute linear combination
vsmooth = evecs(:,1:Nmax)*c';
figure;
trisurf(t,vsmooth(:,1),vsmooth(:,2),vsmooth(:,3),'EdgeColor','none','FaceColor','flat'); 
axis equal;
light('Position',[10 0 0],'Style','infinite')
light('Position',[-10 0 0],'Style','infinite') 
lighting flat
axis tight equal off
zoom(1) 
material dull
fprintf('Press Enter to continue...\n');
pause;

%Exercise 14
fprintf('Exercise 14\n');
% filter geometry with heat kernel
Nmax = 800;
time = 0.5;
figure;
plot((1:Nmax),exp(- evals(1:Nmax)*time),'k-') 
hkw = (exp(- evals(1:Nmax)*time)*ones(1,3))'; 
% represent vnoise as lin comb of eigenvectors 
c = (vnoise'*B)*evecs(:,1:Nmax);
% compute filtered linear combination
vsmooth = evecs(:,1:Nmax)*(c.*hkw)';
figure;
trisurf(t,vsmooth(:,1),vsmooth(:,2),vsmooth(:,3),'EdgeColor','none','FaceColor','flat'); 
axis equal
light('Position',[10 0 0],'Style','infinite') 
light('Position',[-10 0 0],'Style','infinite') 
lighting flat
axis tight equal off
zoom(1) 
material dull
fprintf('Press Enter to continue...\n');
pause;

% Exercise 15
fprintf('Exercise 15\n');
% now let's do hot stuff, eg, plot heat kernel at specific locations
time = 400;
head = 932;
rfoot = 6719;
lfoot = 7570;
vfix = lfoot;
h = heatkernel(time,vfix,evecs,evals,Nmax);
figure;
trisurf(t,v(:,1),v(:,2),v(:,3),h(:,1),'EdgeColor','none','FaceColor','interp'); 
axis equal
light('Position',[10 0 0],'Style','infinite') 
light('Position',[-10 0 0],'Style','infinite') 
lighting phong
axis tight equal off
zoom(1) 
material dull
fprintf('Press Enter to continue...\n');
pause;

% Exercise 16
fprintf('Exercise 16\n');
time = [0.5:0.5:100];
x = [1,head,rfoot,lfoot];
hdiag = heatdiag(time,x,evecs,evals,Nmax);
figure;
plot (time,hdiag(:,:))
xlabel('time');
ylabel('K(t,x,x)');
title('Heat Kernel Diagonal at different t');
lcells = strread(sprintf('x = %d, ',x), '%s','delimiter',','); 
legend(lcells);
fprintf('Press Enter to continue...\n');
pause;

% Exercise 17
fprintf('Exercise 17\n');
%heat trace
time = 0.5:0.5:100;
z = sum(exp(-evals(1:Nmax)*time));
figure;
plot(time,z);
xlabel('time');
ylabel('Heat Trace Z(t)');
title('Heat Trace at different t');
fprintf('Press Enter to continue...\n');
pause;

% Exercise 18
fprintf('Exercise 18\n');
% compute 4pi s^2 sum exp(-lambda_i s^2)
s = sqrt(time);
y = 4*pi*s.^2.*sum(exp(-evals(1:Nmax)*s.^2));
plot(s,y);
xlabel('sqrt(t)');
ylabel('4 pi t Z(t)'); 
title('Modified Heat Trace');
% Regression in linear interval to extrapolate 
interval = find(s>6 & s<8);
X = [s(interval);ones(size(interval))]';
B = regress(y(interval)',X);
B(2)
hold on; 
plot([0,s],[0,s]*B(1)+B(2),':r');
get_area(v,t)

