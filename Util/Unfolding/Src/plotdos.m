% plotdos.m

clear
wdir = '../Examples/Si/Defects/Vacancy/';
fname = [wdir,'unfoldedBandLine1.out'];

% Read DOS
file = fopen(fname);
dat = fscanf(file,'%d %d %f %f',4);
nq = dat(1);
ne = dat(2);
emin = dat(3);
emax = dat(4);
q = zeros(3,nq);
dos = zeros(nq,ne);
for iq = 1:nq
    q(:,iq) = fscanf(file,'%f %f %f',3);
    dos(iq,:) = fscanf(file,'%f',ne);
end
fclose(file)

% dos=load([wdir,'unfoldedBandLine1.out']);
% q = dos(:,1:3)';
% dos = dos(:,4:end);
% [nq,ne] = size(dos);

if nq==1
    qpath = 0;
else
    for ix=1:3
        dq(ix,:) = diff(q(ix,:));
    end
    dq = sqrt(sum(dq.^2));
    qpath = [0,cumsum(dq)];
end
de = (emax-emin)/(ne-1);
e = emin:de:emax;
[qpath,e] = ndgrid(qpath,e);
if nq==1
    plot(e,dos)
%    plot(e,cumsum(dos))
    xlabel('energy')
    ylabel('dos')
    grid on
else
    surf(qpath',e',dos',min(dos',1.0))
%    surf(qpath',e',cumsum(dos,2)',min(cumsum(dos,2)',2))
%    surf(qpath',e',log10(dos+1e-3)')
    xlabel('qpath')
    ylabel('energy')
    zlabel('dos')
    axis([0,qpath(end),emin,30,0,2])
end
colorbar
