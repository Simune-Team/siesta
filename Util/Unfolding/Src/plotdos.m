% plotdos.m
% Plot qDOS writen by unfold
% SGM & JMS, Oct.2018

clear
wdir = '../Examples/Si/Bulk/Si8/';
syslabel = 'si8';
fname = [wdir,syslabel,'.unfoldedBands'];

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
    iline(iq) = fscanf(file,'%i',1);
    dos(iq,:) = fscanf(file,'%f',ne);
end
fclose(file);

% Plot qDOS
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
%    s = surf(qpath',e',dos')
%    s = surf(qpath',e',cumsum(dos,2)')
    s = surf(qpath',e',max(log10(dos),-2)');
    s.EdgeColor = 'none';
    xlabel('qpath')
    ylabel('energy')
    zlabel('dos')
%    axis([0,qpath(end),emin,emax,0,1/de])
%    axis([0,qpath(end),emin,emax,-2,0.5])
end

% Set viewpoint and colormap
view(2);
colormap default
% cmap = colormap(bone);    % white on black
% cmap = cmap(end:-1:1,:);  % black on white
% colormap(cmap);
colorbar

% Plot line edges
lineEnd = find(iline(1:nq-1)<iline(2:nq));
hold on
for jl = 1:numel(lineEnd)
    il = lineEnd(jl);
    plot(qpath(il)*[1,1],[emin,emax])
end
hold off
