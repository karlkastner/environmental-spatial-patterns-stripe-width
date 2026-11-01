% Wed 29 Apr 13:46:05 CEST 2026
% Karl Kastner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% schematically plot a transect through a striped pattern with rectangular
%% wave shape and the corresponding surface water balance
%
function sw_noisy_oscillator_schematic(fid,pflag)
	if (nargin()<1)
		fid = 1;
	end
	if (nargin()<2)
		pflag = 0;
	end
	fflag = pflag;
	init_plot();
	colour={'r','b'};
	meta = sw_metadata();
	% plot scale
	plotscale = 3.8;
	
	% characteristic wavelength
	lc = 1; 
	% relative stripe width
	pw = 0.4;
	% characteristic stripe width
	wc = pw*lc;
	% wave number
	kc = 2*pi/lc;
	% phase, at which stripes terminate
	phi_c = kc*wc;
	% standard deviation of exogenous heterogneity
	s  = [0,0.1];
	% spatial resolution
	dx = lc/100;
	% spatial extent
	L=(3+(1-pw))*lc;
	% number of grid cells
	nx = L/dx;
	% distance along pattern
	x = linspace(0,L,nx)';
	
	% offset to start with a gap instead of a stripe
	phi0 = 2*pi*(pw);
	
	% deterministic phase along pattern
	phi_det = phi0 + kc*x;
	
	xlim_ = [0,3+(1-pw)];
	
	% initialize random number generator
	% this seed gives a good example
	rng(105)
	
	% spatial heterogeneity
	e = -randn(nx,1);
	% Wiender process (Brownian motion) with unit standard deviation
	W = sqrt(dx)*cumsum(e);
	
	phi = zeros(nx,2);
	for sdx=1:2
	
	% phase deviation
	phi_e = 2*pi*s(sdx)*W;
	% phase, including deviation
	phi(:,sdx) = phi_det + phi_e;
	
	% thresholded pattern
	bt = (mod((phi(:,sdx))/kc,1)<=pw);
	% surface water
	water = interp1([0,lc*pw,lc],[1,0,1],mod((phi(:,sdx))/kc,lc),'linear');
	
	splitfigure([2,3],[fid,sdx],fflag);
	yyaxis left
	cla
	a=area(x-0*phi0/(2*pi),bt,'facecolor',meta.colour.lightgreen,'edgecolor','none');
	hold on
	xlim(limits(x))
	ylim([0,1.2])
	set(gca,'ytick',[])
	xlabel('Position $x/\lambda_c$');
	ylabel('Biomass density $b$');
	set(gca,'ycolor',[0,0.5,0])
	yyaxis right
	cla
	plot(x+0*phi0/(2*pi),water,'-b','linewidth',meta.linewidth);
	ylabel('Flow depth $\int(R-I)\mathrm{d}x$');
	set(gca,'ycolor',meta.colour.blue)
	ylim([0,1.4*1.2])
	set(gca,'ytick',[])
	
	splitfigure([2,3],[fid,3+sdx],fflag);
	cla();
	% plot diagonal stripes in phase space
	for idx=1:5
		patch([-1,-1+pw,1,1]+(idx-2)+(1-pw),[1,1,-1+pw,-1],[0.5,0.8,0.5]);
		hold on;
	end
	hline(0);
	% proability cone
	y=s(sdx)*sqrt(x);
	area(x,y,'facecolor',[0.5,0.5,1],'facealpha',0.5);
	area(x,-y,'facecolor',[0.5,0.5,1],'facealpha',0.5);
	% plot phase deviation
	if (1==sdx)
	plot(x,phi_e/(2*pi),'b','linewidth',1.2);
	else
	plot(x,phi_e/(2*pi),'b','linewidth',meta.linewidth);
	end
	%$plot(x,phi_e/(2*pi),colour{sdx},'linewidth',meta.linewidth);
	hline(0);
	hold on;
	
	% plot four stripes, the start and end has to be found as it is random
	for jdx=1:4
		x0 = -(1-pw);
		fdx1  = find(phi(:,sdx)/(2*pi)>0+(jdx-1),1,'first');
		fdx2 = find(phi(:,sdx)/(2*pi)>pw+(jdx-1),1,'first');
		if (~isempty(fdx2))
			patch([x(fdx1:fdx2);x(fdx2);x(fdx1)],[phi_e(fdx1:fdx2)/(2*pi);1;1],[0.5,0.5,0.5],'facealpha',0.5,'edgecolor','none');
		end
	end
	ylabel('Phase deviation $\frac{\varphi_e}{2 \pi}$');
	xlabel('Position $x/\lambda_c$');
	 xlim(xlim_);
	 ylim([-1,1]*0.35);
	if (pflag)
		daspect([3,1,1])
		pdfprint(13+sdx,sprintf('img/phase-error-schematic-s-%g.pdf',s(sdx)),plotscale);
	end
	end % for sdx
	
	splitfigure([2,3],[fid,3],fflag);
	y = phi/(2*pi);
	%[x,(x+2*W)];
	for pdx=1:2
	if (2==pdx)
	yyaxis left
	end % if 2 == pdx
	cla();
	plot(x,y(:,1:pdx),'color',meta.colour.blue,'linewidth',meta.linewidth)
	set(gca,'ycolor',meta.colour.blue)
	ylabel('Phase $\varphi/(2 \pi)$')
	xlabel('Position $x/\lambda_c$');
	if (2 == pdx)
		yyaxis right
		plot(x,e,'k.')
		ylim(limits(x))
		xlim(limits(x))
		%grid on
		set(gca,'ycolor','k')
		ylabel('Phase noise $\varphi_e$');
	end % of if 2 == pdx
	
	if (pflag)
		aspect = 3/4;
		pdfprint(13,sprintf('img/pattern-schematic-phase-%d.pdf',pdx),plotscale,aspect);
	end % of if pflag
	
	end % for pdx
	
	if (pflag)
		aspect = 0.83;
		figure(10*fid+1); 
		pdfprint(10*fid+1,'img/pattern-schematic-biomass-periodic.pdf',plotscale,aspect);
		figure(10*fid+1);
		set(gca,'xtick',[],'xlabel',[]);
		pdfprint(10*fid+1,'img/pattern-schematic-biomass-periodic-without-x.pdf',plotscale,aspect);
		figure(10*fid+2);
		pdfprint(10*fid+2,'img/pattern-schematic-biomass-perturbed.pdf',plotscale,aspect);
		figure(10*fid+2);
		set(gca,'xtick',[],'xlabel',[]);
		pdfprint(10*fid+2,'img/pattern-schematic-biomass-perturbed-without-x.pdf',plotscale,aspect);
	end % if pflag
end % sw_noisy_oscillator_schematic

