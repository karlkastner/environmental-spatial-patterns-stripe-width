% Apr 15:07:02 CEST 2026
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
function sw_noisy_oscillator_simulated_with_correlation(fid,pflag)
	if (nargin()<1)
		fid =1;
	end
	if (nargin()<2)
		pflag =0;
	end
	fflag = pflag;
	plotscale = 5;
	init_plot();
	meta = sw_metadata();
	
	% characteristic wavelength (determines drift rate)
	lc = 100;
	% stripe width
	wc = 0.3*lc;
	% standard deviation of exogenous heterogeneity
	s  = 2.5;
	% correlation length
	lr =  [0.01/2,0.01,0.1*2]*lc*2.5;

	% threshold biomass density
	b_thresh = 0.5;
	% hysteresis in x
	x_hyst = 0;
	% plot offset of pattern
	x0_ = 0;
	
	% for process simulation
	% spatial resolution
	dx = 0.1;
	% spatial extent
	Lx = 1e4*lc;
	% number of grid points
	nx = Lx/dx;
	
	% simulate the process
	f_h = [];
	f_ksd = [];
	% histogram (density) bin edgces
	xh = (0:80)';
	% centres
	xc = mid(xh);
	b = [];
	clear out;
	% simulate the process
	for idx=1:length(lr)
		rng(1);
		%[stripe{idx}, gap, out(idx)] = stripe_width_2_pdf_simulate(Lx/lc,nx,wc/lc,s/sqrt(lc),x_hyst/lc,lr(idx)/lc);
		[f, out(idx)] = stripe_width_2_pdf_simulate(xc,wc/lc,s/sqrt(lc),lr(idx)/lc,x_hyst/lc,Lx/lc,nx);
		f_h(:,idx)   = histcounts(lc*out(idx).stripe.length,xh,'normalization','pdf');
		f_ksd(:,idx) = ksdensity(lc*out(idx).stripe.length,xc ...
				,'support','positive' ...
				,'boundarycorrection','reflection'...
				,'bandwidth',meta.bandwidth ...
				);
		[b(:,idx),p] = pattern_simple(out(idx).phi,wc/lc);
	end
	fdx = (out(1).x < 4);
	bt = (b>b_thresh);
	
	% normalize
	f_h = f_h./(sum(f_h)*(xh(2)-xh(1)));
	
	splitfigure([2,2],[fid,1],fflag);
	cla();
	a=area(out(1).x(fdx)-x0_/lc,[bt(fdx,idx)],'facecolor',meta.colour.lightgreen,'edgecolor','none');
	hold on
	plot(out(1).x(fdx)-x0_/lc,b(fdx,idx),'-','linewidth',meta.linewidth);
	set(gca,'colororder',[0,0,0.6; 0.9,0,0; 0,0.6,0]);
	hline(0.5);
	ylabel('Biomass density $b/(2 \, b_t)$')
	xlim([0,4])
	set(gca,'xtick',0:4);
	if (0)
	plot(xs-x0_/lc,[1.1,1.1],'-k','linewidth',meta.linewidth);
	ylim([0,1.2])
	else
	ylim([0,1.1])
	end
	xlabel('Position $x/\lambda_c$');
	axis square
	
	% biomass legend
	figure(fid+1)
	cla();
	plot(NaN,NaN,'-','color',meta.colour.blue,'linewidth',meta.linewidth)
	hold on
	plot(NaN,NaN,'color',meta.colour.lightgreen,'linewidth',10);
	axis off
	legend('Plain','Thresholded');
	if (pflag)
		pdfprint(fid+1,'img/biomass-legend.pdf',plotscale,[],'pdf');
	end
	
	% individual biomass plots
	for idx=1:3
		splitfigure([2,3],[fid+1,idx],fflag);
		cla();
		a=area(out(1).x(fdx)-x0_/lc,[bt(fdx,idx)],'facecolor',meta.colour.lightgreen,'edgecolor','none');
		hold on
		set(gca,'colororder',[0,0,0.6; 0.9,0,0; 0,0.6,0]);
		set(gca,'colororderindex',1);
		plot(out(1).x(fdx)-x0_/lc,b(fdx,idx),'-','linewidth',meta.linewidth);
		%hold on
		hline(0.5);
		ylabel('$\hspace{-1em}$Biomass density $b/(2 \, b_t)$');
		xlim([0,4]-x0_/lc);
		set(gca,'xtick',0:4);
		ylim([0,1.1])
		xlabel('Position $x/\lambda_c$')
		axis square
		if (pflag)
			pdfprint(10*(fid+1)+idx,sprintf('img/noisy-oscillator-correlated-biomass-se-%0.2g-lr-%0.2e.pdf',s,lr(idx)),plotscale,[],'pdf');
		end % if pflag
	end % for idx
	
	% SWD
	splitfigure([2,2],[fid,3],fflag);
	cla();
	plot(mid(xh)/lc,f_h,'linewidth',1.5)
	xlabel('Stripe width $w/\lambda_c$');
	ylabel('$\hspace{-1em}$Probability density $f \cdot \lambda_c$')
	set(gca,'colororder',[0,0,0.6; 0.9,0,0; 0,0.6,0]);
	ylim([0,0.2]);
	axis square
	
	% individual plots
	for idx=1:3
		splitfigure([2,3],[fid+1,3+idx],fflag);
		cla();
		set(gca,'colororderindex',1);
		switch (meta.pdf_method)
		case {'histogram'}
			a=histplot(xh/lc,lc*f_h(:,idx),'facecolor',meta.colour.lightgrey,'edgecolor','none')
		case {'ksdensity'}
			area(xc/lc,lc*f_ksd(:,idx),'facecolor',meta.colour.lightgrey,'edgecolor','none');
		otherwise
			error('here')
		end
		xlabel('Stripe width $w/\lambda_c$');
		ylabel('$\hspace{-1em}$Probability density $f \cdot \lambda_c$')
		set(gca,'colororder',[0,0,0.6; 0.9,0,0; 0,0.6,0]);
		xlim([0,0.6]);
		ylim([0,1.1*max(f_ksd(:,idx))*lc]);
		axis square
		if (pflag)
			pdfprint(10*(fid+1)+3+idx,sprintf('img/noisy-oscillator-correlated-stripe-width-distribution-se-%0.2g-lr-%0.2e.pdf',s,lr(idx)),plotscale,[],'pdf');
		end
	end
	
	for idx=1:length(lr)
		splitfigure([2,length(lr)],[11,idx],fflag);
		plot(out(idx).x(fdx),out(idx).phi(fdx)/(2*pi)-out(idx).x(fdx)/lc);
		hline(0);
	end % for idx
	
	if (pflag)
		figure(fid*10+1);
		xlabel('');
		ylabel('');
		xr=diff(xlim); yr=diff(ylim);
		d = 0.04;
		xlim([-d,d]);
		ylim(0.5+d*yr/xr*[-1,1]);
		pdfprint(fid*10+1,'img/noisy-oscillator-biomass-correlated-inset.pdf',plotscale,[],'pdf');
		pdfprint(fid*10+3,'img/noisy-oscillator-correlated-stripe-width-distribution.pdf',plotscale);
	end % if pflag
end % sw_noisy_oscillator_simulated_with_correlation

