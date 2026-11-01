% 2026-03-26 16:47:37.249433457 +0100
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
%% simulate a striped pattern generate by a noisy oscillator
%% plot the pattern and stripe width distribution (SWD) after tresholding with
%% hysteresis
%% note that the fpe-function can only determine the asymptotic SWD, but not
%% a pattern with individual stripes
function sw_noisy_oscillator_simulated_two_lobes(fid,pflag)
	if (nargin()<1)
		fid = 1;
	end
	if (nargin()<2)
		pflag = 0;
	end
	fflag = pflag;
	plotscale = 3;
	init_plot();
	meta = sw_metadata();
	
	% characteristic wavelength (determines drift rate)
	lc = 100;
	% stripe width
	wc = 0.5*lc;
	% standard deviation of exogenous heterogeneity
	s = 2.75;
	% correlation length
	lr=0;
	
	% threshold biomass density
	b_thresh    = 0.5;
	% hysteresis in x
	x_hyst = [0.05,0.1,0.2]*wc*2.2;
	
	% hystersis in b
	% b_hyst = x_hyst*db/dx
	% db/dx ~ 0.5 * 2 pi/lc (for sine at transition)
	b_hyst = pi/lc*x_hyst;
	
	% for process simulation
	% spatial resolution
	dx = 0.1;
	% spatial extent
	Lx = 1e4*lc;
	% number of grid points
	nx = Lx/dx;
	
	% plot offset (skip first period since first band is skipped by the thres alg)
	x0_ = 0.75*lc;
	% offset for pattern_simple and swd_2_pdf_simulte to aligin
	% TODO this should be determined by pattern_simple depending on wc/lc
	phi_0 = 0.5*pi;
	
	% histogram (density) bin edgces
	xh = (0:80)';
	% centres
	xc = mid(xh);
	
	% declare output variables
	f_h  = [];
	f_ksd = [];
	b  = [];
	bt = [];
	clear out;
	% simulate the process
	for idx=1:length(x_hyst)
		rng(1);
		%[stripe{idx}, gap, out(idx)] = stripe_width_2_pdf_simulate(Lx/lc,nx,wc/lc,s/sqrt(lc),x_hyst(idx)/lc);
		[f, out(idx)] = stripe_width_2_pdf_simulate(xc,wc/lc,s/sqrt(lc),lr/lc,x_hyst(idx)/lc,Lx/lc,nx);
		f_h(:,idx)     = histcounts(lc*out(idx).stripe.length,xh,'normalization','pdf');
		f_ksd(:,idx) = ksdensity(lc*out(idx).stripe.length,xc ...
				,'support','positive' ...
				,'boundarycorrection','reflection' ...
				,'bandwidth',meta.bandwidth ...
				);
		b(:,idx)  = pattern_simple(out(idx).phi+phi_0,wc/lc);
		bt(:,idx) = out(idx).stripe.in;
	end % for idx
	fdx = (out(1).x < 5);
	%[bt,is,ie] = threshold_with_hysteresis(b,b_thresh,b_hyst(idx));
	
	% normalize
	f_h = f_h./(sum(f_h)*(xh(2)-xh(1)));
	
	for idx=1:3
		splitfigure([2,3],[fid,idx],fflag);
		cla();
		a=area(out(1).x(fdx)-x0_/lc,[out(idx).stripe.in(fdx)],'facecolor',[0.5,0.8,0.5],'edgecolor','none');
		hold on
		plot(out(1).x(fdx)-x0_/lc,b(fdx,idx),'color',meta.colour.blue,'linewidth',meta.linewidth);
		set(gca,'xtick',0:4);
		xlabel('Position $x/\lambda_c$','interpreter','latex')
		ylabel('Biomass $b/(2 \, b_t)$','interpreter','latex')
		hline(b_thresh-b_hyst(idx),'color',meta.colour.red,'linewidth',meta.linewidth)
		hline(b_thresh,'color','k','linewidth',meta.linewidth)
		hline(b_thresh+b_hyst(idx),'color',[0,0.6,0],'linewidth',meta.linewidth)
		xlim([0,4]);
		ylim([0,1.25]);
		axis square
		text(+0.01*Lx,0.55,'Entry');
		text(+0.01*Lx,0.45,'Exit');
		drawnow();
	if (pflag)
		pdfprint(10*fid+idx,sprintf('img/noisy-oscillator-biomass-two-sided-entry-cover-%g-x-hyst-%0.2e.pdf',wc/lc,x_hyst(idx)),plotscale);
	end
	end % for idx
	
	splitfigure([2,2],[fid,3],fflag);
	cla();
	area(out(1).x(fdx),2*pi*bt(fdx,idx),'facecolor',meta.colour.lightgrey,'edgecolor','none');
	hold on
	plot(out(1).x(fdx),(2*pi*out(1).x(fdx)),':','color',meta.colour.blue,'linewidth',meta.linewidth);
	sc=2;
	hline(2*pi*0.25+sc*b_hyst(idx),'color',[0,0.5,0],'linewidth',1);
	hline(2*pi*0.25-sc*b_hyst(idx),'color',[0.9,0,0],'linewidth',1);
	hline(2*pi*0.75-sc*b_hyst(idx),'color',[0,0.5,0],'linewidth',1);
	hline(2*pi*0.75+sc*b_hyst(idx),'color',[0.9,0,0],'linewidth',1);
	ylim(2*pi*[0,1]);
	xlim([0,1]);
	axis square
	xlabel('Position $x/\lambda_c$','interpreter','latex');
	ylabel('Phase $\varphi$','interpreter','latex');
	drawnow
	if (pflag)
		pdfprint(10*fid+2,'img/noisy-oscillator-phase-two-sided-entry.pdf',plotscale);
	end
	
	splitfigure([2,2],[fid,4],fflag);
	cla
	plot(mid(xh)/lc,f_h,'linewidth',meta.linewidth)
	xlabel('Stripe width $w/\lambda_c$');
	ylabel('Probability density $f \cdot \lambda_c$')
	set(gca,'colororder',[meta.colour.blue; meta.colour.red; 0,0.6,0]);
	ylim([0,0.2]);
	axis square
	
	if (pflag)
		pdfprint(13,'img/noisy-oscillator-stripe-width-two-sided-distribution.pdf',plotscale);
	end
	
	
	% individual plots
	for idx=1:3
	splitfigure([2,3],[fid+1,3+idx],fflag);
	cla();
	set(gca,'colororderindex',1);
	switch (meta.pdf_method)
	case {'histogram'}
		a=histplot(xh/lc,lc*f_h(:,idx),'facecolor',[1,1,1]*2/3,'edgecolor','none');
	case {'ksdensity'}
		area(xc/lc,lc*f_ksd(:,idx),'facecolor',0.5*[1,1,1],'edgecolor','none');
	otherwise
		error('here');
	end
	hold on
	xlabel('Stripe width $w/\lambda_c$');
	ylabel('Probability density $f \cdot \lambda_c$')
	set(gca,'colororder',[0,0,0.6; 0.9,0,0; 0,0.6,0]);
	xlim([0.0, 0.6]);
	ylim([0, 1.1*lc*max(f_ksd(:,idx))]);
	axis square
	if (pflag)
		pdfprint(10*(fid+1)+3+idx,sprintf('img/noisy-oscillator-two-sided-stripe-width-distribution-cover-%g-x-hyst-%0.2e.pdf',wc/lc,x_hyst(idx)),plotscale,[],'pdf');
	end
	meta = sw_metadata();
	Lfit = 2*lc;
	hh = lc*stripe_width_2_pdf_fpe(xh/lc,wc/lc,s/sqrt(lc),x_hyst(idx)/lc,Lfit/lc,meta.n_fpe);
	hh = hh/sum(hh);
	plot(xh/lc,hh*lc,'r','linewidth',meta.linewidth)
	if (pflag)
		pdfprint(10*(fid+1)+3+idx,sprintf('img/noisy-oscillator-two-sided-stripe-width-distribution-cover-%g-x-hyst-%0.2e-with-fit.pdf',wc/lc,x_hyst(idx)),plotscale,[],'pdf');
	end
	end % for idx
end % sw_noisy_oscillator_simulated_two_lobes

