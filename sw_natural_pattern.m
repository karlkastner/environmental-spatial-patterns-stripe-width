% Wed 29 Apr 18:19:40 CEST 2026
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
%% process the natural example pattern
%%
function sw_natural_pattern(fid,pflag)
	if (nargin()<1)
		fid = 1;
	end
	if (nargin()<2)
		pflag = 0;
	end
	fflag = pflag;
	init_plot();
	plotscale = 3.2;
	meta = sw_metadata();
	
	% spatial filter window size
	nf_spatial = [2,4,8];
	% spectral filter window size
	nf_spectral = 2.5;
	
	% load satellite image
	g=GeoImg();
	g.read('input/pattern-sw-utm.png');
	% coordinate axes
	x=g.x;
	y=g.y;
	dx = diff(x(1:2));
	% zero axes
	x=x-x(1);
	y=y-y(end);
	% spatial extent
	L = min(max(x),max(y));
	x=x-max(x)+max(y);
	% spatial extent
	Lx = x(end);
	% spectral axis
	fx=fourier_axis(x);
	% spectral resolution
	df = 1/Lx;
	% biomass proxy
	b=mean(g.img,3);

	% periodogram 
	Shat2d = abs(fft2(b)).^2;
	n=size(Shat2d);
	% zero mean
	Shat2d(1)=0;
	if (1)
		Sbar2d=gaussfilt2(Shat2d,nf_spectral);
	end
	Sbarxp = mean(Sbar2d,1);
	% noramlize of the positive half-axis
	Sbarxp = 2*Sbarxp/sum(Sbarxp*df);
	[Sxpc,mdx] = max(Sbarxp);
	% characteristic frequency
	fc=abs(fx(mdx));
	% characteristic wavelength
	lc = 1./fc;
	
	% plot pattern
	splitfigure([2,3],[fid,1],fflag);
	imagesc(x,y,g.img);
	axis xy
	axis square
	axis tight
	axis([0,L,0,L])
	xlabel('Position $x$/m');
	ylabel('Position $y$/m');
	
	% plot pattern adjust contrast for display
	splitfigure([2,3],[fid,2],fflag);
	o=0.05;
	is=1+o;
	img_ = single(g.img)/255;
	lim = limits(flat(mean(img_,3)));
	img_ = max(0,min(1,is*(img_-o)/(lim(2)-o)));
	imagesc(x,y,img_);
	axis xy
	axis square
	axis tight
	axis([0,L,0,L])
	xlabel('Position $x$/m');
	ylabel('Position $y$/m');

	% plot spectral density
	splitfigure([2,3],[fid,3],fflag);
	cla
	fdx = (fx>=0 & fx<=3.5*fc);
	plot(fx(fdx)/fc,Sbarxp(fdx)*fc,'color',meta.colour.blue,'linewidth',meta.linewidth);
	xlim([0,3.5])
	xlabel('Wavenumber $k/k_c$');
	ylabel('Spectral density $S_x^+/\lambda_c$');
	text(1.75,1*Sxpc*fc,sprintf('$\\lambda_c$ = %g m',round(1./fc)));
	ylim([0,1.1*Sxpc*fc]);
	axis square

	if (pflag)
		pdfprint(fid*10+1,'img/stripe-width-natural-pattern-satellite-image.pdf',plotscale);
		pdfprint(fid*10+2,'img/stripe-width-natural-pattern-satellite-image-scaled.pdf',plotscale);
		pdfprint(fid*10+3,'img/stripe-width-natural-pattern-spectral-density.pdf',plotscale);
	end

	ncol = length(nf_spatial);
	for idx=1:ncol

	% extract stripes
	[s,gw,id,bt] = stripe_width(-b',L,nf_spatial(idx));
	bt = bt';

	% plot thresholded pattern
	splitfigure([2,ncol],[fid+1,idx],fflag);
	bt = double(bt);
	dw = (-6:6);
	b0 = bt(end/2+dw,:);
	bt(end/2+dw,:) = 2*b0 + 3*(~b0);
	imagesc(x,y,bt);
	axis xy
	axis square
	axis tight
	axis([0,L,0,L])
	colormap([1,1,1;
                  0,0.3,0;
                  1,0,0;
                  0.5,0.5,0.5])
	xlabel('Position $x$/m');
	ylabel('Position $y$/m');

	% plot stripe width distribution
	splitfigure([2,ncol],[fid+1,ncol+idx],fflag);
	cla();
	switch (meta.pdf_method)
	case {'histogram'}
		dx=(x(2)-x(1));
		id = 0:2:300;
		xh = dx*id;
		dxh = (xh(2)-xh(1));
		f_h = histcounts(s,xh);
		% normalize
		f_h = f_h/(dxh*sum(f_h));
		histplot(xh,f_h,'edgecolor','none','facecolor',meta.colour.lightgrey)
		hold on
	case {'ksdensity'}
		xc = mid([0.1,0:300]);
		dxc = xc(2)-xc(1);
		f_ksd  = ksdensity(s,xc,'support','positive', ...
					'boundarycorrection','reflection', ...
					'bandwidth',meta.bandwidth);
		area(xc/lc,f_ksd*lc,'facecolor',meta.colour.lightgrey,'edgecolor','none')
	end
	area(xc*fc,f_ksd/fc,'facecolor',meta.colour.lightgrey,'edgecolor','none')
	xlabel('Stripe width $w/\lambda_c$');
	ylabel('Probability density $f \cdot \lambda_c$');
	xlim([0,fc*50]);
	% TODO nf
	wmu = mean(s);
	wsd = std(s);
	wsk = skewness(s);
	%text(0.2,6,sprintf('$\\bar{w}$ = %g',wmu),'interpreter','latex');
	x0 = [0.225,0.005*[1,1]];
	text(x0(idx),5.15,sprintf('$L_s$ \\hspace{1.35em} = %0.2g m\n$\\bar{w}$ \\hspace{1.5em} = %0.3g m\n{cv($w$) = %0.2g}\n{sk($w$) = %0.2g}',dx*nf_spatial(idx),wmu,wsd/wmu,wsk),'interpreter','latex');

	%ylim([0,1.1*max(ksd)/fc])
	ylim([0,6.25])
	axis square

	if (pflag)
		base = sprintf('img/stripe-width-natural-pattern-nf-%g',nf_spatial(idx));
		pdfprint((fid+1)*10+idx,[base,'-thresholded.pdf'],plotscale);
		pdfprint((fid+1)*10+idx+ncol,[base,'-distribution.pdf'],plotscale);
	end % if pflag
	
	end % for idx

	figure(30)
	clf();
	nf=[0,3,10];
	for idx=1:length(nf);
		s = stripe_width(-b',L,nf(idx));
		x = 2:0.1:100;
		h = ksdensity(s,x,'Bandwidth',1);
		loglog(fc*x,h/fc,'linewidth',1);
		hold on;
	end;
	leg = arrayfun(@(x) sprintf('%3.1f',x),(dx*nf),'uniformoutput',false);
	lh = legend(leg);
	title(lh,'$L_s$');
	xlim(fc*[2,1e2]);
	ylim([1e-1,1e2]);
	%ylim([1e-3,1]/fc);
	xlabel('Stripe width $w/\lambda_c$');
	ylabel('Probability density $\lambda_c \cdot f$'); 
	set(gca,'colororder',[0,0,0; meta.colour.red; meta.colour.blue])
	axis square
	
	if (pflag)
		pdfprint(30,['img/stripe-width-natural-pattern-loglog.pdf'],plotscale);
	end

end % sw_natural_pattern

