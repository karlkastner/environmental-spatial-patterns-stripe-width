% Sat 17 Oct 10:34:58 CEST 2026
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
function out = sw_noisy_oscillator_special_cases(fid,pflag)

if (nargin()<1||isempty(fid))
	fid = 1;
end

if (nargin<2)
	pflag = 0;
end
fflag = pflag;
set(groot,'DefaultTextInterpreter','latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

cmap = [0,0,0;
        0.8,0.2,0; 
        0,0.2,0.8;
        0.5,0,0.5];

lc   = 100;
% spatial extent
L    = lc;
% spatial resolution
dx   = 1;
% number of points along SWD
nx   = 200;
% points where SWD is plotted
xi   = logspace(log10(dx),log10(4*L),nx)';
% number of grid points to solve FPE
n    = [1e3,1e3];
m    = 2;
% characteristic stripe width
wc     = [0.2,0.2, 4,       4]*lc;
% degree of random exogenous heterogeneity
s      = [0.2,  0.2,   0.2, 1]*sqrt(lc);
% hystereis when threholding biomass
x_hyst = [0.1, 0.01,   0.01,         0.01]*lc;

f=[];
F=[];
out = {};
for idx=1:length(wc)
 	[f(:,idx),F(:,idx),out{idx}]=stripe_width_2_pdf_fpe(xi/lc,wc(idx)/lc,s(idx)/sqrt(lc),x_hyst(idx)/lc,m*L/lc,n);
 	%[f(:,idx),F(:,idx),out{idx}]=stripe_width_2_pdf_fpe(xi,wc(idx),s(idx),x_hyst(idx),m*L,n);
	%leg_C{idx} = sprintf('$\\frac{w_c}{L}$ = %0.1f, $\\frac{s_e}{\\sqrt{L}}$ = %0.2g',wc(idx)/L,s(idx)/sqrt(L));
	leg_C{idx} = sprintf('$\\displaystyle\\frac{w_c}{L}$ = %0.1f, $\\displaystyle\\frac{s_e}{\\sqrt{L}}$ = %0.1f $\\displaystyle\\frac{\\varphi_h}{\\varphi_c}$ = %0.2f',wc(idx)/L,s(idx)/sqrt(L),x_hyst(idx)/wc(idx));
end

% normalize
f = f./(sum(mid(f).*diff(xi)));
F = inner2outer(cumsum(mid(f).*diff(xi)));

% truncate to spatial extent to L
fdx = (xi<=L);
xi  = xi(fdx);
f   = f(fdx,:);
F   = F(fdx,:);

splitfigure([2,2],[fid,1],fflag);
loglog(xi/lc,lc*f,'linewidth',1);
xlim([1/lc,1]);
ylim(lc*[0.001,1]);
axis square
xlabel('Stripe width $w/\lambda_c$');
ylabel('Probability density $f \cdot \lambda_c$');
set(gca,'colororder',cmap);

splitfigure([2,2],[fid,2],fflag);
loglog(xi/lc,1-F,'linewidth',1);
xlim([1/lc,1]);
ylim(10.^[-2,0.1])
axis square
xlabel('Stripe width $w/\lambda_c$');
ylabel('Survival function $1 - F$');
set(gca,'colororder',cmap);

splitfigure([2,2],[fid,3],fflag);
plot(NaN(length(wc)),'linewidth',1.5);
set(gca,'colororder',cmap);
axis off
legend(leg_C)


if (pflag)
	pdfprint(10*fid+1,'img/swd-noisy-oscillator-special-cases.pdf',3);
	pdfprint(10*fid+2,'img/swd-noisy-oscillator-special-cases-survival.pdf',3);
	pdfprint(10*fid+3,'img/swd-noisy-oscillator-special-legend.pdf',3);
end

end

