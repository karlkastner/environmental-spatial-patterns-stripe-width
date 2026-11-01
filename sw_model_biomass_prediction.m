% Fri 30 Oct 22:03:13 CET 2026
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

function sw_model_biomass_prediction(fid,pflag)

if (nargin()<1||isempty(fid))
	fid = 1;
end

if (nargin()<2)
	pflag = 0;
end
fflag = pflag;
init_plot();
%get(groot, 'defaultFigureColormap');
set(groot, 'defaultAxesColorOrder', get(groot, 'factoryAxesColorOrder'));                                            

model = {'Rietkerk','Tlidi','Joint'};
legend_C = {'$p_c$','$\bar{w}/\lambda_c$','cv($w$)','sk($w$)','$\lambda_c$','$S_{xc}^+/\lambda_c$'}

R2_individual = 0;
R2_pw = [];
% repeat with linear and logged variables
for fundx=1:2

% regression matrices combining data for both models
yy = [];
AA = [];

if (1==fundx)
	fun = @(x) x;
else
	fun = @log;
end

nn = [];
% for each model
for mdx=1:length(model)

switch (model{mdx})
case {'Rietkerk'}
	load('mat/rk-vx-10-ex-50-Lc-1-Lx-16384-T-73000/summary.mat');
	RR = repmat(rvec(o.R),length(o.cva),1);
	nn(1) = numel(RR);
	w = ones(nn(1),1);
case {'Tlidi'}
	load('mat/tlidi-L-16384-T-73000-dx-1-sl-1/summary.mat');
	RR = repmat(rvec(o.water_availability),length(o.cva),1);
	nn(2) = numel(RR);
	w = ones(nn(2),1);
case {'Joint'}
	% weight each run inverse proportionally to the number with runs by each model,
	% so that both models receive equal weight
	w = [1./nn(1)*ones(nn(1),1);
	     1./nn(2)*ones(nn(2),1)];
otherwise
	error('here');
end % switch


% normalize
w = w/sum(w);
% effective sample size
% ndof = sum(w)^2/sum(w^2)
ndof = 1./sum(w.*w); 

W = diag(sparse(w));

if (mdx<3)
	y = o.b_mean(:);
	y = y/mean(y);
	y = fun(y);

	% regression matrix
	A = [ones(numel(RR),1), ...
	fun(o.cover(:)), ...
	fun(o.w(1).mean(:).*o.fc.bar(:)), ...
	fun(o.w(1).sd(:)./o.w(1).mean(:)), ...
	o.w(1).sk(:), ...
	-fun(o.fc.bar(:)), ...
	fun(o.Sc.bar(:).*o.fc.bar(:))];
	% concatenate joint regression matrix
	AA = [AA;A];
	yy = [yy;y];
else
	A = AA;
	y = yy;
end % mdx


% number of coefficients
np = size(A,2);
% number of runs
n  = length(y);

% variance between runs
sy2 = var(y);
sy2_(mdx,fundx) = sy2;

% regression with one predictor
for idx=1:np-1
	Ai  = A(:,[1,idx+1]);
	[ci,se2_individual(idx),R2_individual(idx+1,mdx,fundx)] = fit(Ai,y,W);
end % idx

% pairwise, excluding constant (column 1) and cover (column 2)
for idx=3:np-1
	for jdx=idx+1:np
		Ai = A(:,[1,idx,jdx]);
		[ci,se2,R2_pw(idx-1,jdx-1,mdx)] = fit(Ai,y,W);
	end % for jdx
end % for idx

end % for mdx


% for each model dataset, find the best pw predictor
pair_C = {};
for mdx=1:length(model)
	%[~,mdx] = max(R2_pw(:,:,mdx)
	[~,mid] =max(R2_pw(:,:,mdx),[],'all');
	[mid,mjd] = ind2sub([size(R2_pw,1),size(R2_pw,2)],mid);
	R2_pw_(:,mdx,fundx) = R2_pw(mid,mjd,:);
	pair_C{mdx,1} = [legend_C{mid},' and ',legend_C{mjd}];
end

splitfigure([2,3],[fid,1+3*(fundx-1)],fflag);
cla();
bar(R2_individual(2:end,:,fundx));
for idx=1:np-1
	for jdx=1:3
		text(idx+(jdx-2)/4,R2_individual(idx+1,jdx,fundx)+0.02 ...
			, sprintf('%0.3f',R2_individual(idx+1,jdx,fundx)) ...
			, 'rot',90);
	end % jdx
end % idx
ax=set(gca,'XTickLabel',legend_C);
set(gca,'TickLabelInterpreter','latex')
ylabel('Goodness of fit $R^2$');
ylim([0,1.19]);
legend(model)
title('Single predictor')

splitfigure([2,3],[fid,2+3*(fundx-1)],fflag);
cla();
bar(R2_pw_(:,:,fundx)');
ylabel('Goodness of fit $R^2$');
set(gca,'xticklabel',pair_C);
set(gca,'TickLabelInterpreter','latex')
title('Pairwise predictors');

% joint only
splitfigure([2,3],[fid,3+3*(fundx-1)],fflag);
cla();
pair_C = {};
for idx=3:np-1
	pair_C{idx-2} = [legend_C{2},' and ',legend_C{idx}];
end
%R2_ = [R2_individual(2,3,fundx),R2_pw(2,3:end,3)];
R2_ = [R2_pw(2,3:end,3)];
cmap=colororder();
bar(R2_,'facecolor',cmap(3,:))
for idx=1:length(R2_)
	text(idx,R2_(idx)+0.02,sprintf('%0.3f',R2_(idx)),'rot',90)
end % jdx
set(gca,'xticklabel',pair_C);
ylabel('Goodness of fit $R^2$');
set(gca,'TickLabelInterpreter','latex')
title('Pairwise predictors');
ylim([0,1.19]);
xlim([0.3,np-2.3])

end % fundx

if (pflag)
	pdfprint(10*fid+1,'img/sw-biomass-indicator-lin.pdf',3);
	pdfprint(10*fid+4,'img/sw-biomass-indicator-log.pdf',3);
	pdfprint(10*fid+3,'img/sw-biomass-indicator-pairwise-ling.pdf',3);
	pdfprint(10*fid+6,'img/sw-biomass-indicator-pairwise-log.pdf',3);
end % if pflag

end % function sw_model_biomass_prediction

function [c,se2,R2] = fit(A,y,W)
	n  = size(A,1);
	np = size(A,2);
	sy2 = var(y);
	c  = (A'*W*A) \ (A'*W)*y;
	res = y - A*c;
	se2 = (n/(n-np))*(res'*W*res);
	R2 = 1 - se2/sy2;

end

