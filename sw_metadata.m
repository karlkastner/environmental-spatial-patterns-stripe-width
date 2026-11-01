% Tue 17 Mar 17:01:41 CET 2026
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
function meta = sw_metadata()
	meta = struct();
	% use a single cpu-core per model run
	% run several matlab instances in parallel for utilizing all cores
	% parallelization within a single matlab instance is considerably less efficient
	meta.maxNumCompThreads = 1;
	% number of grid points for FPE discretization
    	% number of points for solving the FPE when fitting the NO2
	% this setting is very conservative
	meta.n_fpe = 1600*[1,1];
	meta.linewidth = 1;
	meta.pdf_method = 'ksdensity';
	meta.bandwidth = 2;
	% consistent project-wide colours
	meta.colour            = struct();
	meta.colour.blue       = [0.0,0.0,0.7];
	meta.colour.red        = [0.7,0.0,0.0];
	meta.colour.lightgrey  = 0.6*[1,1,1];
	meta.colour.lightgreen = [0.5,0.8,0.5];
end

