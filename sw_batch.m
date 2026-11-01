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
%
%% batch script to reproducing simulations and plots

	% create auxiliar directories
	mkdir('./lib/');
	mkdir('./mat/');
	mkdir('./img/');
	mkdir('./lib/auxiliar/');
	addpath('./lib/auxiliar');

	% activate this code block at the first run to fetch dependencies
	% if you downloaded the code directly from github
	% this is not necessary if you have the prepackaged archive from zenodo
if (0)
	% fetch the script for fetching library files
	url  = 'https://raw.githubusercontent.com/karlkastner/auxiliar/master/dependencies_fetch.m';
	dest = './lib/auxiliar/dependencies_fetch.m';
	urlwrite(url,dest);

	% fetch library files, not required for code packed with prefeteched
	% libraries
	repository_url = 'https://github.com/karlkastner/';
	% dependencies_determine('dependencies.csv','mat/profile.mat',{'sw_batch'})
	dependencies_fetch(repository_url,'dependencies.csv');
end
	% determine dependencies (prior to packing)
	% dependencies_determine(meta.filename.dependencies,meta.filename.profile,{'pattern_periodicity_test_batch'}) 

	% add library subfolders to path
	addpath_recursive('./lib');

	% set to true to write images to disk
	pflag = false;

	% figure #1 in the manuscript
	% fid is figure handle index in matlab, all scripts plot to different figures
	% note that later on matlab figure handles do not match the figure
	% number in the manuscript, as more figures open in matlab, as there
	% are figures displayed in the manuscript
	fid = 1;
	sw_natural_pattern(fid,pflag);
	
	% figure 2
	fid = fid+2;
	sw_noisy_oscillator_schematic(fid,pflag);
	
	% figure 3 and 5
	fid = fid+1;
	sw_noisy_oscillator_fpe(fid,pflag)
	
	% figure 4 in manuscript
	fid = fid+2;
	sw_noisy_oscillator_simulated_two_lobes(fid,pflag)

	% figure 6
	fid = fid+1;
	sw_noisy_oscillator_limit_cases(fid,pflag);
	
	% figure 7a,b
	fid = fid+3;
	sw_noisy_oscillator_simulated_with_correlation(fid,pflag)
	
	% figure 7c
	fid = fid+1;
	sw_experiment_correlated_heterogeneity(fid,pflag);

	% figure 8
	fid  = fid+1;
	sw_rk1d_homogeneous(fid,pflag);

	fid = fid+1;
	sw_tlidi_homogeneous(fid,pflag);

	% run simulations to generate patterns for figure 8,9,10,11,12
	fid = fid+3;
	sw_rk1d_heterogeneous_simulate(fid,[],[],pflag);
	fid = fid+1;
	sw_tlidi_heterogeneous_simulate(fid,[],[],pflag);

	% set to true to generate output
	pflag = false;
	% figures 10, 12
	fid = fid+2;
	sw_model_plot_summary('mat/rk-vx-10-ex-50-Lc-1-Lx-16384-T-73000/',fid,pflag);
	fid = fid+2;
	sw_model_plot_summary('mat/tlidi-L-16384-T-73000-dx-1-sl-1/',fid,pflag);

	% figure 11
	fid = fid+1;
	sw_model_biomass_prediction();

