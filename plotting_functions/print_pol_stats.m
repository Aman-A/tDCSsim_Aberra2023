function print_pol_stats(pol_vals, comp_func,comp_type,plot_region_name,varargin)
in.units = 'uV';
in.cutoff = []; % cutoff fraction to calculate focality below (fraction above this value)
in = sl.in.processVarargin(in,varargin);
if strcmp(in.units,'mV')
   scale_factor = 1;  
elseif strcmp(in.units,'uV')
   scale_factor = 1e3; % scale factor to convert to uV
end

fprintf('%s %s polarization stats in region: %s:\n',comp_func,comp_type,plot_region_name)
for i = 1:size(pol_vals,1)
    fprintf(' L%g of %g\n',i,size(pol_vals,1));
   for j = 1:size(pol_vals,2)    
       p = pol_vals{i,j}(:)*scale_factor;
       fprintf('  E%g: min %.3f, mean %.3f, median %.3f max %.3f %s\n',...
                            j,nanmin(p),nanmean(p),nanmedian(p),nanmax(p),in.units);
       if ~isempty(in.cutoff)
           per_above_cutoff = 100*nansum(abs(p) > max(abs(p))*in.cutoff)/sum(~isnan(p));
           fprintf('       %g %% above %g %% of max\n',per_above_cutoff,100*in.cutoff);
       end
   end
end