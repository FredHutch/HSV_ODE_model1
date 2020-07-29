This directory contains some perl and shell scripts that can be used to guide fitting
if you opt to fit outside of the program (allows for fitting related parameters together)
Like the built-in fitting, these scripts use high/low variable values to set the 
ranges for each "latin-hypercube"-style search.

The scripts are detailed below (though they haven't been run in a number of years...)

    run_pairs.sh - used to run rounds of fitting with pairs of variables
		   (results of each step fed into the next simulation via the 
		    input file)
    top_scores.sh - find the best score in a file of results
    top_summary.pl - prints out the best score found

    Files for finding best parameter values and adjusting high/low range for next search
	top_betae.pl
	top_beta.pl
	top_c.pl
	top_delta.pl
	top_eclipse.pl
	top_inf.pl
	top_p.pl
	top_rho.pl
	top_rinf.pl
	top_r.pl
	top_theta.pl
	top_vburstrate.pl
