
# usage: runall.sh ../benchs | awk -f evalresults.awk | sort

# for each directory, returns:
# directory mean_initial_upperbound 
#           mean_optimum mean_initlb mean_backtracks mean_nodes mean_time(in seconds) 
#           number_of_problems_completely_solved 
#           number_of_problems


!/^[#]/ {
    dir = $1;
    n = split(dir,path,"/");
    name = path[1];
    for (i=2; i<=n; i++) {
	d = name " ";
	dirs[d]++;
	ub[d] += $2;
	opt[d] += $3;
	if ($3 != "-" && $3 >= 0) { nbopt[d]++; }
	initlb[d] += $4;
	bt[d] += $5;
	nodes[d] += $6;
	maxdepth[d] += $7;
	totaldepth[d] += $8;
	time[d] += $9;
	name = name "/" path[i];
    }
}

END {
    for (d in dirs) {
	printf("------------------------------------------------------------\n");
	n=dirs[d];
	if (d in nbopt) {
	    printf("%s optimum: %.2f  initlb: %.2f  backtracks: %.2f  nodes: %.2f  maxdepth: %.2f  totaldepth: %.2f  time: %.2f  optimal: %d of %d\n", d, opt[d]/n, initlb[d]/n, bt[d]/n, nodes[d]/n, maxdepth[d]/n, totaldepth[d]/n, time[d]/n, nbopt[d], n);
	} else {
	    printf("%s optimum: %.2f  initlb: %.2f  backtracks: %.2f  nodes: %.2f  maxdepth: %.2f  totaldepth: %.2f  time: %.2f  optimal: %d of %d\n", d, opt[d]/n, initlb[d]/n, bt[d]/n, nodes[d]/n, maxdepth[d]/n, totaldepth[d]/n, time[d]/n, 0, n);
	}
    }
}
