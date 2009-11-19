
# usage: runall.sh ../benchs | awk -f evalresults.awk | sort

# for each directory, returns:
# directory mean_initial_upperbound 
#           mean_optimum mean_nodes mean_time(in seconds) 
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
	if ($2 == $3) { nbopt[d]++; }
	nodes[d] += $4;
	time[d] += $5;
	name = name "/" path[i];
    }
}

END {
    for (d in dirs) {
	printf("------------------------------------------------------------\n");
	n=dirs[d];
	if (d in nbopt) {
	    printf("%s time: %.2f  nodes:  %.2f  optimal: %d of %d\n", d, time[d]/n, nodes[d]/n, nbopt[d], n);
	} else {
	    printf("%s time: %.2f  nodes:  %.2f  optimal: %d of %d\n", d, time[d]/n, nodes[d]/n, 0, n);
	}
    }
}
