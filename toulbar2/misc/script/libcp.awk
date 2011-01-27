
# Library of awk functions loaded at each evaluation of constraint formulae
# during the translation process from cp format to wcsp format 

# hard constraint operator
function hard(e) {
    if (e == 0) return(-1);
    else return(0);
}

# soft constraint operator
function soft(v, e) {
    if (e == 0) return(v);
    else return(0);
}

# objective function operator
function obj(e) {
    return(e);
}

# integer functions
function abs(x) {
    if (x < 0) return(-x); 
    else return(x);
}

function min(x,y) {
    if (x < y) return(x); 
    else return(y);
}

function max(x,y) {
    if (x > y) return(x); 
    else return(y);
}

# boolean connectors
function implies(c1,c2) {
    if (c1 != 0 && c2 == 0) return(0);
    else return(1);
}

function equiv(c1,c2) {
    if ((c1 == 0 && c2 == 0) || (c1 != 0 && c2 != 0)) return(1);
    else return(0);
}

# global constraint alldifferent
function alldiff(a,b,c,d,e,f,g,h,i,j) {
  delete tmpval;
  if (a == "") return 1;
  if (a in tmpval) return 0;
  else tmpval[a] = 1;
  if (b == "") return 1;
  if (b in tmpval) return 0;
  else tmpval[b] = 1;
  if (c == "") return 1;
  if (c in tmpval) return 0;
  else tmpval[c] = 1;
  if (d == "") return 1;
  if (d in tmpval) return 0;
  else tmpval[d] = 1;
  if (e == "") return 1;
  if (e in tmpval) return 0;
  else tmpval[e] = 1;
  if (f == "") return 1;
  if (f in tmpval) return 0;
  else tmpval[f] = 1;
  if (g == "") return 1;
  if (g in tmpval) return 0;
  else tmpval[g] = 1;
  if (h == "") return 1;
  if (h in tmpval) return 0;
  else tmpval[h] = 1;
  if (i == "") return 1;
  if (i in tmpval) return 0;
  else tmpval[i] = 1;
  if (j == "") return 1;
  if (j in tmpval) return 0;
  else tmpval[j] = 1;
  return 1;
}

# celar trick for dividing the number of variables by two
function celar(x) {
  if (x >= 652 || (x >= 254 && x <= 408) || x == 478) return(x - 238);
  else return(x + 238);
}

# pedigree: returns first lowest allele number
# i,j: awk local variables
function low(x, n, i, j) {
    if (!(x in pedlow)) {
	for (i=1; i<= n; i++) {
	    for (j=i; j <= n; j++) {
		pedlow[i "" j] = i;
		pedhigh[i "" j] = j;
	    }
	}
    }
    return(pedlow[x]);
}

# pedigree: returns first lowest allele number
# i,j: awk local variables
function high(x, n, i, j) {
    if (!(x in pedhigh)) {
	for (i=1; i<= n; i++) {
	    for (j=i; j <= n; j++) {
		pedlow[i "" j] = i;
		pedhigh[i "" j] = j;
	    }
	}
    }
    return(pedhigh[x]);
}
