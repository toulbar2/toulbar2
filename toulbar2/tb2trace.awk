BEGIN {
  first = 1;
  curdepth = 0;
  depth = 0;
}

match($1,"[[][0-9]*[]]") {
  depthstring = $1;
  sub("[[]","",depthstring);
  depth = 0 + depthstring;
  if (first) {
    curdepth = depth;
    first = 0;
  } else if (depth < curdepth) {
    curdepth = depth;
  }
}

match($1,"[[][0-9]*[]]") && !first && (depth <= curdepth) {
  print $0;
  if (match($0,"Try ")) {
    curdepth = depth - 1;
  }
}

!match($1,"[[][0-9]*[]]") {
  print $0;
}
