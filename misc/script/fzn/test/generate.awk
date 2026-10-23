NF==1&&/[(]coef,vars,res[)]/{
name=$1;
sub("[(].*","",name);
print $0;
file = "test_" name ".fzn";
print "" > file;
if(match(name,"bool_")) {
    print "var bool: x :: output_var;" >> file;
    print "var bool: y :: output_var;" >> file;
} else {
    print "array [1..4] of int: coef = [1,-1,2,-2];" >> file;
    print "array [1..4] of var -2..2: vars:: output_array([1..4]);" >> file;
    sub(",res",",3",$1);
}
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}

NF==1&&/[(]x,y,z[)]/ {
name=$1;
sub("[(].*","",name);
print $0;
file = "test_" name ".fzn";
print "" > file;
if(match(name,"bool_")) {
    print "var bool: x :: output_var;" >> file;
    print "var bool: y :: output_var;" >> file;
} else {
    print "var -10..10: x :: output_var;" >> file;
    print "var -10..10: y :: output_var;" >> file;
}
if(match(name,"_reif") || match(name,"bool_")) {
    print "var bool: z :: output_var;" >> file;
} else {
    print "var -10..10: z :: output_var;" >> file;
}
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}

NF==1&&/[(]x,y[)]/ {
name=$1;
sub("[(].*","",name);
print $0;
file = "test_" name ".fzn";
print "" > file;
if(match(name,"bool_") || match(name,"bool2int")) {
    print "var bool: x :: output_var;" >> file;
    if (match(name,"bool2int")) {
        print "var -10..10: y :: output_var;" >> file;
    } else {   
        print "var bool: y :: output_var;" >> file;
    }
} else {
    print "var -10..10: x :: output_var;" >> file;
    print "var -10..10: y :: output_var;" >> file;
}
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}


