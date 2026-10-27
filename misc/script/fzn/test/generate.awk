NF==1&&/[(]coef,vars,res,z[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "array [1..4] of int: coef = [1,-1,2,-2];" >> file;
if(match(name,"bool_")) {
    print "array [1..4] of var bool: vars:: output_array([1..4]);" >> file;
} else {
    print "array [1..4] of var -2..2: vars:: output_array([1..4]);" >> file;
}
sub(",res",",3",$1);
print "var bool: z :: output_var;" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}

NF==1&&/[(]coef,vars,res[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "array [1..4] of int: coef = [1,-1,2,-2];" >> file;
if(match(name,"bool_")) {
    print "array [1..4] of var bool: vars:: output_array([1..4]);" >> file;
} else {
    print "array [1..4] of var -2..2: vars:: output_array([1..4]);" >> file;
}
sub(",res",",3",$1);
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}

NF==1&&/array_/&&/_element/&&/[(]x,y,z[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
if(match(name,"_bool_")) {
    if(match(name,"_var_")) {
        print "array [1..4] of var bool: y:: output_array([1..4]);" >> file;
    } else {
        print "array [1..4] of bool: y = [true,false,true,false];" >> file; # warning! ALL parameters must be printed before variables
    }
    print "var bool: z :: output_var;" >> file;
} else {
    if(match(name,"_var_")) {
        print "array [1..4] of var -2..2: y:: output_array([1..4]);" >> file;
    } else {
        print "array [1..4] of int: y = [1,-1,2,-2];" >> file;
    }
    print "var -10..10: z :: output_var;" >> file;
}
print "var 1..4: x :: output_var;" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
next;
}

NF==1&&/[(]x,y,z[)]/{
name=$1;
sub("[(].*","",name);
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

NF==1&&/array_bool_and/&&/[(]x,y[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "array [1..4] of var bool: x:: output_array([1..4]);" >> file;
print "var bool: y :: output_var;" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
next;
}

NF==1&&/bool_clause/&&/[(]x,y[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "array [1..4] of var bool: x:: output_array([1..4]);" >> file;
print "array [1..4] of var bool: y:: output_array([1..4]);" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
next;
}

NF==1&&/array_bool_or/&&/[(]x,y[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "array [1..4] of var bool: x:: output_array([1..4]);" >> file;
print "var bool: y :: output_var;" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
next;
}

NF==1&&/[(]x,y[)]/{
name=$1;
sub("[(].*","",name);
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

NF==1&&/[(]x,dom,z[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "set of int: dom = {-1,1,3,5,7};" >> file;
print "var -10..10: x :: output_var;" >> file;
print "var bool: z :: output_var;" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}

NF==1&&/[(]x,dom[)]/{
name=$1;
sub("[(].*","",name);
file = "test_" name ".fzn";
print "" > file;
print "set of int: dom = {-1,1,3,5,7};" >> file;
print "var -10..10: x :: output_var;" >> file;
print "constraint",$1,";" >> file;
print "solve satisfy;" >> file;
close(file);
}


