#!/bin/tcsh

echo "*** extract flatzinc built-in functions ***"
grep -E "def .*[(].*[)]:" ../fzn2py.py | sed 's/def //' | sed 's/://' | sed 's/ //g' >! builtin.txt
echo ""

echo "*** generate flatzinc built-in tests ***"
awk -f generate.awk builtin.txt
echo ""

echo "*** run flatzinc built-in tests ***"
foreach f (test_*.fzn)
    echo "** $f **"
    ../fzn2py $f -a | python | awk '\!/^%/&&\!/---/{printf $0} /----/{print ""} END{print ""}' | sed 's/ //g' | sed 's/=UNKNOWN=/=UNSATISFIABLE=/' | sort >! ${f:r}_tb2
    fzn-gecode -f -a  -s $f  | awk '\!/^%/&&\!/---/{printf $0} /----/{print ""} END{print ""}' | sed 's/ //g' | sed 's/false/0/g' | sed 's/true/1/g' | sort >! ${f:r}_gecode
    diff -q ${f:r}_tb2 ${f:r}_gecode
end

# special int_pow builtin
cat int_pow_chuffed | awk '\!/^%/&&\!/---/{printf $0} /----/{print ""} END{print ""}' | sed 's/ //g' | sed 's/false/0/g' | sed 's/true/1/g' | sort >! test_int_pow_chuffed
diff -q test_int_pow_tb2 test_int_pow_chuffed

