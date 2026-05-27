import pytoulbar2 

pb = pytoulbar2.CFN(1)

a = pb.AddVariable('a' , range(1,20))
b = pb.AddVariable('b' , range(1,20))
c = pb.AddVariable('c' , range(1,20))
d = pb.AddVariable('d' , range(1,20))
e = pb.AddVariable('e' , range(1,20))
f = pb.AddVariable('f' , range(1,20))
g = pb.AddVariable('g' , range(1,20))
h = pb.AddVariable('h' , range(1,20))
i = pb.AddVariable('i' , range(1,20))
j = pb.AddVariable('j' , range(1,20))
k = pb.AddVariable('k' , range(1,20))
l = pb.AddVariable('l' , range(1,20))
m = pb.AddVariable('m' , range(1,20))
n = pb.AddVariable('n' , range(1,20))
o = pb.AddVariable('o' , range(1,20))
p = pb.AddVariable('p' , range(1,20))
q = pb.AddVariable('q' , range(1,20))
r = pb.AddVariable('r' , range(1,20))
s = pb.AddVariable('s' , range(1,20))

LD = [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s]

pb.AddAllDifferent( LD , encoding = 'binary' )

pb.AddSumConstraint([  a , b , c ], '==',   38 )
pb.AddSumConstraint([ d , e , f , g ], '==',   38 )
pb.AddSumConstraint([ h , i , j , k , l ], '==',   38 ) 
pb.AddSumConstraint([ m , n , o , p ], '==',   38 ) 
pb.AddSumConstraint([ q , r , s ], '==',   38 ) 
pb.AddSumConstraint([ a , d , h ], '==',   38 ) 
pb.AddSumConstraint([ b , e , i , m ], '==',   38 ) 
pb.AddSumConstraint([ c , f , j , n , q ], '==',   38 ) 
pb.AddSumConstraint([ g , k , o , r ], '==',   38 ) 
pb.AddSumConstraint([ l , p , s ], '==',   38 ) 
pb.AddSumConstraint([ c , g , l ], '==',   38 ) 
pb.AddSumConstraint([ b , f , k , p ], '==',   38 ) 
pb.AddSumConstraint([ a , e , j , o , s ], '==',   38 ) 
pb.AddSumConstraint([ d , i , n , r ], '==',   38 ) 
pb.AddSumConstraint([ h , m , q ], '==',   38 ) 

pb.AddLinearConstraint([1, -1], [ a , c], '<', 0)
pb.AddLinearConstraint([1, -1], [ a , h], '<', 0)
pb.AddLinearConstraint([1, -1], [ a , l], '<', 0)
pb.AddLinearConstraint([1, -1], [ a , q], '<', 0)
pb.AddLinearConstraint([1, -1], [ a , s], '<', 0)
pb.AddLinearConstraint([1, -1], [ c , h], '<', 0)

pb.Dump('magichexagon3.cfn')
print(pb.Solve( showSolutions = 3))

