# problem name
5warehouses_10stores_opencost30

# warehouse open state
warehouse0 0 1
warehouse1 0 1
warehouse2 0 1
warehouse3 0 1
warehouse4 0 1

# warehouse opening costs
soft(30, warehouse0 == 0)
soft(30, warehouse1 == 0)
soft(30, warehouse2 == 0)
soft(30, warehouse3 == 0)
soft(30, warehouse4 == 0)

# store allocation to warehouses
store0 0 1 2 3 4
store1 0 1 2 3 4
store2 0 1 2 3 4
store3 0 1 2 3 4
store4 0 1 2 3 4
store5 0 1 2 3 4
store6 0 1 2 3 4
store7 0 1 2 3 4
store8 0 1 2 3 4
store9 0 1 2 3 4

# channeling constraints between stores and warehouses
hard( implies(store0 == 0, warehouse0 == 1) )
hard( implies(store1 == 0, warehouse0 == 1) )
hard( implies(store2 == 0, warehouse0 == 1) )
hard( implies(store3 == 0, warehouse0 == 1) )
hard( implies(store4 == 0, warehouse0 == 1) )
hard( implies(store5 == 0, warehouse0 == 1) )
hard( implies(store6 == 0, warehouse0 == 1) )
hard( implies(store7 == 0, warehouse0 == 1) )
hard( implies(store8 == 0, warehouse0 == 1) )
hard( implies(store9 == 0, warehouse0 == 1) )
hard( implies(store0 == 1, warehouse1 == 1) )
hard( implies(store1 == 1, warehouse1 == 1) )
hard( implies(store2 == 1, warehouse1 == 1) )
hard( implies(store3 == 1, warehouse1 == 1) )
hard( implies(store4 == 1, warehouse1 == 1) )
hard( implies(store5 == 1, warehouse1 == 1) )
hard( implies(store6 == 1, warehouse1 == 1) )
hard( implies(store7 == 1, warehouse1 == 1) )
hard( implies(store8 == 1, warehouse1 == 1) )
hard( implies(store9 == 1, warehouse1 == 1) )
hard( implies(store0 == 2, warehouse2 == 1) )
hard( implies(store1 == 2, warehouse2 == 1) )
hard( implies(store2 == 2, warehouse2 == 1) )
hard( implies(store3 == 2, warehouse2 == 1) )
hard( implies(store4 == 2, warehouse2 == 1) )
hard( implies(store5 == 2, warehouse2 == 1) )
hard( implies(store6 == 2, warehouse2 == 1) )
hard( implies(store7 == 2, warehouse2 == 1) )
hard( implies(store8 == 2, warehouse2 == 1) )
hard( implies(store9 == 2, warehouse2 == 1) )
hard( implies(store0 == 3, warehouse3 == 1) )
hard( implies(store1 == 3, warehouse3 == 1) )
hard( implies(store2 == 3, warehouse3 == 1) )
hard( implies(store3 == 3, warehouse3 == 1) )
hard( implies(store4 == 3, warehouse3 == 1) )
hard( implies(store5 == 3, warehouse3 == 1) )
hard( implies(store6 == 3, warehouse3 == 1) )
hard( implies(store7 == 3, warehouse3 == 1) )
hard( implies(store8 == 3, warehouse3 == 1) )
hard( implies(store9 == 3, warehouse3 == 1) )
hard( implies(store0 == 4, warehouse4 == 1) )
hard( implies(store1 == 4, warehouse4 == 1) )
hard( implies(store2 == 4, warehouse4 == 1) )
hard( implies(store3 == 4, warehouse4 == 1) )
hard( implies(store4 == 4, warehouse4 == 1) )
hard( implies(store5 == 4, warehouse4 == 1) )
hard( implies(store6 == 4, warehouse4 == 1) )
hard( implies(store7 == 4, warehouse4 == 1) )
hard( implies(store8 == 4, warehouse4 == 1) )
hard( implies(store9 == 4, warehouse4 == 1) )

# transportation cost from warehouses to stores
store0 0
0 20
1 24
2 11
3 25
4 30
store1 0
0 28
1 27
2 82
3 83
4 74
store2 0
0 74
1 97
2 71
3 96
4 70
store3 0
0 2
1 55
2 73
3 69
4 61
store4 0
0 46
1 96
2 59
3 83
4 4
store5 0
0 42
1 22
2 29
3 67
4 59
store6 0
0 1
1 5
2 73
3 59
4 56
store7 0
0 10
1 73
2 13
3 43
4 96
store8 0
0 93
1 35
2 63
3 85
4 46
store9 0
0 47
1 65
2 55
3 71
4 95
