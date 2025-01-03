import pytoulbar2 as tb2

model = tb2.CFN()

for i in range(1,7):
    model.AddVariable(f"x{i}", [1,2,3,4,5,6])

model.AddAllDifferent(['x1', 'x2', 'x3', 'x4', 'x5', 'x6'])
model.Solve()

