import pytoulbar2 as pytb2

# print a solution as a sudoku grid
def print_grid(solution):

    print('-------------------------')

    var_index = 0

    for row in range(9):
        line = ''
        for col in range(9):
            if col % 3 == 0:
                line += '|'
            line += ' ' 
            line += str(solution[var_index]+1)
            if col % 3 == 2:
                line += ' '
            var_index += 1
        line += '|'
        print(line)
        if row % 3 == 2:
            print('-------------------------')
    

cfn = pytb2.CFN()

# variables
for row in range(9):
    for col in range(9):
        cfn.AddVariable('cell_'+str(row)+'_'+str(col), range(9))

# define known values
initial_values = [[5,3,0,0,7,0,0,0,0],
                [6,0,0,1,9,5,0,0,0],
                [0,9,8,0,0,0,0,6,0],
                [8,0,0,0,6,0,0,0,3],
                [4,0,0,8,0,3,0,0,1],
                [7,0,0,0,2,0,0,0,6],
                [0,6,0,0,0,0,2,8,0],
                [0,0,0,4,1,9,0,0,5],
                [0,0,0,0,8,0,0,7,9]]

var_index = 0
for row in range(9):
    for col in range(9):
        if initial_values[row][col] != 0:
            cfn.Assign(var_index, initial_values[row][col]-1)
        var_index += 1

# row constraints
for row_ind in range(9):
    cfn.AddAllDifferent([row_ind*9+col_ind for col_ind in range(9)])

# column constraints
for col_ind in range(9):
    cfn.AddAllDifferent([row_ind*9+col_ind for row_ind in range(9)])

# sub grids constraints
for sub_ind1 in range(3): # row offset
    for sub_ind2 in range(3): # column offset
        cfn.AddAllDifferent([(sub_ind1*3+row_ind)*9+ sub_ind2*3+col_ind for col_ind in range(3) for row_ind in range(3)])
                
result = cfn.Solve(showSolutions = 3, allSolutions=1)
print_grid(result[0])