
import sys
import CFN

f = open(sys.argv[1], 'r').readlines()

tokens = []
for l in f:
    tokens += l.split()

pos = 0

def token():
   global pos, tokens
   if (pos == len(tokens)):
      return None
   s = tokens[pos]
   pos += 1
   return int(float(s)) 

N = token()
token() # skip freeze time

LT = []
PC = []
ST = []

for i in range(N):
   token()  # skip appearance time
# Times per plane: {earliest landing time, target landing time, latest landing time}
   LT.append([token(), token(), token()])

# Penalty cost per unit of time per plane:
# [for landing before target, after target]
   PC.append([token(), token()])

# Separation time required after i lands before j can land
   ST.append([token() for j in range(N)])

top = 99999

Problem = CFN.CFN(top)
for i in range(N):
    Problem.AddVariable('x' + str(i), range(LT[i][0],LT[i][2]+1))

for i in range(N):
    Problem.AddFunction([i], [PC[i][0]*abs(a-LT[i][1]) for a in range(LT[i][0], LT[i][2]+1)])

for i in range(N):
    for j in range(i+1,N):
        Problem.AddFunction([i, j], [top*(a+ST[i][j]>b and b+ST[j][i]>a) for a in range(LT[i][0], LT[i][2]+1) for b in range(LT[j][0], LT[j][2]+1)])

Problem.Dump('airplane.cfn')
Problem.NoPreprocessing()
Problem.Solve()

