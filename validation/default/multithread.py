import glob
import random
from threading import Thread, Lock
import pytoulbar2 as pytb2

def solve_instance(instance, id):
    print('start solving', instance)
    pytb2.tb2.init()
    cfn = pytb2.CFN(verbose = -1)
    cfn.Read(instance)
    res = cfn.Solve()
    lock.acquire()
    results[id] = res
    lock.release()
    print('done solving', instance)

random.seed(14)

path = '../validation/default/'
instances = ['max.cfn', 'magic3.wcsp', 'knapsack.wcsp', 'cap131.wcsp', 'golomb4-salldiff-reverse.wcsp', 'weightedcspconstraint.wcsp', 'samong.cfn', 'samongdp.cfn', 'smaxdp.cfn', '10_1.wcsp']

lock = Lock()
results = [None for _ in range(len(instances))]

# multi threads version
threads = []
for id, inst in enumerate(instances): # Start each thread
    t = Thread(target=solve_instance, args=(path+inst, id))
    t.start()
    threads.append(t)
for t in threads: # Wait for all threads to finish
    t.join()
#-------------------

# single thread version
# for id, inst in enumerate(instances):
#     solve_instance(path+inst, id)
#-------------------

print('\nsolutions:')
for id, sol in enumerate(results):
    print(instances[id], ':', sol)
