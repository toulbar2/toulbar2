# Discrete Lot Sizing Problem (CSPLib 058)
# dataset from MiniZinc Challenge 20219 and 2020 

import sys
import pytoulbar2


import sys
import sys
import os
import re


def parse_dzn(path):
	"""
	Reads a .dzn file (MiniZinc data format) and returns a
	dictionary {variable_name: value}.
	Handles scalars (integers/floats) and arrays (lists of
	numbers between square brackets [ ... ]), possibly spanning
	multiple lines.
	"""
	with open(path, 'r') as f:
		text = f.read()

	# Strip "% ..." comments
	text = re.sub(r'%.*', '', text)

	data = {}
	# Look for assignments of the form:  name = value ;
	# value can be:
	#   - a plain array [ ... ] (possibly spanning multiple lines)
	#   - a 2D array written as array2d(Set1, Set2, [ ... ])
	#   - a numeric scalar
	pattern = r'(\w+)\s*=\s*(array2d\s*\([^,]+,[^,]+,\s*\[[^\]]*\]\s*\)|\[[^\]]*\]|[-+]?\d+\.?\d*)\s*;'
	for match in re.finditer(pattern, text, re.DOTALL):
		name, value = match.group(1), match.group(2).strip()

		if value.startswith('array2d'):
			# Extract only the inner [ ... ] part of array2d(Set1, Set2, [ ... ])
			inner_match = re.search(r'\[[^\]]*\]', value, re.DOTALL)
			inner = inner_match.group(0)[1:-1]
			nums = re.findall(r'[-+]?\d+\.?\d*', inner)
			data[name] = [float(n) if '.' in n else int(n) for n in nums]
		elif value.startswith('['):
			inner = value[1:-1]
			nums = re.findall(r'[-+]?\d+\.?\d*', inner)
			data[name] = [float(n) if '.' in n else int(n) for n in nums]
		else:
			data[name] = float(value) if '.' in value else int(value)

	return data


def load_from_dzn(path):
	data = parse_dzn(path)

	required = ['nb_item_types', 'nb_orders', 'nb_periods', 'inventory_cost',
				'due_period', 'change_cost']
	missing = [r for r in required if r not in data]
	if missing:
		print('Error: missing variables in .dzn file: %s' % ', '.join(missing))
		sys.exit(1)

	return (data['nb_item_types'], data['nb_orders'], data['nb_periods'],
			data['inventory_cost'], data['due_period'], data['change_cost'])


def load_named_instance(instance_name):
	if instance_name == 'pigment15a':
		nb_item_types = 5
		nb_orders = 14
		nb_periods = 15
		inventory_cost = 10

		due_period = [8, 14, 5, 12, 15, 7, 12, 14, 8, 11, 15, 9, 12, 15, ]
		change_cost = [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 105, 105, 105, 154, 154, 154, 130, 130, 130, 100, 100, 100,
			0, 0, 0, 105, 105, 105, 154, 154, 154, 130, 130, 130, 100, 100, 100,
			0, 146, 146, 0, 0, 0, 135, 135, 135, 139, 139, 139, 167, 167, 167,
			0, 146, 146, 0, 0, 0, 135, 135, 135, 139, 139, 139, 167, 167, 167,
			0, 146, 146, 0, 0, 0, 135, 135, 135, 139, 139, 139, 167, 167, 167,
			0, 101, 101, 183, 183, 183, 0, 0, 0, 193, 193, 193, 113, 113, 113,
			0, 101, 101, 183, 183, 183, 0, 0, 0, 193, 193, 193, 113, 113, 113,
			0, 101, 101, 183, 183, 183, 0, 0, 0, 193, 193, 193, 113, 113, 113,
			0, 188, 188, 112, 112, 112, 111, 111, 111, 0, 0, 0, 103, 103, 103,
			0, 188, 188, 112, 112, 112, 111, 111, 111, 0, 0, 0, 103, 103, 103,
			0, 188, 188, 112, 112, 112, 111, 111, 111, 0, 0, 0, 103, 103, 103,
			0, 179, 179, 117, 117, 117, 161, 161, 161, 124, 124, 124, 0, 0, 0,
			0, 179, 179, 117, 117, 117, 161, 161, 161, 124, 124, 124, 0, 0, 0,
			0, 179, 179, 117, 117, 117, 161, 161, 161, 124, 124, 124, 0, 0, 0
		]

	elif instance_name == 'pigment15b':
		nb_item_types = 5
		nb_orders = 13
		nb_periods = 15
		inventory_cost = 10

		due_period = [6, 9, 11, 15, 10, 13, 15, 7, 14, 11, 14, 9, 15, ]
		change_cost = [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 128, 128, 128, 117, 117, 197, 197, 178, 178,
			0, 0, 0, 0, 0, 128, 128, 128, 117, 117, 197, 197, 178, 178,
			0, 0, 0, 0, 0, 128, 128, 128, 117, 117, 197, 197, 178, 178,
			0, 0, 0, 0, 0, 128, 128, 128, 117, 117, 197, 197, 178, 178,
			0, 125, 125, 125, 125, 0, 0, 0, 102, 102, 146, 146, 175, 175,
			0, 125, 125, 125, 125, 0, 0, 0, 102, 102, 146, 146, 175, 175,
			0, 125, 125, 125, 125, 0, 0, 0, 102, 102, 146, 146, 175, 175,
			0, 179, 179, 179, 179, 191, 191, 191, 0, 0, 190, 190, 161, 161,
			0, 179, 179, 179, 179, 191, 191, 191, 0, 0, 190, 190, 161, 161,
			0, 156, 156, 156, 156, 139, 139, 139, 119, 119, 0, 0, 189, 189,
			0, 156, 156, 156, 156, 139, 139, 139, 119, 119, 0, 0, 189, 189,
			0, 109, 109, 109, 109, 121, 121, 121, 111, 111, 187, 187, 0, 0,
			0, 109, 109, 109, 109, 121, 121, 121, 111, 111, 187, 187, 0, 0,
		]

	elif instance_name == 'pigment15c':
		nb_item_types = 8
		nb_orders = 13
		nb_periods = 15
		inventory_cost = 10

		due_period = [15, 10, 11, 14, 9, 14, 9, 15, 11, 13, 14, 10, 14, ]
		change_cost = [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 149, 134, 134, 110, 110, 137, 137, 191, 191, 104, 192, 192,
			0, 131, 0, 139, 139, 192, 192, 118, 118, 117, 117, 161, 132, 132,
			0, 128, 115, 0, 0, 195, 195, 196, 196, 174, 174, 127, 144, 144,
			0, 128, 115, 0, 0, 195, 195, 196, 196, 174, 174, 127, 144, 144,
			0, 134, 154, 119, 119, 0, 0, 165, 165, 118, 118, 163, 109, 109,
			0, 134, 154, 119, 119, 0, 0, 165, 165, 118, 118, 163, 109, 109,
			0, 169, 118, 117, 117, 140, 140, 0, 0, 166, 166, 174, 116, 116,
			0, 169, 118, 117, 117, 140, 140, 0, 0, 166, 166, 174, 116, 116,
			0, 134, 153, 139, 139, 158, 158, 159, 159, 0, 0, 100, 170, 170,
			0, 134, 153, 139, 139, 158, 158, 159, 159, 0, 0, 100, 170, 170,
			0, 170, 199, 185, 185, 193, 193, 143, 143, 146, 146, 0, 186, 186,
			0, 162, 101, 130, 130, 115, 115, 193, 193, 190, 190, 150, 0, 0,
			0, 162, 101, 130, 130, 115, 115, 193, 193, 190, 190, 150, 0, 0,
		]

	elif instance_name == 'pigment20b':
		nb_item_types = 10
		nb_orders = 18
		nb_periods = 20
		inventory_cost = 10

		due_period = [4, 17, 5, 19, 5, 18, 6, 17, 19, 8, 16, 6, 18, 9, 4, 20, 10, 19, ]
		change_cost = [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 114, 114, 160, 160, 127, 127, 127, 172, 172, 176, 176, 132, 199, 199, 197, 183,
			0, 0, 0, 114, 114, 160, 160, 127, 127, 127, 172, 172, 176, 176, 132, 199, 199, 197, 183,
			0, 153, 153, 0, 0, 123, 123, 114, 114, 114, 138, 138, 198, 198, 147, 170, 170, 195, 171,
			0, 153, 153, 0, 0, 123, 123, 114, 114, 114, 138, 138, 198, 198, 147, 170, 170, 195, 171,
			0, 115, 115, 183, 183, 0, 0, 147, 147, 147, 189, 189, 111, 111, 135, 198, 198, 185, 138,
			0, 115, 115, 183, 183, 0, 0, 147, 147, 147, 189, 189, 111, 111, 135, 198, 198, 185, 138,
			0, 106, 106, 159, 159, 185, 185, 0, 0, 0, 134, 134, 164, 164, 122, 166, 166, 151, 104,
			0, 106, 106, 159, 159, 185, 185, 0, 0, 0, 134, 134, 164, 164, 122, 166, 166, 151, 104,
			0, 106, 106, 159, 159, 185, 185, 0, 0, 0, 134, 134, 164, 164, 122, 166, 166, 151, 104,
			0, 135, 135, 122, 122, 166, 166, 110, 110, 110, 0, 0, 100, 100, 118, 130, 130, 166, 199,
			0, 135, 135, 122, 122, 166, 166, 110, 110, 110, 0, 0, 100, 100, 118, 130, 130, 166, 199,
			0, 116, 116, 116, 116, 114, 114, 150, 150, 150, 138, 138, 0, 0, 188, 114, 114, 124, 126,
			0, 116, 116, 116, 116, 114, 114, 150, 150, 150, 138, 138, 0, 0, 188, 114, 114, 124, 126,
			0, 152, 152, 124, 124, 125, 125, 119, 119, 119, 136, 136, 118, 118, 0, 169, 169, 199, 123,
			0, 101, 101, 161, 161, 171, 171, 177, 177, 177, 111, 111, 147, 147, 180, 0, 0, 110, 101,
			0, 101, 101, 161, 161, 171, 171, 177, 177, 177, 111, 111, 147, 147, 180, 0, 0, 110, 101,
			0, 142, 142, 186, 186, 156, 156, 167, 167, 167, 191, 191, 122, 122, 137, 180, 180, 0, 123,
			0, 182, 182, 181, 181, 132, 132, 181, 181, 181, 138, 138, 195, 195, 165, 110, 110, 123, 0,
		]

	elif instance_name == 'pigment30a':
		nb_item_types = 5
		nb_orders = 12
		nb_periods = 30
		inventory_cost = 10

		due_period = [5, 20, 28, 16, 28, 14, 16, 27, 2, 10, 21, 30, ]
		change_cost = [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 184, 184, 194, 194, 194, 130, 112, 112, 112,
			0, 0, 0, 0, 184, 184, 194, 194, 194, 130, 112, 112, 112,
			0, 0, 0, 0, 184, 184, 194, 194, 194, 130, 112, 112, 112,
			0, 114, 114, 114, 0, 0, 172, 172, 172, 100, 185, 185, 185,
			0, 114, 114, 114, 0, 0, 172, 172, 172, 100, 185, 185, 185,
			0, 143, 143, 143, 175, 175, 0, 0, 0, 125, 173, 173, 173,
			0, 143, 143, 143, 175, 175, 0, 0, 0, 125, 173, 173, 173,
			0, 143, 143, 143, 175, 175, 0, 0, 0, 125, 173, 173, 173,
			0, 104, 104, 104, 110, 110, 191, 191, 191, 0, 156, 156, 156,
			0, 127, 127, 127, 195, 195, 172, 172, 172, 139, 0, 0, 0,
			0, 127, 127, 127, 195, 195, 172, 172, 172, 139, 0, 0, 0,
			0, 127, 127, 127, 195, 195, 172, 172, 172, 139, 0, 0, 0,
		]

	elif instance_name == 'pigment30c':
		nb_item_types = 10
		nb_orders = 16
		nb_periods = 30
		inventory_cost = 10

		due_period = [8, 29, 28, 15, 10, 30, 13, 9, 10, 11, 13, 14, 10, 7, 15, 14, ]
		change_cost = [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 161, 186, 172, 172, 194, 193, 193, 193, 193, 193, 126, 111, 153, 149,
			0, 0, 0, 161, 186, 172, 172, 194, 193, 193, 193, 193, 193, 126, 111, 153, 149,
			0, 120, 120, 0, 140, 110, 110, 108, 120, 120, 120, 120, 120, 121, 137, 154, 198,
			0, 131, 131, 109, 0, 142, 142, 144, 116, 116, 116, 116, 116, 129, 183, 154, 135,
			0, 148, 148, 109, 141, 0, 0, 121, 192, 192, 192, 192, 192, 167, 154, 185, 104,
			0, 148, 148, 109, 141, 0, 0, 121, 192, 192, 192, 192, 192, 167, 154, 185, 104,
			0, 155, 155, 196, 103, 186, 186, 0, 105, 105, 105, 105, 105, 189, 164, 128, 124,
			0, 173, 173, 118, 195, 145, 145, 176, 0, 0, 0, 0, 0, 136, 170, 193, 150,
			0, 173, 173, 118, 195, 145, 145, 176, 0, 0, 0, 0, 0, 136, 170, 193, 150,
			0, 173, 173, 118, 195, 145, 145, 176, 0, 0, 0, 0, 0, 136, 170, 193, 150,
			0, 173, 173, 118, 195, 145, 145, 176, 0, 0, 0, 0, 0, 136, 170, 193, 150,
			0, 173, 173, 118, 195, 145, 145, 176, 0, 0, 0, 0, 0, 136, 170, 193, 150,
			0, 186, 186, 117, 155, 154, 154, 159, 157, 157, 157, 157, 157, 0, 197, 104, 102,
			0, 149, 149, 181, 156, 143, 143, 148, 168, 168, 168, 168, 168, 107, 0, 119, 174,
			0, 118, 118, 129, 188, 125, 125, 175, 142, 142, 142, 142, 142, 176, 100, 0, 118,
			0, 171, 171, 184, 167, 157, 157, 136, 126, 126, 126, 126, 126, 194, 158, 167, 0,
		]

	else:
		print('Sorry, instance unknown!')
		sys.exit(1)

	return (nb_item_types, nb_orders, nb_periods, inventory_cost, due_period, change_cost)


def load_instance(arg):
	"""
	Single entry point:
	- if arg ends with '.dzn' or matches an existing file on disk,
	  the data is read from that .dzn file;
	- otherwise, arg is treated as the name of a hardcoded instance
	  (pigment15a, pigment15b, pigment15c, pigment20b, pigment30a, pigment30c).
	"""
	if arg.lower().endswith('.dzn') or os.path.isfile(arg):
		return load_from_dzn(arg)
	else:
		return load_named_instance(arg)


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: python pigment_instances.py <instance_name | path_to_file.dzn>')
		sys.exit(1)

	(nb_item_types, nb_orders, nb_periods, inventory_cost,
	 due_period, change_cost) = load_instance(sys.argv[1])
	
	Problem = pytoulbar2.CFN(verbose = 0)

	#Model with nb_periods variables having domains in [0,nb_orders]
	for p in range(1, nb_periods + 1):
		x = Problem.AddVariable('p' + str(p), range(0, nb_orders + 1))
		# Add inventory costs as unary cost functions
		Problem.AddFunction([x], [inventory_cost * (due_period[i-1] - p) if i != 0 and p <= due_period[i-1] else (Problem.Top if i != 0 else 0) for i in range(0, nb_orders + 1)])

	#Semi-assignment problem modelled with a closed GCC with unit capacity, except for value 0 (meaning no item production)
	Problem.AddGlobalCardinalityConstraint(['p' + str(p) for p in range(1, nb_periods + 1)], [(i, 1, 1) if i != 0 else (0, nb_periods - nb_orders, nb_periods - nb_orders) for i in range(0, nb_orders + 1)])

	#Add change costs between successive periods
	for p in range(1, nb_periods):
		Problem.AddFunction(['p' + str(p), 'p' + str(p+1)], [change_cost[i * (nb_orders + 1) + j]  for i in range(0, nb_orders + 1)  for j in range(0, nb_orders + 1)])

	#Build a weighted automaton for encoding change costs between non successive periods separated by some 0s
	#TODO: reduce the number of states by using item types instead of orders
	#TODO: or merge this automatom in previous binary cost functions by changing domains of variables (adding specific 0 values for each item type) (warning, GCC will be non squared with unit demands)
	parameters = [nb_orders * 2 + 1, 1, (0, 0), nb_orders * 2 + 1] + [(q, 0) for q in range(nb_orders * 2 + 1)]
	transitions = [(0, 0, 0, 0)]
	for i in range(1,nb_orders + 1):
		transitions.append((0,i,i,0))
		transitions.append((i,0,nb_orders + i,0))
		transitions.append((nb_orders + i,0,nb_orders + i,0))
		for j in range(1,nb_orders + 1):
			if i == j:
				transitions.append((i,i,i,0))
				transitions.append((nb_orders + i,i,i,0))
			else:
				transitions.append((i,j,j,0))
				transitions.append((nb_orders + i,j,j,change_cost[i * (nb_orders + 1) + j]))
	parameters += [len(transitions)] + transitions
	Problem.AddGlobalFunction(['p' + str(p) for p in range(1, nb_periods + 1)], 'wregular', parameters)

	#Problem.Dump(f'lotsizing_{sys.argv[1]}.wcsp')
	Problem.Option.singletonConsistency = nb_periods * 20
	Problem.Option.GilmoreLawler = 0
	res = Problem.Solve(showSolutions = 3, timeLimit = 300)
	if res:
		print('Sequence:', res[0])
		print('NbBacktracks:', Problem.GetNbBacktracks())
		print('NbNodes:', Problem.GetNbNodes())
		print('Objective:', int(res[1]))
	else:
		print('No solution!')
    
