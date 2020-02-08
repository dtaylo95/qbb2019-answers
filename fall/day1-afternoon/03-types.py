#!/usr/bin/env Python3

# Integers
an_int = 42

# Booleans
truthy = True
falsish = False

# Floats
a_float = 42.0

# Strings
a_string = "42"

# Dictionary
a_dict = {
	'a' : 1,
	'b' : 5
}

#Lists (By **convention** lists are homogenous and tuples inhomogenous)
a_list = [5,10,15,20]

# Tuples (These are immutable, while lists are mutable)
a_tuple = (1, True, "Seven")


another_list = list(a_list)
another_list[3] = 4000
print(a_list)

for value in an_int, a_float, a_string, truthy, a_list, a_tuple:
	print(value, type(value))