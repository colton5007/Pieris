from copy import copy

def copy_matrix(A):
	rows = len(A)
	cols = len(A[0])
	AC = [[A[i][j] for j in range(rows)] for i in range(cols)]
	return AC

# Computes the determinant (recursively) using minor reduction of a matrix of Schubert cycles
def gdeterminant(A,total=None):
	if not total:
		total = []
	indices = list(range(len(A)))
	# Recursive base case
	if len(A) == 2 and len(A[0]) == 2:
		a00 = copy(A[0][0])
		a01 = copy(A[0][1])
		a10 = copy(A[1][0])
		a11 = copy(A[1][1])
		a10.scalar *= -1
		return [[a00, a11], [a10,a01]]
	# Computes minors
	for fc in indices:
		As = copy_matrix(A)
		As = As[1:]
		height = len(As)

		for i in range(height):
			As[i] = As[i][0:fc] + As[i][fc+1:]

		sign = (-1) ** (fc % 2)

		sub_det = gdeterminant(As)
		kfc = copy(A[0][fc])
		kfc.scalar *= sign
		
		this_det = [[kfc] + entry for entry in sub_det]
		total += this_det
	return total

# Class for handling sums of Schubert cycles in the cohomology ring
class CSum:
	def __init__(self):
		self.summands = []
		
	# Adds a cycle to the sum
	def add_cycle(self, other):
		existing = False
		for i in range(len(self.summands)):
			summand = self.summands[i]
			# Check if the currect summand in the list iteration is the same as the to be added summand
			if summand == other:
				cycle = SchubertCycle(a=summand.a, k=summand.k, n=summand.n, scalar=summand.scalar + other.scalar)
				# Deletes summand if scalars add to 0
				if cycle.scalar == 0:
					del self.summands[i]
					existing = True
					break
				# Modifies summand in the list with updated scalar
				else:
					self.summands[i] = cycle
					existing = True
					break
		# Adds new cycle if not present in the list of summands
		if not existing and other.scalar != 0:
			self.summands.append(copy(other))
			
	# Adds a cycle or another sum of cycles to the sum
	def add(self, other):
		if type(other) is SchubertCycle:
			self.add_cycle(other)
		elif type(other) is CSum:
			for cycle in other.summands:
				self.add_cycle(cycle)

	# Scales the cycles in the sum by an integer scalar
	def scale(self, alpha):
		cs = copy(self)
		for i in range(len(cs.summands)):
			cs.summands[i].scalar *= alpha
		return cs
				
	# Python method overload for adding a+b
	def __add__(self, other):
		self.add(other)
		return self

	# Python method overload for subtracting a-b
	def __sub__(self, other):
		o = other.scale(-1)
		return self+o

	# Python method overload for negation -a
	def __neg__(self):
		return self.scale(-1)
	
	# Produces nice string representation for sum class
	def __repr__(self):
		return str(self)
		
	# For stringifying Schubert cycles
	def __str__(self):
		if len(self.summands) == 0:
			return "0"
		else:
			str_summands = [str(summand) for summand in self.summands]
			return "+".join(str_summands)

	# Computes the length of the stringified version of the class. To be used in computing padding for prettified output
	def str_length(self):
		str_summands = [str(summand) for summand in self.summands]
		o = "+".join(str_summands)
		return len(o)

	# Python method overload for multiplication a*b, uses multiplication in cycle class
	def __mul__(self, other):
		output = CSum()
		if type(other) is SchubertCycle:
			for summand in self.summands:
				output.add(summand*other)
		elif type(other) is CSum:
			for summand1 in self.summands:
				for summand2 in other.summands:
					output.add(summand1*summand2)
		return output

	# Python method overload for left multiplication by scalars 
	def __rmul__(self, alpha):
		return self.scale(alpha)

class SchubertCycle:
	def __init__(self, a=[],k=0,n=0,scalar=1):
		# Initializes cycle and checks if initial paramaters define a valid cycle
		self.n = n
		self.k = k
		self.zpos = len(a)
		for i in range(len(a)):
			if a[i] == 0:
				self.zpos = i
				break
		
		ordered = True
		if len(a) == 0:
			self.a = [0 for i in range(k)]
			self.scalar = scalar
		elif a[0] <= n+1-k:
			for i in range(len(a)-1):
				if a[i] < a[i+1] or a[i+1] < 0:
					ordered = False
					break
			if ordered:
				self.a = a + [0 for i in range(k-len(a))]
				self.scalar = scalar
			else:
				print(f"Not ordered {str(a)}")
				self.a = [0 for i in range(k)]
				self.scalar = 0
		else:
			self.a = [0 for i in range(k)]
			self.scalar = 0

	# Produces stringified version of a cycle for textual output
	def __str__(self):
		scalar_str = ""
		if self.grading() == 0:
			scalar_str = str(self.scalar)
			return scalar_str
		
		if self.scalar == -1:
			scalar_str = '-'
		elif self.scalar != 1:
			scalar_str = str(self.scalar)
			
		index_str = (str(self.a[:self.zpos])[1:-1]).replace(' ', '')
		
		return f"{scalar_str}σ_({index_str})" 
		
	# Checks if cycles belong to the same cohomology class, and not equality since scalars are disregarded
	def __eq__(self, other):
		if self.k == other.k and self.n == other.n:
			for i in range(self.k):
				if self.a[i] != other.a[i]:
					return False
		else:
			return False
		return True
		
	# Returns the grading of the cycle in the cohomology ring
	def grading(self):
		return 2*(sum(self.a))

	# Produces a scaled cycle using an integer scalar
	def scale(self, alpha):
		c = copy(self)
		c.scalar *= alpha
		return c
		
	# Uses the CSum class to compute sums of cycles
	def __add__(self, other):
		if type(other) is SchubertCycle:
			cs = CSum()
			cs.add_cycle(self)
			cs.add_cycle(other)
			return cs
		elif type(other) is CSum:
			other.add_cycle(self)
			return other
			
	# Uses the CSum class to compute difference of cycles
	def __sub__(self, other):
		c = other.scale(-1)
		return self+c

	# Bootstraps multiplication to perform exponention, not well-optimized
	def __pow__(self, k):
		if k == 0:
			return SchubertCycle(a=[], k=self.k,n=self.n)
		t = self
		for i in range(k-1):
			t = t*self
		return t

	# Negates the Schubert cycle, equivalent to scaling by -1
	def __neg__(self):
		return self.scale(-1)

	# Computes the matrix to be used in the Giambelli's formula
	def giambellis(self):
		A = [[SchubertCycle(a=[self.a[i] + j-i ],k=self.k, n=self.n) for j in range(self.k)] for i in range(self.k)]
		return A

	# Multiplies two cycles a*b using Pieri's formula
	def multiply(self, a, mu, scalar):
		cs = CSum()
		possibilities = []
		# Obtains all of the k-tuples that satisfy the range requirements given in Pieri's formula
		for i in range(len(a)):
			cis = []
			lb = a[i]
			ub = self.n
			if i > 0:
				ub = a[i-1] + 1
			for ci in range(lb, ub):
				cis.append([ci])
			if len(possibilities) == 0:
				possibilities = cis
			else:
				np = []
				for p in possibilities:
					for ci in cis:
						np.append(p + ci)
				possibilities = np

		# Filters the k-tuples from above using the second condition of Pieri's formula, and produces the output
		for p in possibilities:
			if sum(p) == sum(a) + mu:
				cs.add_cycle(SchubertCycle(a=p, k=self.k, n=self.n, scalar=scalar))
		return cs
		
	# Checks if the cycle is primitive, i.e., is a tuple of the form (a1,0,...,0)
	def is_primitive(self):
		return self.zpos == 1

	# Returns the string representation for the cycle for console output
	def __repr__(self):
		return str(self)
	
	# General case manager for cycle multiplication a*b
	def __mul__(self, other):
		# Distributes product across sum
		if type(other) is CSum:
			return other*self
		# Checks if grading of product is greater than the dimension of the cohomology ring
		if self.grading() + other.grading() > 2*(self.k * (self.n + 1 - self.k)):
			return CSum()
		# Checks if cycle is in the 0-th cohomology class
		elif self.grading() == 0:
			other.scalar *= self.scalar
			return copy(other)
		elif other.grading() == 0:
			self.scalar *= other.scalar
			return copy(self)
		# Checks if cycle is primitive, and if so, performs Pieri's formula
		elif self.is_primitive():
			return self.multiply(other.a, self.a[0], other.scalar * self.scalar)
		elif other.is_primitive():
			return self.multiply(self.a, other.a[0], other.scalar * self.scalar)
		# Application of giambell's formula in cases when Pieri's formula cannot be applied directly
		else:
			A = self.giambellis()
			g = gdeterminant(copy(A))
			cs = CSum()
			for entry in g:
				temp = copy(other)
				for cycle in entry:
					temp = temp*cycle
				cs.add(temp)
			return cs.scale(self.scalar)
	
	# Performs scalar multiplication of a cycle by an integer
	def __rmul__(self, alpha):
		c= copy(self)
		if type(alpha) is int:
			c.scalar *= alpha
		else:
			raise NotImplementedError
		return c

# Computes partitions of n-1 as used in finding the genrators for the cohomology classes
def find_partitions(k,n):
	if k == 1:
		return [[i] for i in range(n)]
	else:
		a1s = [i for i in range(n)]
		parts = []
		for a in a1s:
			a2s = find_partitions(k-1, a+1)
			a_parts = [[a] + a2sp for a2sp in a2s]
			parts += a_parts
		return parts
			
# Computes the cohomology generators for the Grassmannian G(k,n+1)
def compute_cohomology_groups(k,n):
	groups = [[] for i in range(2*k*(n+1-k) + 1)]
	partitions = find_partitions(k,n-k+2)
	for a in partitions:
		cycle = SchubertCycle(a=a, k=k, n=n)
		idx = cycle.grading()
		# print(idx)
		groups[idx].append(cycle)
	return groups
			
# Prints the cohomology groups with genrators
def print_cohomology_groups(k,n):
	groups = compute_cohomology_groups(k,n)
	n = groups[0][0].n
	k = groups[0][0].k
	for i in range(len(groups)):
		str_gens = [str(cycle) for cycle in groups[i]]
		zstr_gens = [f"ℤ{cycle}" for cycle in str_gens]
		if len(str_gens) == 0:
			print(f"H^{i} (G({k},{n+1})) = 0")
		else:
			print(f"H^{i} (G({k},{n+1})) = {'⊕'.join(zstr_gens)}")
			
# Prints the multiplication table of a corresponding ring or subring of the Grassmannian, allows for LaTeX output as well as console print
def print_multiplication_table(groups, latex=False):
	gens = []
	output = []
	max_length = 1
	for d in groups:
		for c in d:
			gens.append(c)
	if latex:
		strgens = ["\\cdot"]
	else:
		strgens = ['*']
	for g in gens:
		sg = str(g)
		l = len(str(g))
		if max_length < l:
			max_length = l
		strgens.append(sg)
		
	output.append(strgens)
	for c1 in gens:
		o = [str(c1)]
		for c2 in gens:
			c3s = str(c1*c2)
			l = len(c3s)
			if max_length < l:
				max_length = l
			o.append(c3s)
		output.append(o)
	if not latex:
		spacing = lambda n : ''.join([' ' for i in range(n)])
		padding = spacing(int(max_length/2))
		headers = True
		for o in output:
			line = "|"
			for e in o:
				diff = max_length - len(e)
				tstr = padding + e + spacing(diff) + padding + '|'
				line += tstr
			print(line)
			if headers:
				headers = False
				print(''.join(['-' for i in range(len(line))]))
	else:
		cols = ' | '.join(['c' for i in strgens])
		print (f'\\begin{{tabular}}{{ | {cols} | }}')
		print('\\hline')
		headers = True
		for o in output:
			line = ""
			for e in o:
				tstr = f" ${e}$ & "
				line += tstr
			tline = line[:-2].replace('σ', '\\sigma').replace('(', '{').replace(')', '}').replace('+-','-') + '\\\\'
			print(tline)
			print('\\hline')
		print ('\\end{tabular}')