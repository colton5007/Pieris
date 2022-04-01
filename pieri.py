
class CSum:
	def __init__(self):
		self.summands = []
		
	def add_cycle(self, other):
		existing = False
		for i in range(len(self.summands)):
			summand = self.summands[i]
			if summand == other:
				cycle = SchubertCycle(a=summand.a, b=summand.b, n=summand.n)
				cycle.scalar = summand.scalar + other.scalar
				if cycle.scalar == 0:
					del self.summands[i]
				else:
					self.summands[i] = cycle
					existing = True
				break
		if not existing:
			self.summands.append(other)
			
	def add(self, other):
		if type(other) is SchubertCycle:
			self.add_cycle(other)
		elif type(other) is CSum:
			for cycle in other.summands:
				self.add_cycle(cycle)
				
	def __add__(self, other):
		self.add(other)
		return self
	
	def __repr__(self):
		return str(self)
		
	def __str__(self):
		if len(self.summands) == 0:
			return "0"
		else:
			str_summands = [str(summand) for summand in self.summands]
			return "+".join(str_summands)

	def str_length(self):
		str_summands = [str(summand) for summand in self.summands]
		o = "+".join(str_summands)
		return len(o)

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

class SchubertCycle:
	def __init__(self, a=0,b=0,n=0):
		if n-1 >= a and a >= b and b>= 0:
			self.n = n
			self.a = a
			self.b = b
			self.scalar = 1
		else:
			raise Exception('Not defined: n-1>= a >= b >= 0')
		
	def set_expr(self):
		if self.b != 0:
			return f"{{ L ϵ G(2,{self.n + 1}) | L ∩ P(V_({self.n - self.a}) ≠ Ø, L ⊆ P(V_({self.n + 1 - self.b}) }}"
		else:
			return f"{{ L ϵ G(2,{self.n + 1}) | L ∩ P(V_({self.n - self.a}) ≠ Ø }}"
		
	def __str__(self):
		scalar_str = ""
		if self.scalar != 1 or (self.a == 0 and self.b == 0):
			scalar_str = str(self.scalar)
			
		if self.b != 0:
			return f"{scalar_str}σ_({self.a},{self.b})" 
		elif self.a != 0:
			return f"{scalar_str}σ_{self.a}"
		else:
			return f"{scalar_str}"
		
	def __eq__(self, other):
		return (self.a == other.a and self.b == other.b and self.n == other.n)
		
	def grading(self):
		return 2*(self.a+self.b)
		
	def __add__(self, other):
		if type(other) is SchubertCycle:
			cs = CSum()
			cs.add_cycle(self)
			cs.add_cycle(other)
			return cs
		elif type(other) is CSum:
			other.add_cycle(self)
			return other
			
	def __pow__(self, k):
		if k == 0:
			return SchubertCycle(a=0,b=0,n=self.n)
		elif k%2 == 1:
			return self*(self**(k-1))
		else:
			return ((self**(int(k/2)))*(self**(int(k/2))))

			
	def multiply(self, a,b,c, scalar):
		cs = CSum()
		for d in range(b,self.n):
			for e in range(c,b+1):
				if d+e == a+b+c:
					cycle = SchubertCycle(a=d, b=e, n=self.n)
					cycle.scalar = scalar
					cs.add(cycle)
		return cs
		
	def __repr__(self):
		return str(self)
	
	def __mul__(self, other):
		if type(other) is CSum:
			return other*self
		if self.b != 0:
			if other.b == 0:
				return self.multiply(other.a, self.a, self.b, self.scalar * other.scalar)
		if self.b == 0:
			return self.multiply(self.a, other.a, other.b, self.scalar*other.scalar)
			
		if self.grading() + other.grading() > 4*(self.n-1):
			return CSum()
		elif other.a == self.n - 1 - self.b and other.b == self.n - 1 - self.a:
			cs = CSum()
			cs.add_cycle(SchubertCycle(a=self.n-1,b=self.n-1, n=self.n))
			return cs
		else:
			return CSum()
		
			
def compute_cohomology_groups(n):
	groups = [[] for i in range(4*(n-1) + 1)]
	for b in range(n):
		for a in range(b,n):
			groups[2*(a+b)].append(SchubertCycle(a=a,b=b,n=n))
	return groups
			
def print_cohomology_groups(groups):
	td = len(groups)
	n = int((td-1)/4+1)
	for i in range(len(groups)):
		str_gens = [str(cycle) for cycle in groups[i]]
		zstr_gens = [f"ℤ{cycle}" for cycle in str_gens]
		if len(str_gens) == 0:
			print(f"H^{i} (G(2,{n+1})) = 0")
		else:
			print(f"H^{i} (G(2,{n+1})) = {'⊕'.join(zstr_gens)}")
			
def print_multiplication_table(groups):
	gens = []
	output = []
	max_length = 1
	for d in groups:
		for c in d:
			gens.append(c)
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
			c3 = c1*c2
			c3s = str(c1*c2)
			l = len(c3s)
			if max_length < l:
				max_length = l
			o.append(c3s)
		output.append(o)
	
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
	
		

