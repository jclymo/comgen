from z3 import Solver, sat

class BaseSolver:
	def __init__(self):
		self.solver = Solver()
		self._variable_cache = {}
		self.constraints = []
 
	def _variables(self, var_name: str, var_type, ids: list, init_func=None) -> dict:
		vars = self._variable_cache.get(var_name)
		if vars is None:
			self._variable_cache[var_name] = self.new_variables(var_name, var_type, ids, init_func)
		return self._variable_cache[var_name]

	def new_variables(self, var_name: str, var_type, ids: list, init_func=None) -> dict:
		if not isinstance(ids[0], int): 
			ids = [str(idx) for idx in ids]
		vars = {idx: var_type(f'{var_name}_{idx}') for idx in ids}
		if init_func:
			self.constraints.extend(init_func(vars)) # constraints that define these variables wrt others
		return vars
			
	# def unset_variables(self, var_name: str) -> None:
	# 	local_variables = self._variable_cache.pop(var_name, None)
	# 	if local_variables is None:
	# 		raise ValueError("Attempting to unset variables that do not exist with name={var_name}")
		
	def remap_variables(self, var_name: str, mapping: dict) -> dict:
		"""
		Return dictionary of existing variables with new keys.

		params:
			var_name: references item in _variable_cache
			mapping: {old_key: new_key}
		"""
		local_variables = self._variable_cache[var_name]
		return {str(mapping[idx]): var for idx, var in local_variables.items()}

	# def build(self):
	# 	raise NotImplementedError('Missing required method to build constraint system.')

	def solve(self):
		# raise NotImplementedError('Missing required method to solve constraint system.')
		if self.solver.check() == sat:
			return self.solver.model()
		return self.solver.check()




