from comgen.constraintsystems import BaseSolver
from z3 import Int, Real, Bool, And, Or, Not, Implies, Sum, sat
import pymatgen.core as pg
# from comgen import PolyAtomicSpecies
from functools import partial
from comgen.util import composition_to_pettifor_dict
import fractions

class BalanceCharges:
	@staticmethod
	def charge_constraints(Species_Quantities: dict, ion_oxi_states: dict) -> list:      
		"""
		Enforce weighted sum of atom charges is zero.
		params:
			Species_Quantities: {ion_label: z3.Var}
			ion_oxi_states: {ion_label: int}
		"""  
		return [Sum([Species_Quantities[ion_label]*int(oxi_state) for ion_label, oxi_state in ion_oxi_states.items()]) == 0]

	@staticmethod
	def element_positive_defn_constraints(
		Element_Positive: dict,
		Species_Quantities: dict,
		ions: dict,
		ions_to_elements: dict) -> list:
		"""
		Record whether an element appears with positive charge.
		params:
			Species_Quantities: {ion_label: z3.Var}
			ions: {ion_label: pg.Species}
			ions_to_elements: {ion_label: element_label}
			Element_Positive: {element_label: z3.Var}
		"""
		cons = []

		for ion_label, ion in ions.items():
			if ion.oxi_state > 0:
				cons.append(
					Implies(
						Species_Quantities[ion_label] > 0, 
						Element_Positive[ions_to_elements[ion_label]]
					)
				)
			if ion.oxi_state < 0:
				cons.append(
					Implies(
						Species_Quantities[ion_label] > 0,
						Not(Element_Positive[ions_to_elements[ion_label]])
					)
				)
		
		return cons

	@staticmethod
	def electronegativity_constraints(
		Element_Positive: dict,
		elements: dict
		) -> list:
		"""
		Enforce cannot have more electronegative element is positive 
		while less electronegative element is negative. 
		TODO not sure this handles zero charge ideally. Rare, and usually comes with undefined e-neg. 
		undefined e-neg means el.X returns np.nan, so no constraints are added. 

		Applies to simple (mono) ions only.

		params:
			Element_Positive: {element_label: z3.Var}
			elements: {element_label: pg.Element}
		"""
		cons = []

		for el_label_1, el_1 in elements.items():
			for el_label_2, el_2 in elements.items():
				# X is electronegativity
				if el_1.X > el_2.X:
					cons.append(
						Implies(
							Element_Positive[el_label_1],
							Element_Positive[el_label_2]))

		return cons

class ElMD:
	PETTI_MAX = 103

	@staticmethod
	def emd_setup_constraints(
		Distance,
		Element_Quantities: dict,  
		Local_Diffs: dict,
		Abs_Diffs: dict,
		known: dict) -> list:
		"""
		Enforce ElMD relationship between known compositions and the composition to be determined.		
		TODO would a more efficient encoding help? i.e. ignoring positions we know are zero? 
		Should be easy to propagate constraints, but it's also a lot of real valued variables to cope with.
		"""
		cons = []
		# for i in range(len(Distance)):
		for c, known_petti_dict in known.items():
			cons.append(Local_Diffs[c][0] == 0)
			cons.append(Abs_Diffs[c][0] == 0)
			cons.append(Sum([d for _, d in Abs_Diffs[c].items()]) == Distance[c])

			for p_num in range(1, ElMD.PETTI_MAX+1):
				cons.append(
					Sum(
						Element_Quantities.get(p_num-1, 0),
						Local_Diffs[c][p_num-1],
						- known_petti_dict.get(p_num-1, 0)) == Local_Diffs[c][p_num])
			
				cons.append(
					Implies(
						Local_Diffs[c][p_num] >= 0, 
						Abs_Diffs[c][p_num] == Local_Diffs[c][p_num]))
				cons.append(
					Implies(
						Local_Diffs[c][p_num] < 0,
						Abs_Diffs[c][p_num] ==  -Local_Diffs[c][p_num]))

		return cons

	@staticmethod	
	def lower_bound_emd(Distance, lb) -> list:
		"""
		Enforce minimum value for ElMD between known compositions and the composition to be determined.
		"""
		return [Or(*[d >= lb for d in Distance.values()])]

	@staticmethod	
	def upper_bound_emd(Distance, ub) -> list:
		"""
		Enforce maximum value for ElMD between known compositions and the composition to be determined.		
		"""
		return [Or(*[d <= ub for d in Distance.values()])]

class Elements:
	@staticmethod
	def known_element_quantity_constraints(
		Element_Quantities: dict,
		known_element_quantities: dict) -> list:
		"""
		Enforce any known quantities.
		With mono ions that just means sum of ion quantities equals element quantity.
		With poly ions also consider how many atoms of this element are in the poly ion. 

		params:
			Element_Quantities: {element_label: z3.Var}
			known_element_quantities: {element_label: (lower_bound, upper_bound)} both bounds are int|float
			normed: do element quantities sum to 1? (if not, expect that they are Integer types)
		"""
		cons = []

		for el_label, quant in known_element_quantities.items():
			lb, ub = quant[0], quant[1]
			cons.extend([Element_Quantities[el_label] >= lb, Element_Quantities[el_label] <= ub])
		
		return cons

	@staticmethod
	def element_group_quantity_constraints(
		Element_Quantities: dict,
		bounds: tuple
	) -> list:
		"""
		Total quantity across the set of elements must be between given bounds.

		params:
		Element_Quantities: {element_label: z3.Var}
		bounds: (lower_bound, upper_bound), bounds are int|float
		"""
		cons = []
		cons.append(Sum(list(Element_Quantities.values())) <= bounds[1])
		cons.append(Sum(list(Element_Quantities.values())) >= bounds[0])
		return cons

	@staticmethod
	def non_negativity_constraints(Vars: dict) -> list:
		return [And([q >= 0 for q in Vars.values()])]

	@staticmethod
	def species_quantity_defn_constraints(
		Species_Quantities: dict,
		Element_Quantities: dict,
		element_ion_weights: dict) -> list:
		"""
		Species_Quantities are positive.
		And require that the quantities of each element and of each species agree. 
		i.e. sum of quantities of the various Fe ions equals quantity of Fe element.

		params:
			Species_Quantities: {ion_label}: z3.Var}
			Element_Quantities: {element_label: z3.Var}
			element_ion_weights: {element_label: {ion_label: multiplier}}
		"""
		cons = Elements.non_negativity_constraints(Species_Quantities)

		for el_label, ion_weights in element_ion_weights.items():
			cons.append(
				Sum(
					[Species_Quantities[ion_label]*m for ion_label, m in ion_weights.items()]
				) == Element_Quantities[el_label])

		return cons

	@staticmethod
	def element_present_defn_constraints(
		Element_Present: dict,
		Element_Quantities: dict) -> list:
		"""
		Definition constraints. Element is "present" if and only if its quantity is greater than 0. 
		"""
		cons = [Implies(q > 0, Element_Present[el_label] == 1) for el_label, q in Element_Quantities.items()]
		cons.extend([Implies(q <= 0, Element_Present[el_label] == 0) for el_label, q in Element_Quantities.items()])
		
		return cons

	@staticmethod
	def normed_quantity_constraints(Vars: dict) -> list:
		"""
		This set of variables represents a normed distribution and so sums to 1. 
		"""
		return [Sum(list(Vars.values())) == 1]

	@staticmethod
	def normed_element_quantity_defn_constraints(Element_Quantities: dict) -> list:
		cons = Elements.non_negativity_constraints(Element_Quantities)
		cons.extend(Elements.normed_quantity_constraints(Element_Quantities))
		return cons

	@staticmethod
	def distinct_elements_min_constraints(Element_Present: dict, min_elements: int) -> list:
		"""
		Fix how many elements are permitted. TODO Special case of select_from_set_min_constraints
		"""
		return [Sum(list(Element_Present.values())) >= min_elements]

	@staticmethod
	def distinct_elements_max_constraints(Element_Present: dict, max_elements: int) -> list:
		"""
		Fix how many elements are permitted.
		"""
		return [Sum(list(Element_Present.values())) <= max_elements]

	@staticmethod
	def quantity_integer_multiplier_constraints(
		Element_Quantities: dict, 
		Element_Integer_Quantities: dict,
		ub: int) -> list:
		"""
		Find a multiplier within given bounds 
		s.t. q * multiplier is an integer for all the normed quantities q. 
		equivalent to constraining the sum of quantities but for Real (normed) quantities instead of Integer. 
		"""
		cons = []
		for multiplier in range(1, ub+1):
			and_cons = []
			for el_label, q in Element_Quantities.items():
				and_cons.append(q * int(multiplier) == Element_Integer_Quantities[el_label])
			cons.append(And(*and_cons))

		return [Or(*cons)]

	@staticmethod
	def sum_of_quantities_constraints(
		Element_Quantities: dict, 
	    bounds: tuple) -> list:
		"""
		This should only be used when Element_Quantities variables are Integer types. 
		Otherwise bound the multiplier to transform the Real quantities to Integer - see quantity_integer_multiplier_constraints.
		"""
		return [Sum(Element_Quantities.values()) >= bounds[0], Sum(Element_Quantities.values()) <= bounds[1]]

	@staticmethod
	def select_from_set_min_constraints(
		Element_Present: dict, 
		elements: set,
		num: int) -> list:
		"""
		Include a non-zero quantity of at least num elements from the given sets.
		"""
		var_list = []
		for el_label in elements:
			if (var := Element_Present.get(el_label)) is not None:
				var_list.append(var)

		return [Sum(*var_list) >= num]

	@staticmethod
	def select_from_set_max_constraints(
		Element_Present: dict, 
		elements: set,
		num: int) -> list:
		"""
		Include a non-zero quantity of at most num elements from the given sets.
		"""
		var_list = []
		for el_label in elements:
			if (var := Element_Present.get(el_label)) is not None:
				var_list.append(var)

		return [Sum(*var_list) <= num]


	@staticmethod
	def total_quantity_from_set_constraints(
		Element_Quantities: dict, 
		elements: set,
		bounds: tuple) -> list:
		"""
		The total quantity (or ratio to whole, if normed quantities) of elements in this set.
		
		params:
			Element_Quantities: {el_label: z3.Var}
			elements: {el_label}
			bounds: (lower_bound, upper_bound)
		"""
		incl = []
		for el_label in elements:
			incl.append(Element_Quantities[el_label])

		return [Sum(incl) >= bounds[0], Sum(incl) <= bounds[1]]

	@staticmethod
	def relative_quantity_from_sets_constraints(
		Element_Quantities: dict,
		elements_1: set,
		elements_2: set,
		bounds: tuple) -> list:
		"""
		Require that total quantity of elements in elements_1 / total quantity of elements in elements_2 
		is between bounds[0] and bounds[1]. 
		TODO not used currently!! Need to expose these constraints in API
		"""
		cons = []

		incl_1 = [q for el_label, q in Element_Quantities.items() if el_label in elements_1]
		incl_2 = [q for el_label, q in Element_Quantities.items() if el_label in elements_2]

		# sum(q1) / sum(q2) >= lb
		# sum(q1) - sum(q2)*lb >= 0
		# sum(q1 \cup -q2*lb) >= 0
		cons.append(Sum([q*bounds[0] for q in incl_2]+incl_1) >= 0)
		cons.append(Sum([q*bounds[1] for q in incl_2]+incl_1) <= 0)

		return cons


	@staticmethod
	def starting_material_constraints(
		Element_Quantities: dict,
		Ingredient_Weights: dict,
		ingredient_compositions: dict
	) -> list:
		"""
		Given some input starting materials (ingredients) require that the final composition can be made from some combination of these ingredients.
		i.e. weighted (by Ingredient_Weights) sum of ingredient compositions equals final composition (specified by Element_Quantities)

		params:
			Element_Quantities: {el_label: z3.Var}
			Ingredient_Weights: {composition: z3.Var}
			ingredient_compositions: {composition: {el_label: quantity}}
		"""
		cons = []
		for el, q in Element_Quantities.items():
			weighted_ingredients = [ingredient_compositions[comp][el]*Ingredient_Weights[comp] for comp in Ingredient_Weights.keys()]
			cons.append(Sum(weighted_ingredients) == q)
		return [And(*cons)]
	
	@staticmethod
	def ingredient_definition_constraints(
		Ingredient_Weights: dict
	) -> list:
		"""
		Only positive quantities of each ingredient composition are allowed.
		"""
		return [w >= 0 for w in Ingredient_Weights.values()]

class IonicCompositionGenerator(BaseSolver):
	def __init__(self, ions=None, precision=0.1):
		super().__init__()
		
		if ions is None: # this is a hack to give a more specific error message.
			raise TypeError("Missing input for permitted ions. Please provide a SpeciesCollection.")
		
		self._ions = ions
		self._elements = self._ions.group_by_element_view().keys()

		self.precision = precision # TODO check if this is useful. 

		self.constraints = []
		self.constraints_summary = []
		self._set_basic_constraints()
		
	@property
	def element_labels(self):
		# use pettifor number instead of element symbol
		return [int(el.mendeleev_no) for el in self._elements] # TODO check if any elements are missing a mendeleev_no... what then? 
	
	# TODO consider functools cached_property instead - avoid a little bit of repeated effort.
	@property
	def ion_labels(self):
		return [str(ion) for ion in self._ions.ungrouped_view()]

	def element_ion_weights(self): # TODO not sure I like this here
		out = {}

		for el, ions in self._ions.group_by_element_view().items():
			ion_weights = {}
			for ion in ions: 
				if isinstance(ion, pg.Species):
					ion_weights[str(ion)] = 1 
				# elif isinstance(ion, PolyAtomicSpecies):
				else:
					ion_weights[str(ion)] = ion.multiplier(el)

			out[int(el.mendeleev_no)] = ion_weights

		return out

	def _element_quantity_variables(self, ids=None):
		eq_vars = self._variables(
			'Element_Quantities', 
			Real, 
			self.element_labels,
			Elements.normed_element_quantity_defn_constraints)
		if ids is None:
			return eq_vars
		return {id: eq_vars[id] for id in ids}
	
	def _element_present_variables(self):
		init_func = partial(
			Elements.element_present_defn_constraints, 
			Element_Quantities = self._element_quantity_variables())

		return self._variables(
			'Element_Present', 
			Int, # use an integer rather than Bool to enable constraints on count of "True" values
			self.element_labels, 
			init_func)		

	def _species_quantity_variables(self):
		init_func = partial(
			Elements.species_quantity_defn_constraints,
			Element_Quantities = self._element_quantity_variables(),
			element_ion_weights = self.element_ion_weights())

		return self._variables('Species_Quantities', Real, self.ion_labels, init_func)
	
	def _element_positive_variables(self):
		mono_ions = {str(ion): ion for ion in self._ions.filter_mono_species().ungrouped_view()}
		mono_ions_to_elements = {str(ion): int(ion.element.mendeleev_no) for ion in self._ions.filter_mono_species().ungrouped_view()}
	
		init_func = partial(
			BalanceCharges.element_positive_defn_constraints,
			Species_Quantities = self._species_quantity_variables(),
			ions = mono_ions,
			ions_to_elements = mono_ions_to_elements)

		return self._variables('Positive', Bool, self.element_labels, init_func)

	def _set_basic_constraints(self):
		"""
		The standard constraints that should always be required for ionic materials. 
		- Balanced charges
		- Electronegativity respected by charge assignment
		"""
		Species_Quantities = self._species_quantity_variables()
		Element_Positive = self._element_positive_variables()

		ion_oxi_states = {str(ion): ion.oxi_state for ion in self._ions.ungrouped_view()}
		self.constraints.extend(
			BalanceCharges.charge_constraints(Species_Quantities, ion_oxi_states))
		
		elements = {int(el.mendeleev_no): el for el in self._elements}
		self.constraints.extend(
			BalanceCharges.electronegativity_constraints(
				Element_Positive,
				elements))
		
		self.constraints_summary.append("Balanced charges.")
		self.constraints_summary.append("Charges respect electronegativity.")

	def max_total_atoms(self, ub):
		Element_Integer_Quantities = self.new_variables(
			'Element_Int_Count',
			Int, 
			self.element_labels)
		
		self.constraints.extend(
			Elements.quantity_integer_multiplier_constraints(
				self._element_quantity_variables(),
				Element_Integer_Quantities,
				ub))
		
		self.constraints_summary.append(f"Total atoms at most {ub}.")

	def distinct_elements(self, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')

		if exact is not None: lb, ub = exact, exact

		if lb is not None:
			self.constraints.extend(
				Elements.distinct_elements_min_constraints(
					self._element_present_variables(), 
					lb))
			self.constraints_summary.append(f"Number of distinct elements at least {lb}.")

		if ub is not None:
			self.constraints.extend(
				Elements.distinct_elements_max_constraints(
					self._element_present_variables(), 
					ub))			
			self.constraints_summary.append(f"Number of distinct elements at most {ub}.")

	def include_element_from(self, element_set, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			exact = 1
			# raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		
		if exact is not None: lb, ub = exact, exact
		
		# expect element_set will be pg.Elements or symbols. 
		# But the solver is using pettifor values to refer to elements. 
		element_set_ids = {int(pg.Element(el).mendeleev_no) for el in element_set}
		
		if lb is not None:
			self.constraints.extend(
				Elements.select_from_set_min_constraints(
					self._element_present_variables(), 
					element_set_ids,
					lb))
			self.constraints_summary.append(f"Include at least {lb} element(s) from {element_set}.")			
		if ub is not None:
			self.constraints.extend(
				Elements.select_from_set_max_constraints(
					self._element_present_variables(), 
					element_set_ids,
					ub))					
			self.constraints_summary.append(f"Include at most {ub} element(s) from {element_set}.")			

	# TODO specify bounds on quantity across a set of elts
	def fix_elements_quantity(self, elts: set, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		
		if ub is not None and lb is None: lb = 0
		if lb is not None and ub is None: ub = 1
		if exact is not None: lb, ub = exact, exact

		self.constraints_summary.append(f"Fix quantity of elements {elts} between {lb} and {ub}.")

		if isinstance(elts, str) or isinstance(elts, pg.Element): 
			elts = {elts} # in case a single element symbol is passed in. 
		elts = {int(pg.Element(el).mendeleev_no) for el in elts}

		self.constraints.extend(
			Elements.element_group_quantity_constraints(
				self._element_quantity_variables(elts),
				(lb, ub)))

	def construct_from(self, compositions):
		self.constraints_summary.append(f"Construct from a weighted sum of {compositions}.")

		composition_dicts = {str(comp): composition_to_pettifor_dict(comp) for comp in compositions}

		# create local variables - not for reuse outside of these constraints
		Ingredient_Weights = self.new_variables('Ingredient_Weight', Real, compositions, Elements.ingredient_definition_constraints)

		self.constraints.extend(
			Elements.starting_material_constraints(
				self._element_quantity_variables(),
				Ingredient_Weights,
				composition_dicts))


	def emd_comparison_compositions(self, compositions, *, lb=None, ub=None):
		if lb is None and ub is None:
			raise ValueError("Either a lower or upper bound (or both) are required for the ElMD from the given compositions.")
		if not compositions:
			raise ValueError("Comparison compositions to measure ElMD from are required.")
		
		if lb is not None:
			self.constraints_summary.append(f"At least {lb} distance from one of {compositions}.")
		if ub is not None:
			self.constraints_summary.append(f"At most {ub} distance from one of {compositions}.")

		if not (isinstance(compositions, list) or isinstance(compositions, set)):
			compositions = [compositions]
		
		composition_dicts = {str(comp): composition_to_pettifor_dict(comp) for comp in compositions}

		# set up local variables - not for reuse outside of these constraints
		Distance = self.new_variables('EMD', Real, compositions)
				
		Local_Diffs, Abs_Diffs = {}, {}
		for comp in compositions:
			c = str(comp)
			Local_Diffs[c] = self.new_variables(f'Local_{c}', Real, list(range(ElMD.PETTI_MAX+1)))
			Abs_Diffs[c] = self.new_variables(f'Abs_{c}', Real, list(range(ElMD.PETTI_MAX+1)))

		self.constraints.extend(
			ElMD.emd_setup_constraints(
				Distance,
				self._element_quantity_variables(),
				Local_Diffs,
				Abs_Diffs,
				composition_dicts))

		if lb is not None:
			self.constraints.extend(
				ElMD.lower_bound_emd(Distance, lb))
		if ub is not None:
			self.constraints.extend(
				ElMD.upper_bound_emd(Distance, ub))

	def get_constraints_summary(self):
		return self.constraints_summary

	def exclude_solution(self, solution):
		Element_Quantities = self._element_quantity_variables() # el_label : var
		exclusions = []
		
		for quant_var in Element_Quantities.values():
			quantity = solution[quant_var]
			if self.precision:
				exclusions.append(quant_var <= quantity - self.precision)
				exclusions.append(quant_var >= quantity + self.precision)
			else:
				exclusions.append(Not(quant_var == quantity))

		self.solver.add(Or(*exclusions))

	def exclude_composition(self, composition):
		comp = composition_to_pettifor_dict(composition) # el_label : normed_quantity
		Element_Quantities = self._element_quantity_variables() # el_label : var
		# check comp can be represented
		for el in comp.keys():
			if not el in Element_Quantities.keys():
				return # nothing to do - can't generate this anyway. 
		comp_as_solution = {var: comp.get(el_label, 0.0) for el_label, var in Element_Quantities.items()} # var: normed_quantity
		self.exclude_solution(comp_as_solution)

	def format_solution(self, model, as_frac=False):
		Element_Quantities = self._element_quantity_variables()
		out = {}
		for el in self._elements:
			quant = model[Element_Quantities[int(el.mendeleev_no)]]
			if quant.numerator_as_long() != 0:
				if as_frac:
					out[str(el)] = str(fractions.Fraction(quant.numerator_as_long(), quant.denominator_as_long()))
				else:
					out[str(el)] = round(float(quant.numerator_as_long()) / float(quant.denominator_as_long()), 3)

		# if self.distance_var:
		# 	for i, emd in self.distance_var.items():
		# 		out[f'emd_{i}'] = solution[emd]
		# 		for j, var in self.abs_diff_var[i].items():
		# 			if solution[var].numerator_as_long() < 0:
		# 				out[f'abs_diff_{i}_{j}'] = solution[var]
		# 				out[f'local_diff_{i}_{j}'] = solution[self.local_diff_var[i][j]]


		# 				return out
		return out
		# return {str(el): solution[Element_Quantities[int(el.mendeleev_no)]] for el in self._elements}

	def get_next(self, max_results = 1, *, as_frac=False):
		solutions = []
		for con in self.constraints:
			self.solver.add(con) # pass any new constraints to the solver.
		self.constraints = [] # reset to avoid re-adding same constraints later.
		while self.solver.check() == sat and len(solutions) < max_results:
			solution = self.format_solution(self.solver.model(), as_frac)
			solutions.append(solution)
			self.exclude_solution(self.solver.model())

		return solutions

	def solve(self):
		Quantities = self._species_quantity_variables()

		for con in self.constraints:
			self.solver.add(con) # pass any new constraints to the solver.
		self.constraints = []
		
		while self.solver.check() == sat:
			solution = self.format_solution(self.solver.model(), as_frac=False)
			# print(self.solver.model()[Int('Multiplier_0')])
			for name, var in Quantities.items():
				print(f'{name}: {self.solver.model()[var]}')
			
			# model = self.solver.model()
			# with open('output.txt', 'w') as f:
			# 	for sp, var in Quantities.items():
			# 		r = model[var]
			# 		f.write(str(round(float(r.numerator_as_long())/float(r.denominator_as_long()),3)))
			# 	f.write(str(model['Multiplier'].as_long()))
			
			self.exclude_solution(self.solver.model())
			yield solution
