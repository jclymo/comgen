from comgen.constraintsystems.base import BaseSolver, Abs
from z3 import Int, Real, Bool, If, And, Or, Not, Implies, Sum, sat, Solver, Q
import pymatgen.core as pg
# from comgen import PolyAtomicSpecies
from functools import partial
from comgen.util import composition_to_pettifor_dict
import fractions
from comgen.species import SpeciesCollection

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

				cons.append(Abs_Diffs[c][p_num] == Abs(Local_Diffs[c][p_num]))

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
		lb: float,
		ub: float) -> list:
		"""
		Total quantity across the set of elements must be between given bounds.

		params:
		Element_Quantities: {element_label: z3.Var}
		lb: lower bound
		ub: upper bound
		"""
		cons = []
		cons.append(Sum(list(Element_Quantities.values())) <= ub)
		cons.append(Sum(list(Element_Quantities.values())) >= lb)
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
			Species_Quantities: {ion_label: z3.Var}
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
	def element_quantity_defn_constraints(Element_Quantities: dict) -> list:
		cons = Elements.non_negativity_constraints(Element_Quantities)
		cons.extend(Elements.normed_quantity_constraints(Element_Quantities))
		return cons

	@staticmethod
	def distinct_elements_constraints(
		Element_Present: dict,
		lb: int,
		ub: int) -> list:
		"""
		Fix how many elements are permitted.
		"""
		cons = []
		cons.append(Sum(list(Element_Present.values())) >= lb)
		cons.append(Sum(list(Element_Present.values())) <= ub)
		return cons

	@staticmethod
	def select_from_set_constraints(
		Element_Present: dict, 
		elements: set,
		lb: int,
		ub: int) -> list:
		"""
		Constrain the number of elements included from the given sets. 
		"""
		var_list = []
		for el_label in elements:
			if (var := Element_Present.get(el_label)) is not None:
				var_list.append(var)
		
		cons = []
		cons.append(Sum(var_list) >= lb)
		cons.append(Sum(var_list) <= ub)
		return cons

	@staticmethod
	def relative_quantity_from_sets_constraints(
		Element_Quantities: dict,
		elements_1: set,
		elements_2: set,
		lb: float, 
		ub: float) -> list:
		"""
		Require that total quantity of elements in elements_1 / total quantity of elements in elements_2 
		is between lb and ub.
		TODO not used currently!! Need to expose these constraints in API
		"""
		cons = []

		incl_1 = [q for el_label, q in Element_Quantities.items() if el_label in elements_1]
		incl_2_lower = [-q*lb for el_label, q in Element_Quantities.items() if el_label in elements_2]
		incl_2_upper = [-q*ub for el_label, q in Element_Quantities.items() if el_label in elements_2]

		# sum(q1) / sum(q2) >= lb
		# sum(q1) - sum(q2)*lb >= 0
		# sum(q1 \cup -q2*lb) >= 0
		cons.append(Sum(incl_2_lower + incl_1) >= 0)
		cons.append(Sum(incl_2_upper + incl_1) <= 0)

		return cons


	@staticmethod
	def starting_material_constraints(
		Element_Quantities: dict,
		Ingredient_Quantities: dict,
		ingredient_compositions: dict
	) -> list:
		"""
		Given some input starting materials (ingredients) require that the final composition can be made from some combination of these ingredients.
		i.e. weighted (by Ingredient_Quantities) sum of ingredient compositions equals final composition (specified by Element_Quantities)

		params:
			Element_Quantities: {el_label: z3.Var}
			Ingredient_Quantities: {composition: z3.Var}
			ingredient_compositions: {composition: {el_label: quantity}}
		"""
		cons = []
		for el, q in Element_Quantities.items():
			weighted_ingredients = [ingredient_compositions[comp].get(el, 0)*i_quant for comp, i_quant in Ingredient_Quantities.items()]
			cons.append(Sum(weighted_ingredients) == q)
		return [And(*cons)]

	@ staticmethod
	def starting_material_cost_constraints(
		Ingredient_Quantities: dict,
		ingredient_costs: dict,
		total_cost: float) -> list:
		"""
		Total cost of selected ingredients must be less than specified bound. 
		Ingredient costs are expected to be per mole of stated composition. 
			Ingredient_Quantities: {composition: z3.Var}
			ingredient_costs: {composition: float}
			total_cost: float
		"""
		return [Sum([ingredient_costs[comp]*q for comp, q in Ingredient_Quantities.items()]) <= total_cost]

	
	@staticmethod
	def ingredient_definition_constraints(
		Ingredient_Quantities: dict
	) -> list:
		"""
		Only non-negative quantities of each ingredient composition are allowed.
		"""
		return [w >= 0 for w in Ingredient_Quantities.values()]
	
	@staticmethod
	def total_atom_constraints(
		Total_Atoms,
		Species_Quantities: dict,
		lb: int,
		ub: int
	) -> list:
		"""
		Enforce total number of atoms, i.e. quantity of each species must be rational with denominator matching number of atoms. 
			Total_Atoms: z3.Var
			Species_Quantities: {species_label: z3.Var}
			lb: lower bound on number of atoms int
			ub: upper bound on number of atoms int
		"""
		cons = []

		for num_atoms in range(lb, ub+1):
			for q in Species_Quantities.values():
				select_q_from = Or(*[q == Q(n, num_atoms) for n in range(num_atoms+1)])
			cons.append(Implies(Total_Atoms == num_atoms, select_q_from))

		cons.append(Total_Atoms >= lb)
		cons.append(Total_Atoms <= ub)

		return cons
	
	@staticmethod
	def select_species_pair_constraints(
		Species_Pairs: dict,
		Species_Quantities: dict) -> list:
		"""
		Species_Pairs: {(species_label1, species_label2): z3.Var}
		Species_Quantities: {species_label: z3.Var}
		"""
		cons = []
		cons.append(Or(*[sp == True for sp in Species_Pairs.values()]))
		for pair, sp in Species_Pairs.items():
			a, b = pair
			cons.append(Implies(sp == True, Species_Quantities[a] > 0))
			cons.append(Implies(sp == True, Species_Quantities[b] > 0))
		return cons

class IonicCompositionGenerator(BaseSolver):
	def __init__(self, ions=None, precision=0.1):
		super().__init__()
		
		if not isinstance(ions, SpeciesCollection): 
			raise TypeError("Missing input for permitted ions. Please provide a SpeciesCollection.")
		
		self._ions = ions
		self._elements = self._ions.group_by_element_view().keys()

		self.precision = precision # TODO make this a variable dependent on number of atoms. 

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
		"""
		Retrieve or create variables for element quantities.
		Pass in a set of ids to retrieve only those variables (but will always create all variables).
		Note that if ids contains an id not in the list of element labels, it will be ignored.
		"""
		eq_vars = self._variables(
			'Element_Quantities', 
			Real, 
			self.element_labels,
			Elements.element_quantity_defn_constraints)
		if ids is None:
			return eq_vars
		return {id: eq_vars[id] for id in ids if id in eq_vars.keys()}
	
	def _element_present_variables(self):
		init_func = partial(
			Elements.element_present_defn_constraints, 
			Element_Quantities = self._element_quantity_variables())

		return self._variables(
			'Element_Present', 
			Int, # use an integer rather than Bool to enable constraints on count of "True" values
			self.element_labels, 
			init_func)		

	def _species_quantity_variables(self, ids=None):
		init_func = partial(
			Elements.species_quantity_defn_constraints,
			Element_Quantities = self._element_quantity_variables(),
			element_ion_weights = self.element_ion_weights())
	
		sps_vars = self._variables('Species_Quantities', Real, self.ion_labels, init_func)
	
		if ids is None:
			return sps_vars
		return {id: sps_vars[id] for id in ids if id in sps_vars.keys()}
	
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

	def total_atoms(self, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')

		if exact is not None: lb, ub = exact, exact

		Total_Atoms = self.new_variables(
			'Total_Atoms',
			Int,
			[0] # fake id because can't yet support no ids
		)[0]

		self.constraints.extend(
			Elements.total_atom_constraints(
				Total_Atoms,
				self._species_quantity_variables(),
				lb,
				ub))

		self.constraints_summary.append(f"Total atoms between {lb} and {ub}.")

	def distinct_elements(self, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')

		if exact is not None: lb, ub = exact, exact
		if ub is not None and lb is None: lb = 0
		if lb is not None and ub is None: ub = 120 # i.e. more than the number of known elements

		self.constraints.extend(
			Elements.distinct_elements_constraints(
				self._element_present_variables(),
				lb,
				ub))
		
		self.constraints_summary.append(f"Number of distinct elements between {lb} and {ub}.")

	def include_element_from(self, element_set, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			exact = 1
			# raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		
		# if exact is not None: lb, ub = exact, exact
		
		if exact is not None: lb, ub = exact, exact	
		if ub is not None and lb is None: lb = 0
		if lb is not None and ub is None: ub = len(element_set)
		
		# expect element_set will be pg.Elements or symbols. 
		# But the solver is using pettifor values to refer to elements. 
		element_set_ids = {int(pg.Element(el).mendeleev_no) for el in element_set}
		
		self.constraints.extend(
			Elements.select_from_set_constraints(
				self._element_present_variables(),
				element_set_ids,
				lb,
				ub
			)
		)
		self.constraints_summary.append(f"Include between {lb} and {ub} elements from {element_set}.")	

	# TODO specify bounds on quantity across a set of elts
	def fix_elements_quantity(self, elts: set, exact=None, *, lb=None, ub=None):
		if exact is not None and (lb is not None or ub is not None):
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		if exact is None and lb is None and ub is None:
			raise ValueError('Please provide exactly one of: a) exact quantity b) lower (lb) and / or upper (ub) bounds on quantity')
		
		if ub is not None and lb is None: lb = 0
		if lb is not None and ub is None: ub = 1
		if exact is not None: lb, ub = exact, exact

		self.constraints_summary.append(f"Fix total quantity of elements from {elts} between {lb} and {ub}.")

		if isinstance(elts, str) or isinstance(elts, pg.Element): 
			elts = {elts} # in case a single element symbol is passed in. 
		elts = {int(pg.Element(el).mendeleev_no) for el in elts}

		self.constraints.extend(
			Elements.element_group_quantity_constraints(
				self._element_quantity_variables(elts),
				lb,
				ub))

	def construct_from(self, compositions, costs=None, total_cost=None):
		"""
		Specify ingredient compositions from which the new composition can be made. 
		Optionally, specify costs of each ingredient and a bound on the total cost. 
			compositions: list of pg.Composition
			costs: list of floats; if specified must have same length and order as compositions.
			total_cost: float
		"""
		self.constraints_summary.append(f"Construct from a weighted sum of {compositions}.")
		composition_dicts = {str(comp): composition_to_pettifor_dict(comp) for comp in compositions}

		# create local variables - not for reuse outside of these constraints
		Ingredient_Quantities = self.new_variables('Ingredient_Quantities', Real, compositions, Elements.ingredient_definition_constraints)

		self.constraints.extend(
			Elements.starting_material_constraints(
				self._element_quantity_variables(),
				Ingredient_Quantities,
				composition_dicts))
		
		if costs is not None:
			assert len(costs) == len(compositions)
			assert total_cost is not None
			costs = {str(comp): cost for comp, cost in zip(compositions, costs)}
			self.constraints_summary.append(f"Total cost of ingredients cannot exceed {total_cost}.")

		self.constraints.extend(
			Elements.starting_material_cost_constraints(
				Ingredient_Quantities,
				costs,
				total_cost))
		

	def radius_ratio(self, sps_1, sps_2, lb=None, ub=None):
		"""
		sps_1 and sps_2 are {sp_label: radius}

		require that something is selected from sps_1, something is selected from sps_2, 
		and the ratio radius_1 / radius_2 is between lb and ub.
		"""
		if lb is None and ub is None:
			raise ValueError('Please provide at least one of lower (lb) and upper (ub) bounds on radius ratios.')

		pairs = []
		for sp1, r1 in sps_1.items():
			for sp2, r2 in sps_2.items():
				if lb is None or r1 / r2 >= lb:
					if ub is None or r1 / r2 <= ub:
						pairs.append((sp1, sp2))
		
		# TODO check what happens if multiple of this type of constraint is used -- will vars conflict? 
		Species_Pairs = self.new_variables('Species_Pairs', Bool, pairs) 

		self.constraints.extend(
			Elements.select_species_pair_constraints(
				Species_Pairs,
				self._species_quantity_variables()))

	def emd_comparison_compositions(self, compositions, *, lb=None, ub=None):
		if lb is None and ub is None:
			raise ValueError("Either a lower or upper bound (or both) are required for the ElMD from the given compositions.")
		if not compositions:
			raise ValueError("Comparison compositions to measure ElMD from are required.")
		
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
			self.constraints_summary.append(f"At least {lb} 1-D pettifor earth movers distance from at least one of {compositions}.")			
		if ub is not None:
			self.constraints.extend(
				ElMD.upper_bound_emd(Distance, ub))
			self.constraints_summary.append(f"At most {ub} 1-D pettifor earth movers distance from at least one of {compositions}.")			
			

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

		self.constraints.append(Or(*exclusions))

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

	def get_next(self, *, as_frac=False):
		solutions = []
		s = Solver() # create a new solver each time because incremental mode makes solving very slow
		for con in self.constraints:
			s.add(con)
		if s.check() == sat:
			solution = self.format_solution(s.model(), as_frac)
			solutions.append(solution)
			self.exclude_solution(s.model()) # doesn't actually exclude them yet, ready for next call to get_next. 

		return solutions
