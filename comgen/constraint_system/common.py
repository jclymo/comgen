from z3 import And, Sum, If

def check_bounds(exact, lb, ub):
    bounds = [exact, lb, ub]
    assert not all([b is None for b in bounds])
    assert not (exact is not None and lb is not None)
    assert not (exact is not None and ub is not None)

def apply_bounds(var, exact=None, *, lb=None, ub=None):
    check_bounds(exact, lb, ub)
    
    if exact is not None:
        return var == exact
    
    constraints = []
    if lb is not None:
        constraints.append(var >= lb)
    if ub is not None:
        constraints.append(var <= ub)
    
    if len(constraints) > 1:
        return And(constraints)
    return constraints[0]

def zero_weighted_sum(vars, weights):
    return weighted_sum(vars, weights, 0)

def weighted_sum(vars, weights, exact=None, *, lb=None, ub=None):
    assert len(vars) == len(weights)
    weighted_vars = [v * w for v, w in zip(vars, weights)]
    return apply_bounds(Sum(weighted_vars), exact, lb=lb, ub=ub)

def bound_weighted_average_value_ratio(weight_vars1, weight_vars2, values1, values2, exact=None, *, lb=None, ub=None):
    """Assumes all vars to be non-negative.
    """
    ratio = weighted_sum(weight_vars1, values1) * Sum(weight_vars2) / weighted_sum(weight_vars2, values2) / Sum(weight_vars1)
    return apply_bounds(ratio, exact, lb=lb, ub=ub)

def Abs(x):
	return If(x >= 0, x, -x)

def ReLU(x):
	return If(x >= 0, x, 0)