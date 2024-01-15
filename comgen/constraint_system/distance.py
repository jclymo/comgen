from comgen.constraint_system.common import apply_bounds, Abs
from z3 import Sum, Real

class Distance:
    def _setup_distance_calculation(self, object_vars):
        raise NotImplementedError
    
    def _distance_var(self, id_1, id_2):
        raise NotImplementedError
    
    def bound_distance(self, object_vars, lb=None, ub=None):
        id_1, id_2 = object_vars.keys()
        if (var := self._distance_var(id_1, id_2)) is None:
            self._setup_distance_calculation(object_vars)
            var = self._distance_var(id_1, id_2)
        
        return apply_bounds(var, lb=lb, ub=ub)

class EMD(Distance):
    def __init__(self, id_mapping_func, ordered_metric_ids, constraint_log):
        self.name = f'EMD{id(self)}'
        self._distance_var_collection = {}
        self._difference_var_collection = {}
        self.ordered_metric_ids = ordered_metric_ids
        self.id_mapping_func = id_mapping_func
        self.constraint_log = constraint_log

    def _difference_var(self, id_1, id_2, var_id=None):
        vars = self._difference_var_collection.get(str((id_1, id_2)))
        if vars is None: 
            vars = self._difference_var_collection.get(str((id_2, id_1)))
        if vars is not None:
            if var_id is None: return vars
            if vars.get(var_id) is not None: return vars[var_id]
        return None
    
    def _new_difference_var(self, id_1, id_2, var_id):
        var = Real(f'{self.name}_{str((id_1, id_2))}_{var_id}_EMDdiff')
        if not (vars := (self._difference_var_collection.get(str((id_1, id_2))) or \
            self._difference_var_collection.get(str((id_2, id_1))))):
            self._difference_var_collection[str((id_1, id_2))] = {}
            vars = self._difference_var_collection[str((id_1, id_2))]
        vars[var_id] = var
        return var
    
    def _distance_var(self, id_1, id_2):
        var = self._distance_var_collection.get(str((id_1, id_2)))
        if var is None:
            var = self._distance_var_collection.get(str((id_2, id_1)))
        return var
    
    def _new_distance_var(self, id_1, id_2):
        var = Real(f'{self.name}_{str((id_1, id_2))}_EMDdistance')
        self._distance_var_collection[str((id_1, id_2))] = var
        return var

    def _setup_distance_calculation(self, object_vars):
        """ 
        objects_vars: {object_id: {var_id: var}} 
        expects exactly 2 objects
        """
        assert len(object_vars) == 2

        id_1, id_2 = object_vars.keys()
        object_vars[id_1] = {self.id_mapping_func(var_id): var for var_id, var in object_vars[id_1].items()}
        object_vars[id_2] = {self.id_mapping_func(var_id): var for var_id, var in object_vars[id_2].items()}

        assert all([k in self.ordered_metric_ids for obj in object_vars.values() for k in obj.keys()])

        ob_1, ob_2 = object_vars[id_1], object_vars[id_2]
        prev_m_id = None
        for m_id in self.ordered_metric_ids:
            if not (var := self._difference_var(id_1, id_2, m_id)):
                var = self._new_difference_var(id_1, id_2, m_id)

            prev_difference = 0 if prev_m_id is None else self._difference_var(id_1, id_2, prev_m_id)
            cons = (var == Sum(ob_1.get(m_id, 0), -1*ob_2.get(m_id, 0), prev_difference))
            self.constraint_log.append(cons)
            prev_m_id = m_id

        if not (var := self._distance_var(id_1, id_2)):
            var = self._new_distance_var(id_1, id_2)
        
        cons = var == Sum([Abs(self._difference_var(id_1, id_2, m_id)) for m_id in self.ordered_metric_ids])
        self.constraint_log.append(cons)
