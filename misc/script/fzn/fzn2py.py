
# library of flatzinc predicates translated into numberjack constraints

MAXCOEF = 2147483647

class Var:
    def __init__(self, index):
        self.ind = index
    
def Variable(lb, ub, name):
    return Var(model.AddVariable(name, range(lb, ub+1)))
    
def VariableInDomain(dom, name):
    lb = min(dom)
    ub = max(dom)
    x = Var(model.AddVariable(name, range(lb, ub+1)))
    set_in(x, dom)
    return x
    
def Boolean():
    return Variable(0, 1, 'BOOL__' + str(model.GetNbVars()) + '__')

def VarArrayBoolean(nb, name):
    l = []
    for i in range(nb):
        l.append(Variable(0, 1, name + '_' + str(i) + '_'))
    return l

def VarArray(nb, lb, ub, name):
    l = []
    for i in range(nb):
        l.append(Variable(lb, ub, name + '_' + str(i) + '_'))
    return l

Constants = dict() # warning: fzn2py automatically replaces {} by []
def Constant(v):
    global Constants
    if type(v) is Var:
        return v
    if v in Constants:
        return Constants[v]
    else:
        Constants[v] = Variable(v, v, 'CONST__' + str(v) + '__')
        return Constants[v]

def ConstantNewVariable(v):
    if type(v) is Var:
        return v
    else:
        return Variable(v, v, 'CONST__' + str(v) + '__' + str(model.GetNbVars()) + '__')

def scope(s):
    if type(s) is int:
        return [Constant(s).ind]
    elif type(s) is Var:
        return [s.ind]
    else:
        return [Constant(x).ind if type(x) is int else x.ind for x in s]

def scopeWithDuplicateConstantVariables(s):
    if type(s) is int:
        return [ConstantNewVariable(s).ind]
    elif type(s) is Var:
        return [s.ind]
    else:
        return [ConstantNewVariable(x).ind if type(x) is int else x.ind for x in s]

def get_values(assignment, vars):
    return [assignment[e.ind] for e in vars]
    
def array_bool_and(x,y):
    if ((type(y) is int) and y != 0):
        int_lin_eq([1]*len(x), x, len(x)) # (Sum(x) == len(x))
    elif ((type(y) is int) and y == 0):
        int_lin_ne([1]*len(x), x, len(x)) # (Sum(x) != len(x))
    else:
        int_lin_eq_reif([1]*len(x), x, len(x), y) # (y == (Sum(x) == len(x)))

def array_bool_or(x,y):
    if ((type(y) is int) and y != 0):
        model.AddSumConstraint(scope(x),'>=',1)
    else:
        model.AddLinearConstraint([-MAXCOEF] + [1]*len(x), scope(y) + scope(x),'>=',-MAXCOEF+1)

def array_bool_xor(x):
    y = Variable(0, len(x)-1, 'XOR__' + str(model.GetNbVars()) + '__')
    int_lin_eq([1]*len(x), x, y) # y == Sum(x)
    int_mod(y,2,1) # ((Sum(x) % 2) == 1)

def array_int_element(x, y, z):
    int_le(1,x)
    int_le(x,len(y))
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    u = set([])  
    for e in y:
        u = u | set([e] if type(e) is int else model.Domain(e.ind))
    set_in(z, u)        
    for i, e in enumerate(y):
        if type(e) is int:
            e = Constant(e)
        sizee = model.GetDomainInitSize(e.ind)
        costs = [model.Top]*sizex*sizee*sizez
        for zval in model.Domain(z.ind):
            for xval in model.Domain(x.ind):
                for eval_ in model.Domain(e.ind):
                    if (xval != i+1) or (zval == eval_):
                       costs[model.GetValueIndex(z.ind, zval)*sizex*sizee + model.GetValueIndex(x.ind, xval)*sizee + model.GetValueIndex(e.ind, eval_)] = 0
        model.AddFunction(scope([z, x, e]), costs)
    # [(x >= 1), (x <= len(y)), set_in(z, u)] + [((z == (Variable(e,e,str(e)) if type(e) is int else e)) | (x != i+1)) for i, e in enumerate(y)]

def array_var_int_element(x,y,z):
    array_int_element(x,y,z)

def array_bool_element(x,y,z):
    array_int_element(x,y,z)

def array_var_bool_element(x,y,z):
    array_var_int_element(x,y,z)

def bool2int(x, y):
    int_eq(x,y) # (x == y)

def bool_and(x, y, z):
    if (type(z) is int) and (z != 0):
        int_lin_le([-1,-1], [x,y], -2)
    else:
        int_lin_le_reif([-1,-1], [x,y], -2, z)
    # (And(x, y) if ((type(z) is int) and (z != 0)) else (z == And(x, y)))

def bool_clause(x, y):
    int_lin_le([-1]*len(x) + [1]*len(y), x + y, -1 + len(y))

def bool_le(x, y):
    int_le(x, y) # ((x == 0) | (y != 0))

def bool_le_reif(x, y, z):
    int_le_reif(x, y, z) # [((x != 0) | (z != 0)), ((y != 0) | (z != 0)), ((x == 0) | (y != 0) | (z == 0))]

def bool_lt(x, y):
    int_eq(x, 0)
    int_neq(y, 0) # [(x == 0), (y != 0)]

def bool_lt_reif(x, y, z):
    int_lt_reif(x, y, z) #  [((x == 0) | (z == 0)), ((y != 0) | (z == 0)), ((x != 0) | (y == 0) | (z != 0))]

def bool_not(x, y):
    int_neq(x,y) # [((x == 0) | (y == 0)), ((x != 0) | (y != 0))]

def bool_or(x, y, z):
    int_lin_le_reif([-1,-1], [x,y], -1, z) # (z == (x | y ))

def bool_xor(x, y, z):
    int_neq_reif(x, y, z) # (z == (x != y))

def int_eq(x,y):
    model.AddLinearConstraint([1,-1], scope([x,y]), '==', 0) # (x == y)

def int_eq_reif(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if bool(zval) == (xval == yval):
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) # [((x != y) | (z != 0)), ((x == y) | (z == 0))]

def bool_eq(x, y):
    int_eq(x,y)

def bool_eq_reif(x, y, z):
    int_eq_reif(x, y, z)

def int_le(x,y):
    model.AddLinearConstraint([1,-1], scope([x,y]), '<=', 0) # (x <= y)

def int_le_reif(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if bool(zval) == (xval <= yval):
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) #  [(z == (x <= y))]

def int_lt(x,y):
    model.AddLinearConstraint([1,-1], scope([x,y]), '<', 0) # (x < y)

def int_lt_reif(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if bool(zval) == (xval < yval):
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) #  [(z == (x < y))]

def int_ne(x,y):
    #neq = Boolean()
    #model.AddLinearConstraint([MAXCOEF,1,-1], scope([neq,x,y]), '<', MAXCOEF)
    #model.AddLinearConstraint([MAXCOEF,1,-1], scope([neq,x,y]), '>', 0) # (x != y)
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    costs = [model.Top]*sizex*sizey
    for xval in model.Domain(x.ind):
        for yval in model.Domain(y.ind):
            if xval != yval:
                costs[model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([x, y]), costs) #  [(x != y)]

def int_ne_reif(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if bool(zval) == (xval != yval):
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) #  [(z == (x != y))]

def int_lin_eq(coef,vars,res):
    if type(res) is int:
        model.AddLinearConstraint(coef, scope(vars), '==', res)
    else:
        model.AddLinearConstraint([-1] + coef, scope(res) + scope(vars), '==', 0) # (res == Sum(vars,coef))

def bool_lin_eq(coef,vars,res):
    int_lin_eq(coef,vars,res)

def int_lin_eq_reif(coef,vars,res,z):
    if type(z) is int:
        if z == 0:
            neq = Boolean()
            if type(res) is int:
                model.AddLinearConstraint([MAXCOEF] + coef, scope([neq] + vars), '<', res + MAXCOEF)
                model.AddLinearConstraint([MAXCOEF] + coef, scope([neq] + vars), '>', res)
            else:
                model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([neq, res] + vars), '<', MAXCOEF)
                model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([neq, res] + vars), '>', 0)
        else: # z!=0
            if type(res) is int:
                model.AddLinearConstraint(coef, scope(vars), '>=', res)
                model.AddLinearConstraint(coef, scope(vars), '<=', res)
            else:
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '>=', 0)
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '<=', 0)        
    else:
        neq = Boolean()
        if type(res) is int:
            model.AddLinearConstraint([-MAXCOEF, MAXCOEF] + coef, scope([z, neq] + vars), '<', res + MAXCOEF)
            model.AddLinearConstraint([MAXCOEF, MAXCOEF] + coef, scope([z, neq] + vars), '>', res)
            model.AddLinearConstraint([-MAXCOEF] + coef, scope([z] + vars), '>=', res - MAXCOEF)
            model.AddLinearConstraint([MAXCOEF] + coef, scope([z] + vars), '<=', res + MAXCOEF)
        else:
            model.AddLinearConstraint([-MAXCOEF, MAXCOEF, -1] + coef, scope([z, neq, res] + vars), '<', MAXCOEF)
            model.AddLinearConstraint([MAXCOEF, MAXCOEF, -1] + coef, scope([z, neq, res] + vars), '>', 0)
            model.AddLinearConstraint([-MAXCOEF, -1] + coef, scope([z, res] + vars), '>=', -MAXCOEF)
            model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([z, res] + vars), '<=', MAXCOEF)
    # (z == (res == Sum(vars, coef)))

def int_lin_le(coef,vars,res):
    if type(res) is int:
        model.AddLinearConstraint(coef, scope(vars), '<=', res)
    else:
        model.AddLinearConstraint([-1] + coef, scope(res) + scope(vars), '<=', 0) # (res >= Sum(vars,coef))

def bool_lin_le(coef,vars,res):
    int_lin_le(coef,vars,res)

def int_lin_le_reif(coef,vars,res,z):
    if type(z) is int:
        if z == 0:
            if type(res) is int:
                model.AddLinearConstraint(coef, scope(vars), '>', res)
            else:
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '>', 0)
        else: # z!=0
            if type(res) is int:
                model.AddLinearConstraint(coef, scope(vars), '<=', res)
            else:
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '<=', 0)
    else:
        if type(res) is int:
            model.AddLinearConstraint([MAXCOEF] + coef, scope([z] + vars), '<=', res + MAXCOEF)
            model.AddLinearConstraint([MAXCOEF] + coef, scope([z] + vars), '>', res)
        else:
            model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([z, res] + vars), '<=', MAXCOEF)
            model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([z, res] + vars), '>', 0)
    # (z == (res >= Sum(vars,coef)))

def int_lin_lt(coef,vars,res):
    if type(res) is int:
        model.AddLinearConstraint(coef, scope(vars), '<', res)
    else:
        model.AddLinearConstraint([-1] + coef, scope(res) + scope(vars), '<', 0) # (res > Sum(vars,coef))

def int_lin_lt_reif(coef,vars,res,z):
    if type(z) is int:
        if z == 0:
            if type(res) is int:
                model.AddLinearConstraint(coef, scope(vars), '>=', res)
            else:
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '>=', 0)
        else: # z!=0
            if type(res) is int:
                model.AddLinearConstraint(coef, scope(vars), '<', res)
            else:
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '<', 0)
    else:
        if type(res) is int:
            model.AddLinearConstraint([MAXCOEF] + coef, scope([z] + vars), '<', res + MAXCOEF)
            model.AddLinearConstraint([MAXCOEF] + coef, scope([z] + vars), '>=', res)
        else:
            model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([z, res] + vars), '<', MAXCOEF)
            model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([z, res] + vars), '>=', 0)
    # (z == (res > Sum(vars,coef)))

def int_lin_ne(coef,vars,res):
    neq = Boolean()
    model.AddLinearConstraint([MAXCOEF] + coef, scope([neq] + vars), '<', res + MAXCOEF)
    model.AddLinearConstraint([MAXCOEF] + coef, scope([neq] + vars), '>', res) # (res != Sum(vars,coef))

def int_lin_ne_reif(coef,vars,res,z):
    if type(z) is int:
        if z != 0:
            neq = Boolean()
            if type(res) is int:
                model.AddLinearConstraint([MAXCOEF] + coef, scope([neq] + vars), '<', res + MAXCOEF)
                model.AddLinearConstraint([MAXCOEF] + coef, scope([neq] + vars), '>', res)
            else:
                model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([neq, res] + vars), '<', MAXCOEF)
                model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([neq, res] + vars), '>', 0)
        else: # z==0
            if type(res) is int:
                model.AddLinearConstraint(coef, scope(vars), '>=', res)
                model.AddLinearConstraint(coef, scope(vars), '<=', res)
            else:
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '>=', 0)
                model.AddLinearConstraint([-1] + coef, scope([res] + vars), '<=', 0)        
    else:
        neq = Boolean()
        if type(res) is int:
            model.AddLinearConstraint([MAXCOEF, MAXCOEF] + coef, scope([z, neq] + vars), '<', res + 2*MAXCOEF)
            model.AddLinearConstraint([-MAXCOEF, MAXCOEF] + coef, scope([z, neq] + vars), '>', res - MAXCOEF)
            model.AddLinearConstraint([MAXCOEF] + coef, scope([z] + vars), '>=', res)
            model.AddLinearConstraint([-MAXCOEF] + coef, scope([z] + vars), '<=', res)
        else:
            model.AddLinearConstraint([MAXCOEF, MAXCOEF, -1] + coef, scope([z, neq, res] + vars), '<', 2*MAXCOEF)
            model.AddLinearConstraint([-MAXCOEF, MAXCOEF, -1] + coef, scope([z, neq, res] + vars), '>', -MAXCOEF)
            model.AddLinearConstraint([MAXCOEF, -1] + coef, scope([z, res] + vars), '>=', 0)
            model.AddLinearConstraint([-MAXCOEF, -1] + coef, scope([z, res] + vars), '<=', 0)
    # (z == (res != Sum(vars,coef)))

def int_abs(x,y):
    if type(x) is Var:
        sizex = model.GetDomainInitSize(x.ind)
        if type(y) is int:
            y = Constant(y)
        sizey = model.GetDomainInitSize(y.ind)
        costs = [model.Top]*sizex*sizey
        for yval in model.Domain(y.ind):
            for xval in model.Domain(x.ind):
                if yval == abs(xval):
                    costs[model.GetValueIndex(y.ind, yval)*sizex + model.GetValueIndex(x.ind, xval)] = 0
        model.AddFunction(scope([y, x]), costs)
    else:
        int_eq(y, abs(x)) # (y == Abs(x))

def int_div(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if yval != 0 and zval == xval // yval:
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) # (z == (x / y))

def int_min(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if zval == min(xval, yval):
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) # (z == Min([x, y]))

def int_max(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if zval == max(xval, yval):
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) # (z == Max([x, y]))

def int_mod(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if zval == xval % yval:
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) # (z == (x % y))

def int_plus(x,y,z):
    model.AddLinearConstraint([1,-1,-1], scope([z,x,y]), '==', 0) # (z == (x + y))

def int_times(x,y,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(y) is int:
        y = Constant(y)
    sizey = model.GetDomainInitSize(y.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizey*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            for yval in model.Domain(y.ind):
                if zval == xval * yval:
                    costs[model.GetValueIndex(z.ind, zval)*sizex*sizey + model.GetValueIndex(x.ind, xval)*sizey + model.GetValueIndex(y.ind, yval)] = 0
    model.AddFunction(scope([z, x, y]), costs) # (z == (x * y))

def set_in(x,dom):
    model.AddCompactFunction(scope(x), model.Top, [[v] for v in dom], [0]*len(dom)) # x in dom

def set_in_reif(x,dom,z):
    if type(x) is int:
        x = Constant(x)
    sizex = model.GetDomainInitSize(x.ind)
    if type(z) is int:
        z = Constant(z)
    sizez = model.GetDomainInitSize(z.ind)
    costs = [model.Top]*sizex*sizez
    for zval in model.Domain(z.ind):
        for xval in model.Domain(x.ind):
            if bool(zval) == (xval in dom):
                costs[model.GetValueIndex(z.ind, zval)*sizex + model.GetValueIndex(x.ind, xval)] = 0
    model.AddFunction(scope([z, x]), costs) # (z == Disjunction([(x == v) for v in dom]))
    
def Minimize(x):
    if type(x) is Var:
        model.AddFunction(scope(x), [model.GetValue(x.ind, index) for index in range(model.GetDomainInitSize(x.ind))])
    
def Maximize(x):
    if type(x) is Var:
        model.AddFunction(scope(x), [-model.GetValue(x.ind, index) for index in range(model.GetDomainInitSize(x.ind))])

#-----------------------------------------
# Specific global constraints for toulbar2
#-----------------------------------------

def fzn_all_different_int(x):
    if len(x) >= 2:  # Some models specified alldiff on 1 variable
        model.AddAllDifferent(scope(x))  # [Variable(e,e,str(e)) if type(e) is int else e for e in x])
        
def fzn_global_cardinality(x, values, counts):
    assert(len(values) == len(counts))
    l = []
    for j in range(len(values)):
        if type(counts[j]) is Var:
            if len(model.Domain(counts[j].ind)) > 1:
                model.AddGeneralizedLinearConstraint([[scope(x[i])[0], values[j], 1] for i in range(len(x))] + [[counts[j].ind, v, -v] for v in model.Domain(counts[j].ind)], '==', 0)
            else:
                l.append((values[j], model.Domain(counts[j].ind)[0], model.Domain(counts[j].ind)[0]))
        else:
            l.append((values[j], counts[j], counts[j]))
    if len(l) > 0:
        model.AddGlobalCardinalityConstraint(scope(x), l)

def fzn_global_cardinality_closed(x, values, counts):
    assert(len(values) == len(counts))
    for i in range(len(x)):
        set_in(x[i], values)
    fzn_global_cardinality(x, values, counts)
        
def fzn_global_cardinality_low_up(x, values, lb, ub):
    assert(len(values) == len(lb))
    assert(len(values) == len(ub))
    l = []
    for j in range(len(values)):
        if (type(lb[j]) is Var) or (type(ub[j]) is Var):
            mylb = 0
            myub = len(x)
            if type(lb[j]) is Var:
                if len(model.Domain(lb[j].ind)) > 1:
                    model.AddGeneralizedLinearConstraint([[scope(x[i])[0], values[j], 1] for i in range(len(x))] + [[lb[j].ind, v, -v] for v in model.Domain(counts[j].ind)], '>=', 0)
                else:
                    mylb = model.Domain(lb[j].ind)[0]
            else:
                mylb = lb[j]
            if type(ub[j]) is Var:
                if len(model.Domain(ub[j].ind)) > 1:
                    model.AddGeneralizedLinearConstraint([[scope(x[i])[0], values[j], 1] for i in range(len(x))] + [[ub[j].ind, v, -v] for v in model.Domain(counts[j].ind)], '<=', 0)
                else:
                    myub = model.Domain(ub[j].ind)[0]
            else:
                myub = ub[j]
            if mylb > 0 or myub < len(x):
                l.append((values[j], mylb, myub))
        else:
            l.append((values[j], lb[j], ub[j]))
    if len(l) > 0:
        model.AddGlobalCardinalityConstraint(scope(x), l)
        
def fzn_global_cardinality_low_up_closed(x, values, lb, ub):
    assert(len(values) == len(lb))
    assert(len(values) == len(ub))
    for i in range(len(x)):
        set_in(x[i], values)
    fzn_global_cardinality_low_up(x, values, lb, ub)

def fzn_table_int(x,t):
    n = len(x)
    assert(len(t) % n == 0)
    nbtuples = len(t) // n
    model.AddCompactFunction(scope(x), model.Top, [[t[i+j] for j in range(n)] for i in range(0,nbtuples,n)], [0]*nbtuples) # x in t
    
def fzn_table_bool(x,t):
    fzn_table_int(x, t)

def fzn_regular(X, Q, S, d, q0, F):
    parameters = [Q+1, 1, (q0, 0), len(F)]
    for e in F:
        parameters.append((e,0))
    assert(len(d) == Q * S)
    parameters.append(Q * S)
    for i in range(Q):
        for j in range(S):
            parameters.append((i+1, j+1, d[i*S+j], 0))
    model.AddGlobalFunction(scopeWithDuplicateConstantVariables(X), 'wregular', parameters)

def fzn_sregular(X, Q, S, d, q0, F):
    parameters = ['var', 1, Q+1, 1, q0, len(F)]
    for e in F:
        parameters.append(e)
    assert(len(d) == Q * S)
    parameters.append(Q * S)
    for i in range(Q):
        for j in range(S):
            parameters.append(i+1)
            parameters.append(j+1)
            parameters.append(d[i*S+j])
    model.AddGlobalFunction(scopeWithDuplicateConstantVariables(X), 'sregular', parameters)

