# ODE solving procedure

import sympy as sp
import ODEfunction as odefx

# Define ode Class with appropriate classifications
class Classified_ODE():
    
    def __init__(self, ode, ind_var, dep_var, order, degree, is_forced, is_linear):
        self.ode = ode
        self.ind_var = ind_var
        self.dep_var = dep_var
        self.order = order
        self.degree = degree
        self.is_forced = is_forced
        self.is_linear = is_linear
    
    def printClassification(self):
        print("This ode has order of: ", self.order)
        print("This ode has degree: ", self.degree)
        print("The independent variable of this ode is ", str(ind_var))
        print("The dependent variable of this ode is ", str(dep_var))
        if self.is_forced:
            print("This ode is forced.")
        else:
            print("This ode is homogeneous.")
        if self.is_linear:
            print("This ode is linear.")
        else:
            print("This ode is not linear.")

x, y = sp.symbols('x, y')

# input the ODE
ode_input = input("Enter the ODE as EQ(lhs,rhs) with x as independent and y as dependent variable: ")

ode_rhs = sp.rhs(ode_input)
ode_lhs = sp.lhs(ode_input)
ode = sp.Eq(sp.sympify(ode_lhs), sp.sympify(ode_rhs))

# Get order, degree, independent and dependent variable, isLinear, isForcing
'''
# Define symbolic variables
def extract_variables(ode):
    variables = sp.atoms(symbols)

    for var in variables:
        if isinstance(var, sp.Dummy):
            continue
        if var.is_symbol and var.diff():
            dep_var.append(var)
        else:
            ind_var.append(var)

    return ind_var, dep_var

ind_var, dep_var = extract_variables(ode)
'''


classifiedODE = Classified_ODE(ode, x, y, sp.ode_order(ode, y), sp.degree(ode.lhs - ode.rhs, x))

printClassification(ode)

if odefx.is_linear_ODE(classifiedODE, dep_var, ind_var):
    print("""This ODE is linear as it satissfies the following conditions: 
            All coefficients, a.#, are functions of independant variables ONLY
            All derivative terms are of degree one or zero
            There are no products of derivative terms
            Coefficient a.n <> 0 for any value of indep.variable in solution range
            Last term on LHS - must be the dep. variable NOT a function of the dep. variable""")
else:
    print("This ODE is not linear since one or more linear conditions are not satisfied.")

"""
    If of the form: 
    𝑎.𝑛(𝑥)*(𝑑^𝑛)𝑦/𝑑𝑥^𝑛 + 𝑎.(𝑛−1)(𝑥)*𝑑^(𝑛−1)𝑦/𝑑𝑥^(𝑛−1) + ⋯ + 𝑎.1(𝑥)*𝑑𝑦/𝑑𝑥 + 𝑎.0(𝑥)𝑦 = F(x)
    Homogeneous ODE if F(x) = 0
    Else forced/non-homogeneous
"""
IS_HOMOGENEOUS_ODE = odefx.is_homogeneous(classifiedODE, x, y)

def solve_ode(classified_ODE, y, x):

    if ORDER == 1 and DEGREE == 1 :
        # If it is linear, solve using first order formula
        if Classified_ODE().is_linear:
            solution = odefx.solve_linear_ODE(ode=classified_ODE, is_linear=True, y=y, x=x)
            return solution
        
        # Check if the ODE is separable
        is_separable = (ode.rhs - ode.lhs / y).simplify() == 0
        if is_separable:
            print(""" Since this ode is seperable solve:
                    dy/dx = f(x)*h(y) => f(x)dx = h(y)dy
                    1. Seperate both sides 
                    2. Integrate both sides
                    3. Rearrange""")
            print("Solution: ", odefx.solve_seperable(ode, y, x))
            return
        
        # Homogeneous substitution
        elif IS_HOMOGENEOUS_ODE:
            # Perform the homogeneous substitution: y = x * v
            v = symbols('v')(x)
            homogeneous_ode = ode.subs(y, x*v)
        
            # Check if homogeneous substitution is seperable
            is_seperable = (homogeneous_ode.rhs - homogeneous_ode.lhs / (x*v)).simplify() == 0
            if is_seperable:
                print("""If dy/dx = F(x/y) or F(y/x) or y(x.0) = y.0
                        Change variable: 
                            v = (y/x) => y = vx
                            dy/dx = x*(dv/dx) + v => x*(dv/dx) + v = F(v)
                        Solve like seperable ODE:
                            1. Seperate ∫dx/x = ∫dv/F(v)-v
                            2. Integrate ln|x| = ln|F(v)-v| + c.1
                            3. Replace v = y/x
                            4. => Solve for y(x)""")
                print("Solution: ", odefx.solve_homogeneous_sub(ode, y, x))
                return
            
        # Check if the ODE is exact
        is_exact = odefx.check_exact(ode)
        if is_exact:
            print("For dy/dx = M(x,y)/-N(x/y) => M(x,y)dx + N(x,y)dy = 0)")
            print("If ∂M/∂y = ∂/∂y(∂U/∂x) = ∂N/∂x = ∂/∂y(∂U/∂y) => Exact")
            print("This ODE is exact.")
            print("Solution: ", odefx.solve_exact_ODE(ode, y, x))
            return
        
        # ODE made exact using integration factor
        integration_factored_ode = odefx.integration_factor_ODE(ode, y, x)
        if odefx.check_exact(integration_factored_ode):
            print("""
                ODE made exact using Integration Factor
                ∂M/∂y = ∂/∂y(∂U/∂x) <> ∂N/∂x = ∂/∂y(∂U/∂y)
                u(x,y) factor multiplies through ODE to make it exact:
                ∂(uM)/∂y = ∂(uN)/∂x
                Four methods to find integration factor
                1. Factor M&N into functions of x,y (M->f(x)*f(y) N->g(x)*g(y)) 
                    seperatly and multiply ODE by 1/(f(x)*g(y)) or 1/(f(y)*g(x)) to make equatioin seperable
                    - Directly integrate
                2. Find if (1/N)(∂M/∂y-∂N/∂x)=f(x) only, then u(x)=e^∫f(x)dx
                3. Find if (1/M)(∂N/∂x-∂M/∂y)=f(y) only, then u(y)=e^∫f(y)dy
                4. Assume form u(x,y) = (x^a)*(y^b), multiply through ODE and solve for a, b to make equation exact
                """)
            print("Solution: ", odefx.solve_exact_int_factor(ode, y, x))
            return
            
        if odefx.is_linear_ode(ode, y, x):
            print("""The procedure for solving a linear ODE is as follows:
                If the ODE is in the linear form: dy/dx + p(x)*y = g(x)
                Apply formula: y(x) = e^(−∫𝑝(x)dx) * [∫𝑔(x)*e^(∫𝑝(x)dx) dx + c]""")
            print("Solution: ", odefx.solve_linear_ODE(ode, is_linear, y, x))
            return
        
        # ODE made linear using variable substitution ie. invert variables to check linearity
        ode_invert_var = odefx.invert_dep_var(ode, y, x)
        if_linear = odefx.is_linear_ODE(ode_invert_var)
        if if_linear:
            print("Solution: ", odefx.solve_linear_ODE(ode_invert_var))
            return
        
        # Sub function of dependent variable to make ode linear
        ode_function_dep_var = odefx.function_dep_var(ode, y, x)
        if_linear = odefx.is_linear_ODE(ode_function_dep_var)
        if is_linear:
            print("Solution: ", odefx.solve_linear_ODE(ode_function_dep_var))
            return
        
        # If of bernoulli form (dy/dx + p(x)y = q(x)y^n), sub w = y^(1-n) and solve
        if ode == odefx.bernoulli_ODE(ode, y, x):
            print("Solution: ", odefx.solve_bernoulli_ODE(ode, y, x))
            return
    
    if ORDER >= 2:    
    # Reduction of Order
        
        # Missing term: sub v = dy/dx to reduce order, then solve for v, then y
        # dw/dx = dw/dy * dy/dx = w * dw/dy
        if ode.rhs == 0:
            print("Solution: ", odefx.RoO_homogeneous_ODE(ode, y, x))
            
        # Factor operators: factor operators, sub v = (D+b)y, solve for v, solve for y
        elif ode.coeff(y.diff(x)) == ode.coeff(y):
            print("Solution: ", odefx.RoO_Factor_Operators(ode, y, x))
            
        # D'Albert method: know one solution y.1(x), assume second solution is y.2(x) = v(x)*y.1(x)
        # Sub in ode and get ode in v(x), Solve for v(x) then solution id Y(x) = y.1(x) + y.2(x)
        else:
            print("Solution: ", odefx.RoO_DAlemberts_method(ode, y, x))  
    
    
    if ORDER >= 2 & ode.constant_coeffs() == True & ode.is_linear() == True:      
    # Method of Undetermined Coefficients
        # Find Complementary Solution
            # Solve ODE as homogeneous using homogeneous techniques that apply
                # Real distinct roots, solution is two distinct terms, ex: c_1*e^(m_1*x) + c_2*e^(m_2*x)
                # Real equal roots, solution same term but one multiplied by the independent variable, ex: (c_1 + t*c_2)e^mx
                # Complex Conjugate roots, real part (m) represented with e^mx, imaginary part (n)
                #   represented with sin and cos, ex: e^(mx)*(c_1*cos(nx) + c_2*sin(nx))
                # Complex Repeated Roots, solution is terms twice but one multiplied by independent variable, ex:
                #   e^(mx)*(c_1*cos(nx) + c_2*sin(nx)) + t*e^(mx)*(c_1*cos(nx) + c_2*sin(nx))
        # Find Particular Solution
            # Assume particular solution to the ode based on the forcing function
            # Sub function and its deribatives into the ode and solve for undetermined coefficients
        # Add particular and complementary solution toguether for general solution to the ode
        # Solve integral constants with intital conditions if required
        print("Solution: ", odefx.method_of_undetermied_coefficients(ode, y, x))
        
    if ORDER >= 2 & ode.constant_coeffs() == True & ode.is_linear() == True:    
    # Laplace Transformations
        # Laplace Transform each term in the ode (from time to frequency domain)
        # Algebraically solve for Y(s)
        # Inverse Laplace transform each term in the ode
        # Solution: y(t) = f(t)
        print("Solution: ", odefx.laplace_transformation(ode, y, x))
        
    if ode_count >= 2
    # Elimination: can do but takes a lot of time and is inefficient
    # Simultateous Sets of ODEs: Factoring by Operators & MUC
        # Write the set of ODEs using factor operators and in the form
            # Ax + By = f(t)
            # Cx + Dy = g(t)
        # Set up a matrix expression of the ODE
        # Apply Cramer's Rule
            # (AD-BC)x=Df-Bg
            # (AD-BC)y=Ag-Cf
        # Solve ODEs to find the solutions to the set for x(t) and y(t)
        print("Solution: ", odefx.SSO_Factor Operators_MUC(ode1, ode2, y, x))    
        
    # Laplace Transformations for sets of ODEs
        # Laplace transform both equations
        # Isolate X(s) and Y(s) in both equations
        # Solve like factor operators and muc (matrix -> Cramer's Rule -> Isolate)
        # Laplace Transform both solutions for x(t) and y(t)
        print("Solution: ", odefx.SSO_Laplace_transformation(ode1, ode2, y, x))