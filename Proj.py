import numpy as np
import scipy as sp
import scipy.interpolate as spint
import matplotlib.pyplot as plt

# NON LAB RELATED FUNCTION DEFINITIONS
# Truncater
def trunc(num, dec = 2):
    """
    Note defined lately in the process therefore unconsistent usage
    """
    rest = round(num, dec)%1
    rest = float(str(rest)[0:2+dec])
    truncated = int(num)+rest
    return truncated

#################################
# LAB #
#################################


steps = 100
x = np.arange(steps+1)*(np.pi)/steps

# Init variables

# Lagom stor hat för kräftskivan
R1 = 7.8
a1 = 5.5
b1 = 11

set1 = (R1, a1, b1)

# Stjärngossstrut
R2 = 8.8
a2 = 12
b2 = 36

set2 = (R2, a2, b2)
    
# dAlpha / du
"Alpha är lutningen för kurvan i punkten (epsilon, eta)"
def dAlpha(u, R, a, b):
    upper = R - a * np.cos(u)
    lower = np.sqrt(np.square(b)+np.square(R - a * np.cos(u)))
    dA = upper / lower
    return dA

# dEpsilon / du
"Epsilon är kurvarns x-värden då parametern u löper från 0 till 2pi"
def dEpsilon(R, alpha):
    dE = R * np.cos(alpha)
    return dE

# dEta / du
"Eta är kurvans y-värden då parametern u löper från 0 till 2pi"
def dEta(R, alpha):
    dn = R * np.sin(alpha)
    return dn

# Runge-Kutta
"Källa: http://www.nada.kth.se/kurser/su/NUMFYS/HT05/l3.pdf , p.5"

"""
Noggranhetsordning 4 - trunkeringsfel approx. c*h^4

"""


"""
xprim och yprim
xprim -> dEdu = R cos a
yprim -> dndu = R sin a

da = ändring av x och y


"""
def klippArket(R,a,b):
    P = np.sqrt(np.square(b) + np.square(R-a))
    return 0, P

###ALPHA IS WRONG
def RungeKutta(set_value = set1, N = 40, k = 1): #k = 3, 12
    R,a,b = set_value
    h = 2*np.pi / float(N) # 0 < x < 2pi in N steps

    H = np.empty(k)

    #yett = np.empty(k)
    alpha_ett = np.empty(k)
    epsilon_ett = np.empty(k)
    eta_ett = np.empty(k)

    # H = []
    # yett = []

    "Printout values"
    #If we use np.empty first value is not appended and thus graphs bugg out
    u_print = np.zeros(N*(2**(k-1))+1) #We append N*k value 0...N*k for last iteration
    alpha_print = np.zeros(N*(2**(k-1))+1)
    epsilon_print = np.zeros(N*(2**(k-1))+1)
    eta_print = np.zeros(N*(2**(k-1))+1)

    for i in range(1,k+1):

        #x = 0.0
        #y = 1.0

        # Startvärden i punkten A:
        u = 0
        alpha = 0 # diffen
        epsilon = 0 # x-coord
        eta = 0 # y-coord

        for n in range(1,N+1):
            # k1 = h * f(x,y)
            # k2 = h * f(x + (h/2), y + k1/2)
            # k3 = h * f(x + (h/2), y +k2/2)
            # k4 = h * f(x + h, y + k3);
            # y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
            # x = x + h

            # Appending values for printout
            #if i == k:
            #    #print(N)
            #    u_print[n-1] = u
            #    alpha_print[n-1] = alpha
            #    epsilon_print[n-1] = epsilon
            #    eta_print[n-1] = eta
                

            """"""
            k1_alpha = h * dAlpha(u, R, a, b)
            k2_alpha = h * dAlpha(u + h/2, R, a, b)
            k3_alpha = h * dAlpha(u + h/2, R, a, b)
            k4_alpha = h * dAlpha(u + h, R, a, b)
            
            k1_epsilon = h * dEpsilon(R, alpha)
            k2_epsilon = h * dEpsilon(R, alpha + (h/2))
            k3_epsilon = h * dEpsilon(R, alpha + (h/2))
            k4_epsilon = h * dEpsilon(R, alpha + h)
            
            k1_eta = h * dEta(R, alpha)
            k2_eta = h * dEta(R, alpha + (h/2))
            k3_eta = h * dEta(R, alpha + (h/2))
            k4_eta = h * dEta(R, alpha + h)

            u = u + h
            alpha = alpha + (k1_alpha + 2 * k2_alpha +
                             2 * k3_alpha + k4_alpha) / 6
            epsilon = epsilon + (k1_epsilon + 2 * k2_epsilon +
                                 2 * k3_epsilon + k4_epsilon) / 6
            eta = eta + (k1_eta + 2 * k2_eta +
                         2 * k3_eta + k4_eta) / 6
            # Appending values for print out
            if i == k:
                #print(N)
                u_print[n] = u
                alpha_print[n] = alpha
                epsilon_print[n] = epsilon
                eta_print[n] = eta
        H[i-1] = h
        #yett[i-1] = y
        alpha_ett[i-1] = alpha
        epsilon_ett[i-1] = epsilon
        eta_ett[i-1] = eta
        #H.append(h)
        #yett.append(k)
        N = 2 * N
        h = h / 2
        
        
    #delta_y = np.diff(yett)
    delta_alpha = np.diff(alpha_ett)
    delta_epsilon = np.diff(epsilon_ett)
    delta_eta = np.diff(eta_ett)

    #m = len(delta_y)
    m_alpha = len(delta_alpha)
    m_epsilon = len(delta_epsilon)
    m_eta = len(delta_eta)
    
    #kvot = delta_y[0:m-1] / delta_y[1:m]
    kvot_alpha = delta_alpha[0:m_alpha-1] / delta_epsilon[1:m_alpha]
    kvot_epsilon = delta_epsilon[0:m_epsilon-1] / delta_epsilon[1:m_epsilon]
    kvot_eta = delta_eta[0:m_eta-1] / delta_eta[1:m_eta]
    
    #return H, yett, delta_y, kvot
    return (H,
            alpha_ett, epsilon_ett, eta_ett,
            delta_alpha, delta_epsilon, delta_eta,
            kvot_alpha, kvot_epsilon, kvot_eta,
            u_print, alpha_print, epsilon_print, eta_print)

# Runge-Kutta plot
def RungeKutta_plot(s, n, k):
    H, alpha_ett, epsilon_ett, eta_ett, \
       delta_alpha, delta_epsilon, delta_eta, \
       kvot_alpha, kvot_epsilon, kvot_eta, \
       u_print, alpha_print, epsilon_print, eta_print = RungeKutta(set_value = s, N = n, k = k)

    # Punkten
    xp, yp = klippArket(R1, a1, b1)
    plt.figure()
    
    plt.subplot(221)
    plt.plot(epsilon_print, eta_print, xp, yp, 'r.' )
    plt.title('$\epsilon$ - $\eta$')

    plt.subplot(222)
    plt.plot(u_print, alpha_print)
    plt.title(r'$\alpha$')
    
    plt.subplot(223)
    plt.plot(u_print, epsilon_print)
    plt.title('$\epsilon$')

    plt.subplot(224)
    plt.plot(u_print, eta_print)
    plt.title('$\eta$')

    
    plt.suptitle('Runge-Kutta')
    plt.show()

def Hermite_Interpolation(s = set1, n = 8, n_step = 20):
    # Hermite interpolation
    """
    Here we view our
        x = u
        y = epsilon, eta
        
    in our formulas

    We generate new x_substitutes with "step" increments as our
    interpolation points.
    """
    # Create Runge-Kutta Values to interpolate from
    H, alpha_ett, epsilon_ett, eta_ett, \
           delta_alpha, delta_epsilon, delta_eta, \
           kvot_alpha, kvot_epsilon, kvot_eta, \
           u_print, alpha_print, epsilon_print, eta_print = RungeKutta(set_value = s, N = n)

    # Create interpolation points
    step = n_step
    x = np.arange(step+1)*(2*np.pi)/step

    # Interpolation functions
    f_alpha = spint.PchipInterpolator(u_print, alpha_print)
    f_epsilon = spint.PchipInterpolator(u_print, epsilon_print)
    f_eta = spint.PchipInterpolator(u_print, eta_print)

    # Mid point
    xp, yp = klippArket(R1, a1, b1)

    # Plotting
    plt.subplot(221)
    plt.plot(epsilon_print, eta_print, 'ro', f_epsilon(x), f_eta(x), '--', xp, yp, 'x')
    plt.title('$\epsilon$ - $\eta$')

    plt.subplot(222)
    plt.plot(u_print, alpha_print, 'ro', x, f_alpha(x), '--')
    plt.title(r'$\alpha$')
        
    plt.subplot(223)
    plt.plot(u_print, epsilon_print, 'ro', x, f_epsilon(x), '--')
    plt.title('$\epsilon$')

    plt.subplot(224)
    plt.plot(u_print, eta_print, 'ro', x, f_eta(x), '--')
    plt.title('$\eta$')

    plt.suptitle('Hermite interpolation')
    plt.show()

def Compare_RungeKutta_Hermite(s = set1, n = 4, n_step = 20):
    # Hermite interpolation
    """
    Here we view our
        x = u
        y = epsilon, eta
        
    in our formulas

    We generate new x_substitutes with "step" increments as our
    interpolation points.
    """
    # Create Runge-Kutta Values to interpolate from
    H, alpha_ett, epsilon_ett, eta_ett, \
           delta_alpha, delta_epsilon, delta_eta, \
           kvot_alpha, kvot_epsilon, kvot_eta, \
           u_print, alpha_print, epsilon_print, eta_print = RungeKutta(set_value = s, N = n)

    # Create interpolation points
    step = n_step
    x = np.arange(step+1)*(2*np.pi)/step

    # Interpolation functions
    f_alpha = spint.PchipInterpolator(u_print, alpha_print)
    f_epsilon = spint.PchipInterpolator(u_print, epsilon_print)
    f_eta = spint.PchipInterpolator(u_print, eta_print)

    # Mid point
    R_s, a_s, b_s = s
    xp, yp = klippArket(R_s, a_s, b_s)

    # Plotting
    plt.subplot(221)
    plt.plot(epsilon_print, eta_print, f_epsilon(x), f_eta(x), '--', xp, yp, 'x')
    plt.title('$\epsilon$ - $\eta$')

    plt.subplot(222)
    plt.plot(u_print, alpha_print, x, f_alpha(x), '--')
    plt.title(r'$\alpha$')
        
    plt.subplot(223)
    plt.plot(u_print, epsilon_print, x, f_epsilon(x), '--')
    plt.title('$\epsilon$')

    plt.subplot(224)
    plt.plot(u_print, eta_print, x, f_eta(x), '--')
    plt.title('$\eta$')

    plt.suptitle('Compare Runga-Kutta and Hermite')
    plt.show()


# Uppgift 1
##RungeKutta_plot(s = set1, n = 8, k = 1)

# Uppgift 2
##Hermite_Interpolation()
##Compare_RungeKutta_Hermite()

# Uppgift 3

"""
alpha(2*pi) = (pi, pi/2) 

R = radie
a = sida
b = höjd

Surface area = pi * R * sqrt(R^2 + b^2)
Volume = 1/3 * pi * R^2 * b

Rektangulär klippark
Surface area = s1 * s2

s1 = lång sida
s2 = kort sida

"""
import matplotlib.image as mpimg
img = mpimg.imread('cone_1.gif')
plt.subplot(311)
plt.imshow(img)
img = mpimg.imread('frustrum_1.gif')
plt.subplot(312)
plt.imshow(img)
img = mpimg.imread('conical_frustum_1.gif')
plt.subplot(313)
plt.imshow(img)
plt.show() #SHOOOOOOOOOOW IMAGES

# Kräfthatt (truncated cone)
R3 = 7.8
a3 = 5.5
set3 = (R3, a3)
"""
Beräkna b så att tillverkningskravet
alpha(2*pi) = pi
"""

# Lussestrut (cone)
R4 = 8.8
a4 = 12
set4 = (R4, a4)
"""
Beräkna b så att tillverkningskravet
alpha(2*pi) = pi/2
"""
# Create Runge-Kutta Values to interpolate from
H, alpha_ett, epsilon_ett, eta_ett, \
   delta_alpha, delta_epsilon, delta_eta, \
   kvot_alpha, kvot_epsilon, kvot_eta, \
   u_print, alpha_print, epsilon_print, eta_print = RungeKutta(set_value = set1, N = 4)

# Create interpolation points
step = 20
x = np.arange(step+1)*(2*np.pi)/step

start = 0
end = np.pi


########################

#
def alpha(R, a, b, u):
    """
    Inputting in Wolfram:
    (R - a cos(x))/sqrt(b^2 + (R-a cos(x))^2

    Gives us following alpha
    """
    alp = (a*np.square(b)*np.sin(u))/((np.square(R-a*np.cos(u))+np.square(b))**(3/2))
    return alp

#Find h
def h_cone(R, a):
    """
    We utilize pythagorian theorem
    sqrt(r^2 + h^2) = sqrt(a^2)

    We know a and r thus
    r^2 + h^2 = a^2
    h^2 = a^2 - r^2
    h = sqrt(a^2 - r^2)
    """

    h = np.sqrt(np.square(a) - np.square(R))
    return h

def h_trunCone(R,r,a):
    h = np.sqrt(np.square(a) - np.square(R-r))
    return h

b2 = h_cone(R4, a4)

def surfaceArea_cone(R, b):
    area = np.pi * R * np.sqrt(np.square(R) + np.square(b))
    return area

def trapets(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, b_set = 11):
#(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, b_set):
    """
    Method from
    http://www.csc.kth.se/utbildning/kth/kurser/DN1240/numi09/forelasningar/f6.pdf

    Excluded delta calculation
    Excluded kvot
    Excluded Ttak

    Maybe insert?

    Behöver trunkeringsfel?

    """

    """
    "TEST SET"
    a = 0
    b = 2*np.pi
    h = 0.5
    R_set = 7.8
    a_set = 5.5
    b_set = 10
    #
    """

    n = (b-a)/h # Antal trapets
    antal = 5; # Antal halveringar
        
    # def dAlpha(u, R, a, b):
    T0 = (dAlpha(a, R_set, a_set, b_set) + dAlpha(b, R_set, a_set, b_set))/2
    T = np.zeros(antal+1)
        
    T[0] = T0
    for i in range(1, antal+1):
        x = np.arange(a, b, h)
        T[i] = (np.sum(dAlpha(x, R_set, a_set, b_set)) - T0)*h
        n = 2*n # fördubblar antal trapetser
        h = h/2 # Steglängd halveras

    return T

def find_b(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, conv_crit = np.pi):
    """
    Own function for iterating to find b until convergence criteria is reached
    """
    i = 0
    b_start = 0
    step = 0.01 #Found for set1 by trial and error
    while True:
        T = trapets(a = a, b = b, h = h, R_set = R_set, a_set = a_set, b_set = b_start)[-1]
        #print(T)
        if (float(str(T)[0:4]) == float(str(conv_crit)[0:4])) or (round(T, 2) == float(str(conv_crit)[0:4])):
            break
        b_start += step
        i += 1 # iter
    return b_start, i

# Finding b for set 1 and 2
# b - höjd
# a - sida
# R - Radie

b_1, i_1 = find_b(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, conv_crit = np.pi)
b_2, i_2 = find_b(a = 0, b = 2*np.pi, h = 0.5, R_set = 8.8, a_set = 12, conv_crit = np.pi/2)
    
    
def surf_area(r,h):
    # R and b
    area = np.pi*r*np.sqrt(np.square(r) + np.square(h))
    return area

def rec_area(a,b):
    # Theoretical
    area = a*b
    return area


sa1 = surf_area(7.8, b_1)
sa2 = surf_area(8.8, b_2)

def find_least_rect(R, b):
    """
    Calculate area of surface
    Calculate area of rectangle
    Assume same height
    Numerically calculate breadth until we find least difference
    """
    sa = surf_area(R, b)
    a = 0
    i = 0
    while True:
        rect = rec_area(a, b)
        if (trunc(sa, 0) == trunc(rect,0) and a*b >= sa):
            break
        a += 0.001 # Trial error
        i += 1 # iter
        #print(trunc(rect,1))
        #print(sa - a*b)
    return a

    
def surf_rect_diff(R, b):
    a = find_least_rect(R, b)
    rect = rec_area(a, b)
    surf = surf_area(R, b)
    diff = rect-surf
    return rect, surf, a, b, diff
    
def answ():
    set1 = surf_rect_diff(7.8, b_1)
    set2 = surf_rect_diff(8.8, b_2)

    spill1 = (set1[0]-set1[1]) / set1[0]
    spill2 = (set2[0]-set2[1]) / set2[0]

    print("SET 1:")
    print("Rectangle area: "+str(set1[0]))
    print("Surface area: "+str(set1[1]))
    print("a: "+str(set1[2]))
    print("b: "+str(set1[3]))
    print("diff: "+str(set1[4]))
    print("spill: "+str(spill1))
    print()
    print("SET 2:")
    print("Rectangle area: "+str(set2[0]))
    print("Surface area: "+str(set2[1]))
    print("a: "+str(set2[2]))
    print("b: "+str(set2[3]))
    print("diff: "+str(set2[4]))
    print("spill: "+str(spill2))
    
    
## UPPG 3 attempt 2
"""
Change find b in ex 3. 

	f(b) = pi, pi/2

We use ekvationslösningsmetoder to find b, 
then we calculate spill by approximating integral under using trapetsmethod

now we take the diff which is side (because pi and pi/2 gives 
straight lines up-down or horsiontal-vertical)
times furthest area away in x-axis ==> rectangle

"""

# f(x) = e^x - sin(x) - 2
# or e^x - sin(x) -2 = 0
# => fp(x) = e^x - cos(x)
x=1
t=1
i = 0
while np.absolute(t) > 5e-5:
    f = np.exp(x)- np.sin(x)-2
    fp = np.exp(x) - np.cos(x)
    g=t
    t=f/fp
    kvad = t/g**2 #check
    linj = t/g
    print("ITER "+str(i))
    print("x "+str(x))
    print("f(x) "+ str(f))
    print("fp(x) " + str(fp))
    print("korr " + str(t))
    print("kvad " + str(kvad))
    print("linj " + str(linj))
    print()
    i+=1
    x = x - t
rot = x
    
# R = given
# a = given
# if f(u) = pi, what is b

# f(x) = 0
# or f(b) = 0
# slutvinkel 2pi = pi/2

R = 8.8
a = 12
u = 2*np.pi
i = 0
b = 0
t = 1
while np.absolute(t) > 5e-5:
    f = ((R - a * np.cos(np.pi)) / 
         np.sqrt(b**2 + (a * np.cos (u))**2)) - np.pi/2
    fp = (a*(b**2)*np.sin(u)/
          ((a**2) * ((np.cos(u))**2) + b**2)**(3/2))
    g=t
    t=f/fp
    kvad = t/g**2 #check
    linj = t/g
    print("ITER "+str(i))
    print("b "+str(b))
    print("f(x) "+ str(f))
    print("fp(x) " + str(fp))
    print("korr " + str(t))
    print("kvad " + str(kvad))
    print("linj " + str(linj))
    print()
    i+=1
    b = b - t
rot = b

# Sekand methoden
x0 = 1
x1 = 1.4
i = 0
while np.absolute(x1-x0) > 5e-6:
    f0 = np.exp(x0) - 2 - np.sin(x0)
    f1 = np.exp(x1) - 2 - np.sin(x1)
    t = f1*(x1-x0) / (f1 - f0)
    print("ITER "+str(i))
    print("x0 "+str(x0))
    print("x1 "+str(x1))
    print("f0 "+str(f0))
    print("f1 "+str(f1))
    print("t "+str(t))
    print()
    x0 = x1
    x1 = x1 - t
    i += 1
rot = x1




R = 8.8
a = 12
u = 2*np.pi

i = 0
b0 = 30
b1 = 40

while np.absolute(b1-b0) > 5e-5:
    f0 = ((R - a * np.cos(np.pi))/
         (np.sqrt(b0**2 + 
                 (a * np.cos(u))**2))) - np.pi/2
    f1 = ((R - a * np.cos(np.pi))/
         (np.sqrt(b1**2 + 
                 (a * np.cos(u))**2))) - np.pi/2
    t = f1*(x1-x0) / (f1 - f0)
    print("ITER "+str(i))
    print("x0 "+str(b0))
    print("x1 "+str(b1))
    print("f0 "+str(f0))
    print("f1 "+str(f1))
    print("t "+str(t))
    print()
    b0 = b1
    b1 = b1 - t
    i += 1
rot = b1

"""
COMPLEMENTARY CODING
"""

def find_b2(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, conv_crit = np.pi, b_start = 0):
    """
    Own function for iterating to find b until convergence criteria is reached
    """
    i = 0
    #b_start = 0
    step = 0.01 #Found for set1 by trial and error
    while True:
        T = trapets(a = a, b = b, h = h, R_set = R_set, a_set = a_set, b_set = b_start)[-1]
        #print(T)
        if (float(str(T)[0:4]) == float(str(conv_crit)[0:4])) or (round(T, 2) == float(str(conv_crit)[0:4])):
            break
        b_start += step
        i += 1 # iter
    return b_start, i


def find_b2(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, conv_crit = np.pi, b_start = 0):
    """
    Own function for iterating to find b until convergence criteria is reached
    """
    i = 0
    #b_start = 0
    step = 0.01 #Found for set1 by trial and error
    while True:
        T = trapets(a = a, b = b, h = h, R_set = R_set, a_set = a_set, b_set = b_start)[-1]
        #print(T)
        if (float(str(T)[0:4]) == float(str(conv_crit)[0:4])) or (round(T, 2) == float(str(conv_crit)[0:4])):
            break
        b_start += step
        i += 1 # iter
    return b_start, i


def find_b3(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, conv_crit = np.pi, b_start = 0):
    """
    Own function for iterating to find b until convergence criteria is reached
    """
    i = 0
    #b_start = 0
    #step = 0.01 #Found for set1 by trial and error
    step = 1
    step_h = 1
    while True:
        T = trapets(a = a, b = b, h = h, R_set = R_set, a_set = a_set, b_set = b_start)[-1]
        #print(T)
        if (float(str(T)[0:4]) == float(str(conv_crit)[0:4])) or (round(T, 2) == float(str(conv_crit)[0:4])):
            break
        if T < conv_crit:
            step_h = step_h / 2
            b_start -= 1
        b_start += 1*step_h
        if i == 0: print(b_start)
        
        i += 1 # iter
    return b_start, i


trapets(a = 0, b = 2*np.pi, h = 0.5, R_set = 7.8, a_set = 5.5, b_set = 10)[-1]
#set 1
find_b3()

#set 2
find_b3(a=0, b=2*np.pi, h=0.5, R_set = 8.8, a_set = 12, conv_crit = np.pi/2, b_start = 0)


#SHOELACE FORMULA

#SHOELACE FORMULA FOR CALCULATING AREA IN POLYGON
def Shoelace(x,y):
    x_shift = np.append(x[-1], x)[:-1]
    y_shift = np.append(y[-1], y)[:-1]
    res =  0.5 * np.abs(np.dot(x,y_shift) - np.dot(y,x_shift))
    return res
print (Shoelace(x,y))

set3 = (7.8, 5.5, 12.125) #b from find_b3()
set4 = (8.8, 12, 31) #b from find_b3()

# set3 spill
# Create Runge-Kutta Values to interpolate from
H, alpha_ett, epsilon_ett, eta_ett, \
    delta_alpha, delta_epsilon, delta_eta, \
    kvot_alpha, kvot_epsilon, kvot_eta, \
    u_print, alpha_print, epsilon_print, eta_print = RungeKutta(set_value = set3, N = 100)
    
print(Shoelace(epsilon_print, eta_print)) #area of polygon, epsilon = x, eta = y

#add middle point to get polygon to go through it
xp, yp = klippArket(set3[0], set3[1], set3[2])
epsilon_print = np.append(epsilon_print, xp)
eta_print = np.append(eta_print, yp)  

Surf_area = Shoelace(epsilon_print, eta_print)
x_dist = epsilon_print.max() - epsilon_print.min()
y_dist = eta_print.max() - eta_print.min()
Rect_area = x_dist*y_dist
diff = Rect_area - Surf_area
spill = (Rect_area - Surf_area) / Rect_area

# set4 spill
# Create Runge-Kutta Values to interpolate from
H, alpha_ett, epsilon_ett, eta_ett, \
    delta_alpha, delta_epsilon, delta_eta, \
    kvot_alpha, kvot_epsilon, kvot_eta, \
    u_print, alpha_print, epsilon_print, eta_print = RungeKutta(set_value = set4, N = 100)

#Add middle point to get polygon to go through it    
print(Shoelace(epsilon_print, eta_print)) #area of polygon, epsilon = x, eta = y
xp, yp = klippArket(set4[0], set4[1], set4[2])
epsilon_print = np.append(epsilon_print, xp)
eta_print = np.append(eta_print, yp)  


Surf_area_2 = Shoelace(epsilon_print, eta_print)
x_dist_2 = epsilon_print.max() - epsilon_print.min()
y_dist_2 = eta_print.max() - eta_print.min()
Rect_area_2 = x_dist_2*y_dist_2
diff_2 = (Rect_area_2 - Surf_area_2)
spill_2 = (Rect_area_2 - Surf_area_2) / Rect_area_2

print("SET 1:")
print("Rectangle area: "+str(Rect_area))
print("Surface area: "+str(Surf_area))
print("a: "+str(x_dist))
print("b: "+str(y_dist))
print("diff: "+str(diff))
print("spill: "+str(spill))
print()
print("SET 2:")
print("Rectangle area: "+str(Rect_area_2))
print("Surface area: "+str(Surf_area_2))
print("a: "+str(x_dist_2))
print("b: "+str(y_dist_2))
print("diff: "+str(diff_2))
print("spill: "+str(spill_2))
    