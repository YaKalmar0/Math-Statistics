import numpy as np
from math import sqrt, cos, sin
import matplotlib.pyplot as plt
from scipy.stats import f
from scipy.stats import t

#____input___

alpha = 0.05

x_i = [i for i in range (1,26)]
file = open("regres_input.txt", "r")
y_i = [float(i) for i in file.read().splitlines()]
file.close()

#___Building_the_Model___

def model(x, beta):
	y = beta.item(0) + x*beta.item(1)
	return y

F_t= np.matrix([np.ones(len(x_i)), x_i]) #Design Matrix
F = F_t.transpose()
#print (F)
A = np.matmul(F_t, F) #Informative Matrix
print("This is A:", A)
#print(np.linalg.det(A))
A_inv = np.linalg.inv(A)
#print(A_inv)

eta = np.matrix(y_i).transpose()

beta = np.matmul(np.matmul(A_inv, F_t), eta)
print("This is beta\n", beta)
print("Our model is:\ty* = %f + %fx" % (beta[0], beta[1]))


#___Check__Model___

eta_mean = sum(eta)/len(eta)
#print("This is eta-mean:", eta_mean)
Dispers1, Dispers2 = 0, 0
for i in eta:
	Dispers1 = Dispers1 + (i - eta_mean)**2
Dispers1 = float(Dispers1/(len(eta)-1))
#print("This is numerator:", Dispers1)
denom_gamma = eta - np.matmul(F, beta)
norm_denom = np.linalg.norm(denom_gamma)
sigma = norm_denom**2/(len(eta)-len(beta))
print("this is denom:", sigma)
#print("sigma^2 = ", sigma)

gamma = Dispers1/sigma
t_cr = f.ppf(0.95, len(eta)-1, len(eta)-len(beta))
print("\ngamma = ", gamma, "\tt_cr = ", t_cr)
print("Is gamma > t_cr?", gamma>t_cr)
if gamma > t_cr:
	print("Model is Adequate")


#____Significance_min_Check___
print("\nChecking significance of beta_minimum...")

j = [abs(float(i)) for i in beta].index(min([abs(float(i)) for i in beta])) #finding min of abs in beta
print("This is beta_min:", beta.item(j))
print("beta.item(j) =", beta.item(j))
print("sigma^2 =", sigma)
print("A_inv.item(%d,%d) =" %(j, j), A_inv.item(j,j))

if (beta.item(j) > 0):
	print("\nHypothesis_0: ß = 0\nHypothesis_1: ß > 0")
	gamma2 = beta.item(j)/sqrt(sigma*A_inv.item(j,j))
	t_cr2 = t.ppf(1-alpha, len(eta)-len(beta))
	print("\ngamma_2=", gamma2, "\tt_cr2=", t_cr2)
	print("Is gamma2 > t_cr2?", gamma2>t_cr2)
	if gamma2 > t_cr2:
		print ("beta_%d is significant" %j)
	else:
		print ("beta_%d is NOT significant" %j)
else:
	print("\nHypothesis_0: ß = 0\nHypothesis_1: ß < 0")
	gamma2 = beta.item(j)/sqrt(sigma*A_inv.item(j,j))
	t_cr2 = t.ppf(alpha, len(eta)-len(beta))
	print("\ngamma_2=", gamma2, "\tt_cr2=", t_cr2)
	print("Is gamma2 < t_cr2?", gamma2<t_cr2)
	if gamma2 < t_cr2:
		print ("beta_%d is significant" %j)
	else:
		print ("beta_%d is NOT significant" %j)


#____Trustworthy_Response_Interval_s____

x0 = 16
x0_m_t = np.matrix([1, x0])
x0_m = np.matrix([1, x0]).transpose()
print("x0_m =", x0_m)
func_x0 = model(x0, beta)
print ("This is f(x0):", func_x0)
t_cr3 = t.ppf((1+0.95)/2, len(eta)-len(beta))
print ("t_cr3 =", t_cr3)

eps = t_cr3*sqrt(sigma*float(np.matmul(x0_m_t.dot(A_inv), x0_m)))
#print("matrix mult = ", float(np.matmul(x0_m_t.dot(A_inv), x0_m)))
#print("eps1 =", eps)
interv = (func_x0 - eps, func_x0+eps)
print("Trustworthy interval for response mean in x0:", interv)

eps2 = t_cr3*sqrt(sigma*(1+float(np.matmul(x0_m_t.dot(A_inv), x0_m))))
print("eps2 =", eps2)
interv2 = (func_x0 - eps2, func_x0+eps2)
print("Trustworthy interval for response value in x0:", interv2)


#______GRAPHICS______
y_star = [model(i,beta) for i in np.linspace(-1, 30, 10)] #y(x) in predicted function
plt.grid(which='major', color = 'k')
plt.grid(which='minor', color = 'grey', linestyle = ':')
graph1 = plt.scatter(x_i, y_i, label = 'Sample')
graph2 = plt.plot(np.linspace(-1, 30, 10), y_star, label = 'Linearization')
plt.minorticks_on()
plt.tight_layout()
plt.axis([0., 26., 1, 84.])
#plt.show()