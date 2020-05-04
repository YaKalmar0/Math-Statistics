import numpy as np
from math import sqrt, sin ,cos
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import f
from scipy.stats import t

#___input___

alpha = 0.05
x1_i, x2_i, y_i = [], [], []
k = 0
with open ("regres_input2.txt", "r") as file:
	for i in file.read().splitlines():
		if (k<10):
			x1_i.append(int(i))
		elif (k>=10 and k<20):
			x2_i.append(int(i))
		else:
			y_i.append(int(i))
		k += 1
k = 7; r = 2
#print(x1_i, x2_i, y_i)


#___Building_the_Model___

def model(x1, x2, beta):
	y = beta.item(0) + x1*beta.item(1) + x2*beta.item(2) + cos(4*x1)*float(beta[3]) + cos(4*x2)*float(beta[4])
	return y
def new_model(x1, x2, beta):
	y = beta.item(0) + x1*beta.item(1) + cos(4*x1)*beta.item(2) + cos(4*x2)*beta.item(3)
	return y
def old_model(x1, x2, beta):
	y = float(beta[0]) + x1*float(beta[1]) + x2*float(beta[2])

F_t = np.matrix([np.ones(len(x1_i)), x1_i, [np.cos(4*i) for i in x1_i], [np.cos(4*i) for i in x2_i], [np.sin(3*i) for i in x1_i]]) #Design Matrix (for expressed x1)
F = F_t.transpose()
print (F)
A = np.matmul(F_t, F) #Informative Matrix
print ("A =", A)
A_inv = np.linalg.inv(A)
#print (A_inv)

#eta = np.matrix(y_i).transpose()
eta = np.matrix(y_i).transpose()

beta = np.matmul(np.matmul(A_inv, F_t), eta)
print("This is beta:\n", beta)
#print("Our model is:\ty* = %f %fx1 %fx2\n" % (-beta[0]/beta[2], -beta[1]/beta[2], 1/beta[2]))


#___Check__Model___
print("\nChecking Adequacy....")

eta_mean = sum(eta)/len(eta)
#print("This is eta-mean:", eta_mean)
Dispers1, Dispers2 = 0, 0
for i in eta:
	Dispers1 = Dispers1 + (i - eta_mean)**2
Dispers1 = float(Dispers1/(len(eta)-1))
print("This is numerator:", Dispers1)
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
print("sigma^2 =", sigma)
print("A_inv.item(%d,%d) =" %(j, j), A_inv.item(j,j))
print("This is beta_min:", beta.item(j))

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

print("n-m=", len(eta)- len(beta))
#____Trustworthy_Response_Interval_s____
print ("\n\nBuilding Trustworthy Intervals....")

x0 = [5, 6]
x0_m_t = np.matrix([1, x0[0], cos(4*x0[0]), cos(4*x0[1]), sin(3*x0[0])])
x0_m = x0_m_t.transpose()
print("x0_m =\n", x0_m, "\nx0_m_t =\n", x0_m_t)
func_x0 = new_model(x0[0], x0[1], beta)
print ("f(x_0) =", func_x0)
t_cr3 = t.ppf((1+0.95)/2, len(eta)-len(beta))
print ("t_cr3 =", t_cr3)

eps = t_cr3*sqrt(sigma*float(np.matmul(x0_m_t.dot(A_inv), x0_m)))
print("eps1 =", eps)
interv = (func_x0 - eps, func_x0+eps)
print("Trustworthy interval for response mean in x0:", interv)

eps2 = t_cr3*sqrt(sigma*(1+float(np.matmul(x0_m_t.dot(A_inv), x0_m))))
print("eps2 =", eps2)
interv2 = (func_x0 - eps2, func_x0+eps2)
print("Trustworthy interval for response value in x0:", interv2)


#______GRAPHICS______

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter3D(x1_i, x2_i, y_i, c=y_i, cmap="Reds")
x1 = np.linspace(-8, 8, 700)
x2 = np.linspace(-5, 20, 700)
x1, x2 = np.meshgrid(x1, x2)
y_star = beta.item(0) + x1*beta.item(1) + np.cos(4*x1)*beta.item(2) + np.cos(4*x2)*beta.item(3)
#y_star = beta.item(0) + x1*beta.item(1)+x2*beta.item(2)
#ax.plot_surface(x1, x2, y_star, color = "blue")
ax.set_title("Sample")

plt.show()