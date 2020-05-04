import matplotlib.pyplot as plt
#from numpy import *
from math import *
#from scipy.stats import geom

def pascal_prob(k, a):
    return (a ** k)/((1 + a)**(k + 1))

sample = [5, 3, 5, 1, 13, 4, 1, 1, 0, 1, 14, 8, 0, 2, 8, 6, 2, 3, 4, 0, 3, 0, 0, 8, 7,
 0, 0, 5, 2, 0, 1, 2, 1, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0, 0, 2, 1, 1, 3, 2, 1, 
 0, 1, 3, 3, 9, 0, 7, 0, 2, 1, 0, 0, 7, 1, 4, 13, 6, 0, 0, 1, 11, 1, 0, 3, 5, 
 0, 1, 3, 3, 0, 1, 0, 2, 0, 0, 1, 0, 0, 0, 6, 6, 4, 0, 2, 6, 0, 5, 4, 8, 3]
sample.sort()
print(sample[49], sample[50])
numb = []
freq = []
for i in range(max(sample)+1):
    if i in sample:
        numb.append(i)
        freq.append(sample.count(i)/len(sample))
cumul_fun = []
for i in range(len(freq)):
    cumul_fun.append(round(sum(freq[0:(i+1)]), 3))
cumul_fun.insert(0,0)

#print(list(zip(numb, freq)))
print(list(zip(numb, cumul_fun)))

#________Mean, Dispersion, Median, Mode, Assymetry______
sample_mean = 0
dispersion = 0
As = 0
for i in numb:
	sample_mean = sample_mean + i*sample.count(i)
sample_mean = sample_mean/len(sample)

for i in sample:
    dispersion = dispersion + (i
     - sample_mean)**2
dispersion = dispersion/len(sample)
dispersion_cor = dispersion*len(sample)/(len(sample)-1)

for i in sample:
	As = As + (i - sample_mean)**3
print("Centred 3rd moment:", As)
As = As/(len(sample)*dispersion**(3/2))

Me = numb[round(len(numb)/2)]
Mo = numb[freq.index(max(freq))]
#___________

print("This is Mode:", Mo)
print("This is Median:", Me)
print("This is Dispersion:", dispersion)
print("This is sample_mean, dispersion_cor:", sample_mean ,dispersion_cor)
print("This is Assymetry:", As)

probab = []
for k in numb:
    probab.append((sample_mean**k)/(sample_mean+1)**(k+1))
#print("These are probabilities", list(zip(numb, probab)))

np_i = [27.7, 20, 14.5, 10.5, 13.1, 12]
n_i = [32, 19, 12, 10, 10, 17]
summ = 0
for i in range(len(np_i)):
    summ = summ + ((n_i[i]-np_i[i])**2)/np_i[i]
print("This is sum of eta:", summ)
gamma = 0.95
t_gamma = 1.96
root1 = (2*sample_mean*len(sample)+t_gamma**2+sqrt((2*sample_mean*len(sample)+t_gamma**2)**2-4*sample_mean**2*len(sample)*(len(sample)-t_gamma**2)))/(2*len(sample)-2*t_gamma**2)
root2 = (2*sample_mean*len(sample)+t_gamma**2-sqrt((2*sample_mean*len(sample)+t_gamma**2)**2-4*sample_mean**2*len(sample)*(len(sample)-t_gamma**2)))/(2*len(sample)-2*t_gamma**2)
print("These are roots for confidence interval:", root1, root2)



#________Graphics_________

grid1 = plt.grid(color='black', lw = 0.5, which = 'major')
plt.minorticks_on()
#plt.axis([0., 14.5, 0., 0.33])
plt.axis([-1., 17, -0.05, 1.05])
grid2 = plt.grid(color='#999999', lw = 0.2, which = 'minor', linestyle='-')
#plt.plot(numb, freq, "-o", label = "Poligon")
#plt.plot([i for i in range(27)], [pascal_prob(i, 5) for i in range(27)], '-.s', color = '0.25', label = r'$\xi \sim Pas(a), P\left\{\xi = 5\right\}$')
plt.savefig('Pictures/pol.png', dpi = 300)

plt.hlines(cumul_fun, [numb[0]-1]+numb, [x+1 for x in [numb[0]-1]+numb], color='tab:blue')
plt.hlines(0, numb[0]-2, numb[0]-1, color='tab:blue')
plt.hlines(1, numb[-1]+1, numb[-1]+2, color='tab:blue');
plt.scatter([x+0.1 for x in numb], cumul_fun[-(len(cumul_fun)-1):], marker='<', label = 'Cumulative Function')

#fig, ax = plt.subplots(1,1)
#ax.hlines(cumul_fun, numb, numb, label = "Cumulative Function")
#plt.step(numb, cumul_fun, label = "Cumulative Function")
#plt.plot(numb, cumul_fun, 'C0o', alpha=0.7)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',   
           ncol=2, mode="expand", borderaxespad=0.)


#plt.show()