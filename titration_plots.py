import sys
from math import log10
from matplotlib import pyplot as plt

L_COLUMN  = 0
RL_COLUMN = 2

for fname in sys.argv[1:]:
    x = [(float(l.split()[L_COLUMN])) for l in open(fname)]
    for i in range(len(x)): x[i] = log10(x[i])
    y = [(float(l.split()[RL_COLUMN])) for l in open(fname)]
    plt.plot(x,y,'o-')
plt.xlabel("$log([L]/c_R $)")
plt.ylabel("$[RL]/c_R$")
plt.savefig("rl-log.png")

plt.clf()
for fname in sys.argv[1:]:
    x = [(float(l.split()[L_COLUMN])) for l in open(fname)]
    y = [(float(l.split()[RL_COLUMN])) for l in open(fname)]
    plt.plot(x,y,'o-')
plt.xlabel("$[L]/c_R $")
plt.ylabel("$[RL]/c_R$")
plt.savefig("rl.png")

