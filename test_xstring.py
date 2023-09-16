import numpy as np

# set the random seed
np.random.seed(42)

import strawberryfields as sf
from strawberryfields.ops import *
from itertools import permutations
from sympy.combinatorics.permutations import Permutation

cutoff=5

print("testing identity")
prog1 = sf.Program(2)
prog2 = sf.Program(2)
s=1.0
with prog1.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
with prog2.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    xstring1(s,0) | q[0]
eng = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results1 = eng.run(prog1)
results2 = eng.run(prog2)
probs1=results1.state.all_fock_probs()
probs2=results2.state.all_fock_probs()
print("sum probs1",np.sum(probs1))
print("sum probs2",np.sum(probs2))
for i in range(min(cutoff,10)):
    for j in range(min(cutoff,10)):
        if abs(probs1[i,j])>0 or abs(probs2[i,j])>0:
            print(i,j,probs1[i,j],probs2[i,j])

print("testing x^2")
prog1 = sf.Program(2)
prog2 = sf.Program(2)
s=0.4
with prog1.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    Pgate(s) | q[0]
with prog2.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    xstring1(s*0.5,2) | q[0]
eng = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results1 = eng.run(prog1)
results2 = eng.run(prog2)
probs1=results1.state.all_fock_probs()
probs2=results2.state.all_fock_probs()
print("sum probs1",np.sum(probs1))
print("sum probs2",np.sum(probs2))
for i in range(min(cutoff,10)):
    for j in range(min(cutoff,10)):
        if abs(probs1[i,j])>0 or abs(probs2[i,j])>0:
            print(i,j,probs1[i,j],probs2[i,j])


print("testing x^3")
prog1 = sf.Program(2)
prog2 = sf.Program(2)
s=0.4
with prog1.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    Vgate(s) | q[0]
with prog2.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    xstring1(s/3,3) | q[0]
eng = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results1 = eng.run(prog1)
results2 = eng.run(prog2)
probs1=results1.state.all_fock_probs()
probs2=results2.state.all_fock_probs()
print("sum probs1",np.sum(probs1))
print("sum probs2",np.sum(probs2))
for i in range(min(cutoff,10)):
    for j in range(min(cutoff,10)):
        if abs(probs1[i,j])>0 or abs(probs2[i,j])>0:
            print(i,j,probs1[i,j],probs2[i,j])

#required transpositions
p=Permutation([0,2,1,3])
print((~p).list())
print((~p).transpositions())

print("testing x_1x_2")
prog1 = sf.Program(2)
prog2 = sf.Program(2)
s=0.001
with prog1.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    CZgate(s) | (q[0],q[1])
with prog2.context as q:
    Fock(1) | q[0]
    Fock(3) | q[1]
    #CZ2gate(s) | (q[0],q[1])
    xstring2(s,1,1) | (q[0],q[1])
eng1 = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results1 = eng1.run(prog1)
eng2 = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results2 = eng2.run(prog2)
probs1=results1.state.all_fock_probs()
probs2=results2.state.all_fock_probs()
print("sum probs1",np.sum(probs1))
print("sum probs2",np.sum(probs2))
for i in range(min(cutoff,10)):
    for j in range(min(cutoff,10)):
        if abs(probs1[i,j])>0 or abs(probs2[i,j])>0:
            print(i,j,probs1[i,j],probs2[i,j])


#required transpositions
p=Permutation([0,2,4,1,3,5])
print((~p).list())
print((~p).transpositions())

print("testing x_0x_1x_2")
prog = sf.Program(3)
s=0.01
with prog.context as q:
    Fock(1) | q[0]
    Fock(2) | q[1]
    Fock(3) | q[2]
    xstring3(s,1,2,3) | (q[0],q[1],q[2])
eng = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results = eng.run(prog)
probs=results.state.all_fock_probs()
print("sum probs",np.sum(probs))
for i in range(min(cutoff,10)):
    for j in range(min(cutoff,10)):
        for k in range(min(cutoff,10)):
            if abs(probs[i,j,k])>1e-8:
                print(i,j,k,probs[i,j,k])

#required transpositions
p=Permutation([0,2,4,6,1,3,5,7])
print((~p).list())
print((~p).transpositions())

print("testing x_0x_1x_2x_3")
prog = sf.Program(4)
s=0.01
with prog.context as q:
  Fock(0) | q[0]
  Fock(1) | q[1]
  Fock(2) | q[2]
  Fock(3) | q[3]
  xstring4(s,1,1,1,1) | (q[0],q[1],q[2],q[3])
eng = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results = eng.run(prog)
probs=results.state.all_fock_probs()
print("sum probs",np.sum(probs))
for i in range(min(cutoff,10)):
  for j in range(min(cutoff,10)):
      for k in range(min(cutoff,10)):
          for l in range(min(cutoff,10)):
              if abs(probs[i,j,k,l])>1e-8:
                  print(i,j,k,l,probs[i,j,k,l])

print("testing x^2")
prog1 = sf.Program(4)
prog2 = sf.Program(4)
s=0.4
with prog1.context as q:
    Fock(0) | q[0]
    Fock(1) | q[1]
    Fock(2) | q[2]
    Fock(3) | q[3]
    Pgate(s) | q[2]
with prog2.context as q:
    Fock(0) | q[0]
    Fock(1) | q[1]
    Fock(2) | q[2]
    Fock(3) | q[3]
    xstring4(s/2,0,0,2,0) | (q[0],q[1],q[2],q[3])
eng = sf.Engine("fock", backend_options={"cutoff_dim": cutoff})
results1 = eng.run(prog1)
results2 = eng.run(prog2)
probs1=results1.state.all_fock_probs()
probs2=results2.state.all_fock_probs()
print("sum probs1",np.sum(probs1))
print("sum probs2",np.sum(probs2))
for i in range(min(cutoff,10)):
    for j in range(min(cutoff,10)):
        for k in range(min(cutoff,10)):
            for l in range(min(cutoff,10)):
                if abs(probs1[i,j,k,l])>0 or abs(probs2[i,j,k,l])>0:
                    print(i,j,k,l,probs1[i,j,k,l],probs2[i,j,k,l])
