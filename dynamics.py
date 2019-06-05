import numpy as np
import re
import subprocess as sp

dyn = sp.run(["maxima", "--very-quiet", "-b", "dynamics.mc"], stdout=sp.PIPE)
out = dyn.stdout.decode("utf-8")
#print(out)

A = re.findall(r"Acode:.+?done", out, re.DOTALL)[0]
B = re.findall(r"Bcode:.+?done", out, re.DOTALL)[0]
xn = re.findall(r"xncode:.+?done", out, re.DOTALL)[0]

def factor(S):
  S = re.sub(r"[A-Za-z_]+code:fortran\([A-Za-z_]+\)\s*\n", r"", S)
  S = re.sub(r"\n\s+done", r"", S)
  S = re.sub(r"\n\s{5}.\s+", r"", S)
  S = re.sub(r"([0-9]+)\.([0-9.]*)[dE]([-+][0-9]+)", r"\1.\2e\3", S)
  S = re.sub(r"([^\(\*e ])(-|\+|\*|\/)([^\* ])", r"\1 \2 \3", S)
  S = re.sub(r"([^\(\*e ])(-|\+|\*|\/)([^\* ])", r"\1 \2 \3", S)
  S = re.sub(r"\*\*", r"^", S)
  S = re.sub(r"([a-zA-Z_]+)\(([0-9,]+)\)", r"\1[\2]", S)
  S = re.sub(r",([0-9])", r", \1", S)
  S = re.sub(r"(^|\n)\s+", r"\1", S, re.MULTILINE)

  return S

print("xn code: #############################################")
print(factor(xn))
print("A code: #############################################")
print(factor(A))
print("B code: #############################################")
print(factor(B))
