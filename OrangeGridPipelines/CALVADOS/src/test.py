from simtk import openmm
system = openmm.XmlSerializer.deserialize(open("system_pro.xml").read())
for fi, f in enumerate(system.getForces()):
    if isinstance(f, openmm.HarmonicBondForce):
        for b in range(f.getNumBonds()):
            i, j, r, k = f.getBondParameters(b)
            if {i, j} & {298, 299}:
                print(f"force#{fi} nbonds={f.getNumBonds()}  {i}-{j}  r={r}  k={k}")
