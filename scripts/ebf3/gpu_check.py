"""Stage-0 GPU check: OpenMM CUDA platform + polysim import + tiny CUDA context."""
import socket, subprocess, sys
print("host", socket.gethostname())
print(subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
                     capture_output=True, text=True).stdout.strip())
smi = subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout
print([l for l in smi.splitlines() if "CUDA Version" in l][0].strip())
import openmm, openmm.unit as u
print("openmm", openmm.__version__, "platforms:",
      [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])
# Build a 2-particle system and actually create a CUDA context (compiles kernels via nvrtc).
s = openmm.System(); [s.addParticle(1.0) for _ in range(2)]
f = openmm.HarmonicBondForce(); f.addBond(0, 1, 0.1, 100.0); s.addForce(f)
ctx = openmm.Context(s, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CUDA"))
ctx.setPositions([[0, 0, 0], [0.2, 0, 0]] * u.nanometer)
ctx.getIntegrator().step(100)
print("CUDA context OK, E =", ctx.getState(getEnergy=True).getPotentialEnergy())
import polysim.extrusion, polysim.sim3d
print("polysim import OK:", polysim.sim3d.__file__)
