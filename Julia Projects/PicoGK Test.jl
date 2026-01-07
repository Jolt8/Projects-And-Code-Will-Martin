using DotNET
using ILUZero

picogk_path = raw"C:\Users\wille\Desktop\C# Repositories\PicoGK\bin\Debug\net10.0\PicoGK.dll"
PicoGK_Assembly = T"System.Reflection.Assembly".LoadFrom(picogk_path)

LatticeType = PicoGK_Assembly.GetType("PicoGK.Lattice")
VoxelsType  = PicoGK_Assembly.GetType("PicoGK.Voxels")

Vector3Type = T"System.Numerics.Vector3" 

lat = LatticeType.new()