from .pdb_datastructures import Atom, Residue, Molecule

def read_pdb_file(filename, realname):
    """
    Reads a PDB file and returns a Molecule object.
    """
    molecules = Molecule(realname)

    with open(filename, 'r') as inp:
        prevresnum = -666
        locnum = 0
        res = None

        for line in inp:
            if line.startswith("END") or line.startswith("TER"):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[16] != ' ' and line[16] != 'A':
                    continue

                resnum = int(line[22:26])
                resname = line[17:20].strip()
                atmname = line[12:16].strip()

                if resnum != prevresnum:
                    prevresnum = resnum
                    if res:
                        # Finalize previous residue
                        pass

                    locnum += 1
                    res = Residue(
                        num=resnum,
                        locnum=locnum,
                        natoms=0,
                        type=0, # Will be set later
                        pdbsg=False,
                        protein=False,
                        name=resname,
                        chain=line[21]
                    )
                    molecules.residues.append(res)
                    molecules.nres += 1

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                atom = Atom(
                    x=x,
                    y=y,
                    z=z,
                    name=atmname,
                    num=int(line[6:11]),
                    locnum=0, # will be set later
                    flag=0, # will be set later
                    cispro=False,
                    res=res
                )
                res.atoms.append(atom)
                res.natoms += 1

    return molecules
