import abc
import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional

from .utils import CmdWrapperBase, ensure_path, local_cwd


class DockBase(object, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def dock(self) -> Path:
        raise NotImplementedError

    @abc.abstractmethod
    def create_complex(self):
        raise NotImplementedError


class HDock(DockBase):
    """
    HDock wrapper. See
    more details at http://hdock.phys.hust.edu.cn/.
    """
    def __init__(self) -> None:
        self.hdock_bin = CmdWrapperBase("hdock")
        self.hdock_create_bin = CmdWrapperBase("createpl")
        self.local_dir = tempfile.TemporaryDirectory(prefix='hdock-')

    def __del__(self):
        if hasattr(self, "local_dir"):
            self.local_dir.cleanup()

    def _filter_dockout(
            self, 
            in_dockout: Path, 
            out_dockout: Path, 
            new_receptor_path: str, 
            new_ligand_path: str, 
            copy: bool = False):
        with open(in_dockout, "r") as f:
            lines = f.readlines()

        pattern = re.compile(r"^(\S+)")

        if copy:
            old_receptor_path = pattern.search(lines[3]).group(0)
            old_ligand_path = pattern.search(lines[4]).group(0)
            shutil.copy(old_receptor_path, new_receptor_path)
            shutil.copy(old_ligand_path, new_ligand_path)

        lines[3] = pattern.sub(str(new_receptor_path), lines[3])
        lines[4] = pattern.sub(str(new_ligand_path), lines[4])

        with open(out_dockout, "w") as f:
            f.writelines(lines)

    def dock(self,
             ligand_pdb: Path,
             receptor_pdb: Path,
             output: Optional[Path] = None,
             itscore: bool = True,
             rsite: Optional[Path] = None,
             lsite: Optional[Path] = None,
             angle: int = 15) -> Path:
        ligand_pdb = ensure_path(ligand_pdb)
        receptor_pdb = ensure_path(receptor_pdb)
        rsite = ensure_path(rsite) if rsite else None
        lsite = ensure_path(lsite) if lsite else None
        output = ensure_path(output) if output else None

        with local_cwd(self.local_dir.name) as cwd:
            ligand_temp_pdb = tempfile.NamedTemporaryFile(
                dir=cwd, prefix="Hdock-ligand-", suffix=".pdb")
            ligand_temp_pdb = Path(ligand_temp_pdb.name)
            shutil.copy(ligand_pdb, ligand_temp_pdb)
            receptor_temp_pdb = tempfile.NamedTemporaryFile(
                dir=cwd, prefix="Hdock-receptor-", suffix=".pdb")
            receptor_temp_pdb = Path(receptor_temp_pdb.name)
            shutil.copy(receptor_pdb, receptor_temp_pdb)
            temp_out = tempfile.NamedTemporaryFile(
                dir=cwd, delete=False, prefix="Hdock-", suffix=".out")
            temp_out = Path(temp_out.name)
            cmds = [
                    receptor_temp_pdb.name,
                    ligand_temp_pdb.name,
                    "-angle", angle,
                    "-out", temp_out.name]
            if not itscore:
                cmds.extend(["-itscore", "false"])
            if rsite:
                temp_rsite = tempfile.NamedTemporaryFile(
                dir=cwd, delete=False, prefix="Hdock-", suffix=".rsite")
                temp_rsite = ensure_path(temp_rsite.name)
                shutil.copy(rsite, temp_rsite)
                cmds.extend(['-rsite', temp_rsite.name])
            if lsite:
                temp_lsite = tempfile.NamedTemporaryFile(
                dir=cwd, delete=False, prefix="Hdock-", suffix=".lsite")
                temp_lsite = ensure_path(temp_lsite.name)
                shutil.copy(lsite, temp_lsite)
                cmds.extend(['-lsite', temp_lsite.name])

            process = self.hdock_bin(*cmds)
            if process.returncode != 0:
                raise RuntimeError(f"HDock failed with exit code {process.returncode}.")

            if output is not None:
                output.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(temp_out.name, output)
            else:
                output = temp_out.name

            self._filter_dockout(
                output,
                output,
                os.path.relpath(str(receptor_pdb), str(output.parent)),
                os.path.relpath(str(ligand_pdb), str(output.parent)))

            return output

    def create_complex(self,
                       dock_result: Path,
                       pdb_name: str,
                       model_num: int = 100,
                       rmsd: float = 5.0,
                       complex: bool = False,
                       models: bool = False,
                       chid: bool = False):
        dock_result = ensure_path(dock_result)
        outpath = dock_result.parent
        outpath.mkdir(parents=True, exist_ok=True)

        with local_cwd(outpath) as cwd:
            # copy to current directory as tmp files
            dock_result_tmp = tempfile.NamedTemporaryFile(
                dir=cwd, prefix="Hdock-dockout-", suffix=".out")
            receptor_tmp = tempfile.NamedTemporaryFile(
                dir=cwd, prefix="Hdock-receptor-", suffix=".pdb")
            ligand_tmp = tempfile.NamedTemporaryFile(
                dir=cwd, prefix="Hdock-ligand-", suffix=".pdb")
            complex_tmp = tempfile.NamedTemporaryFile(
                dir=cwd,
                prefix="Hdock-complex-",
                suffix=".pdb",
                delete=False
            )
            self._filter_dockout(
                in_dockout=dock_result.name,
                out_dockout=dock_result_tmp.name,
                new_receptor_path=os.path.basename(receptor_tmp.name),
                new_ligand_path=os.path.basename(ligand_tmp.name),
                copy=True)

            cmds = [os.path.basename(dock_result_tmp.name),
                    os.path.basename(complex_tmp.name),
                    "-nmax", str(model_num),
                    "-rmsd", str(rmsd)]
            if complex:
                cmds.append("-complex")
            if models:
                cmds.append("-models")
            if chid:
                cmds.append("-chid")

            status = self.hdock_create_bin(*cmds).returncode
            if status != 0:
                raise RuntimeError(
                    f"HDock create complex failed with exit code {status}.")

            os.rename(complex_tmp.name, pdb_name)


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    hdock_parser = subparsers.add_parser('hdock')
    hdock_parser.add_argument("--receptor", "-r", type=Path, required=True)
    hdock_parser.add_argument("--ligand", "-l", type=Path, required=True)
    hdock_parser.add_argument("--output", "-o", type=Path, required=True)
    hdock_parser.add_argument("--model-num", "-n", type=int, default=1)

    args = parser.parse_args()
    if args.cmd == "hdock":
        hdock = HDock()
        dockout = hdock.dock(
            ligand_pdb=args.ligand,
            receptor_pdb=args.receptor,
            output=args.output.with_suffix(".dockout")
        )
        hdock.create_complex(
            dockout,
            pdb_name=args.output.name,
            model_num=args.model_num,
            complex=True
        )
    else:
        raise RuntimeError(f"Unknown subcommand: {args.cmd}")
