"""Bundle only this owned repo's code and prepared real data, without raw/Git."""
from pathlib import Path
import tarfile

root = Path(__file__).resolve().parents[1]
destination = root / "data" / "remote_bundle.tar.gz"
with tarfile.open(destination, "w:gz", compresslevel=1) as out:
    for name in ["src", "tests", "scripts", "pyproject.toml", "README.md", "LICENSE", ".gitignore", ".gitattributes", "data/prepared", "data/source_manifest.json"]:
        path = root / name
        if not path.exists():
            raise FileNotFoundError(path)
        out.add(path, arcname=name, filter=lambda info: None if "__pycache__" in info.name else info)
print(destination, destination.stat().st_size)
