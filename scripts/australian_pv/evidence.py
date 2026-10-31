"""Content-addressed evidence bundles reproducible from immutable Git commits.

No source hashes are changed. Bundles live in an ignored local cache; the
manifest and the pinned Git tree are sufficient to reproduce every byte.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile

REPO = Path(__file__).resolve().parents[2]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def safe_path(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"unsafe archive path: {name}")
    return path


def git_files(repo, revision, roots):
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("evidence revision must be a full Git commit SHA")
    for root in roots:
        safe_path(root)
    raw = subprocess.check_output(["git", "-C", str(repo), "archive", "--format=tar", revision, "--", *roots])
    files = {}
    with tarfile.open(fileobj=io.BytesIO(raw)) as archive:
        for member in archive:
            safe_path(member.name)
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError("evidence must contain regular files only")
            files[member.name] = archive.extractfile(member).read()
    if not files:
        raise ValueError("empty evidence selection")
    return files


def bundle_bytes(files):
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, filename="", mode="wb", mtime=0, compresslevel=9) as compressed:
        with tarfile.open(fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT) as archive:
            for name, data in sorted(files.items()):
                safe_path(name)
                member = tarfile.TarInfo(name)
                member.size = len(data)
                member.mode = 0o644
                archive.addfile(member, io.BytesIO(data))
    return output.getvalue()


def pack(repo, revision, roots, manifest_path):
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", revision+"^{commit}"], text=True).strip()
    files = git_files(repo, revision, roots)
    data = bundle_bytes(files)
    manifest = dict(format=1, source_revision=revision, roots=roots,
                    bundle_sha256=digest(data), bundle_bytes=len(data),
                    files={name: dict(sha256=digest(value), bytes=len(value)) for name, value in sorted(files.items())})
    cache = repo/".evidence-cache"/(manifest["bundle_sha256"]+".tar.gz")
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists() and cache.read_bytes() != data:
        raise ValueError("conflicting content-addressed bundle")
    cache.write_bytes(data)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(manifest, indent=2, sort_keys=True)+"\n"
    if manifest_path.exists() and manifest_path.read_text() != encoded:
        raise ValueError("manifest is immutable; choose a new name")
    manifest_path.write_text(encoded)
    return manifest


def restore(repo, manifest_path, output):
    manifest = json.loads(manifest_path.read_text())
    if manifest["format"] != 1:
        raise ValueError("unsupported evidence format")
    if not re.fullmatch(r"[0-9a-f]{64}", manifest["bundle_sha256"]):
        raise ValueError("invalid bundle digest")
    if not re.fullmatch(r"[0-9a-f]{40}", manifest["source_revision"]):
        raise ValueError("evidence revision must be a full Git commit SHA")
    cache = repo/".evidence-cache"/(manifest["bundle_sha256"]+".tar.gz")
    if cache.exists():
        data = cache.read_bytes()
    else:
        data = bundle_bytes(git_files(repo, manifest["source_revision"], manifest["roots"]))
    if digest(data) != manifest["bundle_sha256"] or len(data) != manifest["bundle_bytes"]:
        raise ValueError("bundle checksum mismatch")
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        for member in archive:
            safe_path(member.name)
            if not member.isfile() or member.name in files:
                raise ValueError("invalid or duplicate evidence member")
            files[member.name] = archive.extractfile(member).read()
    if set(files) != set(manifest["files"]):
        raise ValueError("bundle file coverage mismatch")
    # Validate the complete bundle and destination before writing any member.
    output = output.resolve()
    for name, data in files.items():
        expected = manifest["files"][name]
        if digest(data) != expected["sha256"] or len(data) != expected["bytes"]:
            raise ValueError(f"evidence checksum mismatch: {name}")
        path = output/name
        if not path.resolve().is_relative_to(output):
            raise ValueError("destination escapes output directory")
        if path.exists() and (not path.is_file() or path.read_bytes() != data):
            raise ValueError(f"destination differs; restore to an empty directory: {name}")
    for name, data in files.items():
        path = output/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return len(files)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("pack")
    p.add_argument("revision")
    p.add_argument("manifest", type=Path)
    p.add_argument("roots", nargs="+")
    p = commands.add_parser("restore")
    p.add_argument("manifest", type=Path)
    p.add_argument("--output", type=Path, default=REPO)
    args = parser.parse_args()
    if args.command == "pack":
        result = pack(REPO, args.revision, args.roots, args.manifest)
        print(f"Packed {len(result['files'])} files; SHA256 {result['bundle_sha256']}")
    else:
        print(f"Restored {restore(REPO, args.manifest, args.output)} verified files")
