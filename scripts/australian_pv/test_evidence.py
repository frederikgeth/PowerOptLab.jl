"""Archive restoration must preserve bytes and reject corruption before writes."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from evidence import pack, restore


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name)/"repo"
        self.repo.mkdir()
        subprocess.run(["git", "init", "-q", str(self.repo)], check=True)
        self.run_git("config", "user.name", "Fixture")
        self.run_git("config", "user.email", "fixture@example.invalid")
        (self.repo/"results").mkdir()
        (self.repo/"results/point.json").write_text('{"accepted": false}\n')
        self.run_git("add", "results")
        self.run_git("commit", "-qm", "immutable fixture")
        self.manifest = self.repo/"manifest.json"
        self.meta = pack(self.repo, "HEAD", ["results"], self.manifest)
        self.output = Path(self.temp.name)/"restored"

    def run_git(self, *args):
        subprocess.run(["git", "-C", str(self.repo), *args], check=True)

    def test_cached_restore(self):
        self.assertEqual(restore(self.repo, self.manifest, self.output), 1)
        self.assertEqual((self.output/"results/point.json").read_bytes(), (self.repo/"results/point.json").read_bytes())

    def test_rebuild_from_pinned_git_without_cache(self):
        for file in (self.repo/".evidence-cache").iterdir():
            file.unlink()
        (self.repo/"results/point.json").write_text("changed working copy")
        self.assertEqual(restore(self.repo, self.manifest, self.output), 1)
        self.assertIn('"accepted": false', (self.output/"results/point.json").read_text())

    def test_bundle_corrupted(self):
        next((self.repo/".evidence-cache").iterdir()).write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "bundle checksum"):
            restore(self.repo, self.manifest, self.output)
        self.assertFalse(self.output.exists())

    def test_manifest_digest_corrupted(self):
        self.meta["files"]["results/point.json"]["sha256"] = "0"*64
        self.manifest.write_text(json.dumps(self.meta))
        with self.assertRaisesRegex(ValueError, "evidence checksum"):
            restore(self.repo, self.manifest, self.output)
        self.assertFalse(self.output.exists())

    def test_destination_conflict(self):
        (self.output/"results").mkdir(parents=True)
        (self.output/"results/point.json").write_text("keep this")
        with self.assertRaisesRegex(ValueError, "destination differs"):
            restore(self.repo, self.manifest, self.output)
        self.assertEqual((self.output/"results/point.json").read_text(), "keep this")

    def test_manifest_is_immutable(self):
        (self.repo/"results/point.json").write_text("new generation")
        self.run_git("add", "results")
        self.run_git("commit", "-qm", "new evidence")
        with self.assertRaisesRegex(ValueError, "manifest is immutable"):
            pack(self.repo, "HEAD", ["results"], self.manifest)

    def test_path_traversal_rejected(self):
        self.meta["roots"] = ["../escape"]
        self.meta["bundle_sha256"] = "0"*64
        self.manifest.write_text(json.dumps(self.meta))
        with self.assertRaisesRegex(ValueError, "unsafe archive path"):
            restore(self.repo, self.manifest, self.output)

    def test_revision_cannot_inject_git_options(self):
        self.meta["source_revision"] = "--remote=unexpected"
        self.manifest.write_text(json.dumps(self.meta))
        with self.assertRaisesRegex(ValueError, "full Git commit"):
            restore(self.repo, self.manifest, self.output)

    def test_cache_digest_cannot_escape_cache(self):
        self.meta["bundle_sha256"] = "../unexpected"
        self.manifest.write_text(json.dumps(self.meta))
        with self.assertRaisesRegex(ValueError, "invalid bundle digest"):
            restore(self.repo, self.manifest, self.output)


if __name__ == "__main__":
    unittest.main()
