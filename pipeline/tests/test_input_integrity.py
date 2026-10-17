"""A moved bundle must retain integrity and reject changed or escaping files."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from pipeline.bundle import digest, verify_inputs


class InputIntegrity(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="OSDR inputs with spaces ")
        self.root = Path(self.temp.name)
        self.file = self.root / "counts.csv"
        self.file.write_bytes(b"gene,mouse1,mouse2\nPycard,12,8\n")
        self.row = {"path": "counts.csv", "bytes": self.file.stat().st_size, "sha256": digest(self.file)}

    def tearDown(self):
        self.temp.cleanup()

    def manifest(self, rows):
        (self.root / "input_manifest.json").write_text(json.dumps({"schema_version": 1, "files": rows}))

    def test_valid_bundle_in_path_with_spaces(self):
        self.manifest([self.row])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(len(verify_inputs(self.root)["files"]), 1)

    def test_same_length_count_edit_is_rejected(self):
        self.manifest([self.row])
        self.file.write_bytes(self.file.read_bytes().replace(b"12,8", b"13,8"))
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            verify_inputs(self.root)

    def test_manifest_cannot_escape_bundle(self):
        self.manifest([dict(self.row, path="../unrelated.csv")])
        with self.assertRaisesRegex(ValueError, "Invalid or duplicate"):
            verify_inputs(self.root)

    def test_duplicate_records_rejected(self):
        self.manifest([self.row, self.row])
        with self.assertRaisesRegex(ValueError, "Invalid or duplicate"):
            verify_inputs(self.root)


if __name__ == "__main__":
    unittest.main()
