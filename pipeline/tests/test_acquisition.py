"""Acquisition safeguards: changed sources, ignored ranges and corrupt members."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zlib
from pipeline import acquire


class AcquisitionIntegrity(unittest.TestCase):
    def test_changed_download_never_becomes_an_input(self):
        with tempfile.TemporaryDirectory() as d:
            repo=Path(d)/'code'; (repo/'data').mkdir(parents=True)
            destination=Path(d)/'inputs'
            row={'path':'counts.csv','bytes':4,'sha256':hashlib.sha256(b'1234').hexdigest()}
            (repo/'data/input_manifest.json').write_text(json.dumps({'schema_version':1,'files':[row]}))
            source={**row,'method':'download','source_url':'https://example.org/counts.csv'}
            (repo/'data/download_manifest.json').write_text(json.dumps({'schema_version':1,'files':[source]}))
            def changed(url,expected_bytes,target=None,byte_range=None): target.write_bytes(b'1235')
            with patch.object(acquire,'REPOSITORY',repo),patch.object(acquire,'fetch',changed):
                with self.assertRaisesRegex(ValueError,'SHA256 mismatch'): acquire.acquire(destination)
            self.assertFalse((destination/'counts.csv').exists())
            self.assertFalse((destination/'counts.csv.partial').exists())
            self.assertFalse((destination/'input_manifest.json').exists())
            self.assertFalse((destination/'.acquisition.lock').exists())

    def test_ignored_range_is_rejected_before_reading_archive(self):
        class Response:
            status=200
            headers={}
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def read(self,*args):raise AssertionError('Full archive must not be read')
        with patch.object(acquire.urllib.request,'urlopen',return_value=Response()),patch.object(acquire.time,'sleep'):
            with self.assertRaisesRegex(ValueError,'honor the exact ZIP byte range'):
                acquire.fetch('https://example.org/archive.zip',4,byte_range=(0,3))

    def test_corrupt_zip_member_crc_is_rejected(self):
        raw=b'gene,count\nA,1\n';compressor=zlib.compressobj(wbits=-15)
        compressed=compressor.compress(raw)+compressor.flush()
        record={'method':'deflate_range','range':[0,len(compressed)-1],
                'source_url':'https://example.org/archive.zip','bytes':len(raw),
                'crc32':(zlib.crc32(raw)+1)&0xffffffff}
        with patch.object(acquire,'fetch',return_value=compressed):
            with self.assertRaisesRegex(ValueError,'CRC32'):acquire.zip_bytes(record)

    def test_source_path_cannot_escape_destination(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'Unsafe manifest path'):
                acquire.within(Path(d).resolve(),'../outside/counts.csv')


if __name__=='__main__':unittest.main()
