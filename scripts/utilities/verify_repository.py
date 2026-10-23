"""Verify the immutable exported file manifest and optional full-tree checksums."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]


def verify():
    entries = json.loads((ROOT/'archive/manifests/export_manifest.json').read_text(encoding='utf-8'))['records']
    for entry in entries:
        path = (ROOT/entry['export_path']).resolve()
        if not path.is_relative_to(ROOT): raise ValueError('Invalid manifest path')
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != entry['export_sha256']:
            raise ValueError('Export changed: '+entry['export_path'])
        if not entry['modified'] and actual != entry['source_sha256']:
            raise ValueError('Unchanged source differs: '+entry['export_path'])
    sums = ROOT/'archive/checksums/SHA256SUMS'
    n = 0
    if sums.exists():
        for line in sums.read_text(encoding='utf8').splitlines():
            expected, rel = line.split('  ',1)
            path = (ROOT/rel).resolve()
            if not path.is_relative_to(ROOT): raise ValueError('Invalid checksum path')
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('Checksum mismatch: '+rel)
            n += 1
    print(f'PASS: {len(entries)} export provenance records; {n} package checksums.')


if __name__ == '__main__':
    verify()
