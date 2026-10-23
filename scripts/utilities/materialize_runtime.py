"""Copy allowlisted, checksum-verified derived inputs into the ignored runtime tree."""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[2]


def materialize():
    entries = json.loads((ROOT / 'config/runtime_artifacts.json').read_text(encoding='utf-8'))
    for entry in entries:
        src = (ROOT / entry['export_path']).resolve()
        dest = (ROOT / entry['runtime_path']).resolve()
        kind = entry.get('kind', 'derived_result')
        valid = (kind == 'derived_result' and src.is_relative_to(ROOT/'source_data') and dest.is_relative_to(ROOT/'data'))
        valid |= (kind == 'frozen_rule' and src.is_relative_to(ROOT/'config/frozen_parameters') and dest.parent == ROOT and dest.suffix == '.md')
        valid |= (kind == 'extension_code' and src.is_relative_to(ROOT/'scripts/exploratory') and dest.parent == ROOT/'scripts' and dest.suffix == '.py')
        if not valid:
            raise ValueError('Artifact path leaves its allowed directory')
        digest = hashlib.sha256(src.read_bytes()).hexdigest()
        if digest != entry['sha256']:
            raise ValueError('Source checksum mismatch: ' + entry['export_path'])
        if dest.exists():
            if hashlib.sha256(dest.read_bytes()).hexdigest() != digest:
                raise FileExistsError('Refusing to overwrite a different runtime artifact: ' + entry['runtime_path'])
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
    print(f'Materialized {len(entries)} derived artifacts; archived source files were not changed.')


if __name__ == '__main__':
    materialize()
