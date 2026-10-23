"""Download one explicitly selected GEO input into the ignored data directory.

Provide a public HTTPS GEO supplementary URL or a GEO sample SOFT URL. This
utility never enumerates/downloads an entire study and never overwrites a file.
"""
from pathlib import Path
import argparse
import hashlib
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--url', required=True)
    parser.add_argument('--output', required=True, help='Path relative to ignored data/')
    parser.add_argument('--sha256', help='Expected original file hash, when available')
    args = parser.parse_args()
    parsed = urllib.parse.urlparse(args.url)
    allowed = {'ftp.ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov'}
    if parsed.scheme != 'https' or parsed.hostname not in allowed or parsed.username or parsed.password:
        raise ValueError('Use a public HTTPS NCBI GEO input URL')
    dest = (ROOT/'data'/args.output).resolve()
    if not dest.is_relative_to(ROOT/'data'):
        raise ValueError('Output must remain in data/')
    partial = dest.with_name(dest.name+'.partial')
    if dest.exists() or partial.exists():
        raise FileExistsError('Output or partial download exists; inspect it before retrying')
    dest.parent.mkdir(parents=True,exist_ok=True)
    h = hashlib.sha256()
    with urllib.request.urlopen(args.url, timeout=120) as response, partial.open('xb') as handle:
        final_host=urllib.parse.urlparse(response.geturl()).hostname
        if final_host not in allowed:
            raise ValueError('Redirect left the allowed NCBI hosts')
        for chunk in iter(lambda: response.read(1024*1024), b''):
            handle.write(chunk)
            h.update(chunk)
    digest=h.hexdigest()
    if args.sha256 and digest.lower() != args.sha256.lower():
        raise ValueError('Checksum mismatch; partial file retained for inspection')
    partial.rename(dest)
    print(f'{digest}  {dest.relative_to(ROOT).as_posix()}')
    print('Check deposited value semantics and source/target eligibility before analysis.')


if __name__ == '__main__':
    main()
