"""Small repository-owned alignment gate for explicitly comparable plot rectangles.

This checks physical plot-area edges, repeated gutters, and equal-width rows.
It deliberately requires explicit groups for these heterogeneous figures.
Scientific comparability, annotation collisions and font sizes require separate QA.
"""
from pathlib import Path
import json


def require_matplotlib_panel_alignment(fig, axes=None, panel_ids=None,
                                      row_groups=None, column_groups=None,
                                      exemptions=None, json_out=None, strict=True,
                                      tolerance_pt=1.5, **kwargs):
    fig.canvas.draw()
    axes = list(fig.axes if axes is None else axes)
    if isinstance(panel_ids, dict):
        names = [str(panel_ids[a]) for a in axes]
    else:
        names = list(map(str, panel_ids)) if panel_ids is not None else [str(i) for i in range(len(axes))]
    if len(names) != len(axes) or len(set(names)) != len(names):
        raise ValueError('Panel IDs must be unique and cover the supplied axes')
    w, h = fig.get_size_inches() * 72
    boxes = {}
    for name, ax in zip(names, axes):
        x, y, dx, dy = ax.get_position().bounds
        boxes[name] = {'left': x*w, 'right': (x+dx)*w, 'bottom': y*h,
                       'top': (y+dy)*h, 'width': dx*w, 'height': dy*h}
    checks = []
    exemptions = exemptions or []

    def exempt(ids, kind):
        return any(set(ids) <= set(map(str, e['panels'])) and kind in e.get('checks', [])
                   and e.get('reason') for e in exemptions)

    for kind, groups in [('row', row_groups or []), ('column', column_groups or [])]:
        for i, group in enumerate(groups):
            ids = list(map(str, group['panels'] if isinstance(group, dict) else group))
            unknown = set(ids)-boxes.keys()
            if unknown: raise ValueError('Unknown panel IDs: '+str(sorted(unknown)))
            metrics = ['bottom', 'top'] if kind == 'row' else ['left', 'right']
            if kind == 'row' and exempt(ids, 'row'):
                # Existing asymmetric ED panels share only their top boundary.
                metrics = ['top']
            for metric in metrics:
                span = max(boxes[n][metric] for n in ids)-min(boxes[n][metric] for n in ids)
                checks.append(dict(group=f'{kind}-{i}', metric=metric,
                                   deviation_pt=float(span), passed=bool(span <= tolerance_pt)))
            if kind == 'row' and len(ids) >= 3:
                for metric in ['width']:
                    if exempt(ids, 'panel-width'): continue
                    span = max(boxes[n][metric] for n in ids)-min(boxes[n][metric] for n in ids)
                    checks.append(dict(group=f'{kind}-{i}', metric=metric,
                                       deviation_pt=float(span), passed=bool(span <= tolerance_pt)))
                ordered = sorted(ids, key=lambda n: boxes[n]['left'])
                gaps = [boxes[b]['left']-boxes[a]['right'] for a,b in zip(ordered,ordered[1:])]
                if gaps and not exempt(ids, 'gutter'):
                    span = max(gaps)-min(gaps)
                    checks.append(dict(group=f'{kind}-{i}', metric='gutter',
                                       deviation_pt=float(span), passed=bool(span <= tolerance_pt)))
    failed = [c for c in checks if not c['passed']]
    report = dict(status='FIX BEFORE DELIVERY' if failed else 'PASS' if checks else 'NOT APPLICABLE',
                  tolerance_pt=tolerance_pt, panels=boxes, checks=checks, exemptions=exemptions,
                  scope='Explicit comparable rectangles only; no collision or text audit.')
    if json_out:
        Path(json_out).parent.mkdir(parents=True,exist_ok=True)
        Path(json_out).write_text(json.dumps(report,indent=2),encoding='utf-8')
    if failed and strict: raise ValueError('Panel alignment failed: '+str(failed))
    return report
