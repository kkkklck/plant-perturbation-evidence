from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse, hashlib, json, sys
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RULE_FILE = ROOT / '任务结构优先证据选择规则_v0.1_冻结_2026-10-21.md'
VERSION = 'TaskStructureEvidenceRanker-v0.1-2026-10-21'
SUPPORTED = {'future_time', 'missing_state', 'cross_experiment', 'combination', 'recovery'}
SHARED_THRESHOLDS = {
    'rank1_energy_min': 0.50,
    'median_pairwise_rho_min': 0.35,
    'effective_state_fraction_min': 0.50,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def resolve_path(v: str | Path) -> Path:
    p = Path(v)
    if not p.is_absolute():
        p = ROOT / p
    return p.resolve()


def read_table(path: Path) -> pd.DataFrame:
    n = path.name.lower()
    if n.endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')):
        return pd.read_csv(path, sep='\t')
    if n.endswith(('.csv', '.csv.gz')):
        return pd.read_csv(path)
    raise ValueError(f'不支持的输入格式：{path.name}')


def parse_boolish(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.astype(bool)
    x = s.astype(str).str.strip().str.lower()
    t = {'true', '1', 'yes', 'y', 't'}
    f = {'false', '0', 'no', 'n', 'f'}
    bad = ~(x.isin(t | f))
    if bad.any():
        raise ValueError(f'筛选列存在无法解释为布尔值的内容：{sorted(x[bad].unique())[:8]}')
    return x.isin(t)


def eligible_mask(df: pd.DataFrame, cfg: dict[str, Any]) -> np.ndarray:
    filt = cfg.get('filter')
    if not filt:
        return np.ones(len(df), dtype=bool)
    col = filt.get('column')
    if col not in df.columns:
        raise ValueError(f'筛选列不存在：{col}')
    target = filt.get('equals', True)
    if isinstance(target, bool):
        return (parse_boolish(df[col]).to_numpy() == target)
    return (df[col].astype(str).to_numpy() == str(target))


def direction(v: np.ndarray) -> np.ndarray:
    out = np.full(len(v), '不变', dtype=object)
    out[v > 0] = '上调'
    out[v < 0] = '下调'
    return out


def rank_output(gene: pd.Series, pred: np.ndarray, eligible: np.ndarray, extra: dict[str, Any]) -> pd.DataFrame:
    out = pd.DataFrame({
        'gene': gene.astype(str),
        'eligible': eligible,
        'predicted_response': pred,
        'predicted_direction': direction(np.nan_to_num(pred, nan=0.0)),
        'priority_score': np.abs(pred),
    })
    for k, v in extra.items():
        out[k] = v
    q = out[out.eligible & np.isfinite(out.priority_score)].copy()
    q = q.sort_values(['priority_score', 'gene'], ascending=[False, True], kind='mergesort')
    q.insert(0, 'priority_rank', np.arange(1, len(q) + 1, dtype=int))
    return q.reset_index(drop=True)


def shared_geometry(E: np.ndarray) -> dict[str, float]:
    # E: genes x evidence states. Geometry is computed on genes finite in every state.
    complete = np.all(np.isfinite(E), axis=1)
    M = E[complete].T
    if M.shape[0] < 2 or M.shape[1] < 20:
        raise ValueError('共享结构诊断需要至少2个来源状态和足够的共同基因')
    u, s, vt = np.linalg.svd(M, full_matrices=False)
    e = s * s / np.sum(s * s)
    p = u[:, 0] ** 2
    p = p / p.sum()
    prs, loo = [], []
    for i in range(M.shape[0]):
        for j in range(i + 1, M.shape[0]):
            if np.std(M[i]) > 0 and np.std(M[j]) > 0:
                prs.append(float(spearmanr(M[i], M[j]).statistic))
        oth = np.mean(np.delete(M, i, axis=0), axis=0)
        if np.std(M[i]) > 0 and np.std(oth) > 0:
            loo.append(float(spearmanr(M[i], oth).statistic))
    eff = float(1.0 / np.sum(p * p))
    return {
        'n_complete_genes': int(M.shape[1]),
        'n_source_states': int(M.shape[0]),
        'rank1_energy': float(e[0]),
        'rank1_3_energy': float(e[:min(3, len(e))].sum()),
        'median_pairwise_rho': float(np.nanmedian(prs)),
        'median_loo_shared_rho': float(np.nanmedian(loo)),
        'max_rank1_loading_share': float(p.max()),
        'effective_states': eff,
        'effective_state_fraction': eff / M.shape[0],
    }


def shared_ok(g: dict[str, float]) -> bool:
    return (
        g['rank1_energy'] >= SHARED_THRESHOLDS['rank1_energy_min'] and
        g['median_pairwise_rho'] >= SHARED_THRESHOLDS['median_pairwise_rho_min'] and
        g['effective_state_fraction'] >= SHARED_THRESHOLDS['effective_state_fraction_min']
    )


def signed_mean(E: np.ndarray) -> dict[str, np.ndarray]:
    present = np.isfinite(E)
    n = present.sum(axis=1)
    filled = np.where(present, E, 0.0)
    pred = np.full(E.shape[0], np.nan)
    ok = n > 0
    pred[ok] = filled.sum(axis=1)[ok] / n[ok]
    denom = np.abs(filled).sum(axis=1)
    coh = np.full(E.shape[0], np.nan)
    nz = denom > 0
    coh[nz] = np.abs(filled.sum(axis=1)[nz]) / denom[nz]
    return {'pred': pred, 'support_count': n, 'support_fraction': n / E.shape[1], 'coherence': coh}


def write_abstain(outdir: Path, cfg: dict[str, Any], config_path: Path, reason: str, diagnostics: dict[str, Any] | None = None) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    man = {
        'version': VERSION,
        'status': '拒判',
        'task_type': cfg.get('task_type'),
        'reason': reason,
        'diagnostics': diagnostics,
        'config_sha256': sha256_file(config_path),
        'rule_sha256': sha256_file(RULE_FILE),
        'generated_utc': datetime.now(timezone.utc).isoformat(),
    }
    (outdir / 'manifest.json').write_text(json.dumps(man, ensure_ascii=False, indent=2), encoding='utf-8')
    print('任务结构证据排序：拒判')
    print(json.dumps(man, ensure_ascii=False, indent=2))
    return 0


def run(config_path: Path) -> int:
    cfg = json.loads(config_path.read_text(encoding='utf-8'))
    task = str(cfg.get('task_type', '')).strip()
    outdir = resolve_path(cfg.get('output_dir', 'data/task_structure_output'))
    if task not in SUPPORTED:
        return write_abstain(outdir, cfg, config_path, f'任务类型 {task!r} 尚无冻结的基因级规则')

    inp = resolve_path(cfg['input'])
    df = read_table(inp)
    gene_col = cfg.get('gene_column', 'gene')
    if gene_col not in df.columns:
        raise ValueError(f'基因列不存在：{gene_col}')
    elig = eligible_mask(df, cfg)
    diagnostics: dict[str, Any] = {'task_type': task}
    extra: dict[str, Any] = {'task_type': task}

    if task == 'future_time':
        latest = cfg.get('latest_column')
        if not latest or latest not in df.columns:
            return write_abstain(outdir, cfg, config_path, '未来时间预测缺少同状态最近历史列')
        pred = pd.to_numeric(df[latest], errors='coerce').to_numpy(float)
        elig = elig & np.isfinite(pred)
        diagnostics.update({'primary_evidence': latest, 'decision': '最近同状态历史', 'sharedness_does_not_override': True})
        extra['evidence_role'] = '最近同状态历史'

    elif task in {'missing_state', 'cross_experiment', 'recovery'}:
        donors = list(cfg.get('donor_columns', []))
        if len(donors) < 2 or any(c not in df.columns for c in donors):
            return write_abstain(outdir, cfg, config_path, '跨状态/跨实验迁移至少需要2个可用来源响应列')
        E = np.column_stack([pd.to_numeric(df[c], errors='coerce').to_numpy(float) for c in donors])
        geom = shared_geometry(E[elig])
        diagnostics.update({'source_columns': donors, 'shared_geometry': geom, 'shared_thresholds': SHARED_THRESHOLDS, 'shared_sufficient': shared_ok(geom)})
        if not shared_ok(geom):
            return write_abstain(outdir, cfg, config_path, '来源共享结构不足，不输出细粒度基因候选', diagnostics)
        agg = signed_mean(E)
        pred = agg['pred']
        elig = elig & np.isfinite(pred)
        extra.update({
            'evidence_role': '跨状态共同响应' if task != 'cross_experiment' else '跨实验共同响应',
            'evidence_support_count': agg['support_count'],
            'evidence_support_fraction': agg['support_fraction'],
            'evidence_direction_coherence': agg['coherence'],
        })
        diagnostics['decision'] = extra['evidence_role']

    elif task == 'combination':
        parts = list(cfg.get('component_columns', []))
        if len(parts) < 2 or any(c not in df.columns for c in parts):
            return write_abstain(outdir, cfg, config_path, '组合任务至少需要2个组成单扰动响应列')
        E = np.column_stack([pd.to_numeric(df[c], errors='coerce').to_numpy(float) for c in parts])
        present = np.all(np.isfinite(E), axis=1)
        pred = np.nansum(E, axis=1)
        pred[~present] = np.nan
        elig = elig & present
        strongest_idx = np.argmax(np.abs(np.where(np.isfinite(E), E, 0.0)), axis=1)
        strongest = E[np.arange(len(E)), strongest_idx]
        denom = np.abs(E).sum(axis=1)
        dominance = np.full(len(E), np.nan)
        nz = denom > 0
        dominance[nz] = np.max(np.abs(E[nz]), axis=1) / denom[nz]
        coherence = np.full(len(E), np.nan)
        coherence[nz] = np.abs(pred[nz]) / denom[nz]
        extra.update({'evidence_role': '组成单扰动有符号相加', 'strongest_component_response': strongest, 'component_dominance_share': dominance, 'component_direction_coherence': coherence})
        diagnostics.update({'component_columns': parts, 'decision': '组成单扰动有符号相加', 'important_note': '组合预测与最强单成员应并列解释；相加不保证产生额外增益'})

    ranked = rank_output(df[gene_col], pred, elig, extra)
    outdir.mkdir(parents=True, exist_ok=True)
    ranking = outdir / 'ranking.tsv.gz'
    ranked.to_csv(ranking, sep='\t', index=False, compression='gzip')
    top_ks = [int(x) for x in cfg.get('top_k', [100, 250, 500, 1000]) if int(x) > 0]
    maxk = min(max(top_ks) if top_ks else 100, len(ranked))
    ranked.head(maxk).to_csv(outdir / 'top_candidates.tsv', sep='\t', index=False)
    summary = []
    for k in top_ks:
        kk = min(k, len(ranked))
        if kk <= 0: continue
        q = ranked.head(kk)
        row = {'k': kk, 'median_priority_score': float(q.priority_score.median())}
        for c in ['evidence_support_fraction', 'evidence_direction_coherence', 'component_dominance_share', 'component_direction_coherence']:
            if c in q.columns:
                row[f'median_{c}'] = float(q[c].median(skipna=True))
        summary.append(row)
    pd.DataFrame(summary).to_csv(outdir / 'top_k_summary.tsv', sep='\t', index=False)
    (outdir / 'source_diagnostics.json').write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding='utf-8')
    man = {
        'version': VERSION,
        'status': '已排序',
        'task_name': cfg.get('task_name'),
        'task_type': task,
        'input': str(inp),
        'input_sha256': sha256_file(inp),
        'config_sha256': sha256_file(config_path),
        'rule_sha256': sha256_file(RULE_FILE),
        'ranking_sha256': sha256_file(ranking),
        'n_input_rows': int(len(df)),
        'n_ranked_genes': int(len(ranked)),
        'diagnostics': diagnostics,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
    }
    (outdir / 'manifest.json').write_text(json.dumps(man, ensure_ascii=False, indent=2), encoding='utf-8')
    print('任务结构证据排序：成功')
    print(json.dumps(man, ensure_ascii=False, indent=2))
    print('\n前20候选')
    cols = ['priority_rank', 'gene', 'priority_score', 'predicted_direction', 'predicted_response', 'evidence_role']
    print(ranked[[c for c in cols if c in ranked.columns]].head(20).to_string(index=False))
    return 0


def main() -> None:
    p = argparse.ArgumentParser(description='任务结构优先的植物扰动证据排序入口')
    p.add_argument('--config', required=True)
    args = p.parse_args()
    try:
        rc = run(resolve_path(args.config))
    except Exception as e:
        print(f'运行失败：{type(e).__name__}: {e}', file=sys.stderr)
        raise
    raise SystemExit(rc)

if __name__ == '__main__':
    main()
