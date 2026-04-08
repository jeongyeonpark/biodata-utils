import logging
import os
import argparse
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)


def normalize_pdb_name(pdb_input: str) -> str:
    return os.path.splitext(os.path.basename(pdb_input.strip()))[0]


def compress_to_ranges(positions: list[int]) -> str:
    if not positions:
        return ""
    # [BUG-2] 중복 위치값 제거 후 정렬
    positions = sorted(set(positions))
    ranges = []
    start = end = positions[0]
    for p in positions[1:]:
        if p == end + 1:
            end = p
        else:
            ranges.append(f"{start}-{end}" if start != end else str(start))
            start = end = p
    ranges.append(f"{start}-{end}" if start != end else str(start))
    return ",".join(ranges)


def select_best_altloc(residue) -> str:
    """Select the conformer with highest average atom occupancy from a disordered residue.

    Side effect: calls residue.disordered_select(), modifying the residue object in-place.
    """
    best_altloc = None
    best_occ = -1
    for altloc in residue.disordered_get_id_list():
        child = residue.child_dict[altloc]
        atoms = list(child.get_atoms())
        if atoms:
            occ = sum(a.get_occupancy() or 0 for a in atoms) / len(atoms)
            if occ > best_occ:
                best_occ = occ
                best_altloc = altloc
    if best_altloc:
        residue.disordered_select(best_altloc)
    else:
        # [REC-1] 모든 conformer에 atoms가 없는 경우 경고 기록
        logger.warning(
            "Disordered residue %s has no atoms in any conformer; "
            "returning current resname without selection.",
            residue.id,
        )
    return residue.get_resname()


def is_histag(icode_entries: list[tuple]) -> bool:
    """Return True if insertion-code residues contain >= 4 consecutive HIS residues.

    [BUG-1] 기존 코드는 총 HIS 개수만 세었으나, His-tag는 연속된 HIS 잔기를
    의미하므로 연속 런(run) 길이로 판별하도록 수정.
    """
    max_run = cur = 0
    for _, _, aa, _ in icode_entries:
        if aa == "H":
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 0
    return max_run >= 4


def extract_chain_sequences(pdb_file: str, inscode_mode: int = 1) -> list[SeqRecord]:
    """
    Extract per-chain amino acid sequences from a PDB file.

    Alternate conformations: the conformer with highest average occupancy is selected.

    inscode_mode:
        1 - include insertion code residues in sequence order (default)
        2 - remove all insertion code residues
        3 - remove insertion code residues only if they form a His-tag (>= 4 consecutive HIS)

    Non-consecutive residue numbers are represented as '-' in the sequence.
    Annotations (altconf, inscodes, gaps) are reported in the FASTA header and in
    record.annotations for downstream use.

    Returns:
        list of Bio.SeqRecord
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]
    pdb_name = normalize_pdb_name(pdb_file)
    chain_records = []

    for chain in model:
        chain_id = chain.id
        raw_entries = []
        altconf_resnums = []

        for residue in chain:
            if not is_aa(residue, standard=True):
                continue

            resnum = residue.id[1]
            icode = residue.id[2].strip()
            has_altconf = residue.is_disordered()

            if has_altconf:
                altconf_resnums.append(resnum)
                resname = select_best_altloc(residue)
            else:
                resname = residue.get_resname()

            try:
                aa = seq1(resname)
            except Exception:
                # [REC-2] seq1()은 비표준 잔기명에 대해 KeyError 외 다른 예외도 발생 가능
                aa = "X"

            raw_entries.append((resnum, icode, aa, has_altconf))

        if not raw_entries:
            continue

        icode_entries = [(rn, ic, aa, hac) for rn, ic, aa, hac in raw_entries if ic]
        all_icode_labels = [f"{rn}{ic}" for rn, ic, _, _ in icode_entries]

        if inscode_mode == 2:
            entries = [(rn, ic, aa, hac) for rn, ic, aa, hac in raw_entries if not ic]
            removed_inscodes = all_icode_labels
            kept_inscodes = []
        elif inscode_mode == 3:
            if icode_entries and is_histag(icode_entries):
                entries = [(rn, ic, aa, hac) for rn, ic, aa, hac in raw_entries if not ic]
                removed_inscodes = all_icode_labels
                kept_inscodes = []
            else:
                entries = raw_entries
                removed_inscodes = []
                kept_inscodes = all_icode_labels
        else:
            entries = raw_entries
            removed_inscodes = []
            kept_inscodes = all_icode_labels

        if not entries:
            continue

        entries.sort(key=lambda x: (x[0], x[1]))

        resnum_groups: dict[int, list[tuple[str, str]]] = defaultdict(list)
        for rn, ic, aa, _ in entries:
            resnum_groups[rn].append((ic, aa))

        start = entries[0][0]
        end = entries[-1][0]
        sequence = ""
        gap_positions = []

        for i in range(start, end + 1):
            if i in resnum_groups:
                group = sorted(resnum_groups[i], key=lambda x: (x[0] != "", x[0]))
                for _, aa in group:
                    sequence += aa
            else:
                sequence += "-"
                gap_positions.append(i)

        res_count = len(sequence.replace("-", ""))
        desc_parts = [f"chain={chain_id}", f"residues={start}-{end}", f"res_count={res_count}"]
        if altconf_resnums:
            desc_parts.append(f"altconf={','.join(str(r) for r in sorted(set(altconf_resnums)))}")
        if kept_inscodes:
            desc_parts.append(f"inscodes={','.join(kept_inscodes)}")
        if removed_inscodes:
            desc_parts.append(f"removed_inscodes={','.join(removed_inscodes)}")
        if gap_positions:
            desc_parts.append(f"gaps={compress_to_ranges(gap_positions)}")

        record = SeqRecord(
            Seq(sequence),
            id=f"{pdb_name}_{chain_id}",
            description=" | ".join(desc_parts),
        )
        record.annotations["altconf_resnums"] = sorted(set(altconf_resnums))
        record.annotations["kept_inscodes"] = kept_inscodes
        record.annotations["removed_inscodes"] = removed_inscodes
        record.annotations["gap_positions"] = gap_positions
        chain_records.append(record)

    return chain_records


def write_fasta(pdb_name: str, records: list[SeqRecord], output_dir: str = ".") -> str:
    output_file = os.path.join(output_dir, f"{pdb_name}.fasta")
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fasta")
    return output_file


def process_single_pdb(args_tuple: tuple) -> dict:
    pdb_input, output_dir, inscode_mode = args_tuple
    pdb_name = normalize_pdb_name(pdb_input)

    if os.path.isfile(pdb_input):
        pdb_file = pdb_input
    elif os.path.isfile(f"{pdb_name}.pdb"):
        pdb_file = f"{pdb_name}.pdb"
    else:
        return {
            "pdb": pdb_name,
            "status": "error",
            "message": f"File not found: {pdb_input}",
        }

    try:
        records = extract_chain_sequences(pdb_file, inscode_mode=inscode_mode)
        output_file = write_fasta(pdb_name, records, output_dir)

        chain_info = []
        for r in records:
            # [BUG-3] pdb_name에 '_'가 포함될 경우 split("_")[-1]이 잘못된 값을 반환할 수 있음
            # rsplit("_", 1)로 오른쪽에서 한 번만 분리하여 chain_id를 안전하게 추출
            chain_id = r.id.rsplit("_", 1)[-1]
            seq_str = str(r.seq)
            chain_info.append(
                {
                    "chain": chain_id,
                    "seq_length": len(seq_str),
                    "gap_count": seq_str.count("-"),
                    "has_gap": "-" in seq_str,
                    "altconf_count": len(r.annotations.get("altconf_resnums", [])),
                    "kept_inscodes": ",".join(r.annotations.get("kept_inscodes", [])),
                    "removed_inscodes": ",".join(r.annotations.get("removed_inscodes", [])),
                    "gap_positions": compress_to_ranges(r.annotations.get("gap_positions", [])),
                }
            )

        return {
            "pdb": pdb_name,
            "status": "success",
            "output": output_file,
            "chains": len(records),
            "chain_info": chain_info,
            "has_gap": any(ci["has_gap"] for ci in chain_info),
        }
    except Exception as e:
        return {"pdb": pdb_name, "status": "error", "message": str(e)}


def load_pdb_list(input_file: str, pdb_col: str | None = None) -> list[str]:
    if pdb_col:
        # [REC-3] sep=None 자동 감지는 신뢰성이 낮으므로 탭 우선, 실패 시 쉼표로 재시도
        try:
            df = pd.read_csv(input_file, sep="\t")
            if pdb_col not in df.columns:
                df = pd.read_csv(input_file, sep=",")
        except Exception:
            df = pd.read_csv(input_file, sep=None, engine="python")
        return df[pdb_col].dropna().tolist()
    with open(input_file) as f:
        return [line.strip() for line in f if line.strip()]


def write_report(results: list[dict], output_file: str) -> pd.DataFrame:
    rows = []
    for r in results:
        if r["status"] == "success":
            for ci in r["chain_info"]:
                rows.append(
                    {
                        "pdb": r["pdb"],
                        "chain": ci["chain"],
                        "seq_length": ci["seq_length"],
                        "gap_count": ci["gap_count"],
                        "has_gap": ci["has_gap"],
                        "altconf_count": ci["altconf_count"],
                        "kept_inscodes": ci["kept_inscodes"],
                        "removed_inscodes": ci["removed_inscodes"],
                        "gap_positions": ci["gap_positions"],
                        "output_file": r["output"],
                        "status": "success",
                    }
                )
        else:
            rows.append(
                {
                    "pdb": r["pdb"],
                    "chain": pd.NA,
                    "seq_length": pd.NA,
                    "gap_count": pd.NA,
                    "has_gap": pd.NA,
                    "altconf_count": pd.NA,
                    "kept_inscodes": pd.NA,
                    "removed_inscodes": pd.NA,
                    "gap_positions": pd.NA,
                    "output_file": pd.NA,
                    "status": f"error: {r.get('message', '')}",
                }
            )

    columns = [
        "pdb", "chain", "seq_length", "gap_count", "has_gap",
        "altconf_count", "kept_inscodes", "removed_inscodes", "gap_positions",
        "output_file", "status",
    ]
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(output_file, index=False)
    return df


def run(
    input_file: str,
    pdb_col: str | None = None,
    output_dir: str = ".",
    report: str = "pdb_seq_report.csv",
    parallel: bool = False,
    workers: int = 4,
    inscode_mode: int = 1,
) -> tuple[list[dict], pd.DataFrame]:
    os.makedirs(output_dir, exist_ok=True)
    pdb_list = load_pdb_list(input_file, pdb_col)
    job_args = [(pdb, output_dir, inscode_mode) for pdb in pdb_list]

    if parallel:
        # [OPT-2] chunksize 지정으로 IPC 오버헤드 감소
        chunksize = max(1, len(job_args) // (workers * 4))
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(process_single_pdb, job_args, chunksize=chunksize))
    else:
        results = [process_single_pdb(a) for a in job_args]

    df = write_report(results, report)
    return results, df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract sequences from PDB files and write per-chain FASTA with gap notation (-) for non-consecutive residues."
    )
    parser.add_argument(
        "input_file",
        help="File containing PDB names (one per line) or tabular file with --pdb_col",
    )
    parser.add_argument(
        "--pdb_col",
        default=None,
        help="Column name for PDB identifiers in tabular input",
    )
    parser.add_argument(
        "--output_dir",
        default=".",
        help="Output directory for FASTA files (default: current directory)",
    )
    parser.add_argument(
        "--report",
        default="pdb_seq_report.csv",
        help="Report output file path (default: pdb_seq_report.csv)",
    )
    parser.add_argument(
        "--parallel",
        type=lambda x: x.lower() in ("true", "1", "yes"),
        default=False,
        metavar="True|False",
        help="Enable parallel processing (default: False)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel workers (default: 4)",
    )
    parser.add_argument(
        "--inscode_mode",
        type=int,
        choices=[1, 2, 3],
        default=1,
        help=(
            "Insertion code handling: "
            "1=include in order (default), "
            "2=remove all, "
            "3=remove only if His-tag (>=4 consecutive HIS)"
        ),
    )
    args = parser.parse_args()

    results, df = run(
        input_file=args.input_file,
        pdb_col=args.pdb_col,
        output_dir=args.output_dir,
        report=args.report,
        parallel=args.parallel,
        workers=args.workers,
        inscode_mode=args.inscode_mode,
    )

    success = sum(1 for r in results if r["status"] == "success")
    failed = len(results) - success
    logger.info("Done: %d succeeded, %d failed. Report: %s", success, failed, args.report)


if __name__ == "__main__":
    main()

