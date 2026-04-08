# biodata-utils

Utility scripts for biological data retrieval and processing.

---

## pdb_fasta_extractor.py

PDB 파일에서 체인별 아미노산 서열을 추출하여 FASTA 파일로 저장하는 도구입니다.
비연속 잔기 번호는 `-`로 표기하며, alternate conformation 및 insertion code 처리를 지원합니다.

### 요구사항

- Python 3.12+
- Biopython
- pandas

```bash
pip install biopython pandas
```

### 사용법

```bash
python pdb_fasta_extractor.py <input_file> [options]
```

#### 인수

| 인수 | 필수 여부 | 설명 |
|------|----------|------|
| `input_file` | 필수 | PDB 이름 목록 파일 (한 줄에 하나) 또는 `--pdb_col` 지정 시 테이블 파일 |
| `--pdb_col` | 선택 | 테이블 입력 시 PDB ID가 있는 컬럼명 |
| `--output_dir` | 선택 | FASTA 출력 디렉토리 (기본값: 현재 디렉토리) |
| `--report` | 선택 | 리포트 CSV 파일 경로 (기본값: `pdb_seq_report.csv`) |
| `--parallel` | 선택 | 병렬 처리 활성화 `True\|False` (기본값: `False`) |
| `--workers` | 선택 | 병렬 처리 워커 수 (기본값: `4`) |
| `--inscode_mode` | 선택 | Insertion code 처리 방식 `1\|2\|3` (기본값: `1`) |

#### `--inscode_mode` 옵션

| 값 | 동작 |
|----|------|
| `1` | insertion code 잔기를 서열 순서대로 포함 (기본값) |
| `2` | 모든 insertion code 잔기 제거 |
| `3` | insertion code 잔기가 His-tag (연속 HIS ≥ 4개)를 구성하는 경우에만 제거 |

### 실행 예시

```bash
# 기본 실행 (텍스트 목록 파일)
python pdb_fasta_extractor.py pdb_list.txt

# 출력 디렉토리 지정
python pdb_fasta_extractor.py pdb_list.txt --output_dir ./fasta_out

# 테이블 파일에서 특정 컬럼 사용
python pdb_fasta_extractor.py metadata.tsv --pdb_col pdb_id

# 병렬 처리
python pdb_fasta_extractor.py pdb_list.txt --parallel True --workers 8

# Insertion code 잔기 제거
python pdb_fasta_extractor.py pdb_list.txt --inscode_mode 2

# 리포트 경로 지정
python pdb_fasta_extractor.py pdb_list.txt --report results/report.csv
```

### 입력 파일 형식

**텍스트 목록 파일** (`pdb_list.txt`)
```
1ABC
2XYZ.pdb
/path/to/3DEF.pdb
```

**테이블 파일** (TSV/CSV, `--pdb_col` 사용 시)
```
pdb_id  sample
1ABC    sample1
2XYZ    sample2
```

### 출력

#### FASTA 파일 (`{pdb_name}.fasta`)

각 PDB 파일당 하나의 FASTA 파일이 생성되며, 체인별로 레코드가 구성됩니다.

```
>1ABC_A chain=A | residues=1-150 | res_count=148 | gaps=10-11
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
```

**헤더 필드 설명**

| 필드 | 설명 |
|------|------|
| `chain` | 체인 ID |
| `residues` | 시작-끝 잔기 번호 |
| `res_count` | 실제 아미노산 수 (갭 제외) |
| `altconf` | alternate conformation이 있는 잔기 번호 목록 |
| `inscodes` | 유지된 insertion code 잔기 목록 |
| `removed_inscodes` | 제거된 insertion code 잔기 목록 |
| `gaps` | 비연속 구간 번호 (범위 표기) |

#### 리포트 CSV (`pdb_seq_report.csv`)

| 컬럼 | 설명 |
|------|------|
| `pdb` | PDB 이름 |
| `chain` | 체인 ID |
| `seq_length` | 서열 길이 (갭 포함) |
| `gap_count` | 갭(`-`) 수 |
| `has_gap` | 갭 존재 여부 |
| `altconf_count` | alternate conformation 잔기 수 |
| `kept_inscodes` | 유지된 insertion code 잔기 |
| `removed_inscodes` | 제거된 insertion code 잔기 |
| `gap_positions` | 갭 위치 (범위 표기) |
| `output_file` | 출력 FASTA 파일 경로 |
| `status` | `success` 또는 `error: <메시지>` |
